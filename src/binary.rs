// binary.rs

use colored::*;
use core::panic;
use csv::ReaderBuilder;
use flate2::write::{DeflateDecoder, DeflateEncoder};
use flate2::Compression;
use glob::glob;
use indicatif::{ProgressBar, ProgressStyle};
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::time::Instant;
use std::{collections::HashSet, path::Path};
use tracing::{info, warn};

use hgindex::{GenomicCoordinates, GenomicDataStore};

use crate::io::{InputFile, MafReader, OutputFile, Strand};
use crate::statistics::{
    calc_alignment_block_statistics, compare_bases, is_gap, AlignmentStatistics, RegionStats,
};

// Species dictionary to avoid storing repetitive strings
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpeciesDictionary {
    // Maps species name to its index
    pub(crate) species_to_index: HashMap<String, u32>,
    // Maps index back to species name
    index_to_species: Vec<String>,
}

impl Default for SpeciesDictionary {
    fn default() -> Self {
        Self::new()
    }
}

impl SpeciesDictionary {
    pub fn new() -> Self {
        Self {
            species_to_index: HashMap::new(),
            index_to_species: Vec::new(),
        }
    }

    pub fn get_or_insert(&mut self, species: String) -> u32 {
        if let Some(&idx) = self.species_to_index.get(&species) {
            idx
        } else {
            let idx = self.index_to_species.len() as u32;
            self.index_to_species.push(species.clone());
            self.species_to_index.insert(species, idx);
            idx
        }
    }

    pub fn get_species(&self, idx: u32) -> Option<&str> {
        self.index_to_species.get(idx as usize).map(|s| s.as_str())
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CompressedText(Vec<u8>);

impl CompressedText {
    pub fn new(text: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(text.as_bytes())?;
        Ok(CompressedText(encoder.finish()?))
    }

    pub fn decompress(&self) -> Result<String, Box<dyn std::error::Error>> {
        let mut decoder = DeflateDecoder::new(Vec::new());
        decoder.write_all(&self.0)?;
        let bytes = decoder.finish()?;
        String::from_utf8(bytes).map_err(|e| e.into())
    }

    pub fn compression_ratio(&self) -> f64 {
        self.0.len() as f64 / self.decompress().map(|s| s.len()).unwrap_or(1) as f64
    }
}

/// Compact metadata for each sequence in the alignment
///
/// Note: this is somewhat like Sequence in io.rs, expect it is designed to be a type more
/// optimized for binary serialization. In particular, the species (src in the lower structure) is
/// now just a u32 index.
#[derive(Debug, Serialize, Deserialize)]
pub struct AlignedSequence {
    pub species_idx: u32, // Index into species dictionary
    pub start: u32,       // Start position in source sequence
    pub size: u32,        // Size of aligning region
    pub strand: bool,     // false for forward (+), true for reverse (-)
    pub src_size: u32,    // Total size of source sequence
}

// Modified MafBlock to include metadata
#[derive(Debug, Serialize, Deserialize)]
pub struct MafBlock {
    pub score: f32,
    pub sequence_metadata: Vec<AlignedSequence>,
    compressed_text: CompressedText,
    sequence_length: usize,
}

impl GenomicCoordinates for MafBlock {
    fn start(&self) -> u32 {
        let ref_seq = &self.sequence_metadata[0];
        ref_seq.start
    }

    fn end(&self) -> u32 {
        let ref_seq = &self.sequence_metadata[0];
        ref_seq.start + ref_seq.size
    }
}

impl MafBlock {
    pub fn new(
        score: f32,
        sequences: Vec<(AlignedSequence, String)>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        if sequences.is_empty() {
            return Err("Cannot create MafBlock with no sequences".into());
        }

        let sequence_length = sequences[0].1.len();

        // Combine all sequences into one string with newlines
        let combined_text = sequences
            .iter()
            .map(|(_, text)| text.as_str())
            .collect::<Vec<_>>()
            .join("\n");

        // Extract just the metadata
        let sequence_metadata = sequences.into_iter().map(|(meta, _)| meta).collect();

        Ok(Self {
            score,
            sequence_metadata,
            compressed_text: CompressedText::new(&combined_text)?,
            sequence_length,
        })
    }

    /// Get the reference aligned sequence (first
    /// in the sequences attribute).
    pub fn get_reference_region(&self, species_dict: &SpeciesDictionary) -> (u32, u32) {
        let ref_seq = self.sequence_metadata.first().unwrap();
        let _ref_species = species_dict.get_species(ref_seq.species_idx);
        (ref_seq.start, ref_seq.start + ref_seq.size)
    }

    // Get a specific sequence by index
    pub fn get_sequence(&self, idx: usize) -> Result<String, Box<dyn std::error::Error>> {
        let text = self.compressed_text.decompress()?;
        let lines: Vec<&str> = text.split('\n').collect();
        if idx >= lines.len() {
            return Err("Sequence index out of bounds".into());
        }
        Ok(lines[idx].to_string())
    }

    // Get all sequences
    pub fn get_sequences(&self) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let text = self.compressed_text.decompress()?;
        Ok(text.split('\n').map(|s| s.to_string()).collect())
    }

    pub fn calc_stats(
        &self,
        start: Option<u32>,
        end: Option<u32>,
        species_indices: Option<&HashSet<u32>>,
    ) -> Option<RegionStats> {
        calc_alignment_block_statistics(self, species_indices, start, end)
    }

    pub fn pretty_print_alignments(&self, species_dict: &SpeciesDictionary, color: bool) {
        // Get all sequences at once
        let sequences = match self.get_sequences() {
            Ok(seqs) => seqs,
            Err(_) => return,
        };

        let alignments: Vec<(String, String, u32, bool)> = self
            .sequence_metadata
            .iter()
            .zip(sequences.iter())
            .enumerate()
            .map(|(idx, (meta, text))| {
                let label = &format!("seq_{}", idx);
                let species_name = species_dict.get_species(meta.species_idx).unwrap_or(label);
                (
                    species_name.to_string(),
                    text.to_string(),
                    meta.start,
                    meta.strand,
                )
            })
            .collect();

        if alignments.is_empty() {
            return;
        }

        // Get reference sequence (first sequence)
        let ref_alignment = &alignments[0];
        let reference = &ref_alignment.1;

        // Calculate the maximum identifier length for padding
        let max_id_len = alignments
            .iter()
            .map(|(id, _, _, _)| id.len())
            .max()
            .unwrap_or(0);

        // Print header with sequence positions
        print!("{:width$} ", "", width = max_id_len + 15); // Increased padding for coordinates
        for (i, _) in reference.chars().enumerate() {
            if i % 10 == 0 {
                print!("{:<10}", i);
            }
        }
        println!();

        // Print each sequence with optional color coding
        for (id, sequence, start_pos, is_reverse) in &alignments {
            // Calculate and format coordinates
            let strand = if *is_reverse { "-" } else { "+" };
            let coord_label = format!(
                "{:<10} {:>10} {}",
                id.chars().collect::<String>(),
                start_pos,
                strand
            );

            // Print identifier and coordinates with padding
            print!("{:width$} ", coord_label, width = max_id_len + 15);

            // For the reference sequence, just print it normally
            if id == &alignments[0].0 {
                println!("{}", sequence);
                continue;
            }

            if color {
                // Compare with reference sequence and color-code
                for (ref_byte, query_byte) in
                    reference.as_bytes().iter().zip(sequence.as_bytes().iter())
                {
                    let colored_base = if is_gap(*ref_byte) && is_gap(*query_byte) {
                        (*query_byte as char).to_string().green()
                    } else if is_gap(*ref_byte) || is_gap(*query_byte) {
                        (*query_byte as char).to_string().yellow()
                    } else if compare_bases(*ref_byte, *query_byte) {
                        (*query_byte as char).to_string().green()
                    } else {
                        (*query_byte as char).to_string().red()
                    };
                    print!("{}", colored_base);
                }
            } else {
                print!("{}", sequence);
            }
            println!();
        }
    }

    pub fn format_fasta(
        &self,
        species_dict: &SpeciesDictionary,
        chrom: &str,
    ) -> Result<String, Box<dyn std::error::Error>> {
        let sequences = self.get_sequences()?;
        let mut fasta_entries = Vec::new();

        for (meta, sequence) in self.sequence_metadata.iter().zip(sequences.iter()) {
            let species = species_dict
                .get_species(meta.species_idx)
                .unwrap_or("unknown");
            let strand = if meta.strand { "-" } else { "+" };

            fasta_entries.push(format!(
                ">{}.{} {}:{}-{} {}\n{}",
                species,
                chrom,
                chrom,
                meta.start,
                meta.start + meta.size,
                strand,
                sequence
            ));
        }

        Ok(fasta_entries.join("\n") + "\n")
    }

    /// Write the alignment blocks to a FASTA file.
    pub fn write_fasta(
        &self,
        path: &Path,
        species_dict: &SpeciesDictionary,
        chrom: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = OutputFile::new(path).writer()?;

        // Write the formatted FASTA string to the file
        let fasta_string = self.format_fasta(species_dict, chrom)?;
        writeln!(writer, "{}", fasta_string)?;

        Ok(())
    }
}

pub fn query_alignments(
    data_dir: &Path,
    chrom: &str,
    start: u32,
    end: u32,
    intersect_only: bool,
) -> Result<(Vec<MafBlock>, SpeciesDictionary), Box<dyn std::error::Error>> {
    let total_time = Instant::now();
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::open(data_dir, None)?;
    let mut blocks: Vec<MafBlock> = store.get_overlapping(chrom, start, end)?;
    let query_time = total_time.elapsed();

    // Get the store's metadata for species names
    let species_dict = store.metadata().unwrap();

    if intersect_only {
        // Filter blocks to only keep sequences that overlap in all blocks
        let mut common_species: Option<HashSet<u32>> = None;

        for block in &blocks {
            let block_species: HashSet<u32> = block
                .sequence_metadata
                .iter()
                .map(|seq| seq.species_idx)
                .collect();

            common_species = match common_species {
                None => Some(block_species),
                Some(species) => Some(&species & &block_species),
            };
        }

        if let Some(common_species) = common_species {
            blocks = blocks
                .into_iter()
                .map(|mut block| {
                    block
                        .sequence_metadata
                        .retain(|seq| common_species.contains(&seq.species_idx));
                    block
                })
                .collect();
        }
    }

    info!(
        "Query elapsed {:?}, {} overlapping alignment block{} found.",
        query_time,
        blocks.len(),
        if blocks.len() > 1 { "s" } else { "" }
    );
    Ok((blocks, species_dict.clone()))
}

pub fn calc_block_statistics(
    stats: &RegionStats,
    species_dict: &SpeciesDictionary,
    focal_species: Option<&str>,
) -> Result<DataFrame, Box<dyn std::error::Error>> {
    // Create vectors to store data for DataFrame
    let mut species1 = Vec::new();
    let mut species2 = Vec::new();
    let mut subst_rates = Vec::new();
    let mut gap_rates = Vec::new();
    let mut num_subst = Vec::new();
    let mut num_single_gaps = Vec::new();
    let mut num_double_gaps = Vec::new();
    let mut valid_positions = Vec::new();
    let mut total_positions = Vec::new();

    // Collect data for each species pair
    for ((sp1_idx, sp2_idx), pair_stats) in &stats.pairwise_stats {
        let sp1_name = species_dict.get_species(*sp1_idx).unwrap_or("unknown");
        let sp2_name = species_dict.get_species(*sp2_idx).unwrap_or("unknown");

        if let Some(focal) = focal_species {
            if sp1_name != focal && sp2_name != focal {
                continue;
            }
        }

        species1.push(sp1_name.to_string());
        species2.push(sp2_name.to_string());
        subst_rates.push(pair_stats.substitution_rate() * 100.0);
        gap_rates.push(pair_stats.gap_rate() * 100.0);
        num_subst.push(pair_stats.substitutions);
        num_single_gaps.push(pair_stats.single_gaps);
        num_double_gaps.push(pair_stats.double_gaps);
        valid_positions.push(pair_stats.valid_positions());
        total_positions.push(pair_stats.total_positions);
    }

    // Create DataFrame using df! macro
    let mut df = df! {
        "Species 1" => &species1,
        "Species 2" => &species2,
        "Substitution Rate %" => &subst_rates,
        "Gap Rate %" => &gap_rates,
        "Num Substitutions" => &num_subst,
        "Num Single Gaps" => &num_single_gaps,
        "Num Double Gaps" => &num_double_gaps,
        "Valid Positions" => &valid_positions,
        "Total Positions" => &total_positions,
    }?;

    // Sort by substitution rate (descending)
    df = df
        .lazy()
        .sort_by_exprs(
            vec![col("Substitution Rate %")],
            SortMultipleOptions::default(),
        )
        .collect()?;

    Ok(df)
}

pub fn print_block_statistics(df: &DataFrame, block_idx: usize) {
    println!("Block {} Statistics:", block_idx);
    unsafe {
        if std::env::var("POLARS_FMT_MAX_COLS").is_err() {
            std::env::set_var("POLARS_FMT_MAX_COLS", "-1");
        }
        if std::env::var("POLARS_FMT_MAX_ROWS").is_err() {
            std::env::set_var("POLARS_FMT_MAX_ROWS", "-1");
        }
    }
    // Print formatted table
    println!("\nAlignment Statistics (sorted by substitution rate):");
    println!("{}", df);
}

pub fn print_alignments(blocks: Vec<MafBlock>, species_dict: &SpeciesDictionary, color: bool) {
    let has_blocks = !blocks.is_empty();
    for (i, block) in blocks.into_iter().enumerate() {
        println!("[block {} | score {}]", i, block.score);
        block.pretty_print_alignments(species_dict, color);
    }
    if has_blocks {
        // Print legend
        println!("\nLegend:");
        println!("{} Match", "A".green());
        println!("{} Mismatch", "A".red());
        println!("{} Indel", "-".yellow());
        println!("Coordinates shown as: species start_position [strand]");
    }
}

/// Estimates number of records in a TSV/CSV file
fn estimate_record_count(path: &Path) -> Result<u64, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let file_size = file.metadata()?.len();

    // Sample first few lines to get average line length
    let input_file = InputFile::new(path);
    let mut reader = input_file.reader()?;
    let mut buffer = String::new();
    let mut total_bytes = 0;
    let mut line_count = 0;
    const SAMPLE_LINES: usize = 100;

    while line_count < SAMPLE_LINES {
        buffer.clear();
        let bytes_read = reader.read_line(&mut buffer)?;
        if bytes_read == 0 {
            break;
        } // EOF

        // Skip comments
        if buffer.starts_with('#') {
            continue;
        }

        total_bytes += bytes_read;
        line_count += 1;
    }

    if line_count == 0 {
        return Ok(0);
    }

    // Calculate average bytes per record
    let avg_bytes_per_record = total_bytes as f64 / line_count as f64;
    let estimated_records = (file_size as f64 / avg_bytes_per_record).ceil() as u64;

    Ok(estimated_records)
}

pub fn stats_command(
    regions: &Path,
    output: &Path,
    species: HashSet<String>,
    data_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let total_time = Instant::now();

    // Setup progress bar
    let estimated_records = estimate_record_count(regions)?;
    let progress = ProgressBar::new(estimated_records);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} records ({eta})")?
        .progress_chars("=>-"));

    // Setup input/output
    let input_file = InputFile::new(regions);
    let buf_reader = input_file.reader()?;
    let output_file = OutputFile::new(output);
    let mut stats_writer = AlignmentStatistics::new(output_file, species.clone())?;

    // Create CSV reader for BED file
    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(buf_reader);

    // Open store and get species dictionary
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::open(data_dir, None)?;
    let species_dict = store
        .metadata()
        .expect("Missing species dictionary")
        .clone();

    // Convert species names to indices
    let species_indices: HashSet<u32> = species
        .iter()
        .filter_map(|name| species_dict.species_to_index.get(name).copied())
        .collect();

    // Process each region
    for result in reader.records() {
        let record = result?;
        let chrom = record[0].to_string();
        let start: u32 = record[1].parse()?;
        let end: u32 = record[2].parse()?;

        // Get overlapping blocks
        let blocks: Vec<MafBlock> = store.get_overlapping(&chrom, start, end)?;

        if blocks.is_empty() {
            continue;
        }

        // Process each block independently
        for block in blocks {
            if let Some(mut block_stats) =
                block.calc_stats(Some(start), Some(end), Some(&species_indices))
            {
                block_stats.chrom = chrom.clone(); // Set chromosome
                stats_writer.write_stats(&block_stats, &species_dict, &block)?;
            } else {
                panic!("No statistics returned due to no overlaps: internal error.");
            };
        }

        progress.inc(1);
    }

    progress.finish_with_message("Processing complete");
    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}

/// Helper function to convert a block from MAF format to binary format
fn convert_maf_block_to_binary(
    block: &crate::io::AlignmentBlock,
    store: &mut GenomicDataStore<MafBlock, SpeciesDictionary>,
    min_length: u64,
    compression_stats: &mut Option<&mut CompressionStats>,
) -> Result<(), Box<dyn Error>> {
    if block.sequences[0].size < min_length {
        return Ok(());
    }

    let ref_seq = &block.sequences[0];
    let chr = ref_seq.src.split('.').nth(1).ok_or("Invalid ref name")?;

    // Convert sequences to metadata + alignment strings pairs
    let mut sequences = Vec::new();

    for seq in &block.sequences {
        let species = seq.src.split('.').next().unwrap_or("unknown").to_string();
        let species_idx = if let Some(metadata) = store.metadata_mut() {
            metadata.get_or_insert(species)
        } else {
            return Err("Missing metadata".into());
        };

        sequences.push((
            AlignedSequence {
                species_idx,
                start: seq.start as u32,
                size: seq.size as u32,
                strand: matches!(seq.strand, Strand::Reverse),
                src_size: seq.src_size as u32,
            },
            seq.text.clone(),
        ));
    }

    let score = block.score.unwrap_or(0.0) as f32;
    let record = MafBlock::new(score, sequences)?;

    // Sample compression ratio if requested
    if let Some(stats) = compression_stats {
        stats.add_ratio(record.compressed_text.compression_ratio());
    }

    store.add_record(chr, &record)?;

    Ok(())
}

/// Stores compression statistics for reporting
struct CompressionStats {
    ratios: Vec<f64>,
    sample_size: usize,
    max_samples: usize,
}

impl CompressionStats {
    fn new(max_samples: usize) -> Self {
        Self {
            ratios: Vec::with_capacity(max_samples),
            sample_size: 0,
            max_samples,
        }
    }

    fn add_ratio(&mut self, ratio: f64) {
        if self.sample_size < self.max_samples {
            self.ratios.push(ratio);
            self.sample_size += 1;
        }
    }

    fn report(&self) {
        if self.ratios.is_empty() {
            info!("No compression ratios collected");
            return;
        }

        let sum: f64 = self.ratios.iter().sum();
        let avg = sum / self.ratios.len() as f64;
        let min = self.ratios.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max = self.ratios.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        info!(
            "Compression Statistics (sample size: {}):",
            self.sample_size
        );
        info!("  Average ratio: {:.2}%", avg * 100.0);
        info!("  Min ratio: {:.2}%", min * 100.0);
        info!("  Max ratio: {:.2}%", max * 100.0);
    }
}

/// Helper function to initialize the store and metadata
fn init_store(
    output_dir: &Path,
) -> Result<GenomicDataStore<MafBlock, SpeciesDictionary>, Box<dyn Error>> {
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::create(output_dir, None)?;
    store.set_metadata(SpeciesDictionary::new());
    Ok(store)
}

/// Helper function to process a single MAF file
fn process_maf_file(
    input: &Path,
    store: &mut GenomicDataStore<MafBlock, SpeciesDictionary>,
    min_length: u64,
    progress: Option<&ProgressBar>,
    compression_stats: &mut Option<&mut CompressionStats>,
) -> Result<(), Box<dyn Error>> {
    let mut maf = MafReader::from_file(input)?;
    let _header = maf.read_header()?;

    while let Some(block) = maf.next_block()? {
        convert_maf_block_to_binary(&block, store, min_length, compression_stats)?;
    }

    if let Some(pb) = progress {
        pb.inc(1);
    }

    Ok(())
}

/// Convert a single MAF file to binary format
pub fn convert_to_binary(
    input: &Path,
    output_dir: &Path,
    min_length: u64,
) -> Result<(), Box<dyn Error>> {
    let total_time = Instant::now();

    let mut store = init_store(output_dir)?;
    let mut compression_stats = CompressionStats::new(1000); // Sample up to 1000 blocks
    process_maf_file(
        input,
        &mut store,
        min_length,
        None,
        &mut Some(&mut compression_stats),
    )?;
    store.finalize()?;

    compression_stats.report();
    compression_stats.report();
    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}

/// Convert multiple MAF files matching a glob pattern to binary format
pub fn convert_to_binary_glob(
    input_pattern: &str,
    output_dir: &Path,
    min_length: u64,
) -> Result<(), Box<dyn Error>> {
    let total_time = Instant::now();

    // Collect and sort input files
    let mut input_files: Vec<PathBuf> = glob(input_pattern)?.filter_map(Result::ok).collect();

    if input_files.is_empty() {
        warn!("No files found matching pattern: {}", input_pattern);
        return Ok(());
    }

    input_files.sort();
    info!("Found {} files to process", input_files.len());

    // Setup progress tracking
    let progress = ProgressBar::new(input_files.len() as u64);
    progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({msg})",
            )?
            .progress_chars("=>-"),
    );

    // Initialize store
    let mut store = init_store(output_dir)?;
    let mut compression_stats = CompressionStats::new(1000); // Sample up to 1000 blocks

    // Process each file
    for input_file in input_files {
        let file_name = input_file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        progress.set_message(format!("Processing {}", file_name));

        process_maf_file(
            &input_file,
            &mut store,
            min_length,
            Some(&progress),
            &mut Some(&mut compression_stats),
        )?;
    }

    // Finalize and cleanup
    progress.finish_with_message("Processing complete");
    store.finalize()?;

    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}
