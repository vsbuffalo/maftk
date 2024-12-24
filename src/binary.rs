// binary.rs
use colored::*;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::time::Instant;
use std::{collections::HashSet, path::Path};
use tracing::info;

use hgindex::GenomicDataStore;

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

// Compact metadata for each sequence in the alignment
#[derive(Debug, Serialize, Deserialize)]
pub struct AlignedSequence {
    pub species_idx: u32, // Index into species dictionary
    pub start: u32,       // Start position in source sequence
    pub size: u32,        // Size of aligning region
    pub strand: bool,     // false for forward (+), true for reverse (-)
    pub src_size: u32,    // Total size of source sequence
    pub text: String,     // Raw sequence
}

// Modified MafBlock to include metadata
#[derive(Debug, Serialize, Deserialize)]
pub struct MafBlock {
    pub score: f32,
    pub sequences: Vec<AlignedSequence>,
}

impl MafBlock {
    pub fn calc_stats(
        &self,
        start: Option<u32>,
        end: Option<u32>,
        species_indices: Option<&HashSet<u32>>,
    ) -> Option<RegionStats> {
        calc_alignment_block_statistics(self, species_indices, start, end)
    }
    pub fn pretty_print_alignments(
        &self,
        species_dict: &SpeciesDictionary,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let alignments: Vec<(String, String, u32, bool)> = self
            .sequences
            .iter()
            .enumerate()
            .map(|(idx, aligned)| {
                let label = &format!("seq_{}", idx);
                let species_name = species_dict
                    .get_species(aligned.species_idx)
                    .unwrap_or(label);
                (
                    // species name
                    species_name.to_string(),
                    // aligned sequence
                    aligned.text.clone(),
                    // start position
                    aligned.start,
                    // end position
                    aligned.strand,
                )
            })
            .collect();

        if alignments.is_empty() {
            return Ok(());
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

        // Print each sequence with color coding
        for (id, sequence, start_pos, is_reverse) in &alignments {
            // Calculate and format coordinates
            let strand = if *is_reverse { "-" } else { "+" };
            let coord_label = format!(
                "{:<9} {:>10} {}",
                id.chars().take(8).collect::<String>(),
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
            println!();
        }

        Ok(())
    }
}

pub fn convert_to_binary(
    input: &Path,
    output_dir: &Path,
    min_length: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let total_time = Instant::now();
    let mut maf = MafReader::from_file(input)?;
    let _header = maf.read_header()?;

    // Create store with species dictionary
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::create(output_dir, None)?;

    store.set_metadata(SpeciesDictionary::new());

    while let Some(block) = maf.next_block()? {
        if block.sequences[0].size < min_length {
            continue;
        }

        let ref_seq = &block.sequences[0];
        let chr = ref_seq.src.split('.').nth(1).ok_or("Invalid ref name")?;

        // Convert sequences to metadata + alignment strings
        let mut alignment_strings = Vec::new();

        for seq in &block.sequences {
            let species = seq.src.split('.').next().unwrap_or("unknown").to_string();
            // Get mutable reference to metadata and update species dictionary
            let species_idx = if let Some(metadata) = store.metadata_mut() {
                metadata.get_or_insert(species)
            } else {
                // This shouldn't happen since we initialized with Some metadata
                return Err("Missing metadata".into());
            };

            alignment_strings.push(AlignedSequence {
                species_idx,
                start: seq.start as u32,
                size: seq.size as u32,
                strand: matches!(seq.strand, Strand::Reverse),
                src_size: seq.src_size as u32,
                text: seq.text.clone(),
            });
        }

        let score = block.score.unwrap_or(0.0) as f32;

        let record = MafBlock {
            score,
            sequences: alignment_strings,
        };

        store.add_record(
            chr,
            ref_seq.start as u32,
            (ref_seq.start + ref_seq.size) as u32,
            &record,
        )?;
    }

    // Finalize the store
    store.finalize()?;

    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}

// Query overlapping alignments
pub fn query_alignments(
    data_dir: &Path,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<Vec<MafBlock>, Box<dyn std::error::Error>> {
    let total_time = Instant::now();
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::open(data_dir, None)?;
    let blocks: Vec<MafBlock> = store.get_overlapping(chrom, start - 1, end);
    let query_time = total_time.elapsed();

    // Use the store's metadata to get species names
    let species_dict = store.metadata().unwrap();

    let mut alignments = Vec::new();
    for (i, block) in blocks.into_iter().enumerate() {
        println!("[block {} | score {}]", i, block.score);
        block.pretty_print_alignments(species_dict)?;
        //if let Some(stats) = block.calc_stats(None, None, None) {
        //    println!("{:}", stats);
        //}
        alignments.push(block);
    }
    // Print legend
    println!("\nLegend:");
    println!("{} Match", "A".green());
    println!("{} Mismatch", "A".red());
    println!("{} Indel", "-".yellow());
    println!("Coordinates shown as: species start_position [strand]");

    info!("Query elapsed {:?}", query_time);
    Ok(alignments)
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
        let blocks: Vec<MafBlock> = store.get_overlapping(&chrom, start, end);

        if blocks.is_empty() {
            continue;
        }

        // Process each block independently
        for block in blocks {
            if let Some(mut block_stats) =
                block.calc_stats(Some(start), Some(end), Some(&species_indices))
            {
                block_stats.chrom = chrom.clone(); // Set chromosome
                stats_writer.write_stats(&block_stats, &species_dict)?;
            };
        }

        progress.inc(1);
    }

    progress.finish_with_message("Processing complete");
    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}
