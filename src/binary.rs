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
use tracing::{info, trace};

use hgindex::GenomicDataStore;

use crate::io::{InputFile, MafReader, OutputFile, Strand};
use crate::statistics::{calc_alignment_statistics, AlignmentStatistics};

// Species dictionary to avoid storing repetitive strings
#[derive(Debug, Serialize, Deserialize)]
pub struct SpeciesDictionary {
    // Maps species name to its index
    species_to_index: HashMap<String, u32>,
    // Maps index back to species name
    index_to_species: Vec<String>,
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
pub struct SequenceMetadata {
    pub species_idx: u32, // Index into species dictionary
    pub start: u32,       // Start position in source sequence
    pub size: u32,        // Size of aligning region
    pub strand: bool,     // false for forward (+), true for reverse (-)
    pub src_size: u32,    // Total size of source sequence
}

// Complete metadata for an alignment block
#[derive(Debug, Serialize, Deserialize)]
pub struct AlignmentMetadata {
    pub score: f32,                       // Alignment score
    pub sequences: Vec<SequenceMetadata>, // Metadata for each sequence
}

// Modified MafBlock to include metadata
#[derive(Debug, Serialize, Deserialize)]
pub struct MafBlock {
    pub metadata: AlignmentMetadata,
    pub sequences: Vec<String>,
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
        let mut sequence_metadata = Vec::new();
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

            sequence_metadata.push(SequenceMetadata {
                species_idx,
                start: seq.start as u32,
                size: seq.size as u32,
                strand: matches!(seq.strand, Strand::Reverse),
                src_size: seq.src_size as u32,
            });
            alignment_strings.push(seq.text.clone());
        }

        let metadata = AlignmentMetadata {
            score: block.score.unwrap_or(0.0) as f32,
            sequences: sequence_metadata,
        };

        let record = MafBlock {
            metadata,
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

/// Pretty prints multiple sequence alignments with color-coded bases
pub fn pretty_print_alignments(
    alignments: &[(String, String)],
) -> Result<(), Box<dyn std::error::Error>> {
    if alignments.is_empty() {
        return Ok(());
    }

    // Get reference sequence (first sequence)
    let reference = &alignments[0].1;

    // Calculate the maximum identifier length for padding
    let max_id_len = alignments.iter().map(|(id, _)| id.len()).max().unwrap_or(0);

    // Print header with sequence positions
    print!("{:width$} ", "", width = max_id_len + 2);
    for (i, _) in reference.chars().enumerate() {
        if i % 10 == 0 {
            print!("{:<10}", i);
        }
    }
    println!();

    // Print each sequence with color coding
    for (id, sequence) in alignments {
        // Print identifier with padding
        print!("{:width$} ", id, width = max_id_len + 2);

        // For the reference sequence, just print it normally
        if id == &alignments[0].0 {
            println!("{}", sequence);
            continue;
        }

        // Compare with reference sequence and color-code
        for (ref_base, query_base) in reference.chars().zip(sequence.chars()) {
            let colored_base = match (ref_base, query_base) {
                (r, q) if r == q => query_base.to_string().green(),
                ('-', q) | (q, '-') => q.to_string().yellow(),
                _ => query_base.to_string().red(),
            };
            print!("{}", colored_base);
        }
        println!();
    }

    // Print legend
    println!("\nLegend:");
    println!("{} Match", "A".green());
    println!("{} Mismatch", "A".red());
    println!("{} Indel", "-".yellow());

    Ok(())
}

// Query overlapping alignments
pub fn query_alignments(
    data_dir: &Path,
    chrom: &str,
    start: u32,
    end: u32,
) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::open(data_dir, None)?;
    let blocks = store.get_overlapping(chrom, start - 1, end);

    let mut alignments = Vec::new();
    for block in blocks {
        for sequence in block.sequences {
            alignments.push((format!("seq_{}", alignments.len()), sequence));
        }
    }

    // Pretty print the alignments
    pretty_print_alignments(&alignments)?;

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

    // Estimate total records for progress bar
    let estimated_records = estimate_record_count(regions)?;
    let progress = ProgressBar::new(estimated_records);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} records ({eta})")?
        .progress_chars("#>-"));

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

    // Open genomic data store
    let mut store = GenomicDataStore::<MafBlock, SpeciesDictionary>::open(data_dir, None)?;

    // Get species dictionary from metadata
    let mut species_dict = SpeciesDictionary::new();

    // Process each region
    for result in reader.records() {
        let record = result?;
        let chrom = record[0].to_string();
        let start: u32 = record[1].parse()?;
        let end: u32 = record[2].parse()?;

        // Query overlapping alignments
        let find_start = Instant::now();
        let blocks = store.get_overlapping(&chrom, start - 1, end); // Convert to 0-based
        trace!(elapsed = ?find_start.elapsed(), "Finding blocks");

        if !blocks.is_empty() {
            // Convert blocks to alignments with proper species names
            let alignments: Vec<(String, String)> = blocks
                .into_iter()
                .flat_map(|block| {
                    block
                        .metadata
                        .sequences
                        .iter()
                        .zip(block.sequences.iter())
                        .filter_map(|(meta, seq)| {
                            // Get species name from dictionary
                            species_dict
                                .get_species(meta.species_idx)
                                .map(|species| (species.to_string(), seq.clone()))
                        })
                        .collect::<Vec<_>>()
                })
                .collect();

            let region_stats = calc_alignment_statistics(alignments, &species, chrom, start, end);
            stats_writer.write_stats(&region_stats)?;
        }
        progress.inc(1);
    }

    progress.finish_with_message("Processing complete");
    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
}
