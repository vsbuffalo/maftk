// statistics.rs
use csv::{Writer, WriterBuilder};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use thiserror::Error;

use crate::binary::SpeciesDictionary;
use crate::io::OutputFile;

#[derive(Error, Debug)]
pub enum StatsError {
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("CSV error: {0}")]
    CsvError(#[from] csv::Error),
}

#[derive(Clone, Default, Debug)]
pub struct PairwiseStats {
    substitutions: u32,
    gaps: u32,
    valid_positions: u32,
}

impl PairwiseStats {
    fn substitution_rate(&self) -> f64 {
        if self.valid_positions == 0 {
            0.0
        } else {
            self.substitutions as f64 / self.valid_positions as f64
        }
    }

    fn gap_rate(&self) -> f64 {
        if self.valid_positions == 0 && self.gaps == 0 {
            0.0
        } else {
            self.gaps as f64 / (self.valid_positions + self.gaps) as f64
        }
    }
}

#[derive(Clone, Default, Debug)]
pub struct RegionStats {
    chrom: String,
    start: u32,
    end: u32,
    pairwise_stats: HashMap<(u32, u32), PairwiseStats>,
}

pub struct AlignmentStatistics {
    writer: Writer<Box<dyn Write>>,
    species: HashSet<String>,
    // Store the species indices in the order they appear in headers
    species_indices: Vec<u32>,
}

impl AlignmentStatistics {
    pub fn new(output: OutputFile, species: HashSet<String>) -> Result<Self, StatsError> {
        let writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(output.writer()?);

        // We'll fill species_indices when we write the header
        let species_indices = Vec::new();

        let mut stats = Self {
            writer,
            species,
            species_indices,
        };
        stats.write_header()?;
        Ok(stats)
    }

    fn write_header(&mut self) -> Result<(), StatsError> {
        let mut headers: Vec<String> =
            vec!["chrom".to_string(), "start".to_string(), "end".to_string()];

        // Sort species names for consistent output
        let mut species_vec: Vec<_> = self.species.iter().collect();
        species_vec.sort();

        for i in 0..species_vec.len() {
            for j in (i + 1)..species_vec.len() {
                headers.push(format!("{}_{}_subst_rate", species_vec[i], species_vec[j]));
                headers.push(format!("{}_{}_gap_rate", species_vec[i], species_vec[j]));
            }
        }

        self.writer.write_record(&headers)?;
        Ok(())
    }

    pub fn write_stats(
        &mut self,
        stats: &RegionStats,
        species_dict: &SpeciesDictionary,
    ) -> Result<(), StatsError> {
        let mut record = vec![
            stats.chrom.clone(),
            stats.start.to_string(),
            stats.end.to_string(),
        ];

        // Sort species for consistent ordering
        let mut species_vec: Vec<_> = self.species.iter().collect();
        species_vec.sort();

        // Convert species names to indices
        let species_indices: Vec<_> = species_vec
            .iter()
            .filter_map(|name| species_dict.species_to_index.get(*name).copied())
            .collect();

        // Write stats for each species pair
        for i in 0..species_indices.len() {
            for j in (i + 1)..species_indices.len() {
                let idx1 = species_indices[i];
                let idx2 = species_indices[j];

                // Create the species pair in sorted order
                let pair = if idx1 < idx2 {
                    (idx1, idx2)
                } else {
                    (idx2, idx1)
                };

                match stats.pairwise_stats.get(&pair) {
                    Some(pair_stats) => {
                        record.push(pair_stats.substitution_rate().to_string());
                        record.push(pair_stats.gap_rate().to_string());
                    }
                    None => {
                        record.push("NA".to_string());
                        record.push("NA".to_string());
                    }
                }
            }
        }

        self.writer.write_record(&record)?;
        Ok(())
    }
}

pub fn calc_alignment_statistics(
    alignments: Vec<(u32, String)>,
    species_indices: &HashSet<u32>,
    chrom: String,
    start: u32,
    end: u32,
) -> RegionStats {
    let mut stats = RegionStats {
        chrom,
        start,
        end,
        pairwise_stats: HashMap::new(),
    };

    for i in 0..alignments.len() {
        let (sp1_idx, seq1) = &alignments[i];
        if !species_indices.contains(sp1_idx) {
            continue;
        }

        for alignment in alignments.iter().skip(i + 1) {
            let (sp2_idx, seq2) = alignment;
            if !species_indices.contains(sp2_idx) {
                continue;
            }

            let mut pair_stats = PairwiseStats::default();

            // Rest of sequence comparison logic stays the same...
            for (c1, c2) in seq1.chars().zip(seq2.chars()) {
                match (c1, c2) {
                    ('.', '.') | ('-', '-') => continue,
                    ('.', _) | ('-', _) | (_, '.') | (_, '-') => pair_stats.gaps += 1,
                    (b1, b2) if b1 != b2 => {
                        pair_stats.substitutions += 1;
                        pair_stats.valid_positions += 1;
                    }
                    _ => pair_stats.valid_positions += 1,
                }
            }

            let species_pair = if sp1_idx < sp2_idx {
                (*sp1_idx, *sp2_idx)
            } else {
                (*sp2_idx, *sp1_idx)
            };
            stats.pairwise_stats.insert(species_pair, pair_stats);
        }
    }

    stats
}

pub fn calc_alignment_statistics_bytes(
    alignments: Vec<(u32, String)>,
    species_indices: &HashSet<u32>,
    chrom: String,
    start: u32,
    end: u32,
) -> RegionStats {
    let mut stats = RegionStats {
        chrom,
        start,
        end,
        pairwise_stats: HashMap::new(),
    };

    for i in 0..alignments.len() {
        let (sp1_idx, seq1) = &alignments[i];
        if !species_indices.contains(sp1_idx) {
            continue;
        }
        // Convert to bytes once
        let seq1_bytes = seq1.as_bytes();

        for alignment in alignments.iter().skip(i + 1) {
            let (sp2_idx, seq2) = alignment;
            if !species_indices.contains(sp2_idx) {
                continue;
            }
            let seq2_bytes = seq2.as_bytes();

            let mut pair_stats = PairwiseStats::default();

            // Use byte comparison instead of char comparison
            for (b1, b2) in seq1_bytes.iter().zip(seq2_bytes.iter()) {
                match (b1, b2) {
                    (b'.' | b'-', b'.' | b'-') => continue,
                    (b'.' | b'-', _) | (_, b'.' | b'-') => pair_stats.gaps += 1,
                    (b1, b2) => {
                        // Simple ASCII uppercase comparison
                        let b1_upper = if *b1 >= b'a' && *b1 <= b'z' {
                            b1 - 32
                        } else {
                            *b1
                        };
                        let b2_upper = if *b2 >= b'a' && *b2 <= b'z' {
                            b2 - 32
                        } else {
                            *b2
                        };
                        if b1_upper != b2_upper {
                            pair_stats.substitutions += 1;
                            pair_stats.valid_positions += 1;
                        } else {
                            pair_stats.valid_positions += 1;
                        }
                    }
                }
            }

            let species_pair = if sp1_idx < sp2_idx {
                (*sp1_idx, *sp2_idx)
            } else {
                (*sp2_idx, *sp1_idx)
            };
            stats.pairwise_stats.insert(species_pair, pair_stats);
        }
    }
    stats
}
