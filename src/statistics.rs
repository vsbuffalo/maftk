// statistics.rs
use csv::{Writer, WriterBuilder};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use thiserror::Error;

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
    pairwise_stats: HashMap<(String, String), PairwiseStats>,
}

pub struct AlignmentStatistics {
    writer: Writer<Box<dyn Write>>,
    species: HashSet<String>,
}

impl AlignmentStatistics {
    pub fn new(output: OutputFile, species: HashSet<String>) -> Result<Self, StatsError> {
        let writer = WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(output.writer()?);

        let mut stats = Self { writer, species };
        stats.write_header()?;
        Ok(stats)
    }

    fn write_header(&mut self) -> Result<(), StatsError> {
        let mut headers: Vec<String> =
            vec!["chrom".to_string(), "start".to_string(), "end".to_string()];

        let species_vec: Vec<_> = self.species.iter().collect();
        for i in 0..species_vec.len() {
            for j in (i + 1)..species_vec.len() {
                headers.push(format!("{}_{}_subst_rate", species_vec[i], species_vec[j]));
                headers.push(format!("{}_{}_gap_rate", species_vec[i], species_vec[j]));
            }
        }

        self.writer.write_record(&headers)?;
        Ok(())
    }

    pub fn write_stats(&mut self, stats: &RegionStats) -> Result<(), StatsError> {
        let mut record = vec![
            stats.chrom.clone(),
            stats.start.to_string(),
            stats.end.to_string(),
        ];

        let species_vec: Vec<_> = self.species.iter().collect();
        for i in 0..species_vec.len() {
            for j in (i + 1)..species_vec.len() {
                let pair = (species_vec[i].to_string(), species_vec[j].to_string());
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
    alignments: Vec<(String, String)>,
    species: &HashSet<String>,
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
        let (sp1, seq1) = &alignments[i];
        if !species.contains(sp1) {
            continue;
        }

        for j in (i + 1)..alignments.len() {
            let (sp2, seq2) = &alignments[j];
            if !species.contains(sp2) {
                continue;
            }

            let mut pair_stats = PairwiseStats::default();

            for (c1, c2) in seq1.chars().zip(seq2.chars()) {
                match (c1, c2) {
                    ('.', '.') | ('-', '-') => continue,
                    ('.', _) | ('-', _) | (_, '.') | (_, '-') => pair_stats.gaps += 1,
                    (b1, b2) if b1.to_ascii_uppercase() != b2.to_ascii_uppercase() => {
                        pair_stats.substitutions += 1;
                        pair_stats.valid_positions += 1;
                    }
                    _ => pair_stats.valid_positions += 1,
                }
            }

            let species_pair = if sp1 < sp2 {
                (sp1.clone(), sp2.clone())
            } else {
                (sp2.clone(), sp1.clone())
            };
            stats.pairwise_stats.insert(species_pair, pair_stats);
        }
    }

    stats
}
