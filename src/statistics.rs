// statistics.rs
use csv::{Writer, WriterBuilder};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use thiserror::Error;

use crate::binary::{MafBlock, SpeciesDictionary};
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
    pub(crate) substitutions: u32,
    pub(crate) gaps: u32,
    pub(crate) valid_positions: u32,
}

impl std::fmt::Display for PairwiseStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:} subs, {:} gaps in {:} positions",
            self.substitutions, self.gaps, self.valid_positions
        )
    }
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
    pub(crate) chrom: String,
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) pairwise_stats: HashMap<(u32, u32), PairwiseStats>,
}

impl std::fmt::Display for RegionStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Region {}:{:}-{:}", self.chrom, self.start, self.end)?;

        for ((sp1, sp2), stats) in &self.pairwise_stats {
            writeln!(f, "  Species {}-{}: {}", sp1, sp2, stats)?;
        }

        Ok(())
    }
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

use bio_seq::prelude::*;

pub fn is_gap(base: Iupac) -> bool {
    base == Iupac::X // X is used for gaps in bio-seq
}

pub fn compare_bases(b1: Iupac, b2: Iupac) -> bool {
    if is_gap(b1) || is_gap(b2) {
        return false;
    }
    b1 == b2
}

pub fn calc_alignment_block_statistics(
    block: &MafBlock,
    species_indices: Option<&HashSet<u32>>,
    region_start: Option<u32>,
    region_end: Option<u32>,
) -> Option<RegionStats> {
    let ref_seq = &block.sequences[0];
    let block_start = ref_seq.start;
    let block_end = block_start + ref_seq.size;

    let (start_offset, length) = if let (Some(start), Some(end)) = (region_start, region_end) {
        if end <= block_start || start >= block_end {
            return None;
        }
        let overlap_start = start.max(block_start);
        let overlap_end = end.min(block_end);
        (
            (overlap_start - block_start) as usize,
            (overlap_end - overlap_start) as usize,
        )
    } else {
        (0, ref_seq.size as usize)
    };

    let mut stats = RegionStats {
        chrom: String::new(),
        start: block_start + start_offset as u32,
        end: block_start + start_offset as u32 + length as u32,
        pairwise_stats: HashMap::new(),
    };

    for i in 0..block.sequences.len() {
        let seq1 = &block.sequences[i];
        if species_indices
            .as_ref()
            .map_or(false, |indices| !indices.contains(&seq1.species_idx))
        {
            continue;
        }

        // Convert to vector so we can iterate multiple times
        let seq1_bases: Vec<Iupac> = seq1
            .text
            .as_ref()
            .iter()
            .skip(start_offset)
            .take(length)
            .collect();

        for j in (i + 1)..block.sequences.len() {
            let seq2 = &block.sequences[j];
            if species_indices
                .as_ref()
                .map_or(false, |indices| !indices.contains(&seq2.species_idx))
            {
                continue;
            }

            let seq2_bases: Vec<Iupac> = seq2
                .text
                .as_ref()
                .iter()
                .skip(start_offset)
                .take(length)
                .collect();

            let mut pair_stats = PairwiseStats::default();

            for (b1, b2) in seq1_bases.iter().zip(seq2_bases.iter()) {
                if is_gap(*b1) && is_gap(*b2) {
                    continue; // Matching gaps
                } else if is_gap(*b1) || is_gap(*b2) {
                    pair_stats.gaps += 1; // One gap
                } else {
                    pair_stats.valid_positions += 1;
                    if !compare_bases(*b1, *b2) {
                        pair_stats.substitutions += 1;
                    }
                }
            }

            let species_pair = if seq1.species_idx < seq2.species_idx {
                (seq1.species_idx, seq2.species_idx)
            } else {
                (seq2.species_idx, seq1.species_idx)
            };
            stats.pairwise_stats.insert(species_pair, pair_stats);
        }
    }

    Some(stats)
}

#[cfg(test)]
mod tests {
    use crate::{binary::AlignedSequence, io::standardize_to_iupac};

    use super::*;
    use std::collections::HashSet;

    // Helper function to create a MafBlock with two sequences
    fn create_test_block(seq1: &str, seq2: &str) -> MafBlock {
        let seq1 = standardize_to_iupac(seq1).unwrap();
        let seq2 = standardize_to_iupac(seq2).unwrap();
        MafBlock {
            score: 0.0,
            sequences: vec![
                AlignedSequence {
                    species_idx: 0,
                    start: 0,
                    size: seq1.len() as u32,
                    strand: false,
                    src_size: seq1.len() as u32,
                    text: seq1.try_into().expect("Invalid IUPAC sequence"),
                },
                AlignedSequence {
                    species_idx: 1,
                    start: 0,
                    size: seq2.len() as u32,
                    strand: false,
                    src_size: seq2.len() as u32,
                    text: seq2.try_into().expect("Invalid IUPAC sequence"),
                },
            ],
        }
    }

    #[test]
    fn test_perfect_match() {
        let block = create_test_block("ATCG", "ATCG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 0);
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_all_mismatches() {
        let block = create_test_block("ATCG", "TAGC");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 4);
        assert_eq!(pair_stats.gaps, 0);
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 1.0);
        assert_eq!(pair_stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_single_gap() {
        let block = create_test_block("A-TCG", "AATCG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 1);
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.2); // 1 gap out of 5 total positions
    }

    #[test]
    fn test_multiple_gaps() {
        let block = create_test_block("A-TC-G", "A-TC-G");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 0); // Matching gaps don't count
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_mixed_case() {
        let block = create_test_block("AtCg", "aTcG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0); // Should handle case-insensitive comparison
        assert_eq!(pair_stats.gaps, 0);
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_gaps_and_mismatches() {
        let block = create_test_block("A-TCG", "AT-CG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 2); // Two non-matching gaps
        assert_eq!(pair_stats.valid_positions, 3);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.4); // 2 gaps out of 5 total positions
    }

    #[test]
    fn test_dot_gaps() {
        let block = create_test_block("A.TCG", "AATCG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 1); // Dots count as gaps
        assert_eq!(pair_stats.valid_positions, 4);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.2);
    }

    #[test]
    fn test_empty_sequences() {
        let block = create_test_block("", "");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats =
            calc_alignment_block_statistics(&block, Some(&species_indices), None, None).unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 0);
        assert_eq!(pair_stats.gaps, 0);
        assert_eq!(pair_stats.valid_positions, 0);
        assert_eq!(pair_stats.substitution_rate(), 0.0);
        assert_eq!(pair_stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_mixed_case_stats() {
        let block = create_test_block("ATCGatcg", "atcgATCG");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats = block
            .calc_stats(None, None, Some(&species_indices))
            .unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(
            pair_stats.substitutions, 0,
            "Case difference should not count as substitution"
        );
        assert_eq!(pair_stats.valid_positions, 8);
    }

    #[test]
    fn test_mixed_case_with_gaps_stats() {
        let block = create_test_block("A-ccG", "AtC-g");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats = block
            .calc_stats(None, None, Some(&species_indices))
            .unwrap();

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(
            pair_stats.substitutions, 0,
            "Case difference should not count as substitution"
        );
        assert_eq!(pair_stats.gaps, 2);
        assert_eq!(pair_stats.valid_positions, 3);
    }

    #[test]
    fn test_mixed_dots_and_dashes() {
        let block = create_test_block("A.TcG", "AtC-g");
        let species_indices: HashSet<u32> = vec![0, 1].into_iter().collect();
        let stats = block
            .calc_stats(None, None, Some(&species_indices))
            .unwrap();
        dbg!(&stats);

        let pair_stats = stats.pairwise_stats.get(&(0, 1)).unwrap();
        assert_eq!(pair_stats.substitutions, 1);
        assert_eq!(pair_stats.gaps, 2);
        assert_eq!(pair_stats.valid_positions, 3);
    }
}
