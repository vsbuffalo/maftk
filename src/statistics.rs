// statistics.rs
use csv::{Writer, WriterBuilder};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use thiserror::Error;
use tracing::info;

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
    pub substitutions: u32,
    pub gaps: u32,
    pub valid_positions: u32,
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
    pub fn substitution_rate(&self) -> f64 {
        if self.valid_positions == 0 {
            0.0
        } else {
            self.substitutions as f64 / self.valid_positions as f64
        }
    }

    pub fn gap_rate(&self) -> f64 {
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
                headers.push(format!("{}_{}_num_valid", species_vec[i], species_vec[j]));
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
                        record.push(pair_stats.valid_positions.to_string());
                    }
                    None => {
                        record.push("NA".to_string());
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

pub fn is_gap(b: u8) -> bool {
    b == b'-' || b == b'.'
}

pub fn compare_bases(b1: u8, b2: u8) -> bool {
    // First check if either is a gap
    if is_gap(b1) || is_gap(b2) {
        return false; // gaps are handled separately
    }

    // Convert both to uppercase for comparison
    let b1_upper = b1.to_ascii_uppercase();
    let b2_upper = b2.to_ascii_uppercase();

    b1_upper == b2_upper
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

    // If region is specified, calculate overlap
    let (start_offset, length) = if let (Some(start), Some(end)) = (region_start, region_end) {
        // Check if there's any overlap
        if end <= block_start || start >= block_end {
            info!(
                "No overlap between block [{}, {}) and region [{}, {})",
                block_start, block_end, start, end
            );
            return None;
        }

        // Calculate overlap region
        let overlap_start = start.max(block_start);
        let overlap_end = end.min(block_end);

        // Calculate offset and length
        (
            (overlap_start - block_start) as usize,
            (overlap_end - overlap_start) as usize,
        )
    } else {
        // Use full sequence if no region specified
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
        // Use slice of the sequence for the region of interest
        let seq1_bytes = &seq1.text.as_bytes()[start_offset..start_offset + length];

        for j in (i + 1)..block.sequences.len() {
            let seq2 = &block.sequences[j];
            if species_indices
                .as_ref()
                .map_or(false, |indices| !indices.contains(&seq2.species_idx))
            {
                continue;
            }
            // Use slice of the sequence for the region of interest
            let seq2_bytes = &seq2.text.as_bytes()[start_offset..start_offset + length];

            let mut pair_stats = PairwiseStats::default();

            for (b1, b2) in seq1_bytes.iter().zip(seq2_bytes.iter()) {
                if is_gap(*b1) && is_gap(*b2) {
                    continue; // Matching gaps
                } else if is_gap(*b1) || is_gap(*b2) {
                    pair_stats.gaps += 1; // One gap
                } else {
                    // Both are bases
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
    use crate::binary::AlignedSequence;

    use super::*;

    // Helper function to create a MafBlock with two sequences
    fn create_test_block(seq1: &str, seq2: &str) -> MafBlock {
        MafBlock {
            score: 0.0,
            sequences: vec![
                AlignedSequence {
                    species_idx: 0,
                    start: 0,
                    size: seq1.len() as u32,
                    strand: false,
                    src_size: seq1.len() as u32,
                    text: seq1.to_string(),
                },
                AlignedSequence {
                    species_idx: 1,
                    start: 0,
                    size: seq2.len() as u32,
                    strand: false,
                    src_size: seq2.len() as u32,
                    text: seq2.to_string(),
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
