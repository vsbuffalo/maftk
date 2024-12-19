///
/// ## Developer Notes
///
/// ### Bit Shifting
///
/// Here are some notes on bit-shifting, since this is a bit subtle.
/// A key trick used by BAM indexing, etc are bin sizes that are powers
/// of 2. A position can then be placed in bins that are powers of 2
/// by doing integer division, e.g. position 1020 would be in first
/// bin [0, 1024), which has index 0. The right bit shifts like
/// position >> 10 is equivalent to position / 2^10 but much faster.
///
/// ```
/// let position = 1020_u32;
/// let level = 10;
/// assert_eq!(position >> level, position / 2_u32.pow(10));
/// ```
///
/// ### Hierarchical Binning
///
/// Hierarchical binning (initially used in the UCSC browser,
/// see http://genomewiki.ucsc.edu/index.php/Bin_indexing_system)
/// divides the genome into multiple bins of *varying* power-of-two
/// widths or *levels*. At each level, bins are non-overlapping.
/// Also, lower level (smaller width) bins are perfectly contained
/// by the higher level (larger width) bins (think about this
/// through the right bit shift operations — the smaller bits are just
/// cut off).
///
/// Here is a figure of the hierarchical binning:
///
/// Level 0 (64Mb):   [       Bin 0       ]   (offset: 0)
///                   |                   |
/// Level 1 (1Mb):    [  1  ][  2  ][  3  ]   (offset: 1)
///                   |      \______  
///                   |             \
/// Level 2 (64kb):   [9][10][11][12]         (offset: 9)
///                   |  \___________
///                   |              \
/// Level 3 (8kb):    [73][74][75][76]        (offset: 73)
///                   |   \______________
///                   |                  \
/// Level 4 (1kb):    [585][586][587][588]    (offset: 585)
///
/// Coordinate Space:
/// 0kb          32kb         64kb         96kb
/// |------------|------------|------------|
///
/// A range is put into the *smallest bin it fits entirely into*.
/// This ensures that each range is uniquely placed in a bin.
/// When we want to do an overlap query operation (i.e. finding
/// all ranges that overlap a query range), we take the query
/// range and find all bins *at all levels* that could plausibly
/// contain it. The overlapping bin indices are returned, and then
/// the subset of ranges in each of these bins can be iterated through,
/// checking for overlaps with the query range. This drastically reduces
/// the search time for overlaps compared to checking every range.
///
/// While a range gets stored in a single bin via region_to_bin(),
/// querying for overlaps via region_to_bins() returns bins at all levels
/// since overlapping features could be stored in any of them. Think of
/// region_to_bin() as what's used to uniquely bin a single feature
/// region, and region_to_bins() as finding all possible bins that could
/// contain a feature.
///
/// ### Offsets
///
/// The offsets are the first bin index for each level. It is
/// essential that the offsets are large enough that the indices
/// in the higher level bins do not clash with the bin indices
/// of bin levels below. The binning levels are set with
/// `LEVELS`, which then determines the offsets *at compile time*.
///
use memmap2::{Mmap, MmapOptions};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Result, Seek, Write};
use std::path::Path;
#[allow(unused_imports)]
use tracing::{debug, error, info, trace, warn};

use crate::binary::BlockLocation;

pub mod binning {
    use tracing::trace;
    // LEVELS defines the bin sizes as powers of 2
    // e.g. 2^10 = 1024 (1kb bins), 2^13 = 8192 (8kb bins), etc.
    pub const LEVELS: [u32; 6] = [
        8,  // 2^8 =  256b bins
        10, // 2^10 = 1kb  bins
        13, // 2^13 = 8kb  bins
        16, // 2^16 = 64kb bins
        20, // 2^20 = 1Mb  bins
        26, // 2^26 = 64Mb bins
    ];

    // Get the maximum level value at compile time, so
    // LEVELS can be changed and nothing downstream is
    // hardcoded.
    const fn get_max_level() -> u32 {
        let mut max = 0;
        let mut i = 0;
        while i < LEVELS.len() {
            if LEVELS[i] > max {
                max = LEVELS[i];
            }
            i += 1;
        }
        max
    }

    // Get the level offset for this particular LEVELS,
    // again to make sure that LEVELS can be changed
    // and nothing has to be hardcoded.
    const fn calc_level_offset(idx: usize) -> u32 {
        let mut offset = 0;
        let mut i = 0;
        let max_level = get_max_level();
        while i < idx {
            // Use max_level instead of hardcoded 26
            offset += 1u32 << (max_level - LEVELS[LEVELS.len() - 1 - i]);
            i += 1;
        }
        offset
    }

    // Helper to create array of specified length at compile time
    const fn make_offsets<const N: usize>() -> [u32; N] {
        let mut offsets = [0; N];
        let mut i = 0;
        while i < N {
            offsets[i] = calc_level_offset(i);
            i += 1;
        }
        offsets
    }

    // OFFSETS length automatically matches LEVELS
    pub const OFFSETS: [u32; LEVELS.len()] = make_offsets();

    const _: () = assert!(OFFSETS[0] == 0); // Will fail at compile time if wrong

    // Use 8kb chunks (2^13) for linear index - matches Level 1 bins
    pub const LINEAR_SHIFT: u32 = LEVELS[1];

    /// Maps a genomic interval to the smallest possible bin that fully contains it.
    /// Checks each level from largest (64Mb) to smallest (1kb) bins until finding
    /// a level where start and end fall in the same bin.
    ///
    /// Returns the unique bin ID that can be used to store/retrieve overlapping
    /// intervals. If no single bin contains the interval, returns bin 0 (root).
    ///
    /// Algorithm:
    /// 1. For each level (largest to smallest bins):
    /// 2. Check if start>>shift == end>>shift (same bin at this level)
    /// 3. If yes, return bin = level_offset + (start>>shift)
    /// 4. If no match found, return 0 (root bin)
    pub fn region_to_bin(start: u32, end: u32) -> u32 {
        let mut offset = 0;
        let mut level = 0u32; // Explicitly u32 level counter

        // Iterate through levels from largest to smallest bins
        for &shift in LEVELS.iter().rev() {
            // Check if start and end fall in same bin at this level
            if start >> shift == end >> shift {
                let bin = offset + (start >> shift);
                trace!(start, end, shift, level, offset, bin, "Found matching bin");
                return bin;
            }

            // Calculate offset for next level using level counter
            offset = OFFSETS[level as usize];
            level += 1;
        }

        trace!(start, end, "Falling back to first bin");
        0
    }

    /// Returns all bins that could contain any portion of the given interval.
    /// This includes bins at all levels that overlap any part of [start,end).
    ///
    /// Different from region_to_bin() which finds a single containing bin,
    /// this returns all bins needed for a complete overlap query.
    ///
    /// Algorithm:
    /// 1. For each level:
    /// 2. Calculate first and last overlapping bins:
    ///    b_start = offset + (start>>shift)
    ///    b_end = offset + (end>>shift)
    /// 3. Add all bins in range [b_start,b_end] at this level
    /// 4. Move to next level and adjust offset
    ///
    /// Returns bins sorted by level (largest to smallest) and position.
    pub fn region_to_bins(start: u32, end: u32) -> Vec<u32> {
        let mut list = Vec::with_capacity(LEVELS.len());
        let mut offset = 0;
        let mut level = 0u32; // Explicit u32 level counter

        for &shift in LEVELS.iter().rev() {
            // Calculate the first and last bin that might contain this region
            // at the current level
            let b_start = offset + (start >> shift);
            let b_end = offset + (end >> shift);

            // Add all bins between start and end at this level
            for bin in b_start..=b_end {
                list.push(bin);
            }

            // Calculate offset for next level
            offset = OFFSETS[level as usize];
            level += 1;
        }
        list
    }
}

#[repr(C, packed)]
struct IndexHeader {
    magic: u32,               // "MAFI" in ASCII
    n_bins: u32,              // Number of bins
    n_intervals: u32,         // Total number of intervals
    linear_index_offset: u64, // Where linear index starts
    species_map_offset: u64,  // Where species map starts
}

#[repr(C, packed)]
struct BinTableEntry {
    bin_id: u32,
    n_blocks: u32,
    blocks_offset: u64,
}

pub struct FastIndex {
    mmap: Mmap,
    header: &'static IndexHeader,
    bin_table: &'static [BinTableEntry],
    linear_index: &'static [u64],
    pub species_map: HashMap<u16, String>,
}

impl FastIndex {
    pub fn write(
        blocks: Vec<BlockLocation>,
        species_map: HashMap<u16, String>,
        path: &Path,
    ) -> Result<()> {
        let mut file = BufWriter::new(File::create(path)?);
        let mut bins: HashMap<u32, Vec<BlockLocation>> = HashMap::new();

        // Create linear index using configured chunk size
        let max_pos = blocks.iter().map(|b| b.end).max().unwrap_or(0);
        let mut linear_index = vec![u64::MAX; (max_pos >> binning::LINEAR_SHIFT) as usize + 1];

        let nblocks = blocks.len();
        // Sort blocks into bins and update linear index
        for block in blocks {
            let bin = binning::region_to_bin(block.start, block.end);
            linear_index[(block.start >> binning::LINEAR_SHIFT) as usize] =
                linear_index[(block.start >> binning::LINEAR_SHIFT) as usize].min(block.offset);
            bins.entry(bin).or_default().push(block);
        }

        // Calculate offsets
        let header_size = std::mem::size_of::<IndexHeader>();
        let bin_table_size = bins.len() * std::mem::size_of::<BinTableEntry>();
        let mut blocks_offset = header_size + bin_table_size;

        // Write placeholder header
        let mut header = IndexHeader {
            magic: 0x4D414649, // "MAFI"
            n_bins: bins.len() as u32,
            n_intervals: nblocks as u32,
            linear_index_offset: 0,
            species_map_offset: 0,
        };
        file.write_all(unsafe { as_bytes(&header) })?;

        // Write bin table and blocks
        for (bin_id, bin_blocks) in &bins {
            let entry = BinTableEntry {
                bin_id: *bin_id,
                n_blocks: bin_blocks.len() as u32,
                blocks_offset: blocks_offset as u64,
            };
            file.write_all(unsafe { as_bytes(&entry) })?;
            blocks_offset += bin_blocks.len() * std::mem::size_of::<BlockLocation>();
        }

        // Write block lists
        for blocks in bins.values() {
            for block in blocks {
                file.write_all(unsafe { as_bytes(block) })?;
            }
        }

        // Write linear index
        header.linear_index_offset = file.stream_position()?;
        for offset in &linear_index {
            file.write_all(unsafe { as_bytes(offset) })?;
        }

        // Write species map
        header.species_map_offset = file.stream_position()?;
        bincode::serialize_into(&mut file, &species_map)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Update header
        file.seek(std::io::SeekFrom::Start(0))?;
        file.write_all(unsafe { as_bytes(&header) })?;

        Ok(())
    }

    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };

        let header = unsafe { &*(mmap.as_ptr() as *const IndexHeader) };
        if header.magic != 0x4D414649 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid index file magic number",
            ));
        }

        let bin_table = unsafe {
            std::slice::from_raw_parts(
                (mmap.as_ptr().add(std::mem::size_of::<IndexHeader>())) as *const BinTableEntry,
                header.n_bins as usize,
            )
        };

        let linear_index = unsafe {
            std::slice::from_raw_parts(
                (mmap.as_ptr().add(header.linear_index_offset as usize)) as *const u64,
                (header.species_map_offset - header.linear_index_offset) as usize / 8,
            )
        };

        let species_map: HashMap<u16, String> =
            bincode::deserialize(&mmap[header.species_map_offset as usize..])
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        Ok(FastIndex {
            mmap,
            header,
            bin_table,
            linear_index,
            species_map,
        })
    }

    /// Find alignment blocks overlapping a genomic interval
    ///
    /// # Arguments
    /// * `start` - Start position (0-based)
    /// * `end` - End position (0-based)
    ///
    /// For a region [start, end), returns blocks that overlap this region, sorted
    /// by their offset in the binary MAF file. Uses a binning scheme to efficiently
    /// find candidate blocks and a linear index to optimize seeking within the file.
    ///
    /// The linear index at the `LINEAR_SHIFT` resolution (8kb) is used to
    /// skip blocks that occur before the query region, reducing seek time.
    pub fn find_blocks(&self, start: u32, end: u32) -> Vec<BlockLocation> {
        let mut result = Vec::new();
        let possible_bins = binning::region_to_bins(start, end);
        // debug!("Searching bins: {:?}", possible_bins);
        let min_offset = self
            .linear_index
            .get(start as usize >> binning::LINEAR_SHIFT)
            .copied()
            .unwrap_or(0);

        for bin in possible_bins {
            if let Ok(idx) = self
                .bin_table
                .binary_search_by_key(&bin, |entry| entry.bin_id)
            {
                let entry = &self.bin_table[idx];

                let blocks = unsafe {
                    std::slice::from_raw_parts(
                        (self.mmap.as_ptr().add(entry.blocks_offset as usize))
                            as *const BlockLocation,
                        entry.n_blocks as usize,
                    )
                };

                for block in blocks {
                    if block.offset >= min_offset && block.end > start && block.start < end {
                        result.push(block.clone());
                    }
                }
            }
        }

        result.sort_by_key(|b| b.offset);
        // debug!("FastIndex.find_blocks() found {} blocks.", result.len());
        result
    }

    pub fn species_map(&self) -> &HashMap<u16, String> {
        &self.species_map
    }

    pub fn debug_print(&self) {
        // Basic header info
        let magic = self.header.magic;
        let n_bins = self.header.n_bins;
        let n_intervals = self.header.n_intervals;

        println!("Index Header:");
        println!(
            "  Magic: 0x{:08X} ({})",
            magic,
            String::from_utf8_lossy(&magic.to_be_bytes())
        );
        println!("  Number of bins: {}", n_bins);
        println!("  Number of intervals: {}", n_intervals);

        // Analyze block sizes and distribution
        let mut total_blocks = 0u32;
        let mut block_counts = vec![0u32; binning::LEVELS.len()];

        for entry in self.bin_table.iter() {
            let blocks = unsafe {
                std::slice::from_raw_parts(
                    (self.mmap.as_ptr().add(entry.blocks_offset as usize)) as *const BlockLocation,
                    entry.n_blocks as usize,
                )
            };

            for block in blocks {
                total_blocks += 1;
                let size = block.end - block.start;

                // Find appropriate size bucket using LEVELS
                let bucket = binning::LEVELS
                    .iter()
                    .position(|&shift| size <= (1 << shift))
                    .unwrap_or(binning::LEVELS.len() - 1);
                block_counts[bucket] += 1;
            }
        }

        println!("\nBlock Size Distribution:");
        println!("Size Range      Count    Percentage");
        println!("---------      -----    ----------");

        // Create size range descriptions using LEVELS
        for (i, &shift) in binning::LEVELS.iter().enumerate() {
            let count = block_counts[i];
            let percentage = (count as f64 / total_blocks as f64) * 100.0;

            let range_desc = if i == 0 {
                format!("≤{:>3}kb", 1 << (shift - 10))
            } else {
                let prev_size = 1 << (binning::LEVELS[i - 1] - 10);
                let this_size = 1 << (shift - 10);
                format!("{:>3}kb-{:>3}kb", prev_size, this_size)
            };

            println!("{:<12}    {:>5}    {:>6.2}%", range_desc, count, percentage);
        }

        // Calculate average block size
        let mut total_size = 0u64;
        let mut max_size = 0u32;

        for entry in self.bin_table.iter() {
            let blocks = unsafe {
                std::slice::from_raw_parts(
                    (self.mmap.as_ptr().add(entry.blocks_offset as usize)) as *const BlockLocation,
                    entry.n_blocks as usize,
                )
            };

            for block in blocks {
                let size = block.end - block.start;
                total_size += size as u64;
                max_size = max_size.max(size);
            }
        }

        let avg_size = total_size as f64 / total_blocks as f64;
        println!("\nBlock Size Statistics:");
        println!("  Average size: {:.1} bp", avg_size);
        println!("  Maximum size: {} bp", max_size);
        println!("  Total blocks: {}", total_blocks);
    }
}

unsafe fn as_bytes<T>(x: &T) -> &[u8] {
    std::slice::from_raw_parts((x as *const T) as *const u8, std::mem::size_of::<T>())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::tempdir;

    #[test]
    fn test_binning_functions() {}

    //#[test]
    //fn test_binning_functions() {
    //    // Test region_to_bin with small interval that fits in a single bin
    //    let bin = binning::region_to_bin(1000, 1500);
    //    assert_eq!(bin, 74, "Small interval should map to bin 74");

    //    // Test region_to_bin with larger interval
    //    let bin = binning::region_to_bin(1000, 2000);
    //    assert_eq!(bin, 9, "Larger interval should map to an 8kb bin");

    //    // Test region_to_bins
    //    let bins = binning::region_to_bins(1000, 3000);
    //    // Check for expected bins at each level
    //    assert!(bins.contains(&0), "Should contain 1Mb bin");
    //    assert!(bins.contains(&1), "Should contain 64kb bin");
    //    assert!(bins.contains(&9), "Should contain first 8kb bin");
    //    assert!(bins.contains(&74), "Should contain second 1kb bin");
    //}

    //#[test]
    //fn test_fastindex_basic() -> Result<()> {
    //    let dir = tempdir()?;
    //    let index_path = dir.path().join("test.index");

    //    // Create some test blocks
    //    let blocks = vec![
    //        BlockLocation {
    //            start: 1000,
    //            end: 1500,
    //            offset: 0,
    //            size: 500,
    //        },
    //        BlockLocation {
    //            start: 2000,
    //            end: 3000,
    //            offset: 1000,
    //            size: 1000,
    //        },
    //    ];

    //    // Create a simple species map
    //    let mut species_map = HashMap::new();
    //    species_map.insert(0, "species1".to_string());
    //    species_map.insert(1, "species2".to_string());

    //    // Write index
    //    FastIndex::write(blocks, species_map.clone(), &index_path)?;

    //    // Read it back
    //    let index = FastIndex::open(&index_path)?;

    //    // Test querying blocks
    //    let found_blocks = index.find_blocks(1200, 2500);
    //    assert_eq!(found_blocks.len(), 2, "Should find both blocks");

    //    // Test species map was preserved
    //    assert_eq!(index.species_map(), &species_map);

    //    Ok(())
    //}

    //#[test]
    //fn test_block_retrieval() -> Result<()> {
    //    let dir = tempdir()?;
    //    let index_path = dir.path().join("test.index");

    //    // Create overlapping blocks
    //    let blocks = vec![
    //        BlockLocation {
    //            start: 1000,
    //            end: 2000,
    //            offset: 0,
    //            size: 1000,
    //        },
    //        BlockLocation {
    //            start: 1500,
    //            end: 2500,
    //            offset: 1000,
    //            size: 1000,
    //        },
    //        BlockLocation {
    //            start: 3000,
    //            end: 4000,
    //            offset: 2000,
    //            size: 1000,
    //        },
    //    ];

    //    let species_map = HashMap::new();
    //    FastIndex::write(blocks, species_map, &index_path)?;
    //    let index = FastIndex::open(&index_path)?;

    //    // Test various queries
    //    let found = index.find_blocks(1200, 1800);
    //    assert_eq!(found.len(), 2, "Should find two overlapping blocks");

    //    let found = index.find_blocks(2600, 2900);
    //    assert_eq!(found.len(), 0, "Should find no blocks in gap");

    //    let found = index.find_blocks(3500, 3800);
    //    assert_eq!(found.len(), 1, "Should find one block");

    //    Ok(())
    //}
}
