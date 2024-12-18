// src/index.rs
/// The index uses a hierarchical binning scheme, like those used in BAM/tabix indexes,
/// to efficiently store and query genomic intervals. The scheme divides the genome into
/// bins of different sizes, with larger bins at higher levels containing smaller bins.
///
/// Position on genome:        0    1kb   2kb   3kb   4kb   5kb   6kb   7kb   8kb
/// Level 3 (1kb bins):     |  0  |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  ...
/// Level 2 (8kb bins):     |--------------- 0 ----------------|  ...
/// Level 1 (64kb bins):    |------------------------- 0 -------------------------| ...
/// Level 0 (1Mb bins):     |-------------------------------- 0 --------------------------------|
///
/// # Binning Scheme
///
/// The index uses 4 levels of bins, from largest to smallest:
/// ```text
/// Level 0: 1Mb   bins (2^20 = 1,048,576 bp)
/// Level 1: 64kb  bins (2^16 = 65,536 bp)
/// Level 2: 8kb   bins (2^13 = 8,192 bp)
/// Level 3: 1kb   bins (2^10 = 1,024 bp)
/// ```
///
/// # Bin Numbering
///
/// Bins are numbered sequentially with a unique ID that encodes both the level and position.
/// The first bins at each level start at:
/// ```text
/// Level 0 (1Mb):   0
/// Level 1 (64kb):  1
/// Level 2 (8kb):   9
/// Level 3 (1kb):   73
/// ```
///
/// The offset for each level is calculated as: `level_offset += 1 << (LEVELS.len() * 3 - level * 3)`
///
/// # Examples
///
/// ```
/// // Example 1: Small interval fits in a single 1kb bin
/// let start = 1000;
/// let end = 1500;
/// let bin = region_to_bin(start, end);
/// assert_eq!(bin, 74);  // Level 3, second 1kb bin
///
/// // Example 2: Larger interval spans multiple bins
/// let start = 1000;
/// let end = 3000;
/// let bins = region_to_bins(start, end);
/// // Returns overlapping bins at all levels:
/// // - Level 0 (1Mb):  bin 0
/// // - Level 1 (64kb): bin 1
/// // - Level 2 (8kb):  bins 9, 10
/// // - Level 3 (1kb):  bins 74, 75, 76
/// ```
///
/// # Query Algorithm
///
/// When querying an interval:
/// 1. Calculate all potentially overlapping bins using `region_to_bins()`
/// 2. Search each bin for overlapping intervals
/// 3. Use linear index to optimize searching within bins
///
/// For example, querying interval [1000-2000]:
/// ```text
/// 1. Calculate bins:
///    - Could be in 1Mb bin:   start>>20 == end>>20  (bin 0)
///    - Could be in 64kb bins: start>>16 == end>>16  (bin 1)
///    - Could be in 8kb bins:  start>>13 == end>>13  (bin 9)
///    - Must check 1kb bins:   bins 74, 75
/// 2. Check each bin for overlapping intervals
/// ```
///
/// # Performance Considerations
///
/// The bin sizes were chosen based on empirical analysis of multiz 100-way alignments:
/// - Most blocks (~99.9%) are under 1kb
/// - Mean block size is ~26bp
/// - Very few blocks (< 0.01%) are over 50kb
///
/// This distribution led to choosing smaller bins than traditional BAM indexing
/// to optimize for the common case of small alignment blocks.
use memmap2::{Mmap, MmapOptions};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Result, Seek, Write};
use std::path::Path;
#[allow(unused_imports)]
use tracing::{debug, error, info, trace, warn};

use crate::binary::BlockLocation;

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

mod binning {
    use tracing::trace;
    // LEVELS defines the bin sizes as powers of 2
    // e.g. 2^10 = 1024 (1kb bins), 2^13 = 8192 (8kb bins), etc.
    pub const LEVELS: [u32; 5] = [
        10, // 2^10 = 1kb bins
        13, // 2^13 = 8kb bins
        16, // 2^16 = 64kb bins
        20, // 2^20 = 1Mb bins
        26, // 2^26 = 64Mb bins
    ];

    // Use 8kb chunks (2^13) for linear index - matches Level 1 bins
    pub const LINEAR_SHIFT: u32 = LEVELS[1];

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
            offset += 1 << (LEVELS.len() as u32 * 3 - (level * 3));
            level += 1;
        }

        trace!(start, end, "Falling back to first bin");
        0
    }

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
            offset += 1 << (LEVELS.len() as u32 * 3 - (level * 3));
            level += 1;
        }
        list
    }
}

pub struct FastIndex {
    mmap: Mmap,
    header: &'static IndexHeader,
    bin_table: &'static [BinTableEntry],
    linear_index: &'static [u64],
    species_map: HashMap<u16, String>,
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
    /// # Example
    /// ```rust
    /// let blocks = index.find_blocks(1000, 2000);
    /// // Returns blocks spanning [1000,2000), sorted by file offset
    /// ```
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


