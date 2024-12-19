// binary.rs
use crate::index::FastIndex;
use crate::io::MafReader;
use memmap2::{Mmap, MmapOptions};
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::{create_dir_all, File},
    io::{BufWriter, Seek, Write},
    path::Path,
};
use tracing::{debug, error, trace, warn};

#[repr(C, packed)]
struct SequenceHeader {
    species_id: u16,
}

#[derive(Clone, Serialize, Deserialize)]
#[repr(C, packed)]
pub struct BlockLocation {
    pub start: u32,
    pub end: u32,
    pub offset: u64,
    pub size: u32,
}

// On-disk format for each block
#[repr(C, packed)]
#[derive(Copy, Clone, Debug)]
struct BlockHeader {
    n_seqs: u16,
    seq_length: u32,
    score: f64,
}

#[derive(Clone, Debug)]
pub struct BlockResult {
    /// File offset to block in binary MAF
    offset: u64,
    /// block metadata (sequences, length, score)
    header: BlockHeader,
    /// starting position of block in genome
    start_pos: u32,
}

pub fn convert_to_binary(
    input: &Path,
    output_dir: &Path,
    min_length: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    create_dir_all(output_dir)?;
    let mut maf = MafReader::from_file(input)?;
    let _header = maf.read_header()?;
    // let mut skipped_min_length = 0;
    // let mut lengths = Histogram::new();

    // Track blocks and species maps per chromosome
    let mut chr_data: HashMap<
        String,
        (
            Vec<BlockLocation>,
            HashMap<u16, String>,
            HashMap<String, u16>,
        ),
    > = HashMap::new();
    let mut current_file: Option<(String, File)> = None;

    while let Some(block) = maf.next_block()? {
        if block.sequences[0].size < min_length {
            // skipped_min_length += 1;
            continue;
        }
        // lengths.add(block.sequences[0].size as usize);

        let ref_seq = &block.sequences[0];
        let ref_seq_len = ref_seq.text.len();
        let chr = ref_seq.src.split('.').nth(1).ok_or("Invalid ref name")?;

        // Create/get file for this chromosome
        if current_file.as_ref().map_or(true, |(c, _)| c != chr) {
            current_file = Some((
                chr.to_string(),
                File::create(output_dir.join(format!("{}.bin", chr)))?,
            ));
        }

        let (_, file) = current_file.as_mut().unwrap();
        let mut writer = BufWriter::new(file);
        let offset = writer.stream_position()?;

        // Get or create data structures for this chromosome
        let (blocks, species_map, reverse_map) = chr_data
            .entry(chr.to_string())
            .or_insert_with(|| (Vec::new(), HashMap::new(), HashMap::new()));

        // Write block header
        let header = BlockHeader {
            n_seqs: block.sequences.len() as u16,
            // NOTE: use the text length instead of size, as the size
            // is the length of the *ungapped* alignment!
            seq_length: ref_seq.text.len() as u32,
            score: block.score.unwrap_or(0.0),
        };
        writer.write_all(unsafe {
            std::slice::from_raw_parts(
                &header as *const _ as *const u8,
                std::mem::size_of::<BlockHeader>(),
            )
        })?;
        // Write sequences
        for seq in &block.sequences {
            let species = seq.src.split('.').next().unwrap_or("unknown");

            // Get or assign species ID
            let species_id = reverse_map.entry(species.to_string()).or_insert_with(|| {
                let new_id = species_map.len() as u16;
                species_map.insert(new_id, species.to_string());
                new_id
            });
            // println!("Writing species {} with ID {}", species, species_id);

            // Write sequence header with species ID
            let seq_header = SequenceHeader {
                species_id: *species_id,
            };
            writer.write_all(unsafe {
                std::slice::from_raw_parts(
                    &seq_header as *const _ as *const u8,
                    std::mem::size_of::<SequenceHeader>(),
                )
            })?;

            // Write sequence data
            // dbg!(&block.sequences);
            // Make sure that all sequences are the same length
            // as the first ref sequence.
            assert_eq!(seq.text.len(), ref_seq_len);
            if seq
                .text
                .to_uppercase()
                .chars()
                .any(|c| !"ATGCN-".contains(c))
            {
                panic!("Invalid sequence data found: {}", seq.text);
            }
            writer.write_all(seq.text.as_bytes())?;
        }

        // Store block location
        let seq_total_size = std::mem::size_of::<SequenceHeader>() + header.seq_length as usize;
        let block_entry = BlockLocation {
            start: ref_seq.start as u32,
            end: (ref_seq.start + ref_seq.size) as u32,
            offset,
            size: std::mem::size_of::<BlockHeader>() as u32
                + (seq_total_size as u32 * header.n_seqs as u32),
        };
        blocks.push(block_entry);
    }

    // Write indices
    for (chr, (blocks, species_map, _)) in chr_data {
        let index_path = output_dir.join(format!("{}.index", chr));
        FastIndex::write(blocks, species_map, &index_path)?;
    }

    // lengths.print("Block Length Distribution (bp)", 1000);
    Ok(())
}

pub struct IndexedMafDb {
    mmap: Mmap,
    index: FastIndex,
}

impl IndexedMafDb {
    pub fn open(chr: &str, data_dir: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let data_path = data_dir.join(format!("{}.bin", chr));
        let index_path = data_dir.join(format!("{}.index", chr));

        debug!(data_path = ?data_path, "Opening data file");
        let data_file = File::open(&data_path)?;

        debug!(index_path = ?index_path, "Opening index file");
        let index = FastIndex::open(&index_path)?;

        debug!(species_map = ?index.species_map(), "Loaded species map");
        debug!("Creating memory map");
        let mmap = unsafe { MmapOptions::new().map(&data_file)? };

        Ok(Self { mmap, index })
    }

    /// Find alignment blocks overlapping a genomic interval
    ///
    /// Note: Positions are 1-based, [start, end] coordinates
    pub fn find_blocks(&self, start: u32, end: u32) -> Vec<BlockResult> {
        let internal_start = start - 1; // Convert 1-indexed to 0-indexed
        let internal_end = end; // end stays same as exclusive end = inclusive end + 1
        debug!(
            "Finding blocks: original start={}, end={}, internal_start={}, internal_end={}",
            start, end, internal_start, internal_end
        );
        let blocks = self.index.find_blocks(internal_start, internal_end);
        debug!(
            "IndexedMafReader.find_blocks() found {} blocks",
            blocks.len()
        );
        self.index
            .find_blocks(internal_start, internal_end)
            .into_iter()
            .map(|loc| {
                let header =
                    unsafe { &*(self.mmap[loc.offset as usize..].as_ptr() as *const BlockHeader) };
                BlockResult {
                    offset: loc.offset,
                    header: *header,
                    start_pos: loc.start,
                }
            })
            .collect()
    }

    /// Get the alignments in a region.
    pub fn get_alignments(
        &self,
        start: u32,
        end: u32,
    ) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
        let seq_header_size = std::mem::size_of::<SequenceHeader>();
        let species_map = self.index.species_map();

        let mut alignments: Vec<(String, String)> = Vec::new();

        let blocks = self.find_blocks(start, end);
        for block in blocks {
            let block_size = block.header.seq_length as usize;
            let n_seqs = block.header.n_seqs as usize;

            // Calculate relative positions within the block
            let start_pos = if start > block.start_pos {
                (start - block.start_pos) as usize
            } else {
                0
            };
            let end_pos = ((end - block.start_pos) as usize).min(block_size);

            trace!(
                ?block.header,
                block.start_pos,
                block_size,
                start,
                end,
                start_pos,
                end_pos,
                "Processing alignment block"
            );

            if start_pos < block_size && start_pos < end_pos {
                // The start of the sequence data (species headers and aligned sequences)
                let data_start = block.offset as usize + std::mem::size_of::<BlockHeader>();

                for i in 0..n_seqs {
                    let seq_start = data_start + (i * (seq_header_size + block_size));

                    // Read sequence header
                    let species_id = unsafe {
                        trace!(seq_start, "Reading sequence header");
                        let data = &self.mmap[seq_start..seq_start + std::mem::size_of::<u16>()];
                        trace!(data = ?data, "Raw bytes");
                        let ptr = data.as_ptr();
                        let id = std::ptr::read_unaligned(ptr as *const u16);
                        trace!(species_id = id);
                        id
                    };

                    // Read sequence data
                    let seq_data = &self.mmap
                        [seq_start + seq_header_size..seq_start + seq_header_size + block_size];
                    trace!(seq_data = ?seq_data, "Raw sequence data");

                    match std::str::from_utf8(seq_data) {
                        Ok(sequence) => {
                            trace!(sequence, "Parsed sequence");
                            if start_pos < sequence.len() {
                                let display_end = end_pos.min(sequence.len());
                                let sequence_slice = &sequence[start_pos..display_end];

                                match species_map.get(&species_id) {
                                    Some(species_name) => {
                                        alignments.push((
                                            species_name.to_string(),
                                            sequence_slice.to_string(),
                                        ));
                                    }
                                    None => {
                                        error!(species_id, "Could not find species name");
                                        return Err(format!(
                                            "Could not find species name for {species_id}"
                                        )
                                        .into());
                                    }
                                }
                            } else {
                                warn!(
                                    start_pos,
                                    sequence_len = sequence.len(),
                                    "Start position exceeds sequence length"
                                );
                            }
                        }
                        Err(e) => error!(error = %e, sequence_number = i, "Error reading sequence"),
                    }
                }
            }
        }
        Ok(alignments)
    }

    pub fn species_map(&self) -> &HashMap<u16, String> {
        self.index.species_map()
    }
}
