/// io.rs
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use thiserror::Error;

use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{BufWriter, Error, Write};

#[derive(Debug, Error)]
pub enum MafError {
    #[error("Invalid MAF format: {0}")]
    ParseError(String),
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Represents a MAF file header
#[derive(Debug, Clone)]
pub struct MafHeader {
    pub version: String,
    pub scoring: Option<String>,
    pub program: Option<String>,
}

/// Represents a sequence within an alignment block
#[derive(Debug, Clone)]
pub struct Sequence {
    /// Source sequence name (e.g. the species' genome version)
    pub src: String,
    /// 0-based start position
    pub start: u64,
    /// size of the aligning region
    pub size: u64,
    /// + or -
    pub strand: Strand,
    /// total size of source sequence
    pub src_size: u64,
    /// alignment text including gaps
    pub text: String,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl FromStr for Strand {
    type Err = MafError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Err(MafError::ParseError("Invalid strand".to_string())),
        }
    }
}

/// Information about sequences before/after alignment block
#[derive(Debug, Clone)]
pub struct Info {
    pub src: String,
    pub left_status: StatusChar,
    pub left_count: u64,
    pub right_status: StatusChar,
    pub right_count: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StatusChar {
    Contiguous,  // C
    Intervening, // I
    New,         // N
    NewBridged,  // n
    Missing,     // M
    Tandem,      // T
}

/// Represents a complete alignment block
#[derive(Debug, Clone)]
pub struct AlignmentBlock {
    pub score: Option<f64>,
    pub pass: Option<u32>,
    pub sequences: Vec<Sequence>,
    pub infos: Vec<Info>,
}

/// Checks if a file is gzipped by looking at magic numbers
pub fn is_gzipped(filepath: &Path) -> io::Result<bool> {
    let mut file = File::open(filepath)?;
    let mut magic_bytes = [0u8; 2];

    // Read first two bytes - gzip magic numbers are 0x1f 0x8b
    file.read_exact(&mut magic_bytes)?;
    Ok(magic_bytes == [0x1f, 0x8b])
}

/// Creates an appropriate BufReader based on whether the file is gzipped or not
pub fn get_reader(filepath: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(filepath)?;

    if is_gzipped(filepath)? {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub type BoxedReader = Box<dyn BufRead>;

pub struct MafReader {
    reader: BoxedReader,
    header: Option<MafHeader>,
    current_line: String,
}

impl MafReader {
    pub fn new(reader: BoxedReader) -> Self {
        Self {
            reader,
            header: None,
            current_line: String::new(),
        }
    }

    /// Creates a new MafReader from a file path, automatically detecting if it's gzipped
    pub fn from_file(filepath: &Path) -> io::Result<Self> {
        let reader = get_reader(filepath)?;
        Ok(Self::new(reader))
    }

    pub fn read_header(&mut self) -> Result<&MafHeader, MafError> {
        if self.header.is_some() {
            return Ok(self.header.as_ref().unwrap());
        }

        self.current_line.clear();
        self.reader.read_line(&mut self.current_line)?;

        if !self.current_line.starts_with("##maf") {
            return Err(MafError::ParseError("Missing ##maf header".to_string()));
        }

        // Parse header key-value pairs
        let mut version = None;
        let mut scoring = None;
        let mut program = None;

        for pair in self.current_line["##maf".len()..].split_whitespace() {
            let mut parts = pair.split('=');
            match (parts.next(), parts.next()) {
                (Some("version"), Some(v)) => version = Some(v.to_string()),
                (Some("scoring"), Some(s)) => scoring = Some(s.to_string()),
                (Some("program"), Some(p)) => program = Some(p.to_string()),
                _ => continue,
            }
        }

        self.header = Some(MafHeader {
            version: version.ok_or_else(|| MafError::ParseError("Missing version".to_string()))?,
            scoring,
            program,
        });

        Ok(self.header.as_ref().unwrap())
    }

    pub fn next_block(&mut self) -> Result<Option<AlignmentBlock>, MafError> {
        let mut block: Option<AlignmentBlock> = None;
        let mut sequences = Vec::new();
        let mut infos = Vec::new();

        self.current_line.clear();
        while self.reader.read_line(&mut self.current_line)? > 0 {
            let line = self.current_line.trim();
            if line.is_empty() {
                // End of block
                if block.is_some() {
                    let mut b: AlignmentBlock = block.take().unwrap();
                    b.sequences = sequences;
                    b.infos = infos;
                    return Ok(Some(b));
                }
                continue;
            }

            match &line[0..1] {
                "a" => {
                    // New alignment block
                    let mut score = None;
                    let mut pass = None;

                    for pair in line[1..].split_whitespace() {
                        let mut parts = pair.split('=');
                        match (parts.next(), parts.next()) {
                            (Some("score"), Some(s)) => {
                                score = Some(s.parse().map_err(|_| {
                                    MafError::ParseError("Invalid score".to_string())
                                })?)
                            }
                            (Some("pass"), Some(p)) => {
                                pass = Some(p.parse().map_err(|_| {
                                    MafError::ParseError("Invalid pass".to_string())
                                })?)
                            }
                            _ => continue,
                        }
                    }

                    block = Some(AlignmentBlock {
                        score,
                        pass,
                        sequences: Vec::new(),
                        infos: Vec::new(),
                    });
                }
                "s" => {
                    // Sequence line
                    let parts: Vec<&str> = line[1..].split_whitespace().collect();
                    if parts.len() != 6 {
                        return Err(MafError::ParseError("Invalid sequence line".to_string()));
                    }

                    sequences.push(Sequence {
                        src: parts[0].to_string(),
                        start: parts[1].parse().map_err(|_| {
                            MafError::ParseError("Invalid start position".to_string())
                        })?,
                        size: parts[2]
                            .parse()
                            .map_err(|_| MafError::ParseError("Invalid size".to_string()))?,
                        strand: parts[3].parse()?,
                        src_size: parts[4]
                            .parse()
                            .map_err(|_| MafError::ParseError("Invalid source size".to_string()))?,
                        text: parts[5].to_string(),
                    });
                }
                "i" => {
                    // Info line
                    let parts: Vec<&str> = line[1..].split_whitespace().collect();
                    if parts.len() != 5 {
                        return Err(MafError::ParseError("Invalid info line".to_string()));
                    }

                    infos.push(Info {
                        src: parts[0].to_string(),
                        left_status: match parts[1] {
                            "C" => StatusChar::Contiguous,
                            "I" => StatusChar::Intervening,
                            "N" => StatusChar::New,
                            "n" => StatusChar::NewBridged,
                            "M" => StatusChar::Missing,
                            "T" => StatusChar::Tandem,
                            _ => {
                                return Err(MafError::ParseError("Invalid left status".to_string()))
                            }
                        },
                        left_count: parts[2]
                            .parse()
                            .map_err(|_| MafError::ParseError("Invalid left count".to_string()))?,
                        right_status: match parts[3] {
                            "C" => StatusChar::Contiguous,
                            "I" => StatusChar::Intervening,
                            "N" => StatusChar::New,
                            "n" => StatusChar::NewBridged,
                            "M" => StatusChar::Missing,
                            "T" => StatusChar::Tandem,
                            _ => {
                                return Err(MafError::ParseError(
                                    "Invalid right status".to_string(),
                                ))
                            }
                        },
                        right_count: parts[4]
                            .parse()
                            .map_err(|_| MafError::ParseError("Invalid right count".to_string()))?,
                    });
                }
                _ => (),
            }
            self.current_line.clear();
        }

        if block.is_some() {
            let mut b: AlignmentBlock = block.take().unwrap();
            b.sequences = sequences;
            b.infos = infos;
            Ok(Some(b))
        } else {
            Ok(None)
        }
    }
}

#[derive(Error, Debug)]
pub enum IoError {
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

pub struct InputFile {
    pub filepath: PathBuf,
}

impl InputFile {
    pub fn new(filepath: &Path) -> Self {
        Self {
            filepath: filepath.into(),
        }
    }

    pub fn reader(&self) -> Result<BufReader<Box<dyn Read>>, IoError> {
        let file = File::open(self.filepath.clone())?;
        let reader: Box<dyn Read> = if self.filepath.ends_with(".gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(BufReader::new(reader))
    }

    pub fn has_header(&self, expect: &str) -> Result<bool, IoError> {
        let mut buf_reader = self.reader()?;
        let mut first_line = String::new();
        buf_reader.read_line(&mut first_line)?;
        let has_header = first_line.starts_with(expect);
        Ok(has_header)
    }
}

pub struct OutputFile {
    pub filepath: PathBuf,
}

impl OutputFile {
    pub fn new(filepath: &Path) -> Self {
        Self {
            filepath: filepath.into(),
        }
    }
    pub fn writer(&self) -> Result<Box<dyn Write>, Error> {
        let outfile = &self.filepath;
        let is_gzip = outfile.ends_with(".gz");
        let bed_writer: Box<dyn Write> = if is_gzip {
            Box::new(BufWriter::new(GzEncoder::new(
                File::create(outfile)?,
                Compression::default(),
            )))
        } else {
            Box::new(BufWriter::new(File::create(outfile)?))
        };
        Ok(bed_writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // Helper function to get path to test data
    fn test_data_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("data")
    }

    fn fetch_chr22_maf_reader() -> MafReader {
        let chr22_file = test_data_dir().join("chr22_chunk.maf");
        MafReader::from_file(&chr22_file).unwrap()
    }

    #[test]
    fn test_sequence_lengths_match() {
        let mut maf = fetch_chr22_maf_reader();
        let block = maf.next_block().unwrap().unwrap();

        // Get the length of the first sequence's text
        let expected_length = block.sequences[0].text.len();

        // All sequences should have the same text length (including gaps)
        for seq in &block.sequences {
            assert_eq!(
                seq.text.len(),
                expected_length,
                "Sequence {} has different length",
                seq.src
            );
        }
    }

    #[test]
    fn test_gap_handling() {
        let mut maf = fetch_chr22_maf_reader();
        let block = maf.next_block().unwrap().unwrap();

        // Check that gaps are preserved in the sequences
        let hg38 = &block.sequences[0].text;
        let pantro4 = &block.sequences[1].text;

        assert!(hg38.contains("----")); // Check for gap presence
        assert!(pantro4.contains("----"));

        // Verify specific gap in panTro4 sequence (the single base deletion)
        assert!(pantro4.contains("TTTTCAA-GC"));
    }

    #[test]
    fn test_real_maf_block() {
        let mut maf = fetch_chr22_maf_reader();

        // Check header first
        let header = maf.read_header().unwrap();
        assert_eq!(header.version, "1");
        assert_eq!(header.scoring, Some("roast.v3.3".to_string()));

        let block = maf.next_block().unwrap().unwrap();

        // Test block score
        assert_eq!(block.score, Some(1058584.0));

        // Test number of sequences and info lines
        assert_eq!(block.sequences.len(), 19); // Updated to correct number
        assert_eq!(block.infos.len(), 18); // One info line per sequence except first

        // Test first sequence (hg38)
        let hg38 = &block.sequences[0];
        assert_eq!(hg38.src, "hg38.chr22");
        assert_eq!(hg38.start, 10510072);
        assert_eq!(hg38.size, 138);
        assert_eq!(hg38.strand, Strand::Forward);
        assert_eq!(hg38.src_size, 50818468);
        assert!(hg38.text.starts_with("TTTTCAAAGC"));

        // Test specific sequences from different parts of the alignment
        let check_species = vec![
            ("hg38.chr22", Strand::Forward),
            ("panTro4.chrUn_GL393523", Strand::Forward),
            ("ponAbe2.chrUn", Strand::Reverse),
            ("myoLuc2.GL429790", Strand::Reverse), // Last sequence
        ];

        for (species, expected_strand) in check_species {
            let seq = block
                .sequences
                .iter()
                .find(|s| s.src == species)
                .unwrap_or_else(|| panic!("Missing sequence for {}", species));
            assert_eq!(seq.strand, expected_strand, "Wrong strand for {}", species);
        }

        // Test alignment consistency
        let alignment_length = block.sequences[0].text.len();
        for seq in &block.sequences {
            assert_eq!(
                seq.text.len(),
                alignment_length,
                "Sequence {} has different length",
                seq.src
            );
        }

        // Check specific gaps that we know should be present
        let pantro4 = block
            .sequences
            .iter()
            .find(|s| s.src == "panTro4.chrUn_GL393523")
            .unwrap();
        assert!(pantro4.text.contains("TTTTCAA-GC"));

        // Test info line for a sequence with non-zero counts
        let equcab_info = block
            .infos
            .iter()
            .find(|i| i.src == "equCab2.chr1")
            .unwrap_or_else(|| panic!("Missing info for equCab2.chr1"));
        assert_eq!(equcab_info.left_count, 242);
    }
}
