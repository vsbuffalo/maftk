/// io.rs
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek};
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
///
/// Note: This is the MAF parsing-level struct, compared to AlignedSequence which intended more for
/// the binary-level representation.
/// Thus these are optimized for different parts of our pipeline.
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

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const DEFAULT_BUFFER_SIZE: usize = 128 * 1024;

#[derive(Error, Debug)]
pub enum IoError {
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Invalid or corrupted gzip header")]
    InvalidGzipHeader,
}

pub struct InputStream {
    filepath: PathBuf,
}

impl InputStream {
    pub fn new(filepath: &Path) -> Self {
        Self {
            filepath: filepath.into(),
        }
    }

    fn is_gzipped(file: &mut File) -> Result<bool, IoError> {
        let mut header = [0u8; 2];
        // Read the first two bytes
        file.read_exact(&mut header)?;
        // Reset the file pointer
        file.rewind()?;
        Ok(header == GZIP_MAGIC)
    }

    pub fn reader(&self) -> Result<BufReader<Box<dyn Read>>, IoError> {
        let mut file = File::open(&self.filepath)?;
        let reader: Box<dyn Read> = if Self::is_gzipped(&mut file)? {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(BufReader::with_capacity(DEFAULT_BUFFER_SIZE, reader))
    }

    pub fn has_header(&self, expect: &str) -> Result<bool, IoError> {
        let mut buf_reader = self.reader()?;
        let mut first_line = String::new();
        buf_reader.read_line(&mut first_line)?;
        Ok(first_line.starts_with(expect))
    }
}

#[derive(Clone)]
pub struct OutputStreamBuilder {
    filepath: Option<PathBuf>,
    buffer_size: usize,
    compression_level: Compression,
}

impl Default for OutputStreamBuilder {
    fn default() -> Self {
        Self {
            filepath: None,
            buffer_size: DEFAULT_BUFFER_SIZE,
            compression_level: Compression::default(),
        }
    }
}

impl OutputStreamBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn filepath(mut self, path: Option<impl AsRef<Path>>) -> Self {
        self.filepath = path.map(|p| p.as_ref().to_path_buf());
        self
    }

    pub fn buffer_size(mut self, size: usize) -> Self {
        self.buffer_size = size;
        self
    }

    pub fn compression_level(mut self, level: Option<Compression>) -> Self {
        if let Some(level) = level {
            self.compression_level = level;
        } else {
            self.compression_level = Compression::best();
        }
        self
    }

    pub fn build(self) -> OutputStream {
        OutputStream {
            filepath: self.filepath,
            buffer_size: self.buffer_size,
            compression_level: self.compression_level,
        }
    }
}

pub struct OutputStream {
    filepath: Option<PathBuf>,
    buffer_size: usize,
    compression_level: Compression,
}

impl OutputStream {
    pub fn new(filepath: Option<impl AsRef<Path>>) -> Self {
        OutputStreamBuilder::new().filepath(filepath).build()
    }

    pub fn builder() -> OutputStreamBuilder {
        OutputStreamBuilder::new()
    }

    fn should_compress(&self) -> bool {
        self.filepath
            .as_ref()
            .map_or(false, |p| p.extension().map_or(false, |ext| ext == "gz"))
    }

    pub fn writer(&self) -> Result<Box<dyn Write>, Error> {
        match &self.filepath {
            Some(path) => {
                let file = File::create(path)?;
                let writer: Box<dyn Write> = if self.should_compress() {
                    Box::new(BufWriter::with_capacity(
                        self.buffer_size,
                        GzEncoder::new(file, self.compression_level),
                    ))
                } else {
                    Box::new(BufWriter::with_capacity(self.buffer_size, file))
                };
                Ok(writer)
            }
            None => Ok(Box::new(BufWriter::with_capacity(
                self.buffer_size,
                io::stdout(),
            ))),
        }
    }
}
