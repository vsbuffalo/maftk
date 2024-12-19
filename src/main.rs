use biomaf::binary::{convert_to_binary, IndexedMafDb};
use biomaf::index::FastIndex;
use biomaf::io::{InputFile, MafReader, OutputFile};
use biomaf::statistics::{calc_alignment_statistics, AlignmentStatistics};
use clap::Parser;
use csv::ReaderBuilder;
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs::create_dir_all;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use tracing::{info, trace};

#[derive(Parser)]
#[command(name = "maftools")]
#[command(about = "Tools for working with Multiple Alignment Format (MAF) files")]
struct Cli {
    /// The verbosity level: ("info", "debug", "warn", or "error")
    #[arg(short, long, default_value = "info")]
    verbosity: String,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Parser)]
enum Commands {
    /// Split a MAF file into chunks by region
    Split {
        /// Input MAF file
        #[arg(value_name = "input.maf")]
        input: PathBuf,

        /// Output directory for split MAF files
        #[arg(short, long, default_value = "split_mafs")]
        output_dir: PathBuf,

        /// Minimum length of alignment block to include (default: 0)
        #[arg(short, long, default_value = "0")]
        min_length: u64,
    },
    /// Convert MAF to binary 3-bit format with index
    Binary {
        /// Input MAF file
        #[arg(value_name = "input.maf")]
        input: PathBuf,

        /// Output directory
        #[arg(short, long, default_value = "maf.mmdb")]
        output_dir: PathBuf,

        /// Minimum length of alignment block to include
        #[arg(short, long, default_value = "0")]
        min_length: u64,
    },
    /// Get alignments at a specific position or range
    Query {
        /// Chromosome name (e.g., chr22)
        #[arg(value_name = "chromosome")]
        chromosome: String,

        /// Start position (1-based)
        #[arg(value_name = "start")]
        start: u32,

        /// End position (optional, defaults to start+1)
        #[arg(value_name = "end")]
        end: Option<u32>,

        /// Directory containing binary MAF files
        #[arg(short, long, default_value = "maf.mmdb")]
        data_dir: PathBuf,
    },
    /// Calculate an alignment statistics table for all alignments
    /// in the supplied BED file.
    Stats {
        /// Input BED3 file
        #[arg(value_name = "regions.bed")]
        regions: PathBuf,

        /// Output TSV file
        #[arg(short, long, value_name = "output.tsv.gz")]
        output: PathBuf,

        /// Comma-separated list of species.
        #[arg(short, long, value_name = "hg38,panTro4,ponAbe2")]
        species: String,

        /// Directory containing binary MAF files
        #[arg(short, long, default_value = "maf.mmdb")]
        data_dir: PathBuf,

        /// Print timing information
        #[arg(short, long)]
        verbose: bool,
    },

    /// Debug an index file's structure
    DebugIndex {
        /// Path to index file
        #[arg(value_name = "index_file")]
        path: PathBuf,
    },
}

fn write_maf_header(file: &mut File, header: &biomaf::io::MafHeader) -> std::io::Result<()> {
    write!(file, "##maf version={}", header.version)?;
    if let Some(scoring) = &header.scoring {
        write!(file, " scoring={}", scoring)?;
    }
    if let Some(program) = &header.program {
        write!(file, " program={}", program)?;
    }
    writeln!(file)?;
    writeln!(file) // Extra newline after header
}

fn write_block(file: &mut File, block: &biomaf::io::AlignmentBlock) -> std::io::Result<()> {
    // Write alignment line
    write!(file, "a")?;
    if let Some(score) = block.score {
        write!(file, " score={}", score)?;
    }
    if let Some(pass) = block.pass {
        write!(file, " pass={}", pass)?;
    }
    writeln!(file)?;

    // Write sequence lines
    for seq in &block.sequences {
        writeln!(
            file,
            "s {} {} {} {} {} {}",
            seq.src,
            seq.start,
            seq.size,
            if seq.strand == biomaf::io::Strand::Forward {
                "+"
            } else {
                "-"
            },
            seq.src_size,
            seq.text
        )?;
    }

    // Write info lines
    for info in &block.infos {
        writeln!(
            file,
            "i {} {} {} {} {}",
            info.src,
            match info.left_status {
                biomaf::io::StatusChar::Contiguous => "C",
                biomaf::io::StatusChar::Intervening => "I",
                biomaf::io::StatusChar::New => "N",
                biomaf::io::StatusChar::NewBridged => "n",
                biomaf::io::StatusChar::Missing => "M",
                biomaf::io::StatusChar::Tandem => "T",
            },
            info.left_count,
            match info.right_status {
                biomaf::io::StatusChar::Contiguous => "C",
                biomaf::io::StatusChar::Intervening => "I",
                biomaf::io::StatusChar::New => "N",
                biomaf::io::StatusChar::NewBridged => "n",
                biomaf::io::StatusChar::Missing => "M",
                biomaf::io::StatusChar::Tandem => "T",
            },
            info.right_count,
        )?;
    }

    writeln!(file) // Extra newline between blocks
}

fn get_block_range(block: &biomaf::io::AlignmentBlock) -> Option<(String, u64, u64)> {
    // Get the first sequence in the block
    let first_seq = block.sequences.first()?;

    // Extract chromosome from source (assumes format like "hg38.chr1")
    let chr = first_seq.src.split('.').nth(1)?;

    // Calculate end position
    let end = first_seq.start + first_seq.size;

    Some((chr.to_string(), first_seq.start, end))
}

struct ProcessStats {
    total_blocks: u64,
    filtered_blocks: u64,
}

fn split_maf(
    input: &Path,
    output_dir: &Path,
    min_length: u64,
) -> Result<ProcessStats, Box<dyn std::error::Error>> {
    // Create output directory if it doesn't exist
    create_dir_all(output_dir)?;

    let mut maf = MafReader::from_file(input)?;
    let header = maf.read_header()?.clone();

    let mut current_file: Option<(File, String)> = None;
    let mut filtered_blocks = 0;
    let mut processed_blocks = 0;

    while let Some(block) = maf.next_block()? {
        processed_blocks += 1;

        // Check if block meets minimum length requirement
        if block.sequences.first().map_or(0, |seq| seq.size) < min_length {
            filtered_blocks += 1;
            continue;
        }

        if let Some((chr, start, end)) = get_block_range(&block) {
            let filename = format!("{}_{}_{}.maf", chr, start, end);
            let filepath = output_dir.join(&filename);

            // If this is a new file or different region
            if current_file
                .as_ref()
                .map(|(_, name)| name != &filename)
                .unwrap_or(true)
            {
                // Create new file
                let mut file = File::create(&filepath)?;
                write_maf_header(&mut file, &header)?;
                current_file = Some((file, filename));
            }

            // Write block to current file
            if let Some((ref mut file, _)) = current_file {
                write_block(file, &block)?;
            }

            if processed_blocks % 100 == 0 {
                eprint!(".");
            }
        }
    }

    Ok(ProcessStats {
        total_blocks: processed_blocks,
        filtered_blocks,
    })
}

pub fn query_command(
    chrom: &str,
    start: u32,
    end: u32,
    data_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let db = IndexedMafDb::open(&chrom, data_dir)?;
    let total_time = Instant::now();
    let alignments = db.get_alignments(start, end)?;
    let mut all_species = HashSet::new();
    for (species, seq) in &alignments {
        println!("{}\t{}", species, seq);
        all_species.insert(species.clone());
    }

    // let region_stats = calc_alignment_statistics(alignments, &all_species, chrom.to_string(), start, end);
    // dbg!(region_stats);
    info!("Total elapsed {:?}", total_time.elapsed());
    Ok(())
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
    verbose: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let total_time = Instant::now();

    // Estimate total records for progress bar
    let estimated_records = estimate_record_count(regions)?;
    let progress = ProgressBar::new(estimated_records);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} records ({eta})")?
        .progress_chars("#>-"));

    // Setup input for reading BED file
    let input_file = InputFile::new(regions);
    let buf_reader = input_file.reader()?;

    let output_file = OutputFile::new(output);
    let mut stats_writer = AlignmentStatistics::new(output_file, species.clone())?;

    let mut reader = ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(buf_reader);

    // Cache indexed MAF databases by chromosome
    let mut index_map: HashMap<String, IndexedMafDb> = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let chrom = record[0].to_string();
        let start: u32 = record[1].parse()?;
        let end: u32 = record[2].parse()?;

        // Get or create index for this chromosome
        let index = match index_map.entry(chrom.clone()) {
            Entry::Occupied(o) => o.into_mut(),
            Entry::Vacant(v) => v.insert(IndexedMafDb::open(&chrom, data_dir)?),
        };

        // Using indexed query to find overlapping blocks
        let find_start = Instant::now();
        let alignments = index.get_alignments(start, end)?;
        if verbose {
            trace!(elapsed = ?find_start.elapsed(), "Finding blocks");
        }

        if !alignments.is_empty() {
            // Calculate statistics for this region
            let region_stats = calc_alignment_statistics(alignments, &species, chrom, start, end);
            // Write statistics to output file
            stats_writer.write_stats(&region_stats)?;
        }

        progress.inc(1);
    }

    progress.finish_with_message("Processing complete");

    if verbose {
        info!("Total elapsed {:?}", total_time.elapsed());
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // Set up tracing level based on verbosity argument
    let level = match cli.verbosity.to_lowercase().as_str() {
        "trace" => tracing::Level::TRACE,
        "debug" => tracing::Level::DEBUG,
        "info" => tracing::Level::INFO,
        "warn" => tracing::Level::WARN,
        "error" => tracing::Level::ERROR,
        _ => {
            eprintln!("Invalid verbosity level: {}", cli.verbosity);
            std::process::exit(1);
        }
    };

    // Initialize logger with default level INFO
    tracing_subscriber::fmt().with_max_level(level).init();

    match cli.command {
        Commands::Binary {
            input,
            output_dir,
            min_length,
        } => {
            convert_to_binary(&input, &output_dir, min_length)?;
        }
        Commands::Split {
            input,
            output_dir,
            min_length,
        } => {
            let stats = split_maf(&input, &output_dir, min_length)?;
            println!("\nProcessing complete:");
            println!("Total blocks processed: {}", stats.total_blocks);
            println!(
                "Blocks filtered (length < {}): {}",
                min_length, stats.filtered_blocks
            );
            println!(
                "Blocks written: {}",
                stats.total_blocks - stats.filtered_blocks
            );
        }
        Commands::Query {
            chromosome,
            start,
            end,
            data_dir,
        } => {
            let end = end.unwrap_or(start + 1);
            if end <= start {
                eprintln!("End position must be greater than start position");
                std::process::exit(1);
            }
            query_command(&chromosome, start, end, &data_dir)?;
        }
        Commands::Stats {
            regions,
            output,
            species,
            data_dir,
            verbose,
        } => {
            let species = HashSet::from_iter(species.split(',').map(String::from));
            stats_command(&regions, &output, species, &data_dir, verbose)?
        }
        Commands::DebugIndex { path } => {
            let index = FastIndex::open(&path)?;
            index.debug_print();
        }
    }

    Ok(())
}
