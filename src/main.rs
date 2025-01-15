// main.rs
use clap::Parser;
use hgindex::BinningIndex;
use maftk::binary::{
    calc_block_statistics, convert_to_binary, convert_to_binary_glob, print_alignments,
    print_block_statistics, query_alignments,
};
use maftk::binary::{stats_command_range, stats_command_ranges, SpeciesDictionary};
use maftk::io::{MafReader, OutputStream};
use polars::io::csv::write::CsvWriter;
use polars::io::SerWriter;
use std::collections::HashSet;
use std::fs::create_dir_all;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use tracing::warn;

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
    /// Convert MAF to binary format with index
    Binary {
        /// Input MAF file, or multiple MAF files with glob (e.g. "mafs/chr*.maf")
        /// Files can be gzipped.
        #[arg(value_name = "input.maf")]
        input: String,

        /// Output directory
        #[arg(short, long, default_value = "maf.mmdb")]
        output_dir: PathBuf,

        /// Minimum length of alignment block to include
        #[arg(short, long, default_value = "0")]
        min_length: u64,
    },
    /// Get alignments at a specific position or region, and either write to a directory (one FASTA
    /// file per alignment) or print alignments.
    Query {
        /// Single genomic region in format chr:start-end (1-based coordinates)
        #[arg(value_name = "chr17:7661779-7687538")]
        region: String,

        /// Directory containing binary MAF files
        #[arg(short, long, default_value = "maf.mmdb")]
        data_dir: PathBuf,

        /// Optional directory to output FASTA files or TSV block statistics files
        /// for all overlapping alignment regions.
        #[arg(short, long)]
        output_dir: Option<PathBuf>,

        /// Don't color
        #[arg(short, long)]
        no_color: bool,

        /// Only consider intersecting region only (i.e. inner join,
        /// rather than the default right outer join).
        #[arg(short, long)]
        intersect_only: bool,

        /// Output (if --output-dir set) or print statistics of alignments only
        /// (not actual alignments).
        #[arg(short, long)]
        stats: bool,

        /// Show all pairwise alignments with --print-stats, not just
        /// the pairwise statistics with the reference species/source.
        #[arg(short, long)]
        all_pairwise: bool,
    },
    /// Calculate an alignment statistics table for all alignments
    /// in the supplied BED file. Note that each entry is from
    /// *one* block; merging blocks will not automatically be combined,
    /// since this can be done downstream on the default output (counts only
    /// data, e.g. without --include-rates). Do not merge rates in post-processsing
    /// with --include-rates unless you know how to do this correctly.
    Stats {
        /// Single genomic region in format chr:start-end (1-based coordinates)
        #[arg(
            value_name = "chr17:7661779-7687538",
            required_unless_present_any = ["regions"],
            conflicts_with_all = ["regions"]
        )]
        region: Option<String>,

        /// Input BED file for batch region queries
        #[arg(
            long,
            value_name = "regions.bed",
            required_unless_present_any = ["region"],
            conflicts_with_all = ["region"]
        )]
        regions: Option<PathBuf>,

        /// Output TSV file
        #[arg(short, long, value_name = "output.tsv.gz")]
        output: Option<PathBuf>,

        /// Comma-separated list of species.
        #[arg(short, long, value_name = "hg38,panTro4,ponAbe2")]
        species: String,

        /// Directory containing binary MAF files
        #[arg(short, long, default_value = "maf.mmdb")]
        data_dir: PathBuf,

        /// Whether to include rates.
        #[arg(long)]
        include_rates: bool,
    },

    /// Debug an index file's structure
    DebugIndex {
        /// Path to index file
        #[arg(value_name = "index_file")]
        path: PathBuf,
    },
}

fn write_maf_header(file: &mut File, header: &maftk::io::MafHeader) -> std::io::Result<()> {
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

fn write_block(file: &mut File, block: &maftk::io::AlignmentBlock) -> std::io::Result<()> {
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
            if seq.strand == maftk::io::Strand::Forward {
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
                maftk::io::StatusChar::Contiguous => "C",
                maftk::io::StatusChar::Intervening => "I",
                maftk::io::StatusChar::New => "N",
                maftk::io::StatusChar::NewBridged => "n",
                maftk::io::StatusChar::Missing => "M",
                maftk::io::StatusChar::Tandem => "T",
            },
            info.left_count,
            match info.right_status {
                maftk::io::StatusChar::Contiguous => "C",
                maftk::io::StatusChar::Intervening => "I",
                maftk::io::StatusChar::New => "N",
                maftk::io::StatusChar::NewBridged => "n",
                maftk::io::StatusChar::Missing => "M",
                maftk::io::StatusChar::Tandem => "T",
            },
            info.right_count,
        )?;
    }

    writeln!(file) // Extra newline between blocks
}

fn get_block_region(block: &maftk::io::AlignmentBlock) -> Option<(String, u64, u64)> {
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

        if let Some((chr, start, end)) = get_block_region(&block) {
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
    tracing_subscriber::fmt()
        .without_time()
        .with_max_level(level)
        .init();

    match cli.command {
        Commands::Binary {
            input,
            output_dir,
            min_length,
        } => {
            if input.contains('*') || input.contains('?') {
                // Use glob version for wildcard patterns
                convert_to_binary_glob(&input, &output_dir, min_length)?;
            } else {
                // Use original version for single files
                convert_to_binary(Path::new(&input), &output_dir, min_length)?;
            }
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
            region,
            data_dir,
            output_dir,
            no_color,
            intersect_only,
            stats,
            all_pairwise,
        } => {
            let (chromosome, start, end) = parse_range(&region)?;

            // Get the overlapping alignment blocks
            let (blocks, species_dict) =
                query_alignments(&data_dir, chromosome, start, end, intersect_only)?;

            if stats {
                for (i, block) in blocks.iter().enumerate() {
                    let focal_species = if all_pairwise {
                        None
                    } else {
                        // the first sequence in the block is always the reference
                        let idx = block.sequence_metadata[0].species_idx;
                        Some(species_dict.get_species(idx).unwrap())
                    };

                    if let Some(stats) = block.calc_stats(Some(start), Some(end), None) {
                        let mut df = calc_block_statistics(&stats, &species_dict, focal_species)?;
                        if let Some(dir) = output_dir.as_ref() {
                            create_dir_all(dir)?;
                            let filepath = format!("block_{}.tsv.gz", i);
                            let path = dir.join(&filepath);
                            let writer = OutputStream::builder()
                                .filepath(Some(path))
                                .compression_level(None)
                                .build();
                            CsvWriter::new(writer.writer()?)
                                .include_header(true)
                                .with_separator(b'\t')
                                .finish(&mut df)?;
                        } else {
                            print_block_statistics(&df, i);
                        }
                    } else {
                        warn!("Block {} didn't have statistics", i);
                    }
                }
            } else if let Some(dir) = output_dir {
                create_dir_all(&dir)?;
                for (i, block) in blocks.iter().enumerate() {
                    let filepath = format!("block_{}.fa.gz", i);
                    let path = dir.join(&filepath);
                    block.write_fasta(Some(&path), &species_dict, chromosome)?;
                }
            } else {
                // We print to screen rather than write to FASTA directory.
                print_alignments(blocks, &species_dict, !no_color);
            }
        }
        Commands::Stats {
            region,
            regions,
            output,
            species,
            data_dir,
            include_rates,
        } => {
            let species = HashSet::from_iter(species.split(',').map(String::from));

            if let Some(region_str) = region {
                // Handle single range case
                let (seqname, start, end) = parse_range(&region_str)?;
                stats_command_range(
                    seqname,
                    start,
                    end,
                    output.as_deref(),
                    species,
                    &data_dir,
                    include_rates,
                )?;
            } else if let Some(ranges_path) = regions {
                // Handle multiple regions from BED file
                stats_command_ranges(
                    &ranges_path,
                    output.as_deref(),
                    species,
                    &data_dir,
                    include_rates,
                )?;
            } else {
                panic!("Invalid options");
            }
        }
        Commands::DebugIndex { path } => {
            let _index = BinningIndex::<SpeciesDictionary>::open(&path)?;
            todo!()
        }
    }

    Ok(())
}

// Also converts from user-space 1-based right inclusive to coder space
// 0-based right exclusive interval format.
fn parse_range(range: &str) -> Result<(&str, u32, u32), Box<dyn std::error::Error>> {
    let parts: Vec<&str> = range.split(':').collect();
    if parts.len() != 2 {
        return Err("Invalid range format. Expected seqname:start-end.".into());
    }

    let seqname = parts[0];
    let coords: Vec<&str> = parts[1].split('-').collect();
    if coords.len() != 2 {
        return Err("Invalid range format. Expected start-end.".into());
    }

    let start: u32 = coords[0].parse().map_err(|_| "Invalid start coordinate.")?;
    let end: u32 = coords[1].parse().map_err(|_| "Invalid end coordinate.")?;

    // Convert to 0-based coordinates like in the query command
    let start = start
        .checked_sub(1)
        .ok_or("Start coordinate must be greater than 0")?;

    Ok((seqname, start, end))
}
