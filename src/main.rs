use clap::Parser;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use seq_io::fasta::{Reader, Record};
use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use std::io::stdout;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Input FASTA file containing multiple sequence alignment
    #[clap(short, long)]
    input: std::path::PathBuf,
    /// Output file to write pairwise SNP distance matrix (optional)
    #[clap(short, long)]
    output: Option<std::path::PathBuf>,
    /// Number of threads to use (defaults to 1)
    #[clap(short = 't', long)]
    nthreads: Option<usize>,
    /// Output in CSV format instead of TSV
    #[clap(short, long)]
    csv: bool,
    /// Sparse output i.e. only non-zero distances and in s1,s2,dist format
    #[clap(short, long)]
    sparse: bool,
    /// Distance threshold for sparse output
    #[clap(short = 'd', long)]
    threshold: Option<u64>,
    /// Output indices instead of sequence IDs
    #[clap(short = 'x', long)]
    indices: bool,
}

fn main() -> Result<(), std::io::Error> {
    // parse command line arguments
    let args = Cli::parse();

    // if threshold is set when sparse is not, raise an error
    if args.threshold.is_some() && !args.sparse {
        eprintln!(
            "Error: threshold cannot be set when sparse output is not requested. Use -s or --sparse to request sparse output."
        );
        std::process::exit(1);
    }

    // configure Rayon global thread pool if user requested a specific count, default to 1
    ThreadPoolBuilder::new()
        .num_threads(args.nthreads.unwrap_or(1))
        .build_global()
        .expect("Failed to configure Rayon thread pool");

    // read input FASTA file
    let mut reader = Reader::from_path(args.input)?;

    // initialize variables
    let mut seq_length = 0;
    let mut nseqs = 0;
    let mut ids = Vec::new();
    let mut a_snps = Vec::new();
    let mut c_snps = Vec::new();
    let mut g_snps = Vec::new();
    let mut t_snps = Vec::new();

    // iterate over records in input FASTA file
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let length = record.seq().len() as u64;
        // check if all sequences have the same length
        if seq_length == 0 {
            seq_length = length;
        } else if length != seq_length {
            panic!("Alignment is not consistent – all sequences must have the same length");
        }
        // store sequence ID
        ids.push(record.id().unwrap().to_string());

        // build nucleotide bitmaps in parallel
        let (a_sites, c_sites, g_sites, t_sites) = build_nucleotide_bitmaps(&record);

        // store nucleotide bitmaps for each sequence
        a_snps.push(a_sites);
        c_snps.push(c_sites);
        g_snps.push(g_sites);
        t_snps.push(t_sites);
        nseqs += 1;
    }

    // calculate pairwise SNP distances in parallel (row‑wise)
    let pair_snps_by_row =
        calculate_pairwise_snp_distances(&a_snps, &c_snps, &g_snps, &t_snps, nseqs, seq_length);

    // write pairwise SNP distance matrix to output file if specified, otherwise write to stdout
    let mut buffer: Box<dyn Write> = if args.output.is_some() {
        Box::new(BufWriter::new(File::create(args.output.unwrap())?))
    } else {
        Box::new(stdout())
    };
    let sep = if args.csv { ',' } else { '\t' };
    let indices = if args.indices {
        // use indices if requested
        (0..nseqs).map(|i| i.to_string()).collect()
    } else {
        // default is to use sequence IDs
        ids
    };
    write_matrix(
        &mut buffer,
        &indices,
        &pair_snps_by_row,
        sep,
        args.sparse,
        args.threshold,
    )?;

    Ok(())
}

/// Build nucleotide bitmaps in parallel
fn build_nucleotide_bitmaps<T: Record>(
    record: &T,
) -> (RoaringBitmap, RoaringBitmap, RoaringBitmap, RoaringBitmap) {
    let (a_sites, c_sites, g_sites, t_sites) = record
        .seq()
        .par_iter()
        .enumerate()
        .fold(
            || {
                (
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                )
            },
            |mut acc, (i, &c)| {
                let i = i as u32;
                match c.to_ascii_uppercase() {
                    b'A' => acc.0.insert(i), // A
                    b'C' => acc.1.insert(i), // C
                    b'G' => acc.2.insert(i), // G
                    b'T' => acc.3.insert(i), // T
                    b'M' => {
                        acc.0.insert(i);
                        acc.1.insert(i)
                    } // A or C
                    b'R' => {
                        acc.0.insert(i);
                        acc.2.insert(i)
                    } // A or G
                    b'W' => {
                        acc.0.insert(i);
                        acc.3.insert(i)
                    } // A or T
                    b'S' => {
                        acc.1.insert(i);
                        acc.2.insert(i)
                    } // C or G
                    b'Y' => {
                        acc.1.insert(i);
                        acc.3.insert(i)
                    } // C or T
                    b'K' => {
                        acc.2.insert(i);
                        acc.3.insert(i)
                    } // G or T
                    b'V' => {
                        acc.0.insert(i);
                        acc.1.insert(i);
                        acc.2.insert(i)
                    } // A/C/G
                    b'H' => {
                        acc.0.insert(i);
                        acc.1.insert(i);
                        acc.3.insert(i)
                    } // A/C/T
                    b'D' => {
                        acc.0.insert(i);
                        acc.2.insert(i);
                        acc.3.insert(i)
                    } // A/G/T
                    b'B' => {
                        acc.1.insert(i);
                        acc.2.insert(i);
                        acc.3.insert(i)
                    } // C/G/T
                    b'N' | b'-' => {
                        acc.0.insert(i);
                        acc.1.insert(i);
                        acc.2.insert(i);
                        acc.3.insert(i)
                    } // any
                    b'\n' => false,          // ignore line feed
                    _ => panic!("Invalid character, {}, in sequence", c), // invalid character
                };
                acc
            },
        )
        .reduce(
            || {
                (
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                )
            },
            |mut acc1, acc2| {
                acc1.0 |= acc2.0;
                acc1.1 |= acc2.1;
                acc1.2 |= acc2.2;
                acc1.3 |= acc2.3;
                acc1
            },
        );
    (a_sites, c_sites, g_sites, t_sites)
}

/// Calculate pairwise SNP distance matrix
fn calculate_pairwise_snp_distances(
    a_snps: &Vec<RoaringBitmap>,
    c_snps: &Vec<RoaringBitmap>,
    g_snps: &Vec<RoaringBitmap>,
    t_snps: &Vec<RoaringBitmap>,
    nseqs: usize,
    seq_length: u64,
) -> Vec<Vec<u64>> {
    (0..nseqs)
        .into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(nseqs - i - 1);
            for j in i + 1..nseqs {
                let mut res = &a_snps[i] & &a_snps[j];
                res |= &c_snps[i] & &c_snps[j];
                res |= &g_snps[i] & &g_snps[j];
                res |= &t_snps[i] & &t_snps[j];
                row.push(seq_length - res.len() as u64);
            }
            row
        })
        .collect()
}

/// Write pairwise SNP distance matrix to output file
fn write_matrix(
    out: &mut dyn Write,
    indices: &Vec<String>,
    pair_snps_by_row: &Vec<Vec<u64>>,
    sep: char,
    sparse: bool,
    threshold: Option<u64>,
) -> Result<(), std::io::Error> {
    let nseqs = indices.len();

    // sparse output is of format s1,s2,dist (<=threshold)
    if sparse {
        for i in 0..nseqs {
            for j in 0..nseqs {
                if i < j {
                    let threshold = threshold.unwrap_or(u64::MAX);
                    if pair_snps_by_row[i][j - i - 1] <= threshold {
                        write!(
                            out,
                            "{},{},{}\n",
                            indices[i],
                            indices[j],
                            pair_snps_by_row[i][j - i - 1]
                        )?;
                    }
                }
            }
        }
    } else {
        for i in 0..nseqs {
            write!(out, "{}", indices[i])?;
            for j in 0..nseqs {
                if i == j {
                    write!(out, "{}{}", sep, 0)?;
                } else if i < j {
                    write!(out, "{}{}", sep, pair_snps_by_row[i][j - i - 1])?;
                } else {
                    write!(out, "{}{}", sep, pair_snps_by_row[j][i - j - 1])?;
                }
            }
            writeln!(out)?;
        }
    }

    Ok(())
}
