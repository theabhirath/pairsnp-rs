use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use clap::Parser;

use seq_io::fasta::{Reader, Record};
use roaring::RoaringBitmap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// input FASTA file containing multiple sequence alignment
    #[clap(short, long)]
    input: std::path::PathBuf,
    /// output file to write pairwise SNP distance matrix
    #[clap(short, long)]
    output: std::path::PathBuf,
    /// number of threads for Rayon (defaults to all available)
    #[clap(short = 't', long)]
    nthreads: Option<usize>,
}

fn main() -> Result<(), std::io::Error> {
    let args = Cli::parse();
    // configure Rayon global thread pool if user requested a specific count
    if let Some(n) = args.nthreads {
        ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to configure Rayon thread pool");
    }
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
        let length = record.seq().len() as u32;
        if seq_length == 0 {
            seq_length = length;
        } else if length != seq_length {
            panic!("Alignment is not consistent – all sequences must have the same length");
        }
        // store sequence ID
        ids.push(record.id().unwrap().to_string());

        // build nucleotide bitmaps in parallel
        let (a_sites, c_sites, g_sites, t_sites) = record
            .seq()
            .par_iter()
            .enumerate()
            .fold(
                || (
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                ),
                |mut acc, (i, &c)| {
                    let i = i as u32;
                    match c.to_ascii_uppercase() {
                        b'A' => acc.0.insert(i), // A
                        b'C' => acc.1.insert(i), // C
                        b'G' => acc.2.insert(i), // G
                        b'T' => acc.3.insert(i), // T
                        b'M' => { acc.0.insert(i); acc.1.insert(i) } // A or C
                        b'R' => { acc.0.insert(i); acc.2.insert(i) } // A or G
                        b'W' => { acc.0.insert(i); acc.3.insert(i) } // A or T
                        b'S' => { acc.1.insert(i); acc.2.insert(i) } // C or G
                        b'Y' => { acc.1.insert(i); acc.3.insert(i) } // C or T
                        b'K' => { acc.2.insert(i); acc.3.insert(i) } // G or T
                        b'V' => { acc.0.insert(i); acc.1.insert(i); acc.2.insert(i) } // A/C/G
                        b'H' => { acc.0.insert(i); acc.1.insert(i); acc.3.insert(i) } // A/C/T
                        b'D' => { acc.0.insert(i); acc.2.insert(i); acc.3.insert(i) } // A/G/T
                        b'B' => { acc.1.insert(i); acc.2.insert(i); acc.3.insert(i) } // C/G/T
                        b'N' | b'-' => { acc.0.insert(i); acc.1.insert(i); acc.2.insert(i); acc.3.insert(i) } // any
                        b'\n' => false, // ignore line feed
                        _ => panic!("Invalid character, {}, in sequence", c), // invalid character
                    };
                    acc
                },
            )
            .reduce(
                || (
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                    RoaringBitmap::new(),
                ),
                |mut acc1, acc2| {
                    acc1.0 |= acc2.0;
                    acc1.1 |= acc2.1;
                    acc1.2 |= acc2.2;
                    acc1.3 |= acc2.3;
                    acc1
                },
            );

        a_snps.push(a_sites);
        c_snps.push(c_sites);
        g_snps.push(g_sites);
        t_snps.push(t_sites);
        nseqs += 1;
    }

    // calculate pairwise SNP distances in parallel (row‑wise)
    let pair_snps_by_row: Vec<Vec<u32>> = (0..nseqs).into_par_iter()
        .map(|i| {
            let mut row = Vec::with_capacity(nseqs - i - 1);
            for j in i + 1..nseqs {
                let mut res = &a_snps[i] & &a_snps[j];
                res |= &c_snps[i] & &c_snps[j];
                res |= &g_snps[i] & &g_snps[j];
                res |= &t_snps[i] & &t_snps[j];
                row.push(seq_length - res.len() as u32);
            }
            row
        })
        .collect();

    // write pairwise SNP distance matrix to output file
    let mut buffer = BufWriter::new(File::create(args.output)?);
    for i in 0..nseqs {
        write!(buffer, "{}", ids[i])?;
        for j in 0..nseqs {
            if i == j {
                write!(buffer, " 0")?;
            } else if i < j {
                write!(buffer, " {}", pair_snps[(2 * nseqs - 3 - i) * i / 2 + j - 1])?;
            } else {
                write!(buffer, " {}", pair_snps[(2 * nseqs - 3 - j) * j / 2 + i - 1])?;
            }
        }
        write!(buffer, "\n")?;
    }

    Ok(())
}
