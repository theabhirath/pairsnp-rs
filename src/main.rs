use std::fs::File;
use std::io::BufWriter;
use std::io::prelude::*;
use clap::Parser;

use seq_io::fasta::{Reader, Record};
use roaring::RoaringBitmap;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// input FASTA file containing multiple sequence alignment
    #[clap(short, long)]
    input: std::path::PathBuf,
    /// output file to write pairwise SNP distance matrix
    #[clap(short, long)]
    output: std::path::PathBuf,
}

fn main() -> Result<(), std::io::Error> {
    let args = Cli::parse();
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
    for record in reader.records() {
        let record = record.unwrap();
        let length = record.seq().len() as u32;
        if seq_length == 0 {
            seq_length = length;
        } else if length != seq_length {
            panic!("Alignment is not consistent – all sequences must have the same length");
        }
        // store sequence ID
        ids.push(record.id().unwrap().to_string());

        // initialize bitmaps for each nucleotide
        let mut a_sites = RoaringBitmap::new();
        let mut c_sites = RoaringBitmap::new();
        let mut g_sites = RoaringBitmap::new();
        let mut t_sites = RoaringBitmap::new();

        // iterate over nucleotides in sequence
        for (i, c) in record.seq().iter().enumerate() {
            let i = i as u32; // convert to u32 for bitmap
            match c.to_ascii_uppercase() {
                b'A' => a_sites.insert(i), // A
                b'C' => c_sites.insert(i), // C
                b'G' => g_sites.insert(i), // G
                b'T' => t_sites.insert(i), // T
                b'M' => { // A or C
                    a_sites.insert(i);
                    c_sites.insert(i)
                }
                b'R' => { // A or G
                    a_sites.insert(i);
                    g_sites.insert(i)
                }
                b'W' => { // A or T
                    a_sites.insert(i);
                    t_sites.insert(i)
                }
                b'S' => { // C or G
                    c_sites.insert(i);
                    g_sites.insert(i)
                }
                b'Y' => { // C or T
                    c_sites.insert(i);
                    t_sites.insert(i)
                }
                b'K' => { // G or T
                    g_sites.insert(i);
                    t_sites.insert(i)
                }
                b'V' => { // A, C, or G
                    a_sites.insert(i);
                    c_sites.insert(i);
                    g_sites.insert(i)
                }
                b'H' => { // A, C, or T
                    a_sites.insert(i);
                    c_sites.insert(i);
                    t_sites.insert(i)
                }
                b'D' => { // A, G, or T
                    a_sites.insert(i);
                    g_sites.insert(i);
                    t_sites.insert(i)
                }
                b'B' => { // C, G, or T
                    c_sites.insert(i);
                    g_sites.insert(i);
                    t_sites.insert(i)
                }
                b'N' | b'-' => { // A, C, G, or T
                    a_sites.insert(i);
                    c_sites.insert(i);
                    g_sites.insert(i);
                    t_sites.insert(i)
                }
                _ => panic!("Invalid character in sequence") // invalid character
            };
        }
        a_snps.push(a_sites);
        c_snps.push(c_sites);
        g_snps.push(g_sites);
        t_snps.push(t_sites);
        nseqs += 1;
    }

    // calculate pairwise SNP distances
    let mut pair_snps = Vec::new();
    for i in 0..nseqs {
        for j in i + 1..nseqs {
            let mut res = &a_snps[i] & &a_snps[j];
            res |= &c_snps[i] & &c_snps[j];
            res |= &g_snps[i] & &g_snps[j];
            res |= &t_snps[i] & &t_snps[j];
            pair_snps.push(seq_length - res.len() as u32);
        }
    }

    // write pairwise SNP distance matrix to output file
    let mut buffer = BufWriter::new(File::create(args.output)?);
    for i in 0..nseqs {
        write!(buffer, "{}", ids[i])?;
        for j in 0..nseqs {
            if i == j {
                write!(buffer, " 0")?;
            } else if i < j {
                write!(buffer, " {}", pair_snps[(nseqs*(nseqs-1)/2) - (nseqs-i)*((nseqs-i)-1)/2 + j - i - 1])?;
            } else {
                write!(buffer, " {}", pair_snps[(nseqs*(nseqs-1)/2) - (nseqs-j)*((nseqs-j)-1)/2 + i - j - 1])?;
            }
        }
        write!(buffer, "\n")?;
    }

    Ok(())
}
