use clap::Parser;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use roaring::RoaringBitmap;
use seq_io::fasta::{Reader, Record};
use std::cmp::Ordering;
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
    a_snps: &[RoaringBitmap],
    c_snps: &[RoaringBitmap],
    g_snps: &[RoaringBitmap],
    t_snps: &[RoaringBitmap],
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
                row.push(seq_length - res.len());
            }
            row
        })
        .collect()
}

/// Write pairwise SNP distance matrix to output file
fn write_matrix(
    out: &mut dyn Write,
    indices: &[String],
    pair_snps_by_row: &[Vec<u64>],
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
                        writeln!(
                            out,
                            "{},{},{}",
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
                match i.cmp(&j) {
                    Ordering::Equal => write!(out, "{}{}", sep, 0)?,
                    Ordering::Less => write!(out, "{}{}", sep, pair_snps_by_row[i][j - i - 1])?,
                    Ordering::Greater => write!(out, "{}{}", sep, pair_snps_by_row[j][i - j - 1])?,
                }
            }
            writeln!(out)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use seq_io::fasta::Record;
    use std::str::Utf8Error;

    // Helper struct to implement Record trait for testing
    struct TestRecord {
        id: String,
        seq: Vec<u8>,
    }

    impl Record for TestRecord {
        fn id(&self) -> Result<&str, Utf8Error> {
            std::str::from_utf8(self.id.as_bytes()).map_err(|e| e)
        }
        fn seq(&self) -> &[u8] {
            &self.seq
        }
        fn desc(&self) -> Option<Result<&str, Utf8Error>> {
            None
        }
        fn head(&self) -> &[u8] {
            self.id.as_bytes()
        }
        fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
            write!(
                writer,
                ">{}\n{}",
                self.id,
                String::from_utf8_lossy(&self.seq)
            )
        }
        fn write_wrap<W: Write>(&self, writer: W, _width: usize) -> std::io::Result<()> {
            self.write(writer)
        }
    }

    #[test]
    fn test_build_nucleotide_bitmaps() {
        let record = TestRecord {
            id: "test1".to_string(),
            seq: b"ACGTMRWSYKVHBDN-".to_vec(),
        };

        let (a_sites, c_sites, g_sites, t_sites) = build_nucleotide_bitmaps(&record);

        // Test individual nucleotides
        assert!(a_sites.contains(0)); // A
        assert!(c_sites.contains(1)); // C
        assert!(g_sites.contains(2)); // G
        assert!(t_sites.contains(3)); // T

        // Test combination nucleotides
        assert!(a_sites.contains(4) && c_sites.contains(4)); // M
        assert!(a_sites.contains(5) && g_sites.contains(5)); // R
        assert!(a_sites.contains(6) && t_sites.contains(6)); // W
        assert!(c_sites.contains(7) && g_sites.contains(7)); // S
        assert!(c_sites.contains(8) && t_sites.contains(8)); // Y
        assert!(g_sites.contains(9) && t_sites.contains(9)); // K
        assert!(a_sites.contains(10) && c_sites.contains(10) && g_sites.contains(10)); // V
        assert!(a_sites.contains(11) && c_sites.contains(11) && t_sites.contains(11)); // H
        assert!(c_sites.contains(12) && g_sites.contains(12) && t_sites.contains(12)); // B
        assert!(a_sites.contains(13) && g_sites.contains(13) && t_sites.contains(13)); // D
        assert!(
            a_sites.contains(14)
                && c_sites.contains(14)
                && g_sites.contains(14)
                && t_sites.contains(14)
        ); // N
        assert!(
            a_sites.contains(15)
                && c_sites.contains(15)
                && g_sites.contains(15)
                && t_sites.contains(15)
        ); // -
    }

    #[test]
    #[should_panic]
    fn test_build_nucleotide_bitmaps_invalid_char() {
        let record = TestRecord {
            id: "test1".to_string(),
            seq: b"ACGTX".to_vec(), // X is invalid
        };
        let _ = build_nucleotide_bitmaps(&record);
    }

    #[test]
    fn test_calculate_pairwise_snp_distances() {
        // Create test sequences
        let seq1 = TestRecord {
            id: "seq1".to_string(),
            seq: b"ACGT".to_vec(),
        };
        let seq2 = TestRecord {
            id: "seq2".to_string(),
            seq: b"ACCT".to_vec(),
        };
        let seq3 = TestRecord {
            id: "seq3".to_string(),
            seq: b"AGGT".to_vec(),
        };

        // Build bitmaps for each sequence
        let (a1, c1, g1, t1) = build_nucleotide_bitmaps(&seq1);
        let (a2, c2, g2, t2) = build_nucleotide_bitmaps(&seq2);
        let (a3, c3, g3, t3) = build_nucleotide_bitmaps(&seq3);

        let a_snps = vec![a1, a2, a3];
        let c_snps = vec![c1, c2, c3];
        let g_snps = vec![g1, g2, g3];
        let t_snps = vec![t1, t2, t3];

        let distances = calculate_pairwise_snp_distances(&a_snps, &c_snps, &g_snps, &t_snps, 3, 4);

        // Verify distances
        assert_eq!(distances[0][0], 1); // seq1 vs seq2 (1 SNP difference)
        assert_eq!(distances[0][1], 1); // seq1 vs seq3 (1 SNP difference)
        assert_eq!(distances[1][0], 2); // seq2 vs seq3 (2 SNP differences)
    }

    #[test]
    fn test_write_matrix_dense() {
        let indices = vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()];
        let distances = vec![vec![1, 1], vec![2]];
        let mut output = Vec::new();

        write_matrix(&mut output, &indices, &distances, '\t', false, None).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        let expected = "seq1\t0\t1\t1\nseq2\t1\t0\t2\nseq3\t1\t2\t0\n";
        assert_eq!(output_str, expected);
    }

    #[test]
    fn test_write_matrix_sparse() {
        let indices = vec!["seq1".to_string(), "seq2".to_string(), "seq3".to_string()];
        let distances = vec![vec![1, 1], vec![2]];
        let mut output = Vec::new();

        write_matrix(&mut output, &indices, &distances, ',', true, Some(1)).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        // Should only output distances <= 1
        let expected = "seq1,seq2,1\nseq1,seq3,1\n";
        assert_eq!(output_str, expected);
    }

    #[test]
    fn test_write_matrix_with_indices() {
        let indices = vec!["0".to_string(), "1".to_string(), "2".to_string()];
        let distances = vec![vec![1, 1], vec![2]];
        let mut output = Vec::new();

        write_matrix(&mut output, &indices, &distances, '\t', false, None).unwrap();
        let output_str = String::from_utf8(output).unwrap();

        let expected = "0\t0\t1\t1\n1\t1\t0\t2\n2\t1\t2\t0\n";
        assert_eq!(output_str, expected);
    }

    #[test]
    #[ignore]
    fn test_parallel_build_nucleotide_bitmaps() {
        let record = TestRecord {
            id: "test1".to_string(),
            seq: b"ACGTMRWSYKVHBDN-".to_vec(),
        };
        ThreadPoolBuilder::new()
            .num_threads(2)
            .build_global()
            .unwrap();
        let (a_sites, c_sites, g_sites, t_sites) = build_nucleotide_bitmaps(&record);
        assert!(a_sites.contains(0)); // A
        assert!(c_sites.contains(1)); // C
        assert!(g_sites.contains(2)); // G
        assert!(t_sites.contains(3)); // T
    }
}
