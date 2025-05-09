use seq_io::fasta::{Reader, Record};
use roaring::RoaringBitmap;

fn main() {
    let fasta_file = "data/good.fasta";
    let mut reader = Reader::from_path(fasta_file).unwrap();

    let mut seq_length = 0;
    let mut nseqs = 0;
    let mut ids = Vec::new();
    let mut a_snps = Vec::new();
    let mut c_snps = Vec::new();
    let mut g_snps = Vec::new();
    let mut t_snps = Vec::new();

    for record in reader.records() {
        let record = record.unwrap();
        let length = record.seq().len() as u32;
        if seq_length == 0 {
            seq_length = length;
        } else if length != seq_length {
            panic!("Alignment is not consistent – all sequences must have the same length");
        }
        ids.push(record.id().unwrap().to_string());

        let mut a_sites = RoaringBitmap::new();
        let mut c_sites = RoaringBitmap::new();
        let mut g_sites = RoaringBitmap::new();
        let mut t_sites = RoaringBitmap::new();

        for (i, c) in record.seq().iter().enumerate() {
            match c {
                b'A' => a_sites.insert(i as u32),
                b'C' => c_sites.insert(i as u32),
                b'G' => g_sites.insert(i as u32),
                b'T' => t_sites.insert(i as u32),
                b'M' => {
                    a_sites.insert(i as u32);
                    c_sites.insert(i as u32)
                }
                b'R' => {
                    a_sites.insert(i as u32);
                    g_sites.insert(i as u32)
                }
                b'W' => {
                    a_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'S' => {
                    c_sites.insert(i as u32);
                    g_sites.insert(i as u32)
                }
                b'Y' => {
                    c_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'K' => {
                    g_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'V' => {
                    a_sites.insert(i as u32);
                    c_sites.insert(i as u32);
                    g_sites.insert(i as u32)
                }
                b'H' => {
                    a_sites.insert(i as u32);
                    c_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'D' => {
                    a_sites.insert(i as u32);
                    g_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'B' => {
                    c_sites.insert(i as u32);
                    g_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                b'N' | b'-' => {
                    a_sites.insert(i as u32);
                    c_sites.insert(i as u32);
                    g_sites.insert(i as u32);
                    t_sites.insert(i as u32)
                }
                _ => panic!("Invalid character in sequence")
            };
        }
        a_snps.push(a_sites);
        c_snps.push(c_sites);
        g_snps.push(g_sites);
        t_snps.push(t_sites);
        nseqs += 1;
    }

    let mut pair_snps = Vec::new();
    for i in 0..nseqs {
        let mut res;
        for j in i + 1..nseqs {
            res = &a_snps[i] & &a_snps[j];
            res |= &c_snps[i] & &c_snps[j];
            res |= &g_snps[i] & &g_snps[j];
            res |= &t_snps[i] & &t_snps[j];
            pair_snps.push(seq_length - res.len() as u32);
        }
    }

    for i in 0..nseqs {
        print!("{}", ids[i]);
        for j in 0..nseqs {
            if i == j {
                print!(" 0");
            } else if i < j {
                print!(" {}", pair_snps[(nseqs*(nseqs-1)/2) - (nseqs-i)*((nseqs-i)-1)/2 + j - i - 1]);
            } else {
                print!(" {}", pair_snps[(nseqs*(nseqs-1)/2) - (nseqs-j)*((nseqs-j)-1)/2 + i - j - 1]);
            }
        }
        println!();
    }
}
