# pairsnp-rs

[![CI](https://github.com/theabhirath/pairsnp-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/theabhirath/pairsnp-rs/actions/workflows/rust.yml)
[![Crates.io Version](https://img.shields.io/crates/v/pairsnp-rs)](https://crates.io/crates/pairsnp-rs)

A Rust implementation for calculating pairwise SNP distance matrices using a multiple sequence alignment. This is heavily inspired by https://github.com/gtonkinhill/pairsnp-cpp and is a project for me to learn Rust.

## Installation

```bash
cargo install pairsnp-rs
```

## Input

The input file should be a multiple sequence alignment in fasta format.

## Output

The output will be a matrix of pairwise SNP distances. By default, this will be output to stdout, but can be piped to a file if the user desires or specified using the `-o` flag.

## Usage

The tool can be run from the command line as:

```console
pairsnp-rs -i input.fasta > output.txt
```

For more information, run `pairsnp-rs --help`.

```console
$ pairsnp-rs --help
Calculate pairwise SNP distances given a multiple sequence alignment.

Usage: pairsnp-rs [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>          Input FASTA file containing multiple sequence alignment
  -o, --output <OUTPUT>        Output file to write pairwise SNP distance matrix (optional)
  -t, --nthreads <NTHREADS>    Number of threads to use (defaults to 1)
  -c, --csv                    Output in CSV format instead of TSV
  -s, --sparse                 Sparse output i.e. only non-zero distances and in s1,s2,dist format
  -d, --threshold <THRESHOLD>  Distance threshold for sparse output
  -x, --indices                Output indices instead of sequence IDs
  -h, --help                   Print help
  -V, --version                Print version
```
