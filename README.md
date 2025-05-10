# pairsnp-rs

A Rust implementation for calculating pairwise SNP distance matrices using a multiple sequence alignment. This is heavily inspired by https://github.com/gtonkinhill/pairsnp-cpp and is a project for me to learn Rust.

> [!WARNING]
> This is a work in progress and not yet complete.

## Installation

```bash
cargo install --git https://github.com/theabhirath/pairsnp-rs.git
```

## Input

The input file should be a multiple sequence alignment in fasta format.

## Output

The output will be a matrix of pairwise SNP distances.

## Usage

```bash
pairsnp-rs -i input.fasta -o output.txt
```

For more information, run `pairsnp-rs --help`.
