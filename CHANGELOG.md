# Changelog

## [0.1.1](https://github.com/theabhirath/pairsnp-rs/releases/tag/v0.1.1) - 2025-05-13

### Changed

- Handle both gzipped and uncompressed FASTA files using the `flate2` crate.

## [0.1.0](https://github.com/theabhirath/pairsnp-rs/releases/tag/v0.1.0) - 2025-05-12

### Description

- Initial release of the project.

### Changed

- Created a CLI tool that calculates pairwise SNP distances given a multiple sequence alignment FASTA file.
- Use `seq_io` for reading and writing FASTA files.
- Uses `rayon` for parallelization and `clap` for CLI parsing.
- Uses `roaring` for efficient bitmaps.
