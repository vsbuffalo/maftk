[![CI](https://github.com/vsbuffalo/biomaf/actions/workflows/ci.yml/badge.svg)](https://github.com/vsbuffalo/biomaf/actions/workflows/ci.yml)

# biomaf — experimental MAF parser in Rust 🦀

![biomaf screenshot](screenshot.png)


A very small Multiple Alignment Format (MAF) parsing library and tool in Rust.

Subcommands:

    - `biomaf split`: Split a chromosome's MAF file into a series of MAF files
      per region alignment.
    - `biomaf binary <in.maf.gz>`: indexes and serializes the alignment blocks using [hgindex](https://github.com/vsbuffalo/hgindex).
    - `biomaf query chr22 10960739 10960750`: queries the binary alignment blocks using fast hierarchical indexing.

## Installation

For now, this isn't on [crates.io](http://crates.io). So, use:

```
$ cargo install --git https://github.com/vsbuffalo/biomaf
```
