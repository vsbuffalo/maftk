[![CI](https://github.com/vsbuffalo/maftk/actions/workflows/ci.yml/badge.svg)](https://github.com/vsbuffalo/maftk/actions/workflows/ci.yml)

# maftk (formerly biomaf) — experimental MAF parser and toolkit in Rust 🦀

![maftk screenshot](screenshot.png)


A very small Multiple Alignment Format (MAF) parsing library and tool in Rust.

Subcommands:

    - `maftk split`: Split a chromosome's MAF file into a series of MAF files
      per region alignment.
    - `maftk binary <in.maf.gz>`: indexes and serializes the alignment blocks using [hgindex](https://github.com/vsbuffalo/hgindex).
    - `maftk query chr22 10960739 10960750`: queries the binary alignment blocks using fast hierarchical indexing.

## Installation

For now, this isn't on [crates.io](http://crates.io). So, use:

```
$ cargo install --git https://github.com/vsbuffalo/maftk
```
