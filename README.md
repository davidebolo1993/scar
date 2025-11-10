# scar
SCAR - Single-Cell ATAC Router



**SCAR** is a high-performance C++ tool for routing single-cell ATAC-seq reads from pooled BAM files into celltype- and donor-specific output files with optional peak-based filtering.

## Features

- **Fast BAM routing**: Processes >500k reads/second with optimized binary search and buffering
- **Celltype-aware peak filtering**: Per-celltype BED file support for peak-based read selection
- **Duplicate removal**: Optional PCR/optical duplicate filtering


## Installation

### Requirements

- C++17 compatible compiler (gcc 7+)
- CMake 3.14+
- Libraries: zlib, bzip2, liblzma

### Build

```bash
git clone https://github.com/davidebolo1993/scar.git
cd scar
mkdir build && cd build
cmake ..
make
```

### Singularity

```bash
singularity pull --dir . docker://davidebolo1993/scar:latest
```

The compiled binary will be at `build/scar`.

## Quick Start

### Step 1: Filter metadata and generate celltype-specific BED files

```bash
scar filter
-i metadata.tsv
-o metadata.filtered.tsv
-c celltypes_LEV1
-m 10
-p peak_matrix.tsv
-b celltype_beds
```

**Input files:**
- `metadata.tsv`: Tab-separated file with columns: `barcode`, `pool_id`, `donor_id`, `celltypes_LEV1`, `celltypes_LEV2`, `bamfile`
- `peak_matrix.tsv`: Peak matrix with celltype columns (header) and binary presence/absence values

**Output:**
- `metadata.filtered.tsv`: Filtered metadata with celltype-specific BED paths
- `celltype_beds/*.bed`: Per-celltype BED files for peak filtering


### Step 2: Split BAM files by celltype and donor

```bash
scar split
-i metadata.filtered.tsv
-o output_dir
-c celltypes_LEV1
-d
```

**Output structure:**

```txt
output_dir/
├── B/
│ ├── donor1.bam
│ ├── donor1.bam.bai
│ └── donor2.bam
├── CD4_T/
│ ├── donor1.bam
│ └── donor3.bam
└── NK/
└── donor2.bam
```

## Input File Formats

### Metadata TSV

Tab-separated file with header:

```tsv
barcode pool_id donor_id celltypes_LEV1 celltypes_LEV2 bamfile
AAACGAAAGACCACGA-1 pool1 donor_A B B_naive /path/to/file.bam
AAACGAAAGAGACGGA-1 pool1 donor_A CD4_T T_CD4_naive /path/to/file.bam
```

### Peak Matrix TSV

Header row with celltype names, followed by peak coordinates and binary values:

```tsv
B_naive CD4_T NK Mono
chr1:804782-805075 1 0 0 0
chr1:816959-817612 1 1 1 1
chr1:819836-820526 0 1 0 0
```

## License

MIT License - see LICENSE file for details

## Contact

- **Author**: Davide Bolognini
- **Issues**: https://github.com/davidebolo1993/scar/issues

## Acknowledgments

Built with [htslib](https://github.com/samtools/htslib) for BAM file processing.


