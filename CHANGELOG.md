# Changelog

## phredsort v0.5.0 - 2024-11-30

Initial release.

`phredsort` is a high-performance tool for sorting FASTQ files based on sequence quality metrics.

### Features

- [new] Multiple quality metrics:
  - `avgphred`: average Phred quality score  
  - `maxee`: maximum expected error, as an absolute number  
  - `meep`: maximum expected error percentage per sequence length  
- [new] Flexible input/output:
  - Supports file-based and stdin/stdout streaming workflows  
  - Handles compressed FASTQ files (`gzip`, `bzip2`, `xz`, `zstd`) and uncompressed FASTQ files  
  - Optionally uses in-memory compression for stdin-based processing to reduce memory footprint  
