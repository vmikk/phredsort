# phredsort

`phredsort` is a command-line tool for sorting sequences in a FASTQ file by their quality scores.

## Usage

Basic usage:
```bash
# Read from `input.fastq.gz` and write to `output.fastq.gz`
phredsort -i input.fastq.gz -o output.fastq.gz

# Read from stdin and write to stdout
zcat input.fastq.gz | phredsort --in - --out - | less -S
```

![phredsort help message](assets/phredsort.png)



