// Performance tests 

// ZSTD benchmark
process zstd_benchmark {
    
    tag "${input.simpleName}___${compression_level}"

    input:
        tuple path(input), val(compression_level)
        // path input              // gz-compressed FASTQ file
        // val compression_level   // ZSTD compression level

    output:
        path "done"

    script:
    """
    echo "stdin-mode with ZSTD compression level ${compression_level}"
    zcat ${input} \
      | phredsort --in - --out - --compress ${compression_level} > /dev/null
    touch done
    """
}

// Snappy benchmark
process snappy_benchmark {
    
    tag "${input.simpleName}___snappy"

    input:
        path input

    output:
        path "done"

    script:
    """
    echo "stdin-mode with snappy"
    # NB. notice the older flag format (`-in` and `-out`)
    zcat ${input} \
      | phredsort_snappy -in - -out - > /dev/null
    touch done
    """
}

// File-based benchmarks
process file_benchmark_compressed {
    
    tag "${input.baseName}___filecompressed"

    input:
        path input

    output:
        path "done"

    script:
    """
    echo "file-mode (compressed)"
    phredsort --in ${input} --out - > /dev/null
    touch done
    """
}

// File-based benchmarks
process file_benchmark_uncompressed {
    
    tag "${input.baseName}___file_uncompressed"

    input:
        path input

    output:
        path "done"

    script:
    """
    echo "file-mode (uncompressed)"
    phredsort --in ${input} --out - > /dev/null
    touch done
    """
}

workflow {

    // Channels with input files (repeated N times)
    ch_fastqgz = Channel.fromList( [ file(params.fastqgz) ] * params.iterations).flatten()
    ch_fastq   = Channel.fromList( [ file(params.fastq) ] * params.iterations).flatten()

    // ZSTD compression levels
    ch_zstd_compression = Channel.fromList(params.compression_levels)

    // Preview channels
    // ch_fastqgz.view()
    // ch_fastq.view()
    // ch_zstd_compression.view()

    // Cartesian product of compressed FASTQ files and ZSTD compression levels
    ch_fastqgz
      .combine(ch_zstd_compression)
      .set { ch_zstd_combinations }

    // ch_zstd_combinations.view()

    // Run stdin-based ZSTD benchmarks
    ch_zstd_combinations | zstd_benchmark

    // Run stdin-based snappy benchmarks
    ch_fastqgz | snappy_benchmark

    // Run file-based benchmarks
    ch_fastqgz | file_benchmark_compressed
    ch_fastq   | file_benchmark_uncompressed
}


