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



