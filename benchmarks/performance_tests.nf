// Performance tests 

// Usage:
// ## `phredsort` (and optionally `phredsort_snappy`) binaries should be in the `bin` directory
// ## currently, performs one test at a time (change `-qs` to run tests in parallel)
//
// nextflow run performance_tests.nf \
//   -resume -with-trace -qs 1 \
//   --fastqgz 'data/*.fastq.gz' \
//   --fastq 'data/*.fastq' \
//   --iterations 3

// Define parameters
params.fastqgz = 'data/*.fastq.gz'     // Set of compressed FASTQ files of different sizes
params.fastq = 'data/*.fastq'          // Set of uncompressed FASTQ files
params.compression_levels = 0..4       // ZSTD compression levels from 0 to 20
// params.compression_levels = [0, 1, 3, 5]
params.iterations = 3                  // Number of times to repeat each execution


// ZSTD benchmark
process phredsorting {
    
    tag "${input.simpleName}___${compression_level}___${version}"

    input:
        tuple path(phredsort_bin), val(version), path(input), val(compression_level)
        // path phredsort_bin      // path to the phredsort binary for this run
        // val  version            // version label derived from the binary name
        // path input              // FASTQ file
        // val  compression_level  // ZSTD compression level

    output:
        path "done"

    script:
    """
    echo -e "Input: ${input}"
    echo -e "ZSTD compression level: ${compression_level}"
    echo -e "Phredsort version: ${version}"

    ${phredsort_bin} \
      --in ${input} \
      --out - \
      --compress ${compression_level} \
    > /dev/null
    
    touch done
    """
}

/*
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
*/

workflow {

    // Channels with input files (repeated N times)
    ch_fastqgz = Channel.fromList( [ file(params.fastqgz) ] * params.iterations).flatten()
    ch_fastq   = Channel.fromList( [ file(params.fastq) ] * params.iterations).flatten()

    // ZSTD compression levels
    ch_zstd_compression = Channel.fromList(params.compression_levels)

    // phredsort binaries (different versions to benchmark)
    // Expected names in the `bin` directory: `phredsort`, `phredsort_<version>`, ...
    ch_phredsort_bins = Channel
        .fromPath("bin/phredsort*")
        .map { bin_path ->
            def name = bin_path.getName()
            def version = (name == 'phredsort') ? 'current' : name.replaceFirst(/^phredsort_?/, '')
            tuple(bin_path, version)
        }

    // Preview channels
    // ch_fastqgz.view()
    // ch_fastq.view()
    // ch_zstd_compression.view()
    // ch_phredsort_bins.view()


    // Cartesian product of:
    //   - phredsort binary (version)
    //   - compressed FASTQ files
    //   - ZSTD compression levels
    ch_task_combinations = ch_fastqgz
      .combine(ch_phredsort_bins)
      .map { combo ->
          /*
           * `combine` flattens tuples, so each emitted item here is:
           *   [ input, phredsort_bin, version ]
           */
          def (input, phredsort_bin, version) = combo
          tuple(phredsort_bin, version, input)
      }
      .combine(ch_zstd_compression)
      .map { combo ->
          /*
           * After combining with `ch_zstd_compression`, each item is:
           *   [ phredsort_bin, version, input, compression_level ]
           */
          def (phredsort_bin, version, input, compression_level) = combo
          tuple(phredsort_bin, version, input, compression_level)
      }

    // ch_task_combinations.view()


    // Run benchmarks for all parameter combinations
    ch_task_combinations | phredsorting

    /*
    // Run stdin-based snappy benchmarks
    // ch_fastqgz | snappy_benchmark

    */
}


