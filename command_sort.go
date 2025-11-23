// Main command (`phredsort sort`) for estimating quality scores and sorting sequences
//  Separate functions for stdin-based and file-based modes

package main

import (
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/klauspost/compress/zstd"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)


// runDefaultCommand is the main entry point for the default sort command.
// It handles flag validation, metric parsing, and routes to either stdin-based
// or file-based sorting depending on the input file parameter
//
// When inFile is "-", it uses sortStdin (which loads all records into memory)
// Otherwise, it uses sortFile (which uses a two-pass approach for better memory efficiency)
func runDefaultCommand(cmd *cobra.Command, args []string) {
	// Check version flag
	if version {
		fmt.Printf("phredsort %s\n", VERSION)
		exitFunc(0)
	}

	// Check required flags
	if inFile == "" || outFile == "" {
		fmt.Fprintln(os.Stderr, red("Error: input and output files are required"))
		fmt.Fprintln(os.Stderr, red("Try 'phredsort --help' for more information"))
		exitFunc(1)
	}

	// Validate metric flag
	qualityMetric, err := validateMetric(metric)
	if err != nil {
		fmt.Fprintln(os.Stderr, red("Error: "+err.Error()))
		exitFunc(1)
	}

	// Validate compression level
	if compLevel < 0 || compLevel > 22 {
		fmt.Fprintln(os.Stderr, red("Error: compression level must be between 0 and 22"))
		exitFunc(1)
	}

	// Parse header metrics
	parsedHeaderMetrics, err := parseHeaderMetrics(headerMetrics)
	if err != nil {
		fmt.Fprintln(os.Stderr, red(err.Error()))
		exitFunc(1)
	}

	// Process the files
	if inFile == "-" {
		sortStdin(outFile, ascending, qualityMetric, compLevel, parsedHeaderMetrics, minPhred, minQualFilter, maxQualFilter)
	} else {
		sortFile(inFile, outFile, ascending, qualityMetric, compLevel, parsedHeaderMetrics, minPhred, minQualFilter, maxQualFilter)
	}
}

type CompressedFastqRecord struct {
	Name    []byte
	Data    []byte
	AvgQual float64
}


// sortStdin reads FASTQ records from stdin, calculates quality metrics, sorts them,
// and writes the sorted output to outFile. This function loads all records into memory,
// making it suitable for smaller datasets
//
// When compLevel > 0, records are compressed using ZSTD to reduce memory usage
// When compLevel == 0, records are stored uncompressed (faster but uses more memory)
//
// Parameters:
//   - outFile: Output file path
//   - ascending: If true, sort in ascending order; if false, sort in descending order
//   - metric: Quality metric to use for sorting
//   - compLevel: Compression level (0-22, 0 = disabled)
//   - headerMetrics: Optional metrics to append to headers
//   - minPhred: Minimum Phred threshold for lqcount/lqpercent calculations
//   - minQualFilter: Minimum quality threshold for filtering
//   - maxQualFilter: Maximum quality threshold for filtering
func sortStdin(outFile string, ascending bool, metric QualityMetric, compLevel int, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
	reader, err := fastx.NewReader(seq.DNAredundant, "-", fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating reader: %v\n"), err)
		exitFunc(1)
	}
	defer reader.Close()

	var name2avgQual []QualityFloat

	// Create output file handle at the beginning
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		exitFunc(1)
	}
	defer outfh.Close()

	if compLevel > 0 {
		encoder, err := zstd.NewWriter(nil, zstd.WithEncoderLevel(zstd.EncoderLevelFromZstd(compLevel)))
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD encoder: %v\n"), err)
			exitFunc(1)
		}
		defer encoder.Close()

		decoder, err := zstd.NewReader(nil)
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD decoder: %v\n"), err)
			exitFunc(1)
		}
		defer decoder.Close()

		sequences := make(map[string]*CompressedFastqRecord)

		// Reading records
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break // Exit loop when we reach end of file
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
				exitFunc(1)
			}

			name := string(record.Name)
			avgQual := calculateQuality(record, metric, minPhred)

			nameCopy := make([]byte, len(record.Name))
			copy(nameCopy, record.Name)

			// Compress sequence and quality scores together
			data := append(record.Seq.Seq, record.Seq.Qual...)
			compressed := encoder.EncodeAll(data, make([]byte, 0, len(data)))

			sequences[name] = &CompressedFastqRecord{
				Name:    nameCopy,
				Data:    compressed,
				AvgQual: avgQual,
			}
			name2avgQual = append(name2avgQual, QualityFloat{
				Name:   name,
				Value:  avgQual,
				Metric: metric,
			})
		}

		// Sort records
		qualityList := NewQualityFloatList(name2avgQual, ascending)
		sort.Sort(qualityList)
		name2avgQual = qualityList.Items()

		// Writing records
		for _, kv := range name2avgQual {
			compRecord := sequences[kv.Name]
			decompressed, err := decoder.DecodeAll(compRecord.Data, nil)
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error decompressing record: %v\n"), err)
				exitFunc(1)
			}

			seqLen := len(decompressed) / 2
			record := &fastx.Record{
				Name: compRecord.Name,
				Seq: &seq.Seq{
					Seq:  decompressed[:seqLen],
					Qual: decompressed[seqLen:],
				},
			}
			writeRecord(outfh, record, kv.Value, headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
		}
	} else {
		sequences := make(map[string]*fastx.Record)

		// Read all records
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break // Exit loop when we reach end of file
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
				exitFunc(1)
			}

			name := string(record.Name)
			avgQual := calculateQuality(record, metric, minPhred)

			// Important: Clone the record to avoid reference issues
			sequences[name] = record.Clone()
			name2avgQual = append(name2avgQual, QualityFloat{
				Name:   name,
				Value:  avgQual,
				Metric: metric,
			})
		}

		// Sort records
		qualityList := NewQualityFloatList(name2avgQual, ascending)
		sort.Sort(qualityList)
		name2avgQual = qualityList.Items()

		// Write sorted records
		outfh, err := xopen.Wopen(outFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
			exitFunc(1)
		}
		defer outfh.Close()

		// Output in sorted order
		for _, kv := range name2avgQual {
			record := sequences[kv.Name]
			writeRecord(outfh, record, kv.Value, headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
		}
	}
}

// sortFile performs a two-pass sort of FASTQ records from a file. The first pass
// calculates quality scores and determines sort order. The second pass reads records
// in sorted order and writes them to the output file. This approach is more memory-efficient
// than sortStdin for large files, as it doesn't require loading all records into memory simultaneously
//
// Parameters:
//   - inFile: Input FASTQ file path
//   - outFile: Output FASTQ file path
//   - ascending: If true, sort in ascending order; if false, sort in descending order
//   - metric: Quality metric to use for sorting
//   - compLevel: Compression level (0-22, 0 = disabled)
//   - headerMetrics: Optional metrics to append to headers
//   - minPhred: Minimum Phred threshold for lqcount/lqpercent calculations
//   - minQualFilter: Minimum quality threshold for filtering
//   - maxQualFilter: Maximum quality threshold for filtering
func sortFile(inFile, outFile string, ascending bool, metric QualityMetric, compLevel int, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
	reader, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating reader: %v\n"), err)
		exitFunc(1)
	}
	defer reader.Close()

	// First pass: collect quality scores and positions
	var qualityScores []QualityFloat
	name2offset := make(map[string]int64)
	var currentOffset int64 = 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			exitFunc(1)
		}

		name := string(record.Name)
		avgQual := calculateQuality(record, metric, minPhred)

		qualityScores = append(qualityScores, QualityFloat{
			Name:   name,
			Value:  avgQual,
			Metric: metric,
		})
		name2offset[name] = currentOffset
		currentOffset++
	}

	// Sort by average quality
	qualityList := NewQualityFloatList(qualityScores, ascending)
	sort.Sort(qualityList)
	qualityScores = qualityList.Items()

	// Second pass: read all records into memory
	reader2, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating second reader: %v\n"), err)
		exitFunc(1)
	}
	defer reader2.Close()

	// Open output file early to be ready for writing
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		exitFunc(1)
	}
	defer outfh.Close()

	if compLevel > 0 {
		// Use ZSTD compression for memory efficiency
		encoder, err := zstd.NewWriter(nil, zstd.WithEncoderLevel(zstd.EncoderLevelFromZstd(compLevel)))
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD encoder: %v\n"), err)
			exitFunc(1)
		}
		defer encoder.Close()

		decoder, err := zstd.NewReader(nil)
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD decoder: %v\n"), err)
			exitFunc(1)
		}
		defer decoder.Close()

		compressedRecords := make(map[int64][]byte)
		var offset int64 = 0
		for {
			record, err := reader2.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
				exitFunc(1)
			}

			// Compress sequence and quality scores together
			data := append(record.Seq.Seq, record.Seq.Qual...)
			compressed := encoder.EncodeAll(data, make([]byte, 0, len(data)))
			compressedRecords[offset] = compressed
			offset++
		}

		// Output in sorted order
		for _, qf := range qualityScores {
			offset := name2offset[qf.Name]
			compData, ok := compressedRecords[offset]
			if !ok {
				fmt.Fprintf(os.Stderr, red("Error: could not find record for %s\n"), qf.Name)
				exitFunc(1)
			}

			decompressed, err := decoder.DecodeAll(compData, nil)
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error decompressing record: %v\n"), err)
				exitFunc(1)
			}

			seqLen := len(decompressed) / 2
			record := &fastx.Record{
				Name: []byte(qf.Name), // Reconstruct name from QualityFloat
				Seq: &seq.Seq{
					Seq:  decompressed[:seqLen],
					Qual: decompressed[seqLen:],
				},
			}
			writeRecord(outfh, record, qf.Value, headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
		}

	} else {
		// Uncompressed mode
		records := make(map[int64]*fastx.Record)
		var offset int64 = 0
		for {
			record, err := reader2.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
				exitFunc(1)
			}
			records[offset] = record.Clone()
			offset++
		}

		// Output in sorted order
		for _, qf := range qualityScores {
			offset := name2offset[qf.Name]
			record, ok := records[offset]
			if !ok {
				fmt.Fprintf(os.Stderr, red("Error: could not find record for %s\n"), qf.Name)
				exitFunc(1)
			}
			writeRecord(outfh, record, qf.Value, headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
		}
	}
}
