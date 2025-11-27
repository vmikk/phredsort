// Main command (`phredsort sort`) for estimating quality scores and sorting sequences
//  Separate functions for stdin-based and file-based modes

package main

import (
	"fmt"
	"io"
	"os"
	"sort"
	"sync"

	"github.com/klauspost/compress/zstd"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// Buffer pools for reducing GC pressure during I/O operations
var (
	// Pool for small buffers (sequence/quality data compression)
	smallBufferPool = sync.Pool{
		New: func() interface{} {
			b := make([]byte, 0, 4096)
			return &b
		},
	}
	// Pool for decompression output buffers
	decompBufferPool = sync.Pool{
		New: func() interface{} {
			b := make([]byte, 0, 8192)
			return &b
		},
	}
)

// getSmallBuffer retrieves a buffer from the small buffer pool
func getSmallBuffer() *[]byte {
	return smallBufferPool.Get().(*[]byte)
}

// putSmallBuffer returns a buffer to the small buffer pool
func putSmallBuffer(b *[]byte) {
	*b = (*b)[:0]
	smallBufferPool.Put(b)
}

// getDecompBuffer retrieves a buffer from the decompression buffer pool
func getDecompBuffer() *[]byte {
	return decompBufferPool.Get().(*[]byte)
}

// putDecompBuffer returns a buffer to the decompression buffer pool
func putDecompBuffer(b *[]byte) {
	*b = (*b)[:0]
	decompBufferPool.Put(b)
}

// CompactStorage provides arena-style storage for compressed FASTQ records.
// Instead of allocating separate []byte slices for each compressed record,
// all compressed data is stored in a single contiguous buffer with offsets.
// This reduces memory overhead from ~24 bytes per slice header to ~8 bytes per record.
type CompactStorage struct {
	data    []byte   // Single contiguous buffer for all compressed data
	offsets []uint32 // Start offset for each record
	lengths []uint32 // Length of each compressed record
}

// NewCompactStorage creates a new CompactStorage with pre-allocated capacity
func NewCompactStorage(estimatedRecords int, estimatedDataSize int) *CompactStorage {
	return &CompactStorage{
		data:    make([]byte, 0, estimatedDataSize),
		offsets: make([]uint32, 0, estimatedRecords),
		lengths: make([]uint32, 0, estimatedRecords),
	}
}

// Append adds compressed data to the storage and returns the index
func (s *CompactStorage) Append(compressed []byte) int {
	idx := len(s.offsets)
	s.offsets = append(s.offsets, uint32(len(s.data)))
	s.lengths = append(s.lengths, uint32(len(compressed)))
	s.data = append(s.data, compressed...)
	return idx
}

// Get retrieves compressed data at the given index
func (s *CompactStorage) Get(idx int) []byte {
	start := s.offsets[idx]
	return s.data[start : start+s.lengths[idx]]
}

// Len returns the number of stored records
func (s *CompactStorage) Len() int {
	return len(s.offsets)
}


// runDefaultCommand is the main entry point for the default sort command.
// It handles flag validation, metric parsing, and delegates to sortRecords
// for the actual sorting work.
//
// Uses a single-pass approach that reads all records into memory, calculates
// quality scores, sorts, and writes output. This unified approach works for
// both stdin ("-") and file inputs, since compressed FASTQ files don't support
// random access anyway.
func runDefaultCommand(cmd *cobra.Command, args []string) {
	// Check version flag
	if version {
		fmt.Printf("phredsort %s\n", VERSION)
		exitFunc(0)
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

	// Process input (unified approach for both stdin and file)
	sortRecords(inFile, outFile, ascending, qualityMetric, compLevel, parsedHeaderMetrics, minPhred, minQualFilter, maxQualFilter)
}

// sortRecords reads FASTQ records from input, calculates quality metrics, sorts them,
// and writes the sorted output. This unified function works for both file and stdin input.
//
// Memory optimization: Uses index-based sorting with compact storage to minimize
// memory overhead. Instead of storing names multiple times (in map keys and sort structs),
// names are stored once in a slice and referenced by integer indices.
//
// When compLevel > 0, records are compressed using ZSTD to reduce memory usage
// When compLevel == 0, records are stored uncompressed (faster but uses more memory)
//
// Parameters:
//   - inFile: Input file path (use "-" for stdin)
//   - outFile: Output file path (use "-" for stdout)
//   - ascending: If true, sort in ascending order; if false, sort in descending order
//   - metric: Quality metric to use for sorting
//   - compLevel: Compression level (0-22, 0 = disabled)
//   - headerMetrics: Optional metrics to append to headers
//   - minPhred: Minimum Phred threshold for lqcount/lqpercent calculations
//   - minQualFilter: Minimum quality threshold for filtering
//   - maxQualFilter: Maximum quality threshold for filtering
func sortRecords(inFile, outFile string, ascending bool, metric QualityMetric, compLevel int, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
	reader, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating reader: %v\n"), err)
		exitFunc(1)
	}
	defer reader.Close()

	// Create output file handle at the beginning
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		exitFunc(1)
	}
	defer outfh.Close()

	if compLevel > 0 {
		sortCompressed(reader, outfh, ascending, metric, compLevel, headerMetrics, minPhred, minQualFilter, maxQualFilter)
	} else {
		sortUncompressed(reader, outfh, ascending, metric, headerMetrics, minPhred, minQualFilter, maxQualFilter)
	}
}

// sortCompressed handles sorting with ZSTD compression enabled.
// Uses compact arena storage to minimize per-record memory overhead.
func sortCompressed(reader *fastx.Reader, outfh *xopen.Writer, ascending bool, metric QualityMetric, compLevel int, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
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

	// Memory-efficient storage using arena and index-based sorting
	// Estimate initial capacity (will grow as needed)
	storage := NewCompactStorage(10000, 1024*1024) // 10K records, 1MB data
	names := make([]string, 0, 10000)
	qualityScores := make([]QualityIndex, 0, 10000)

	// Get a reusable buffer for compression
	compBuf := getSmallBuffer()
	defer putSmallBuffer(compBuf)

	// Reading records
	var idx int32 = 0
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

		// Compress sequence and quality scores together using pooled buffer
		dataLen := len(record.Seq.Seq) + len(record.Seq.Qual)
		if cap(*compBuf) < dataLen {
			*compBuf = make([]byte, 0, dataLen*2)
		}
		*compBuf = (*compBuf)[:0]
		*compBuf = append(*compBuf, record.Seq.Seq...)
		*compBuf = append(*compBuf, record.Seq.Qual...)

		compressed := encoder.EncodeAll(*compBuf, make([]byte, 0, len(*compBuf)/2))
		storage.Append(compressed)

		names = append(names, name)
		qualityScores = append(qualityScores, QualityIndex{
			Index: idx,
			Value: float32(avgQual),
		})
		idx++
	}

	// Sort records using index-based sorting
	qualityList := NewQualityIndexList(qualityScores, names, ascending, metric)
	sort.Sort(qualityList)

	// Get a reusable buffer for decompression
	decompBuf := getDecompBuffer()
	defer putDecompBuffer(decompBuf)

	// Writing records in sorted order
	for _, qi := range qualityList.Items() {
		compData := storage.Get(int(qi.Index))

		// Decompress using pooled buffer
		decompressed, err := decoder.DecodeAll(compData, (*decompBuf)[:0])
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error decompressing record: %v\n"), err)
			exitFunc(1)
		}

		seqLen := len(decompressed) / 2
		record := &fastx.Record{
			Name: []byte(names[qi.Index]),
			Seq: &seq.Seq{
				Seq:  decompressed[:seqLen],
				Qual: decompressed[seqLen:],
			},
		}
		writeRecord(outfh, record, float64(qi.Value), headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
	}
}

// sortUncompressed handles sorting without compression.
// Uses index-based sorting with a slice instead of a map for record storage.
func sortUncompressed(reader *fastx.Reader, outfh *xopen.Writer, ascending bool, metric QualityMetric, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
	// Use slices instead of maps for more efficient memory layout
	records := make([]*fastx.Record, 0, 10000)
	names := make([]string, 0, 10000)
	qualityScores := make([]QualityIndex, 0, 10000)

	// Read all records
	var idx int32 = 0
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

		// Clone and store in slice (indexed access)
		records = append(records, record.Clone())
		names = append(names, name)
		qualityScores = append(qualityScores, QualityIndex{
			Index: idx,
			Value: float32(avgQual),
		})
		idx++
	}

	// Sort records using index-based sorting
	qualityList := NewQualityIndexList(qualityScores, names, ascending, metric)
	sort.Sort(qualityList)

	// Output in sorted order using indices
	for _, qi := range qualityList.Items() {
		record := records[qi.Index]
		writeRecord(outfh, record, float64(qi.Value), headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
	}
}

