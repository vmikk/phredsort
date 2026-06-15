// Main command (`phredsort sort`) for estimating quality scores and sorting sequences
//  Separate functions for stdin-based and file-based modes

package main

import (
	"fmt"
	"io"
	"math/bits"
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

const (
	defaultChunkSize          = 16 * 1024 * 1024 // 16 MiB should be a balance between copies and locality
	warningThresholdBytes     = 4 * 1024 * 1024 * 1024
	computedQualityFastqError = "computed quality sorting requires FASTQ input with quality scores; use headersort for FASTA with precomputed metrics"
)

type storageChunk struct {
	data  []byte
	used  int
	start uint64
}

// ChunkedStorage keeps compressed FASTQ records in fixed-size chunks
// Record offsets are tracked in uint64 to safely support multi-gigabyte datasets without re-copying old data
type ChunkedStorage struct {
	chunks       []*storageChunk
	offsets      []uint64
	lengths      []uint64
	chunkSize    int
	currentChunk *storageChunk
	totalSize    uint64
	warned       bool
}

// NewChunkedStorage reserved metadata for the expected number of records.
func NewChunkedStorage(estimatedRecords int, chunkSize int) *ChunkedStorage {
	if chunkSize <= 0 {
		chunkSize = defaultChunkSize
	}
	return &ChunkedStorage{
		chunks:    make([]*storageChunk, 0, 4),
		offsets:   make([]uint64, 0, estimatedRecords),
		lengths:   make([]uint64, 0, estimatedRecords),
		chunkSize: chunkSize,
	}
}

// Append stores compressed data and returns its index
func (s *ChunkedStorage) Append(compressed []byte) int {
	length := len(compressed)
	reserveLength := length
	if reserveLength == 0 {
		reserveLength = 1
	}
	s.ensureCapacity(reserveLength)
	chunk := s.currentChunk
	copy(chunk.data[chunk.used:chunk.used+length], compressed)

	idx := len(s.offsets)
	start := s.totalSize

	s.offsets = append(s.offsets, start)
	s.lengths = append(s.lengths, uint64(length))
	chunk.used += reserveLength
	s.totalSize += uint64(reserveLength)
	s.maybeWarn()

	return idx
}

// Get returns compressed data view for index
func (s *ChunkedStorage) Get(idx int) []byte {
	offset := s.offsets[idx]
	length := s.lengths[idx]
	chunkIdx := s.findChunk(offset)
	chunk := s.chunks[chunkIdx]
	relative := offset - chunk.start
	return chunk.data[int(relative):int(relative+length)]
}

func (s *ChunkedStorage) Len() int {
	return len(s.offsets)
}

func (s *ChunkedStorage) ensureCapacity(length int) {
	if length <= 0 {
		length = 1
	}
	if s.currentChunk != nil && len(s.currentChunk.data)-s.currentChunk.used >= length {
		return
	}

	newSize := s.chunkSize
	if length > newSize {
		newSize = nextPowerOfTwo(length)
	}

	if s.totalSize+uint64(newSize) < s.totalSize {
		fmt.Fprint(os.Stderr, red("Error: compressed data exceeds supported capacity\n"))
		exitFunc(1)
	}

	chunk := &storageChunk{
		data:  make([]byte, newSize),
		used:  0,
		start: s.totalSize,
	}
	s.chunks = append(s.chunks, chunk)
	s.currentChunk = chunk
}

func (s *ChunkedStorage) findChunk(offset uint64) int {
	i := sort.Search(len(s.chunks), func(i int) bool {
		chunk := s.chunks[i]
		return offset < chunk.start+uint64(chunk.used)
	})
	if i == len(s.chunks) {
		fmt.Fprint(os.Stderr, red("Error: chunk lookup failed\n"))
		exitFunc(1)
	}
	return i
}

func (s *ChunkedStorage) maybeWarn() {
	if s.warned {
		return
	}
	if s.totalSize >= warningThresholdBytes {
		fmt.Fprint(os.Stderr, yellow("Warning: compressed data exceeded 4 GiB; memory usage may grow\n"))
		s.warned = true
	}
}

func nextPowerOfTwo(v int) int {
	if v <= 1 {
		return 1
	}
	return 1 << uint(bits.Len(uint(v-1)))
}

func zstdEncodeCapacity(inputLen int) int {
	if inputLen <= 0 {
		return 64
	}
	return nextPowerOfTwo(inputLen + inputLen/8 + 64)
}

func exitIfNotFastq(reader *fastx.Reader, closeReader *bool) {
	if reader.IsFastq {
		return
	}
	*closeReader = false
	fmt.Fprintln(os.Stderr, red("Error: "+computedQualityFastqError))
	exitFunc(1)
}

// runDefaultCommand is the main entry point for the default sort command
// It handles flag validation, metric parsing, and delegates to sortRecords
// for the actual sorting work
//
// Uses a single-pass approach that reads all records into memory,
// calculates quality scores, sorts, and writes output.
// This unified approach works for both stdin ("-") and file inputs,
// since compressed FASTQ files don't support random access anyway
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
// and writes the sorted output. This unified function works for both file and stdin input
//
// Memory optimization: Uses index-based sorting with compact storage to minimize
// memory overhead. Instead of storing names multiple times (in map keys and sort structs),
// names are stored once in a slice and referenced by integer indices
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
	closeReader := true
	defer func() {
		if closeReader {
			reader.Close()
		}
	}()

	// Create output file handle at the beginning
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		exitFunc(1)
	}
	defer outfh.Close()

	if compLevel > 0 {
		sortCompressed(reader, outfh, ascending, metric, compLevel, headerMetrics, minPhred, minQualFilter, maxQualFilter, &closeReader)
	} else {
		sortUncompressed(reader, outfh, ascending, metric, headerMetrics, minPhred, minQualFilter, maxQualFilter, &closeReader)
	}
}

// sortCompressed handles sorting with ZSTD compression enabled
// Uses chunked storage to avoid monolithic compressed-buffer reallocations
func sortCompressed(reader *fastx.Reader, outfh *xopen.Writer, ascending bool, metric QualityMetric, compLevel int, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64, closeReader *bool) {
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

	// Memory-efficient storage using chunks and index-based sorting
	// Estimate initial metadata capacity (will grow as needed)
	storage := NewChunkedStorage(10000, 0)
	names := make([]string, 0, 10000)
	qualityScores := make([]QualityIndex, 0, 10000)

	// Get a reusable buffer for compression
	compBuf := getSmallBuffer()
	defer putSmallBuffer(compBuf)
	encBuf := getSmallBuffer()
	defer putSmallBuffer(encBuf)

	// Reading records
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			exitFunc(1)
		}
		exitIfNotFastq(reader, closeReader)

		name := string(record.Name)
		avgQual := calculateQuality(record, metric, minPhred)
		if avgQual < minQualFilter || avgQual > maxQualFilter {
			continue
		}

		// Compress sequence and quality scores together using pooled buffer
		dataLen := len(record.Seq.Seq) + len(record.Seq.Qual)
		if cap(*compBuf) < dataLen {
			*compBuf = make([]byte, 0, nextPowerOfTwo(dataLen))
		}
		buf := (*compBuf)[:0]
		buf = append(buf, record.Seq.Seq...)
		buf = append(buf, record.Seq.Qual...)
		*compBuf = buf

		encCap := zstdEncodeCapacity(len(buf))
		if cap(*encBuf) < encCap {
			*encBuf = make([]byte, 0, encCap)
		}
		compressed := encoder.EncodeAll(buf, (*encBuf)[:0])
		storage.Append(compressed)
		*encBuf = compressed

		names = append(names, name)
		qualityScores = append(qualityScores, QualityIndex{
			Index: len(names) - 1,
			Value: avgQual,
		})
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
		*decompBuf = decompressed

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

// sortUncompressed handles sorting without compression
// Uses index-based sorting with a slice instead of a map for record storage
func sortUncompressed(reader *fastx.Reader, outfh *xopen.Writer, ascending bool, metric QualityMetric, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64, closeReader *bool) {
	// Use slices instead of maps for more efficient memory layout
	records := make([]*fastx.Record, 0, 10000)
	names := make([]string, 0, 10000)
	qualityScores := make([]QualityIndex, 0, 10000)

	// Read all records
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			exitFunc(1)
		}
		exitIfNotFastq(reader, closeReader)

		name := string(record.Name)
		avgQual := calculateQuality(record, metric, minPhred)
		if avgQual < minQualFilter || avgQual > maxQualFilter {
			continue
		}

		// Clone and store in slice (indexed access)
		records = append(records, record.Clone())
		names = append(names, name)
		qualityScores = append(qualityScores, QualityIndex{
			Index: len(records) - 1,
			Value: avgQual,
		})
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
