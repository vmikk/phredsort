// Subcommand (`phredsort headersort`) for sorting sequences using pre-computed quality scores from headers

package main

import (
	"fmt"
	"math"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/maruel/natural"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// Regular expressions for header parsing
var (
	spaceMetricRe = regexp.MustCompile(`\s+(\w+)=(\d+\.?\d*)`)
	semiMetricRe  = regexp.MustCompile(`;(\w+)=(\d+\.?\d*)`)
	sizeRe        = regexp.MustCompile(`(?:\s|;)size=(\d+)`)
)

// HeaderSortIndex is a memory-efficient struct for sorting pre-computed headers.
// Uses index-based approach to minimize per-record memory overhead.
type HeaderSortIndex struct {
	Index   int32   // Position in records slice
	Quality float32 // Parsed quality value from header
	Size    int32   // Parsed size value from header (0 if not present)
	HasSize bool    // Whether size was present in header
}

// HeaderSortIndexList implements sort.Interface for memory-efficient header-based sorting
type HeaderSortIndexList struct {
	items     []HeaderSortIndex
	ids       []string // External reference for tie-breaking (sequence IDs)
	ascending bool
	metric    QualityMetric
}

// NewHeaderSortIndexList creates a new HeaderSortIndexList
func NewHeaderSortIndexList(items []HeaderSortIndex, ids []string, ascending bool, metric QualityMetric) *HeaderSortIndexList {
	return &HeaderSortIndexList{
		items:     items,
		ids:       ids,
		ascending: ascending,
		metric:    metric,
	}
}

func (list *HeaderSortIndexList) Len() int { return len(list.items) }
func (list *HeaderSortIndexList) Swap(i, j int) {
	list.items[i], list.items[j] = list.items[j], list.items[i]
}

func (list *HeaderSortIndexList) Less(i, j int) bool {
	qi, qj := list.items[i].Quality, list.items[j].Quality

	// Primary sort by quality
	if qi != qj {
		var result bool
		if list.metric == MaxEE || list.metric == Meep || list.metric == LQCount || list.metric == LQPercent {
			result = qi < qj
		} else {
			result = qi > qj
		}
		if list.ascending {
			return !result
		}
		return result
	}

	// Secondary sort by size (if both have size)
	if list.items[i].HasSize && list.items[j].HasSize {
		si, sj := list.items[i].Size, list.items[j].Size
		if si != sj {
			if list.ascending {
				return si < sj
			}
			return si > sj
		}
	}

	// Tertiary sort by ID using natural ordering
	return natural.Less(list.ids[list.items[i].Index], list.ids[list.items[j].Index])
}

// Items returns the underlying items slice
func (list *HeaderSortIndexList) Items() []HeaderSortIndex {
	return list.items
}

// PreSortRecord represents a sequence record with pre-computed quality
// Kept for backward compatibility with parsePreSortRecord function used in tests
type PreSortRecord struct {
	ID      string
	Name    string
	Quality float64
	Size    int
	Record  *fastx.Record
	HasSize bool
	HasQual bool
}

// parseHeaderInfo extracts quality metric, size, and ID from a record header.
// Returns the parsed values for use with index-based sorting.
// This is a more efficient version that doesn't allocate a full PreSortRecord.
func parseHeaderInfo(header string, metric QualityMetric) (id string, quality float32, size int32, hasQual bool, hasSize bool) {
	parts := strings.SplitN(header, " ", 2)
	id = parts[0]
	if strings.HasPrefix(id, ">") || strings.HasPrefix(id, "@") {
		id = id[1:]
	}

	// Try to find size annotation
	if sizeMatch := sizeRe.FindStringSubmatch(header); sizeMatch != nil {
		if s, err := strconv.Atoi(sizeMatch[1]); err == nil {
			size = int32(s)
			hasSize = true
		}
	}

	// Look for quality metric in both formats
	metricStr := metric.String()
	var qualityStr string
	found := false

	// Check space-separated format
	if matches := spaceMetricRe.FindAllStringSubmatch(header, -1); matches != nil {
		for _, match := range matches {
			if match[1] == metricStr {
				qualityStr = match[2]
				found = true
				break
			}
		}
	}

	// Check semicolon-separated format if not found
	if !found {
		if matches := semiMetricRe.FindAllStringSubmatch(header, -1); matches != nil {
			for _, match := range matches {
				if match[1] == metricStr {
					qualityStr = match[2]
					found = true
					break
				}
			}
		}
	}

	if found {
		if q, err := strconv.ParseFloat(qualityStr, 32); err == nil {
			quality = float32(q)
			hasQual = true
		}
	}

	return
}

// parsePreSortRecord parses a FASTQ/FASTA record header to extract quality metric
// and size information. Supports both space-separated and semicolon-separated formats
//
// The function looks for metric annotations in the format "metric=value" (e.g., "maxee=2.5")
// and size annotations as "size=value" (e.g., "size=100")
//
// Returns a PreSortRecord with parsed information, or an error if the required
// metric is missing from the header
//
// Note: This function is kept for backward compatibility with tests.
// The new index-based sorting uses parseHeaderInfo instead.
func parsePreSortRecord(record *fastx.Record, metric QualityMetric) (*PreSortRecord, error) {
	header := string(record.Name)
	id, quality, size, hasQual, hasSize := parseHeaderInfo(header, metric)

	return &PreSortRecord{
		ID:      id,
		Name:    header,
		Record:  record,
		Quality: float64(quality),
		Size:    int(size),
		HasQual: hasQual,
		HasSize: hasSize,
	}, nil
}

// HeaderSortCommand creates the `headersort` subcommand which sorts sequences
// using pre-computed quality scores stored in sequence headers
//
// This command is useful when quality metrics have already been calculated and
// stored in headers (e.g., by a previous run of phredsort with --header flag)
// It supports both space-separated (">seq1 maxee=2") and semicolon-separated
// (">seq1;maxee=2") header formats
//
// Secondary sorting is performed by size annotation (if present) and sequence ID
func HeaderSortCommand() *cobra.Command {
	var (
		inFile        string
		outFile       string
		metric        string
		ascending     bool
		minQualFilter float64
		maxQualFilter float64
	)

	cmd := &cobra.Command{
		Use:   "headersort",
		Short: "Sort FASTA/FASTQ sequences using pre-computed quality scores from headers",
		Long: `Sort sequences using pre-computed quality scores stored in sequence headers.
Supports both FASTA and FASTQ formats with space-separated (">seq1 maxee=2") or 
semicolon-separated (">seq1;maxee=2") quality annotations. Secondary sorting is done 
by size annotation (if present, e.g., "size=123") and sequence ID.`,
		RunE: func(cmd *cobra.Command, args []string) error {
			// Validate metric flag
			qualityMetric, err := validateMetric(metric)
			if err != nil {
				return err
			}

			return runPresort(inFile, outFile, qualityMetric, ascending, minQualFilter, maxQualFilter)
		},
	}

	// Define flags
	flags := cmd.Flags()
	flags.StringVarP(&inFile, "in", "i", "", "Input sequence file (required)")
	flags.StringVarP(&outFile, "out", "o", "", "Output sequence file (required)")
	flags.StringVarP(&metric, "metric", "s", "avgphred", "Quality metric to use from headers")
	flags.BoolVarP(&ascending, "ascending", "a", false, "Sort in ascending order")
	flags.Float64VarP(&minQualFilter, "minqual", "m", -math.MaxFloat64, "Minimum quality threshold")
	flags.Float64VarP(&maxQualFilter, "maxqual", "M", math.MaxFloat64, "Maximum quality threshold")

	cmd.MarkFlagRequired("in")
	cmd.MarkFlagRequired("out")

	return cmd
}

// runPresort reads FASTQ/FASTA records, extracts quality metrics from headers,
// filters records based on quality thresholds, sorts them, and writes the
// sorted output. This function requires that quality metrics are already present
// in sequence headers
//
// Memory optimization: Uses index-based sorting with slices instead of storing
// full PreSortRecord structs. Records are stored in a slice and referenced by
// integer indices during sorting.
//
// Parameters:
//   - inFile: Input sequence file path
//   - outFile: Output sequence file path
//   - metric: Quality metric to extract from headers and use for sorting
//   - ascending: If true, sort in ascending order; if false, sort in descending order
//   - minQual: Minimum quality threshold for filtering
//   - maxQual: Maximum quality threshold for filtering
//
// Returns an error if file I/O fails or if a record is missing the required metric
func runPresort(inFile, outFile string, metric QualityMetric, ascending bool, minQual, maxQual float64) error {
	// Create reader with automatic format detection
	reader, err := fastx.NewDefaultReader(inFile)
	if err != nil {
		return fmt.Errorf("error creating reader: %v", err)
	}
	defer reader.Close()

	// Create buffered output writer
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		return fmt.Errorf("error creating output file: %v", err)
	}
	defer outfh.Close()

	// Memory-efficient storage using slices and index-based sorting
	records := make([]*fastx.Record, 0, 10000)
	ids := make([]string, 0, 10000)
	sortIndices := make([]HeaderSortIndex, 0, 10000)

	minQual32 := float32(minQual)
	maxQual32 := float32(maxQual)

	// Use ChunkChan for asynchronous reading with reasonable buffer sizes
	bufferSize := 100 // Number of chunks to buffer
	chunkSize := 1000 // Records per chunk

	var idx int32 = 0
	for chunk := range reader.ChunkChan(bufferSize, chunkSize) {
		if chunk.Err != nil {
			return fmt.Errorf("error reading chunk: %v", chunk.Err)
		}

		for _, record := range chunk.Data {
			header := string(record.Name)
			id, quality, size, hasQual, hasSize := parseHeaderInfo(header, metric)

			if !hasQual {
				return fmt.Errorf("record missing required quality metric (%s): %s", metric, header)
			}

			// Apply quality filters
			if quality >= minQual32 && quality <= maxQual32 {
				// Store record and add to sort indices
				records = append(records, record) // ChunkChan already provides copies
				ids = append(ids, id)
				sortIndices = append(sortIndices, HeaderSortIndex{
					Index:   idx,
					Quality: quality,
					Size:    size,
					HasSize: hasSize,
				})
				idx++
			}
		}
	}

	// Sort using index-based sorting
	sortList := NewHeaderSortIndexList(sortIndices, ids, ascending, metric)
	sort.Sort(sortList)

	// Write sorted records using indices
	for _, si := range sortList.Items() {
		records[si.Index].FormatToWriter(outfh, 0)
	}

	return nil
}
