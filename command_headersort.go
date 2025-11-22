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

// PreSortRecord represents a sequence record with pre-computed quality
type PreSortRecord struct {
	ID      string
	Name    string
	Quality float64
	Size    int
	Record  *fastx.Record
	HasSize bool
	HasQual bool
}

// Regular expressions for header parsing
var (
	spaceMetricRe = regexp.MustCompile(`\s+(\w+)=(\d+\.?\d*)`)
	semiMetricRe  = regexp.MustCompile(`;(\w+)=(\d+\.?\d*)`)
	sizeRe        = regexp.MustCompile(`(?:\s|;)size=(\d+)`)
)

// PreSortRecordList for sorting
type PreSortRecordList struct {
	items     []*PreSortRecord
	ascending bool
	metric    QualityMetric
}

func NewPreSortRecordList(items []*PreSortRecord, ascending bool, metric QualityMetric) *PreSortRecordList {
	return &PreSortRecordList{
		items:     items,
		ascending: ascending,
		metric:    metric,
	}
}

func (list *PreSortRecordList) Len() int { return len(list.items) }
func (list *PreSortRecordList) Swap(i, j int) {
	list.items[i], list.items[j] = list.items[j], list.items[i]
}
func (list *PreSortRecordList) Less(i, j int) bool {
	// Primary sort by quality, using the same semantics as QualityFloatList:
	// - For AvgPhred (higher is better), default (ascending=false) orders from
	//   highest to lowest.
	// - For MaxEE, Meep, LQCount, LQPercent (lower is better), default
	//   (ascending=false) orders from lowest to highest.
	// The "ascending" flag flips this default ordering.
	if list.items[i].Quality != list.items[j].Quality {
		var result bool
		if list.metric == MaxEE || list.metric == Meep || list.metric == LQCount || list.metric == LQPercent {
			// For these metrics, lower values indicate better quality.
			result = list.items[i].Quality < list.items[j].Quality
		} else {
			// For AvgPhred and other metrics where higher is better.
			result = list.items[i].Quality > list.items[j].Quality
		}

		if list.ascending {
			return !result
		}
		return result
	}

	// Secondary sort by size (if both have size)
	if list.items[i].HasSize && list.items[j].HasSize {
		if list.items[i].Size != list.items[j].Size {
			if list.ascending {
				return list.items[i].Size < list.items[j].Size
			}
			return list.items[i].Size > list.items[j].Size
		}
	}

	// Tertiary sort by ID
	return natural.Less(list.items[i].ID, list.items[j].ID)
}

// parsePreSortRecord parses a FASTQ/FASTA record header to extract quality metric
// and size information. Supports both space-separated and semicolon-separated formats
//
// The function looks for metric annotations in the format "metric=value" (e.g., "maxee=2.5")
// and size annotations as "size=value" (e.g., "size=100")
//
// Returns a PreSortRecord with parsed information, or an error if the required
// metric is missing from the header
func parsePreSortRecord(record *fastx.Record, metric QualityMetric) (*PreSortRecord, error) {
	header := string(record.Name)
	parts := strings.SplitN(header, " ", 2) // Split at first space
	id := parts[0]
	if strings.HasPrefix(id, ">") || strings.HasPrefix(id, "@") {
		id = id[1:]
	}

	result := &PreSortRecord{
		ID:      id,
		Name:    header,
		Record:  record,
		HasQual: false,
		HasSize: false,
	}

	// Try to find size annotation
	if sizeMatch := sizeRe.FindStringSubmatch(header); sizeMatch != nil {
		if size, err := strconv.Atoi(sizeMatch[1]); err == nil {
			result.Size = size
			result.HasSize = true
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
		if quality, err := strconv.ParseFloat(qualityStr, 64); err == nil {
			result.Quality = quality
			result.HasQual = true
		}
	}

	return result, nil
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

	// Read and parse all records
	var records []*PreSortRecord

	// Use ChunkChan for asynchronous reading with reasonable buffer sizes
	bufferSize := 100 // Number of chunks to buffer
	chunkSize := 1000 // Records per chunk

	for chunk := range reader.ChunkChan(bufferSize, chunkSize) {
		if chunk.Err != nil {
			return fmt.Errorf("error reading chunk: %v", chunk.Err)
		}

		for _, record := range chunk.Data {
			// No need to clone since ChunkChan already provides copies
			presortRecord, err := parsePreSortRecord(record, metric)
			if err != nil {
				return fmt.Errorf("error parsing record: %v", err)
			}

			if !presortRecord.HasQual {
				return fmt.Errorf("record missing required quality metric (%s): %s", metric, presortRecord.Name)
			}

			// Apply quality filters
			if presortRecord.Quality >= minQual && presortRecord.Quality <= maxQual {
				records = append(records, presortRecord)
			}
		}
	}

	// Sort records
	recordList := NewPreSortRecordList(records, ascending, metric)
	sort.Sort(recordList)

	// Write sorted records using buffered writer
	for _, record := range recordList.items {
		record.Record.FormatToWriter(outfh, 0)
	}

	return nil
}
