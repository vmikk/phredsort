// Subcommand (`phredsort nosort`) for sequence quality estimation and optional filtering.
// No sorting is performed; sequences are emitted in the same order as input.

package main

import (
	"fmt"
	"io"
	"math"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// NoSortCommand creates the `nosort` subcommand which estimates quality metrics
// and optionally filters/annotates sequences without changing their order.
//
// This command is useful when you want to add quality annotations to headers or
// filter sequences based on quality thresholds, but don't need to reorder them.
// It processes records in a streaming fashion, making it memory-efficient for
// large input files
func NoSortCommand() *cobra.Command {
	var (
		inFile        string
		outFile       string
		metric        string
		minPhred      int
		minQualFilter float64
		maxQualFilter float64
		headerMetrics string
	)

	cmd := &cobra.Command{
		Use:   "nosort",
		Short: "Estimate sequence quality without sorting (optional filtering and header annotation)",
		Long: `Estimate sequence quality metrics for each FASTQ record and optionally filter or
annotate sequence headers, while preserving the original input order. Supports the
same quality metrics and header annotations as the default sorting mode, but does
not perform any reordering of records.`,
		RunE: func(cmd *cobra.Command, args []string) error {

			// Validate metric flag
			qualityMetric, err := validateMetric(metric)
			if err != nil {
				return err
			}

			// Parse header metrics
			parsedHeaderMetrics, err := parseHeaderMetrics(headerMetrics)
			if err != nil {
				return err
			}

			return runNoSort(
				inFile,
				outFile,
				qualityMetric,
				parsedHeaderMetrics,
				minPhred,
				minQualFilter,
				maxQualFilter,
			)
		},
	}

	flags := cmd.Flags()
	flags.StringVarP(&inFile, "in", "i", "-", "Input FASTQ file (default: stdin)")
	flags.StringVarP(&outFile, "out", "o", "-", "Output FASTQ file (default: stdout)")
	flags.StringVarP(&metric, "metric", "s", "avgphred", "Quality metric (avgphred, maxee, meep, lqcount, lqpercent)")
	flags.IntVarP(&minPhred, "minphred", "p", DEFAULT_MIN_PHRED, "Quality threshold for 'lqcount' and 'lqpercent' metrics")
	flags.Float64VarP(&minQualFilter, "minqual", "m", -math.MaxFloat64, "Minimum quality threshold for filtering")
	flags.Float64VarP(&maxQualFilter, "maxqual", "M", math.MaxFloat64, "Maximum quality threshold for filtering")
	flags.StringVarP(&headerMetrics, "header", "H", "", "Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')")

	return cmd
}

// runNoSort streams records from input to output, computing the requested
// quality metric for each record, applying quality filters, and optionally
// appending additional metrics to the header. Record order is preserved
//
// This function reads records sequentially and writes them immediately after
// processing, making it suitable for very large files that don't fit in memory.
// Records that don't meet the quality thresholds are filtered out and not written
//
// Parameters:
//   - inFile: Input FASTQ file path (use "-" for stdin)
//   - outFile: Output FASTQ file path (use "-" for stdout)
//   - metric: Quality metric to calculate for filtering
//   - headerMetrics: Optional metrics to append to headers
//   - minPhred: Minimum Phred threshold for lqcount/lqpercent calculations
//   - minQualFilter: Minimum quality threshold for filtering
//   - maxQualFilter: Maximum quality threshold for filtering
//
// Returns an error if file I/O operations fail
func runNoSort(
	inFile, outFile string,
	metric QualityMetric,
	headerMetrics []HeaderMetric,
	minPhred int,
	minQualFilter, maxQualFilter float64,
) error {
	reader, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		return fmt.Errorf("error creating reader: %v", err)
	}
	defer reader.Close()

	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		return fmt.Errorf("error creating output file: %v", err)
	}
	defer outfh.Close()

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("error reading record: %v", err)
		}

		quality := calculateQuality(record, metric, minPhred)
		// writeRecord handles header annotation and filtering
		writeRecord(outfh, record, quality, headerMetrics, metric, minPhred, minQualFilter, maxQualFilter)
	}

	return nil
}
