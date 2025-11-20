// Subcommand (`phredsort nosort`) for sequence quality estimation and optional filtering.
// No sorting is performed; sequences are emitted in the same order as input.

package main

import (
	"fmt"
	"io"
	"math"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// NoSortCommand creates the `nosort` subcommand which estimates quality metrics
// and optionally filters/annotates sequences without changing their order.
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
			if inFile == "" || outFile == "" {
				return fmt.Errorf("input and output files are required")
			}

			// Validate metric flag
			var qualityMetric QualityMetric
			switch strings.ToLower(metric) {
			case "avgphred":
				qualityMetric = AvgPhred
			case "maxee":
				qualityMetric = MaxEE
			case "meep":
				qualityMetric = Meep
			case "lqcount":
				qualityMetric = LQCount
			case "lqpercent":
				qualityMetric = LQPercent
			default:
				return fmt.Errorf("invalid metric '%s'", metric)
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
	flags.StringVarP(&inFile, "in", "i", "", "Input FASTQ file (required, use '-' for stdin)")
	flags.StringVarP(&outFile, "out", "o", "", "Output FASTQ file (required, use '-' for stdout)")
	flags.StringVarP(&metric, "metric", "s", "avgphred", "Quality metric (avgphred, maxee, meep, lqcount, lqpercent)")
	flags.IntVarP(&minPhred, "minphred", "p", DEFAULT_MIN_PHRED, "Quality threshold for 'lqcount' and 'lqpercent' metrics")
	flags.Float64VarP(&minQualFilter, "minqual", "m", -math.MaxFloat64, "Minimum quality threshold for filtering")
	flags.Float64VarP(&maxQualFilter, "maxqual", "M", math.MaxFloat64, "Maximum quality threshold for filtering")
	flags.StringVarP(&headerMetrics, "header", "H", "", "Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')")

	_ = cmd.MarkFlagRequired("in")
	_ = cmd.MarkFlagRequired("out")

	return cmd
}

// runNoSort streams records from input to output, computing the requested
// quality metric for each record, applying quality filters, and optionally
// appending additional metrics to the header. Record order is preserved.
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


