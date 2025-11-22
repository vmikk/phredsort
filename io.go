// Phredsort I/O utilities for FASTQ/FASTA record handling

package main

import (
	"fmt"
	"io"
	"strings"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

// HeaderMetric represents additional metrics that can be appended to FASTQ/FASTA
// headers when writing sorted output. The IsLength field indicates whether this
// metric represents sequence length rather than a quality metric
type HeaderMetric struct {
	Name     string
	IsLength bool
}

// parseHeaderMetrics parses a comma-separated string of metric names into a
// slice of HeaderMetric structs. Validates that all metric names are supported.
// Returns an error if any metric name is invalid
//
// Supported metrics: avgphred, maxee, meep, lqcount, lqpercent, length
//
// Example:
//   parseHeaderMetrics("avgphred,maxee,length") // Returns 3 HeaderMetric structs
func parseHeaderMetrics(metrics string) ([]HeaderMetric, error) {
	if metrics == "" {
		return nil, nil
	}

	parts := strings.Split(metrics, ",")
	result := make([]HeaderMetric, 0, len(parts))

	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}

		hm := HeaderMetric{Name: p}
		if p == "length" {
			hm.IsLength = true
		} else {
			switch p {
			case "avgphred", "maxee", "meep", "lqcount", "lqpercent":
				// valid metric
			default:
				return nil, fmt.Errorf("Error: invalid header metric: %s", p)
			}
		}
		result = append(result, hm)
	}

	return result, nil
}

// writeRecord writes a FASTQ/FASTA record to the output writer, applying quality
// filters and optionally appending header annotations. Returns true if the record
// was written (passed filters), false if it was filtered out
//
// The function:
//   - Filters records based on minQualFilter and maxQualFilter thresholds
//   - Optionally appends quality metrics and sequence length to the header
//   - Writes the record in FASTQ/FASTA format
//
// Parameters:
//   - outfh: Output writer (must be *xopen.Writer)
//   - record: The FASTQ/FASTA record to write
//   - quality: The calculated quality value for the record
//   - headerMetrics: List of metrics to append to the header (nil/empty = no annotation)
//   - metric: The quality metric type used (for context, not recalculated)
//   - minPhred: Minimum Phred threshold for lqcount/lqpercent calculations
//   - minQualFilter: Minimum quality threshold for filtering (records below this are skipped)
//   - maxQualFilter: Maximum quality threshold for filtering (records above this are skipped)
func writeRecord(outfh io.Writer, record *fastx.Record, quality float64, headerMetrics []HeaderMetric, metric QualityMetric, minPhred int, minQualFilter float64, maxQualFilter float64) bool {
	// Skip records that don't meet quality thresholds
	if quality < minQualFilter || quality > maxQualFilter {
		return false
	}

	if len(headerMetrics) > 0 {
		var additions []string

		for _, hm := range headerMetrics {
			if hm.IsLength {
				additions = append(additions, fmt.Sprintf("length=%d", len(record.Seq.Seq)))
			} else {
				// Calculate the requested metric
				var metricValue float64
				switch hm.Name {
				case "avgphred":
					metricValue = calculateAvgPhred(record.Seq.Qual)
				case "maxee":
					metricValue = calculateMaxEE(record.Seq.Qual)
				case "meep":
					metricValue = calculateMeep(record.Seq.Qual)
				case "lqcount":
					metricValue = countLowQualityBases(record.Seq.Qual, minPhred)
				case "lqpercent":
					metricValue = calculateLQPercent(record.Seq.Qual, minPhred)
				}
				additions = append(additions, fmt.Sprintf("%s=%.6f", hm.Name, metricValue))
			}
		}

		if len(additions) > 0 {
			record.Name = append(record.Name, " "+strings.Join(additions, " ")...)
		}
	}

	writer := outfh.(*xopen.Writer)
	record.FormatToWriter(writer, 0)
	return true
}

