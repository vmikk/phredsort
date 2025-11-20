package main

import (
	"fmt"
	"math"
	"os"
	"strings"

	"github.com/fatih/color"
	"github.com/maruel/natural"
	"github.com/spf13/cobra"
)

const (
	VERSION           = "1.4.0"
	PHRED_OFFSET      = 33
	DEFAULT_MIN_PHRED = 15 // min Phred score threshold for `lqcount` and `lqpercent` metrics
)

// Mock exit function for testing
var exitFunc = os.Exit

// QualityMetric represents different methods for calculating sequence quality
type QualityMetric int

const (
	AvgPhred QualityMetric = iota
	MaxEE
	Meep
	LQCount
	LQPercent
)

// Add this near the QualityMetric type definition at the top of the file
func (m QualityMetric) String() string {
	switch m {
	case AvgPhred:
		return "avgphred"
	case MaxEE:
		return "maxee"
	case Meep:
		return "meep"
	case LQCount:
		return "lqcount"
	case LQPercent:
		return "lqpercent"
	default:
		return "unknown"
	}
}

// QualityRecord stores just the essential info for sorting
type QualityRecord struct {
	Offset  int64   // File offset for a record
	AvgQual float64 // Average quality score
}

// QualityFloat is for sorting sequences by quality scores
type QualityFloat struct {
	Name   string
	Value  float64
	Metric QualityMetric
}

// QualityFloatList is a slice of QualityFloat with sorting direction
type QualityFloatList struct {
	items     []QualityFloat
	ascending bool // true for ascending, false for descending
}

// New constructor function
func NewQualityFloatList(items []QualityFloat, ascending bool) QualityFloatList {
	return QualityFloatList{items: items, ascending: ascending}
}

func (list QualityFloatList) Len() int { return len(list.items) }
func (list QualityFloatList) Swap(i, j int) {
	list.items[i], list.items[j] = list.items[j], list.items[i]
}
func (list QualityFloatList) Less(i, j int) bool {
	metric := list.getMetricFromValue()

	// Compare values based on metric type
	var result bool
	if metric == MaxEE || metric == Meep || metric == LQCount || metric == LQPercent {
		// For these metrics, higher values indicate lower quality
		if list.items[i].Value != list.items[j].Value {
			result = list.items[i].Value < list.items[j].Value
		} else {
			result = natural.Less(list.items[i].Name, list.items[j].Name)
		}
	} else {
		// For other metrics (e.g., AvgPhred), higher values indicate better quality
		if list.items[i].Value != list.items[j].Value {
			result = list.items[i].Value > list.items[j].Value
		} else {
			result = natural.Less(list.items[i].Name, list.items[j].Name)
		}
	}

	// Flip the result if we want ascending order
	if list.ascending {
		return !result
	}
	return result
}

// Helper function to get metric from QualityFloatList
func (list QualityFloatList) getMetricFromValue() QualityMetric {
	if len(list.items) > 0 {
		return list.items[0].Metric
	}
	return AvgPhred // Default
}

// Helper function to get the underlying items
func (list QualityFloatList) Items() []QualityFloat {
	return list.items
}

// Define color functions
var (
	bold   = color.New(color.Bold).SprintFunc()
	cyan   = color.New(color.FgCyan).SprintFunc()
	yellow = color.New(color.FgYellow).SprintFunc()
	red    = color.New(color.FgRed).SprintFunc()
)

func getColorizedLogo() string {
	symbols := []rune{'⣿', '⣶', '⣦', '⣄', '⣀'}
	var logo strings.Builder

	// Create gradient by assigning specific green intensities to each symbol
	colors := []struct{ r, g, b int }{
		{0, 255, 0}, // Brightest green
		{0, 204, 0},
		{0, 153, 0},
		{0, 102, 0},
		{0, 51, 0}, // Darkest green
	}

	for i, symbol := range symbols {
		logo.WriteString(color.RGB(colors[i].r, colors[i].g, colors[i].b).Sprint(string(symbol)))
	}

	return logo.String()
}

// Variable declarations (at package level)
var (
	// Command-line flags for the default sorting command
	inFile        string
	outFile       string
	metric        string
	minPhred      int
	minQualFilter float64
	maxQualFilter float64
	headerMetrics string
	ascending     bool
	compLevel     int
	version       bool
)

func main() {
	// Create root command
	rootCmd := &cobra.Command{
		Use:   "phredsort",
		Short: bold("Sorts FASTQ files by quality metrics"),
		// When no subcommand is specified, run the default sorting behavior
		Run: runDefaultCommand,
	}

	// The default command = quality estimation and sorting
	defaultCmd := &cobra.Command{
		Use:   "sort",
		Short: "Sort sequences by calculating quality metrics",
		Run:   runDefaultCommand,
	}

	// Define flags for the default sorting behavior
	// The same flags are registered on both the root command and the explicit "sort" subcommand,
	// so that both of the following commands are supported:
	//  phredsort      -i in.fq -o out.fq ...
	//  phredsort sort -i in.fq -o out.fq ...

	rootFlags := rootCmd.Flags()
	rootFlags.StringVarP(&inFile, "in", "i", "", "Input FASTQ file (required, use - for stdin)")
	rootFlags.StringVarP(&outFile, "out", "o", "", "Output FASTQ file (required)")
	rootFlags.StringVarP(&metric, "metric", "s", "avgphred", "Quality metric (avgphred, maxee, meep, lqcount, lqpercent)")
	rootFlags.IntVarP(&minPhred, "minphred", "p", DEFAULT_MIN_PHRED, "Quality threshold for 'lqcount' and 'lqpercent' metrics")
	rootFlags.Float64VarP(&minQualFilter, "minqual", "m", -math.MaxFloat64, "Minimum quality threshold for filtering")
	rootFlags.Float64VarP(&maxQualFilter, "maxqual", "M", math.MaxFloat64, "Maximum quality threshold for filtering")
	rootFlags.StringVarP(&headerMetrics, "header", "H", "", "Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')")
	rootFlags.BoolVarP(&ascending, "ascending", "a", false, "Sort sequences in ascending order of quality (default: descending)")
	rootFlags.IntVarP(&compLevel, "compress", "c", 1, "Memory compression level for stdin-based mode (0=disabled, 1-22; default: 1)")
	rootFlags.BoolVarP(&version, "version", "v", false, "Show version information")

	sortFlags := defaultCmd.Flags()
	sortFlags.StringVarP(&inFile, "in", "i", "", "Input FASTQ file (required, use - for stdin)")
	sortFlags.StringVarP(&outFile, "out", "o", "", "Output FASTQ file (required)")
	sortFlags.StringVarP(&metric, "metric", "s", "avgphred", "Quality metric (avgphred, maxee, meep, lqcount, lqpercent)")
	sortFlags.IntVarP(&minPhred, "minphred", "p", DEFAULT_MIN_PHRED, "Quality threshold for 'lqcount' and 'lqpercent' metrics")
	sortFlags.Float64VarP(&minQualFilter, "minqual", "m", -math.MaxFloat64, "Minimum quality threshold for filtering")
	sortFlags.Float64VarP(&maxQualFilter, "maxqual", "M", math.MaxFloat64, "Maximum quality threshold for filtering")
	sortFlags.StringVarP(&headerMetrics, "header", "H", "", "Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')")
	sortFlags.BoolVarP(&ascending, "ascending", "a", false, "Sort sequences in ascending order of quality (default: descending)")
	sortFlags.IntVarP(&compLevel, "compress", "c", 1, "Memory compression level for stdin-based mode (0=disabled, 1-22; default: 1)")
	sortFlags.BoolVarP(&version, "version", "v", false, "Show version information")

	// Add commands
	rootCmd.AddCommand(defaultCmd)          // sort using quality estimation
	rootCmd.AddCommand(NoSortCommand())     // estimate quality without sorting
	rootCmd.AddCommand(HeaderSortCommand()) // sort using pre-computed quality scores

	// Set help function
	rootCmd.SetHelpFunc(helpFunc)

	// Execute
	if err := rootCmd.Execute(); err != nil {
		fmt.Fprintln(os.Stderr, red(err.Error()))
		fmt.Fprintln(os.Stderr, red("Try 'phredsort --help' for more information"))
		exitFunc(1)
	}
}
