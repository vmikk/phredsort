package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/fatih/color"
	"github.com/maruel/natural"

	"github.com/klauspost/compress/zstd"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
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

// Metrics to add to headers
type HeaderMetric struct {
	Name     string
	IsLength bool
}

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
			// Validate metric name
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

// The main sorting function logic (file-based mode or stdin-based mode)
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
		fmt.Fprintf(os.Stderr, red("Error: invalid metric '%s'. Must be one of: avgphred, maxee, meep, lqcount, lqpercent"), metric)
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
		sortFile(inFile, outFile, ascending, qualityMetric, parsedHeaderMetrics, minPhred, minQualFilter, maxQualFilter)
	}
}

type CompressedFastqRecord struct {
	Name    []byte
	Data    []byte
	AvgQual float64
}

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

func sortFile(inFile, outFile string, ascending bool, metric QualityMetric, headerMetrics []HeaderMetric, minPhred int, minQualFilter float64, maxQualFilter float64) {
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

	// Open output file and write sorted records
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		exitFunc(1)
	}
	defer outfh.Close()

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
