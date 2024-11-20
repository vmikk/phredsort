package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/fatih/color"

	"github.com/klauspost/compress/zstd"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

const (
	PHRED_OFFSET = 33
	VERSION      = "0.5.0"
)

// QualityMetric represents different methods for calculating sequence quality
type QualityMetric int

const (
	AvgPhred QualityMetric = iota
	MaxEE
	Meep
)

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

// QualityFloatList is a slice of QualityFloat
type QualityFloatList []QualityFloat

func (list QualityFloatList) Len() int { return len(list) }
func (list QualityFloatList) Less(i, j int) bool {
	metric := list.getMetricFromValue()
	if metric == MaxEE || metric == Meep {
		if list[i].Value > list[j].Value {
			return true
		}
		if list[i].Value == list[j].Value {
			return list[i].Name < list[j].Name
		}
		return false
	} else {
		if list[i].Value < list[j].Value {
			return true
		}
		if list[i].Value == list[j].Value {
			return list[i].Name < list[j].Name
		}
		return false
	}
}
func (list QualityFloatList) Swap(i, j int) { list[i], list[j] = list[j], list[i] }

// ReversedQualityFloatList is for reverse sorting
type ReversedQualityFloatList struct {
	QualityFloatList
}

func (list ReversedQualityFloatList) Less(i, j int) bool {
	metric := list.QualityFloatList.getMetricFromValue()
	if metric == MaxEE || metric == Meep {
		if list.QualityFloatList[i].Value < list.QualityFloatList[j].Value {
			return true
		}
		if list.QualityFloatList[i].Value == list.QualityFloatList[j].Value {
			return list.QualityFloatList[i].Name < list.QualityFloatList[j].Name
		}
		return false
	} else {
		if list.QualityFloatList[i].Value > list.QualityFloatList[j].Value {
			return true
		}
		if list.QualityFloatList[i].Value == list.QualityFloatList[j].Value {
			return list.QualityFloatList[i].Name < list.QualityFloatList[j].Name
		}
		return false
	}
}

var errorProbs [256]float64

func init() {
	// Pre-compute error probabilities for Phred scores
	for i := range errorProbs {
		errorProbs[i] = math.Pow(10, float64(i-PHRED_OFFSET)/-10)
	}
}

// Sum of error probabilities for quality scores
func sumErrorProbs(qual []byte) float64 {
	var sum float64
	for _, q := range qual {
		sum += errorProbs[q]
	}
	return sum
}

// Average Phred score from quality scores
func calculateAvgPhred(qual []byte) float64 {
	if len(qual) == 0 {
		return 0.0
	}
	meanProb := sumErrorProbs(qual) / float64(len(qual))
	return -10 * math.Log10(meanProb)
}

// Maximum expected error (absolute number)
func calculateMaxEE(qual []byte) float64 {
	if len(qual) == 0 {
		return math.Inf(1) // Return positive infinity for zero-length sequences
	}
	return sumErrorProbs(qual)
}

// Maximum expected error rate (percentage per sequence length)
func calculateMeep(qual []byte) float64 {
	if len(qual) == 0 {
		return math.Inf(1) // Return positive infinity for zero-length sequences
	}
	return (sumErrorProbs(qual) * 100) / float64(len(qual))
}

// Wrapper function to calculate quality scores
func calculateQuality(record *fastx.Record, metric QualityMetric) float64 {
	switch metric {
	case AvgPhred:
		return calculateAvgPhred(record.Seq.Qual)
	case MaxEE:
		return calculateMaxEE(record.Seq.Qual)
	case Meep:
		return calculateMeep(record.Seq.Qual)
	default:
		return 0
	}
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

func main() {
	var (
		inFile    string
		outFile   string
		metric    string
		ascending bool
		compLevel int
		version   bool
	)

	// Create custom help function
	helpFunc := func(cmd *cobra.Command, args []string) {
		fmt.Printf(`
%s

%s
  %s
  %s
  %s

%s
  %s
  %s
  %s
  %s
  %s
  %s
  %s

%s
  # File-based mode (reads from a file, lower memory usage)
  %s

  # Stdin-based mode (reads from stdin, higher memory usage)
  %s

%s
  https://github.com/vmikk/phredsort

`,
			bold(getColorizedLogo()+" phredsort v."+VERSION+" - Sorts FASTQ based on different sequence quality metrics"),
			bold(yellow("Quality metrics:")),
			cyan("avgphred")+": average Phred quality score",
			cyan("maxee")+":    maximum expected error (absolute number)",
			cyan("meep")+":     maximum expected error (percentage per sequence length)",
			bold(yellow("Flags:")),
			cyan("-i, --in")+" <string>      : Input FASTQ file (required, use `-` for stdin)",
			cyan("-o, --out")+" <string>     : Output FASTQ file (required, use `-` for stdout)",
			cyan("-m, --metric")+" <string>  : Quality metric (avgphred, maxee, meep) (default, `avgphred`)",
			cyan("-a, --ascending")+" <bool> : Sort sequences in ascending order of quality (default, false)",
			cyan("-c, --compress")+" <int>   : Memory compression level for stdin-based mode (0=disabled, 1-22; default, 1)",
			cyan("-h, --help")+"             : Show help message",
			cyan("-v, --version")+"          : Show version information",
			bold(yellow("Usage examples:")),
			cyan("phredsort --metric avgphred --in input.fq.gz --out output.fq.gz"),
			cyan("cat input.fq | phredsort --compress 0 -i - -o - > sorted.fq"),
			bold(yellow("More information:")),
		)
	}

	rootCmd := &cobra.Command{
		Use:   "phredsort",
		Short: bold("Sorts FASTQ files by quality metrics"),
		Run: func(cmd *cobra.Command, args []string) {
			// Check version flag
			if version {
				fmt.Printf("phredsort %s\n", VERSION)
				os.Exit(0)
			}

			// Check required flags
			if inFile == "" || outFile == "" {
				fmt.Fprintln(os.Stderr, red("Error: input and output files are required"))
				fmt.Fprintln(os.Stderr, red("Try 'phredsort --help' for more information"))
				os.Exit(1)
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
			default:
				fmt.Fprintf(os.Stderr, red("Error: invalid metric '%s'. Must be one of: avgphred, maxee, meep\n"), metric)
				os.Exit(1)
			}

			// Validate compression level
			if compLevel < 0 || compLevel > 22 {
				fmt.Fprintln(os.Stderr, red("Error: compression level must be between 0 and 22"))
				os.Exit(1)
			}

			// Process the files
			if inFile == "-" {
				sortStdin(outFile, ascending, qualityMetric, compLevel)
			} else {
				sortFile(inFile, outFile, ascending, qualityMetric)
			}
		},
	}

	// Set the help function
	rootCmd.SetHelpFunc(helpFunc)

	// Define flags
	flags := rootCmd.Flags()
	flags.StringVarP(&inFile, "in", "i", "", "Input FASTQ file (required, use - for stdin)")
	flags.StringVarP(&outFile, "out", "o", "", "Output FASTQ file (required)")
	flags.StringVarP(&metric, "metric", "m", "avgphred", "Quality metric (avgphred, maxee, meep)")
	flags.BoolVarP(&ascending, "ascending", "a", false, "Sort sequences in ascending order of quality (default: descending)")
	flags.IntVarP(&compLevel, "compress", "c", 1, "Memory compression level for stdin-based mode (0=disabled, 1-22; default: 1)")
	flags.BoolVarP(&version, "version", "v", false, "Show version information")

	// Custom error handling
	if err := rootCmd.Execute(); err != nil {
		fmt.Fprintln(os.Stderr, red(err.Error()))
		fmt.Fprintln(os.Stderr, red("Try 'phredsort --help' for more information"))
		os.Exit(1)
	}
}

type CompressedFastqRecord struct {
	Name    []byte
	Data    []byte
	AvgQual float64
}

func sortStdin(outFile string, ascending bool, metric QualityMetric, compLevel int) {
	reader, err := fastx.NewReader(seq.DNAredundant, "-", fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating reader: %v\n"), err)
		os.Exit(1)
	}
	defer reader.Close()

	var name2avgQual []QualityFloat

	// Create output file handle at the beginning
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		os.Exit(1)
	}
	defer outfh.Close()

	if compLevel > 0 {
		encoder, err := zstd.NewWriter(nil, zstd.WithEncoderLevel(zstd.EncoderLevelFromZstd(compLevel)))
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD encoder: %v\n"), err)
			os.Exit(1)
		}
		defer encoder.Close()

		decoder, err := zstd.NewReader(nil)
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating ZSTD decoder: %v\n"), err)
			os.Exit(1)
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
				os.Exit(1)
			}

			name := string(record.Name)
			avgQual := calculateQuality(record, metric)

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
		if ascending {
			sort.Sort(QualityFloatList(name2avgQual))
		} else {
			sort.Sort(ReversedQualityFloatList{QualityFloatList(name2avgQual)})
		}

		// Writing records
		for _, kv := range name2avgQual {
			compRecord := sequences[kv.Name]
			decompressed, err := decoder.DecodeAll(compRecord.Data, nil)
			if err != nil {
				fmt.Fprintf(os.Stderr, red("Error decompressing record: %v\n"), err)
				os.Exit(1)
			}

			seqLen := len(decompressed) / 2
			record := &fastx.Record{
				Name: compRecord.Name,
				Seq: &seq.Seq{
					Seq:  decompressed[:seqLen],
					Qual: decompressed[seqLen:],
				},
			}
			record.FormatToWriter(outfh, 0)
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
				os.Exit(1)
			}

			name := string(record.Name)
			avgQual := calculateQuality(record, metric)

			// Important: Clone the record to avoid reference issues
			sequences[name] = record.Clone()
			name2avgQual = append(name2avgQual, QualityFloat{
				Name:   name,
				Value:  avgQual,
				Metric: metric,
			})
		}

		// Sort records
		if ascending {
			sort.Sort(QualityFloatList(name2avgQual))
		} else {
			sort.Sort(ReversedQualityFloatList{QualityFloatList(name2avgQual)})
		}

		// Write sorted records
		outfh, err := xopen.Wopen(outFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
			os.Exit(1)
		}
		defer outfh.Close()

		// Output in sorted order
		for _, kv := range name2avgQual {
			record := sequences[kv.Name]
			record.FormatToWriter(outfh, 0)
		}
	}
}

func sortFile(inFile, outFile string, ascending bool, metric QualityMetric) {
	reader, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating reader: %v\n"), err)
		os.Exit(1)
	}
	defer reader.Close()

	// First pass: collect quality scores and positions
	qualityScores := QualityFloatList{}
	name2offset := make(map[string]int64)
	var currentOffset int64 = 0

	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			os.Exit(1)
		}

		name := string(record.Name)
		avgQual := calculateQuality(record, metric)

		qualityScores = append(qualityScores, QualityFloat{
			Name:   name,
			Value:  avgQual,
			Metric: metric,
		})
		name2offset[name] = currentOffset
		currentOffset++
	}

	// Sort by average quality
	if ascending {
		sort.Sort(QualityFloatList(qualityScores))
	} else {
		sort.Sort(ReversedQualityFloatList{qualityScores})
	}

	// Second pass: write records in sorted order
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		os.Exit(1)
	}
	defer outfh.Close()

	reader2, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating second reader: %v\n"), err)
		os.Exit(1)
	}
	defer reader2.Close()

	// Read all records into a map
	records := make(map[int64]*fastx.Record)
	var offset int64 = 0
	for {
		record, err := reader2.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			os.Exit(1)
		}
		records[offset] = record.Clone()
		offset++
	}

	// Write records in sorted order
	outfh, err = xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, red("Error creating output file: %v\n"), err)
		os.Exit(1)
	}
	defer outfh.Close()

	// Output in sorted order
	for _, qf := range qualityScores {
		offset := name2offset[qf.Name]
		record, ok := records[offset]
		if !ok {
			fmt.Fprintf(os.Stderr, red("Error: could not find record for %s\n"), qf.Name)
			os.Exit(1)
		}
		record.FormatToWriter(outfh, 0)
	}
}

// Helper function to get metric from QualityFloatList
func (list QualityFloatList) getMetricFromValue() QualityMetric {
	if len(list) > 0 {
		return list[0].Metric
	}
	return AvgPhred // Default
}
