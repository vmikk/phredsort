package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"

	"github.com/fatih/color"

	"github.com/klauspost/compress/zstd"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

const (
	PHRED_OFFSET = 33
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
	Offset  int64   // File offset for this record
	AvgQual float64 // Average quality score
}

// QualityFloat is for sorting sequences by quality scores
type QualityFloat struct {
	Name  string
	Value float64
}

// QualityFloatList is a slice of QualityFloat
type QualityFloatList []QualityFloat

func (list QualityFloatList) Len() int { return len(list) }
func (list QualityFloatList) Less(i, j int) bool {
	if list[i].Value < list[j].Value {
		return true
	}
	if list[i].Value == list[j].Value {
		// Secondary sort by name for stable sorting
		return list[i].Name < list[j].Name
	}
	return false
}
func (list QualityFloatList) Swap(i, j int) { list[i], list[j] = list[j], list[i] }

// ReversedQualityFloatList is for reverse sorting
type ReversedQualityFloatList struct {
	QualityFloatList
}

func (list ReversedQualityFloatList) Less(i, j int) bool {
	if list.QualityFloatList[i].Value > list.QualityFloatList[j].Value {
		return true
	}
	if list.QualityFloatList[i].Value == list.QualityFloatList[j].Value {
		// Secondary sort by name for stable sorting
		return list.QualityFloatList[i].Name < list.QualityFloatList[j].Name
	}
	return false
}

var errorProbs [256]float64

func init() {
	// Pre-compute error probabilities for Phred scores
	for i := range errorProbs {
		errorProbs[i] = math.Pow(10, float64(i-PHRED_OFFSET)/-10)
	}
}

func calculateQuality(record *fastx.Record, metric QualityMetric) float64 {
	switch metric {
	case AvgPhred:
		return record.Seq.AvgQual(PHRED_OFFSET)
	case MaxEE:
		return calculateMaxEE(record.Seq.Qual)
	case Meep:
		return calculateMeep(record.Seq.Qual)
	default:
		return 0
	}
}

func calculateMaxEE(qual []byte) float64 {
	var sum float64
	for _, q := range qual {
		sum += errorProbs[q]
	}
	return sum
}

func calculateMeep(qual []byte) float64 {
	maxEE := calculateMaxEE(qual)
	return (maxEE * 100) / float64(len(qual))
}

// Define color functions
var (
	bold   = color.New(color.Bold).SprintFunc()
	cyan   = color.New(color.FgCyan).SprintFunc()
	yellow = color.New(color.FgYellow).SprintFunc()
	red    = color.New(color.FgRed).SprintFunc()
)

func main() {
	var (
		inFile    string
		outFile   string
		ascending bool
		metric    string
		compLevel int
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
  # File-based mode (reads from a file, lower memory usage)
  %s

  # Stdin-based mode (reads from stdin, higher memory usage)
  %s

`,
			bold(cyan("⣿⣶⣦⣄⣀ phredsort")+" - Sort FASTQ based on different sequence quality metrics"),
			bold(yellow("Quality metrics:")),
			cyan("avgphred")+": average Phred quality score",
			cyan("maxee")+":    maximum expected error (absolute number)",
			cyan("meep")+":     maximum expected error (percentage per sequence length)",
			bold(yellow("Flags:")),
			cyan("-i, --in")+" <string>      : Input FASTQ file (required, use `-` for stdin)",
			cyan("-o, --out")+" <string>     : Output FASTQ file (required, use `-` for stdout)",
			cyan("-m, --metric")+" <string>  : Quality metric (avgphred, maxee, meep) (default, `avgphred`)",
			cyan("-a, --ascending")+" <bool> : Sort sequences in ascending order (default, false)",
			cyan("-c, --compress")+" <int>   : Memory compression level for stdin-based mode (0=disabled, 1-22)",
			bold(yellow("Usage examples:")),
			cyan("phredsort --metric avgphred --in input.fq.gz --out output.fq.gz"),
			cyan("cat input.fq | phredsort -i - -o - > sorted.fq"))
	}

	rootCmd := &cobra.Command{
		Use:   "phredsort",
		Short: bold("Sort FASTQ files by quality metrics"),
		Run: func(cmd *cobra.Command, args []string) {
			// If no arguments are provided, show help
			if len(args) == 0 && (inFile == "" || outFile == "") {
				helpFunc(cmd, args)
				return
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
	flags.BoolVarP(&ascending, "ascending", "a", false, "Sort sequences in ascending order (default: descending)")
	flags.IntVarP(&compLevel, "compress", "c", 0, "Memory compression level for stdin-based mode (0=disabled, 1-22)")

	// Mark required flags
	rootCmd.MarkFlagRequired("in")
	rootCmd.MarkFlagRequired("out")

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
			if err != nil && err != io.EOF {
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
			name2avgQual = append(name2avgQual, QualityFloat{Name: name, Value: avgQual})
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
			if err != nil && err != io.EOF {
				fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			os.Exit(1)
		}

		name := string(record.Name)
		avgQual := calculateQuality(record, metric)

			// Important: Clone the record to avoid reference issues
		sequences[name] = record.Clone()
		name2avgQual = append(name2avgQual, QualityFloat{Name: name, Value: avgQual})
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
		if err != nil && err != io.EOF {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			os.Exit(1)
		}

		name := string(record.Name)
		avgQual := calculateQuality(record, metric)

		qualityScores = append(qualityScores, QualityFloat{Name: name, Value: avgQual})
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

	// Read all records into a map first
	records := make(map[int64]*fastx.Record)
	var offset int64 = 0
	for {
		record, err := reader2.Read()
		if err != nil && err != io.EOF {
			fmt.Fprintf(os.Stderr, red("Error reading record: %v\n"), err)
			os.Exit(1)
		}
		records[offset] = record.Clone()
		offset++
	}

	// Write records in sorted order
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
