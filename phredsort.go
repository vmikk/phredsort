package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

const (
	PHRED_OFFSET = 33
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

func main() {
	seq.ValidateSeq = false

	inFile := flag.String("in", "", "Input FASTQ file (required, use - for stdin)")
	outFile := flag.String("out", "", "Output FASTQ file (required)")
	reverse := flag.Bool("reverse", true, "Sort in descending order")
	flag.Parse()

	if *inFile == "" || *outFile == "" {
		flag.Usage()
		os.Exit(1)
	}

	// For stdin, we need to store complete records
	if *inFile == "-" {
		sortStdin(*outFile, *reverse)
		return
	}
	// For files, use the two-pass approach
	sortFile(*inFile, *outFile, *reverse)
}

type FastqRecord struct {
	Name    []byte
	Seq     []byte
	Qual    []byte
	AvgQual float64
}

func sortStdin(outFile string, reverse bool) {
	reader, err := fastx.NewReader(seq.DNAredundant, "-", fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating reader: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	// Store sequences in a map to maintain association with quality scores
	sequences := make(map[string]*fastx.Record)
	name2avgQual := []QualityFloat{}

	// Read all records
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintf(os.Stderr, "Error reading record: %v\n", err)
			os.Exit(1)
		}

		name := string(record.Name)
		avgQual := record.Seq.AvgQual(PHRED_OFFSET)

		// Store complete record
		sequences[name] = record.Clone()
		name2avgQual = append(name2avgQual, QualityFloat{Name: name, Value: avgQual})
	}

	// Sort by average quality
	if reverse {
		sort.Sort(ReversedQualityFloatList{QualityFloatList(name2avgQual)})
	} else {
		sort.Sort(QualityFloatList(name2avgQual))
	}

	// Write sorted records
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer outfh.Close()

	// Output in sorted order
	for _, kv := range name2avgQual {
		record := sequences[kv.Name]
		record.FormatToWriter(outfh, 0)
	}
}

func sortFile(inFile, outFile string, reverse bool) {
	reader, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating reader: %v\n", err)
		os.Exit(1)
	}
	defer reader.Close()

	// First pass: collect quality scores and positions
	qualityScores := QualityFloatList{}
	name2offset := make(map[string]int64)
	var currentOffset int64 = 0

	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintf(os.Stderr, "Error reading record: %v\n", err)
			os.Exit(1)
		}

		name := string(record.Name)
		avgQual := record.Seq.AvgQual(PHRED_OFFSET)

		qualityScores = append(qualityScores, QualityFloat{Name: name, Value: avgQual})
		name2offset[name] = currentOffset
		currentOffset++
	}

	// Sort by average quality
	if reverse {
		sort.Sort(ReversedQualityFloatList{qualityScores})
	} else {
		sort.Sort(qualityScores)
	}

	// Second pass: write records in sorted order
	outfh, err := xopen.Wopen(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
		os.Exit(1)
	}
	defer outfh.Close()

	reader2, err := fastx.NewReader(seq.DNAredundant, inFile, fastx.DefaultIDRegexp)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error creating second reader: %v\n", err)
		os.Exit(1)
	}
	defer reader2.Close()

	// Read all records into a map first
	records := make(map[int64]*fastx.Record)
	var offset int64 = 0
	for {
		record, err := reader2.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			fmt.Fprintf(os.Stderr, "Error reading record: %v\n", err)
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
			fmt.Fprintf(os.Stderr, "Error: could not find record for %s\n", qf.Name)
			os.Exit(1)
		}
		record.FormatToWriter(outfh, 0)
	}
}
