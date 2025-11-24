package main

import (
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"sort"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
	"github.com/spf13/cobra"
)

// Helper function to create test FASTX records
func createTestRecord(name string, sequence string, quality string) *fastx.Record {
	return &fastx.Record{
		Name: []byte(name),
		Seq: &seq.Seq{
			Seq:  []byte(sequence),
			Qual: []byte(quality),
		},
	}
}


// Test sorting functionality
func TestQualityFloatListSorting(t *testing.T) {
	tests := []struct {
		name      string
		items     []QualityFloat
		ascending bool
		metric    QualityMetric
		want      []string // Expected order of names after sorting
	}{
		{
			name: "AvgPhred Descending",
			items: []QualityFloat{
				{Name: "seq1", Value: 30.0, Metric: AvgPhred},
				{Name: "seq2", Value: 40.0, Metric: AvgPhred},
				{Name: "seq3", Value: 20.0, Metric: AvgPhred},
			},
			ascending: false,
			metric:    AvgPhred,
			want:      []string{"seq2", "seq1", "seq3"},
		},
		{
			name: "MaxEE Ascending",
			items: []QualityFloat{
				{Name: "seq1", Value: 0.1, Metric: MaxEE},
				{Name: "seq2", Value: 0.01, Metric: MaxEE},
				{Name: "seq3", Value: 1.0, Metric: MaxEE},
			},
			ascending: true,
			metric:    MaxEE,
			want:      []string{"seq3", "seq1", "seq2"},
		},
		{
			name: "MaxEE - Equal values, natural sort by name",
			items: []QualityFloat{
				{Name: "seq10", Value: 0.1, Metric: MaxEE},
				{Name: "seq2", Value: 0.1, Metric: MaxEE},
				{Name: "seq1", Value: 0.1, Metric: MaxEE},
			},
			ascending: false,
			metric:    MaxEE,
			want:      []string{"seq1", "seq2", "seq10"},
		},
		{
			name: "Meep - Mixed values ascending",
			items: []QualityFloat{
				{Name: "seq1", Value: 5.0, Metric: Meep},
				{Name: "seq2", Value: 2.0, Metric: Meep},
				{Name: "seq3", Value: 10.0, Metric: Meep},
			},
			ascending: true,
			metric:    Meep,
			want:      []string{"seq3", "seq1", "seq2"},
		},
		{
			name: "LQPercent - Zero values",
			items: []QualityFloat{
				{Name: "seq1", Value: 0.0, Metric: LQPercent},
				{Name: "seq2", Value: 0.0, Metric: LQPercent},
			},
			ascending: false,
			metric:    LQPercent,
			want:      []string{"seq1", "seq2"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			list := NewQualityFloatList(tt.items, tt.ascending)
			sort.Sort(list)

			got := make([]string, len(list.items))
			for i, item := range list.items {
				got[i] = item.Name
			}

			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Sort() got %v, want %v", got, tt.want)
			}
		})
	}
}

// Test header metric parsing
func TestParseHeaderMetrics(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		want    []HeaderMetric
		wantErr bool
	}{
		{
			name:    "Empty string",
			input:   "",
			want:    nil,
			wantErr: false,
		},
		{
			name:  "Valid metrics",
			input: "avgphred,maxee,length",
			want: []HeaderMetric{
				{Name: "avgphred", IsLength: false},
				{Name: "maxee", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantErr: false,
		},
		{
			name:    "Invalid metric",
			input:   "avgphred,invalid,length",
			want:    nil,
			wantErr: true,
		},
		{
			name:  "Multiple length metrics",
			input: "length,avgphred,length",
			want: []HeaderMetric{
				{Name: "length", IsLength: true},
				{Name: "avgphred", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantErr: false,
		},
		{
			name:  "Whitespace handling",
			input: " avgphred , maxee , length ",
			want: []HeaderMetric{
				{Name: "avgphred", IsLength: false},
				{Name: "maxee", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantErr: false,
		},
		{
			name:    "Mixed valid and invalid",
			input:   "avgphred,invalid1,maxee,invalid2",
			want:    nil,
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := parseHeaderMetrics(tt.input)
			if (err != nil) != tt.wantErr {
				t.Errorf("parseHeaderMetrics() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("parseHeaderMetrics() = %v, want %v", got, tt.want)
			}
		})
	}
}

// Test record writing with quality filters
func TestWriteRecord(t *testing.T) {
	tests := []struct {
		name          string
		record        *fastx.Record
		quality       float64
		minQualFilter float64
		maxQualFilter float64
		headerMetrics []HeaderMetric
		wantWrite     bool
		wantHeader    string
	}{
		{
			name:          "Quality within bounds",
			record:        createTestRecord("test1", "ACGT", "IIII"),
			quality:       30.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			headerMetrics: []HeaderMetric{
				{Name: "avgphred", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantWrite:  true,
			wantHeader: "test1 avgphred=40.000000 length=4",
		},
		{
			name:          "Quality below minimum",
			record:        createTestRecord("test2", "ACGT", "$$$$"),
			quality:       10.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			wantWrite:     false,
			wantHeader:    "",
		},
		{
			name:          "Quality above maximum",
			record:        createTestRecord("test3", "ACGT", "IIII"),
			quality:       45.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			wantWrite:     false,
			wantHeader:    "",
		},
		{
			name:          "No header metrics",
			record:        createTestRecord("test4", "ACGT", "IIII"),
			quality:       30.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			headerMetrics: nil,
			wantWrite:     true,
			wantHeader:    "test4",
		},
		{
			name:          "Header with maxee metric",
			record:        createTestRecord("test5", "ACGT", "IIII"),
			quality:       30.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			headerMetrics: []HeaderMetric{
				{Name: "maxee", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantWrite:  true,
			wantHeader: "test5 maxee=0.000400 length=4",
		},
		{
			name:          "Header with meep metric",
			record:        createTestRecord("test5", "ACGT", "IIII"),
			quality:       30.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			headerMetrics: []HeaderMetric{
				{Name: "meep", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantWrite:  true,
			wantHeader: "test5 meep=0.010000 length=4",
		},
		{
			name:          "Header with lqpercent metric",
			record:        createTestRecord("test6", "ACGT", "II$$"),
			quality:       30.0,
			minQualFilter: 20.0,
			maxQualFilter: 40.0,
			headerMetrics: []HeaderMetric{
				{Name: "lqpercent", IsLength: false},
				{Name: "length", IsLength: true},
			},
			wantWrite:  true,
			wantHeader: "test6 lqpercent=50.000000 length=4",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Create a temporary file for testing
			tmpfile, err := os.CreateTemp("", "test*.fastq")
			if err != nil {
				t.Fatal(err)
			}
			defer os.Remove(tmpfile.Name())
			defer tmpfile.Close()

			// Create writer using the temp file
			writer, err := xopen.Wopen(tmpfile.Name())
			if err != nil {
				t.Fatal(err)
			}
			defer writer.Close()

			// Test writeRecord
			got := writeRecord(writer, tt.record, tt.quality, tt.headerMetrics, AvgPhred, DEFAULT_MIN_PHRED, tt.minQualFilter, tt.maxQualFilter)

			if got != tt.wantWrite {
				t.Errorf("writeRecord() = %v, want %v", got, tt.wantWrite)
			}

			// If the record should be written, verify the header
			if tt.wantWrite {
				// Close the writer to ensure all data is written
				writer.Close()

				// Read the file content
				content, err := os.ReadFile(tmpfile.Name())
				if err != nil {
					t.Fatal(err)
				}

				// Extract the header from the FASTQ format (first line)
				lines := strings.Split(string(content), "\n")
				if len(lines) > 0 {
					gotHeader := strings.TrimPrefix(lines[0], "@")
					if gotHeader != tt.wantHeader {
						t.Errorf("Header = %q, want %q", gotHeader, tt.wantHeader)
					}
				}
			}
		})
	}
}

// Test quality metric string representation
func TestQualityMetricString(t *testing.T) {
	tests := []struct {
		metric QualityMetric
		want   string
	}{
		{AvgPhred, "avgphred"},
		{MaxEE, "maxee"},
		{Meep, "meep"},
		{LQCount, "lqcount"},
		{LQPercent, "lqpercent"},
		{QualityMetric(999), "unknown"},
	}

	for _, tt := range tests {
		t.Run(tt.want, func(t *testing.T) {
			if got := tt.metric.String(); got != tt.want {
				t.Errorf("QualityMetric.String() = %v, want %v", got, tt.want)
			}
		})
	}
}
// TestSortFile tests the file-based sorting functionality
func TestSortFile(t *testing.T) {
	tests := []struct {
		name          string
		records       []*fastx.Record
		metric        QualityMetric
		ascending     bool
		compLevel     int
		headerMetrics []HeaderMetric
		minPhred      int
		minQual       float64
		maxQual       float64
		wantOrder     []string
		wantErr       bool
	}{
		{
			name: "Basic sorting by AvgPhred descending",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // Phred 40
				createTestRecord("seq2", "ACGT", "$$$$"), // Phred 3
				createTestRecord("seq3", "ACGT", "@@@@"), // Phred 31
			},
			metric:    AvgPhred,
			ascending: false,
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,         // Allow all quality scores
			maxQual:   math.Inf(1), // Allow all quality scores
			wantOrder: []string{"seq1", "seq3", "seq2"},
		},
		{
			name: "MaxEE ascending with quality filters",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // Very low MaxEE
				createTestRecord("seq2", "ACGT", "$$$$"), // Very high MaxEE
				createTestRecord("seq3", "ACGT", "@@@@"), // Medium MaxEE
			},
			metric:    MaxEE,
			ascending: true,
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,
			maxQual:   1.0, // Should filter out seq2
			wantOrder: []string{"seq3", "seq1"},
		},
		{
			name: "LQCount with header metrics",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "III$"),     // 1 low quality out of 4
				createTestRecord("seq2", "ACGTAA", "$$$$$$"), // all 6 low quality
				createTestRecord("seq3", "ACGT", "&&&@"),     // 3 low quality out of 4
			},
			metric:    LQCount,
			ascending: false,
			headerMetrics: []HeaderMetric{
				{Name: "lqcount", IsLength: false},
				{Name: "length", IsLength: true},
			},
			minPhred:  30,
			minQual:   0.0,
			maxQual:   math.Inf(1),
			wantOrder: []string{"seq1", "seq3", "seq2"},
		},
		{
			name: "Natural sort on equal values",
			records: []*fastx.Record{
				createTestRecord("seq2", "ACGT", "IIII"),
				createTestRecord("seq10", "ACGT", "IIII"),
				createTestRecord("seq1", "ACGT", "IIII"),
			},
			metric:    AvgPhred,
			ascending: false,
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,
			maxQual:   math.Inf(1),
			wantOrder: []string{"seq1", "seq2", "seq10"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Create temporary input and output files
			inFile, err := os.CreateTemp("", "test_in_*.fastq")
			if err != nil {
				t.Fatal(err)
			}
			defer os.Remove(inFile.Name())

			outFile, err := os.CreateTemp("", "test_out_*.fastq")
			if err != nil {
				t.Fatal(err)
			}
			defer os.Remove(outFile.Name())

			// Write test records to input file
			writer, err := xopen.Wopen(inFile.Name())
			if err != nil {
				t.Fatal(err)
			}

			// Write records in FASTQ format
			for _, record := range tt.records {
				fmt.Fprintf(writer, "@%s\n%s\n+\n%s\n",
					record.Name,
					record.Seq.Seq,
					record.Seq.Qual)
			}
			writer.Close()

			// Run sortFile
			sortFile(
				inFile.Name(),
				outFile.Name(),
				tt.ascending,
				tt.metric,
				tt.compLevel,
				tt.headerMetrics,
				tt.minPhred,
				tt.minQual,
				tt.maxQual,
			)

			// Read and verify output
			reader, err := fastx.NewReader(seq.DNAredundant, outFile.Name(), fastx.DefaultIDRegexp)
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()

			var gotOrder []string
			for {
				record, err := reader.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatal(err)
				}
				name := strings.Split(string(record.Name), " ")[0]
				gotOrder = append(gotOrder, name)
			}

			// Only show file contents if there's a test failure
			if len(gotOrder) == 0 || !reflect.DeepEqual(gotOrder, tt.wantOrder) {
				// Read and log input file contents
				inContent, err := os.ReadFile(inFile.Name())
				if err != nil {
					t.Logf("Failed to read input file: %v", err)
				} else {
					t.Logf("Input file contents:\n%s", string(inContent))
				}

				// Read and log output file contents
				outContent, err := os.ReadFile(outFile.Name())
				if err != nil {
					t.Logf("Failed to read output file: %v", err)
				} else {
					t.Logf("Output file contents:\n%s", string(outContent))
				}
			}

			if len(gotOrder) == 0 {
				t.Error("No records were read from the output file")
			}

			if !reflect.DeepEqual(gotOrder, tt.wantOrder) {
				t.Errorf("sortFile() got order = %v, want %v", gotOrder, tt.wantOrder)
			}

			// If header metrics were specified, verify they were added correctly
			if len(tt.headerMetrics) > 0 {
				reader, _ = fastx.NewReader(seq.DNAredundant, outFile.Name(), fastx.DefaultIDRegexp)
				record, _ := reader.Read()
				header := string(record.Name)

				// Check that all requested metrics are present
				for _, metric := range tt.headerMetrics {
					if metric.IsLength {
						if !strings.Contains(header, "length=") {
							t.Errorf("Header missing length metric: %s", header)
						}
					} else {
						if !strings.Contains(header, metric.Name+"=") {
							t.Errorf("Header missing metric %s: %s", metric.Name, header)
						}
					}
				}
			}
		})
	}
}

// TestSortStdin tests the stdin-based sorting functionality
func TestSortStdin(t *testing.T) {
	tests := []struct {
		name          string
		records       []*fastx.Record
		metric        QualityMetric
		ascending     bool
		compLevel     int
		headerMetrics []HeaderMetric
		minPhred      int
		minQual       float64
		maxQual       float64
		wantOrder     []string
		wantErr       bool
	}{
		{
			name: "Basic sorting with compression",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // Phred 40
				createTestRecord("seq2", "ACGT", "$$$$"), // Phred 3
				createTestRecord("seq3", "ACGT", "@@@@"), // Phred 31
			},
			metric:    AvgPhred,
			ascending: false,
			compLevel: 1,
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,
			maxQual:   math.Inf(1),
			wantOrder: []string{"seq1", "seq3", "seq2"},
		},
		{
			name: "Sorting with quality filters",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // Phred 40
				createTestRecord("seq2", "ACGT", "$$$$"), // Phred 3
				createTestRecord("seq3", "ACGT", "@@@@"), // Phred 31
			},
			metric:    AvgPhred,
			ascending: false,
			compLevel: 0, // No compression
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   30.0,             // Only keep sequences with AvgPhred >= 30
			maxQual:   35.0,             // Only keep sequences with AvgPhred <= 35
			wantOrder: []string{"seq3"}, // Only seq3 falls within the quality range
		},
		{
			name: "Sorting with header metrics",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"),
				createTestRecord("seq2", "ACGTAA", "$$$$$$"),
			},
			metric:    AvgPhred,
			ascending: false,
			compLevel: 1,
			headerMetrics: []HeaderMetric{
				{Name: "avgphred", IsLength: false},
				{Name: "length", IsLength: true},
			},
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,
			maxQual:   math.Inf(1),
			wantOrder: []string{"seq1", "seq2"},
		},
		{
			name: "Natural sort with equal qualities",
			records: []*fastx.Record{
				createTestRecord("seq2", "ACGT", "IIII"),
				createTestRecord("seq10", "ACGT", "IIII"),
				createTestRecord("seq1", "ACGT", "IIII"),
			},
			metric:    AvgPhred,
			ascending: false,
			compLevel: 1,
			minPhred:  DEFAULT_MIN_PHRED,
			minQual:   0.0,
			maxQual:   math.Inf(1),
			wantOrder: []string{"seq1", "seq2", "seq10"},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Create temporary input and output files
			tmpInFile, err := os.CreateTemp("", "test_stdin_*.fastq")
			if err != nil {
				t.Fatal(err)
			}
			defer os.Remove(tmpInFile.Name())

			tmpOutFile, err := os.CreateTemp("", "test_stdout_*.fastq")
			if err != nil {
				t.Fatal(err)
			}
			defer os.Remove(tmpOutFile.Name())

			// Write test records to input file
			for _, record := range tt.records {
				fmt.Fprintf(tmpInFile, "@%s\n%s\n+\n%s\n",
					record.Name,
					record.Seq.Seq,
					record.Seq.Qual)
			}
			tmpInFile.Close()

			// Redirect stdin to read from the temp file
			oldStdin := os.Stdin
			newStdin, err := os.Open(tmpInFile.Name())
			if err != nil {
				t.Fatal(err)
			}
			os.Stdin = newStdin
			defer func() {
				os.Stdin = oldStdin
				newStdin.Close()
			}()

			// Run sortStdin
			sortStdin(
				tmpOutFile.Name(),
				tt.ascending,
				tt.metric,
				tt.compLevel,
				tt.headerMetrics,
				tt.minPhred,
				tt.minQual,
				tt.maxQual,
			)

			// Read and verify output
			reader, err := fastx.NewReader(seq.DNAredundant, tmpOutFile.Name(), fastx.DefaultIDRegexp)
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()

			var gotOrder []string
			for {
				record, err := reader.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatal(err)
				}
				name := strings.Split(string(record.Name), " ")[0] // Extract base name without metrics
				gotOrder = append(gotOrder, name)
			}

			// Verify the results
			if len(gotOrder) != len(tt.wantOrder) {
				t.Errorf("Got %d records, want %d records", len(gotOrder), len(tt.wantOrder))
			}

			if !reflect.DeepEqual(gotOrder, tt.wantOrder) {
				t.Errorf("sortStdin() got order = %v, want %v", gotOrder, tt.wantOrder)
			}

			// If header metrics were specified, verify they were added correctly
			if len(tt.headerMetrics) > 0 {
				reader, _ = fastx.NewReader(seq.DNAredundant, tmpOutFile.Name(), fastx.DefaultIDRegexp)
				record, _ := reader.Read()
				header := string(record.Name)

				// Check that all requested metrics are present
				for _, metric := range tt.headerMetrics {
					if metric.IsLength {
						if !strings.Contains(header, "length=") {
							t.Errorf("Header missing length metric: %s", header)
						}
					} else {
						if !strings.Contains(header, metric.Name+"=") {
							t.Errorf("Header missing metric %s: %s", metric.Name, header)
						}
					}
				}
			}
		})
	}
}

// TestRunNoSort verifies that runNoSort preserves input order while applying
// quality filters and optional header annotations
func TestRunNoSort(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "nosort_test_*")
	if err != nil {
		t.Fatal(err)
	}
	defer os.RemoveAll(tmpDir)

	createInput := func(records []*fastx.Record) string {
		path := filepath.Join(tmpDir, "input.fq")
		f, err := os.Create(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()

		for _, record := range records {
			fmt.Fprintf(f, "@%s\n%s\n+\n%s\n",
				record.Name,
				record.Seq.Seq,
				record.Seq.Qual)
		}
		return path
	}

	tests := []struct {
		name          string
		records       []*fastx.Record
		metric        QualityMetric
		minPhred      int
		minQual       float64
		maxQual       float64
		headerSpec    string
		wantOrder     []string
		wantFirstHead []string // substrings expected in first header (if any)
	}{
		{
			name: "Order preserved without filtering",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // high quality
				createTestRecord("seq2", "ACGT", "$$$$"), // low quality
				createTestRecord("seq3", "ACGT", "@@@@"), // medium quality
			},
			metric:   AvgPhred,
			minPhred: DEFAULT_MIN_PHRED,
			// No filtering on metric
			minQual:   -math.MaxFloat64,
			maxQual:   math.MaxFloat64,
			headerSpec: "",
			wantOrder: []string{"seq1", "seq2", "seq3"},
		},
		{
			name: "Filtering and header metrics",
			records: []*fastx.Record{
				createTestRecord("seq1", "ACGT", "IIII"), // very high quality
				createTestRecord("seq2", "ACGT", "$$$$"), // very low quality
				createTestRecord("seq3", "ACGT", "@@@@"), // medium quality
			},
			metric:   AvgPhred,
			minPhred: DEFAULT_MIN_PHRED,
			// Reuse thresholds from TestSortStdin ("Sorting with quality filters"):
			// only the medium-quality record should remain.
			minQual:    30.0,
			maxQual:    35.0,
			headerSpec: "avgphred,length",
			wantOrder:  []string{"seq3"},
			wantFirstHead: []string{
				"avgphred=",
				"length=",
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			inPath := createInput(tt.records)
			outPath := filepath.Join(tmpDir, "output.fq")

			var headerMetrics []HeaderMetric
			if tt.headerSpec != "" {
				var err error
				headerMetrics, err = parseHeaderMetrics(tt.headerSpec)
				if err != nil {
					t.Fatalf("parseHeaderMetrics(%q) error: %v", tt.headerSpec, err)
				}
			}

			err := runNoSort(
				inPath,
				outPath,
				tt.metric,
				headerMetrics,
				tt.minPhred,
				tt.minQual,
				tt.maxQual,
			)
			if err != nil {
				t.Fatalf("runNoSort() error: %v", err)
			}

			// Read and verify output
			reader, err := fastx.NewReader(seq.DNAredundant, outPath, fastx.DefaultIDRegexp)
			if err != nil {
				t.Fatal(err)
			}
			defer reader.Close()

			var gotOrder []string
			var firstHeader string

			for {
				record, err := reader.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatal(err)
				}
				if firstHeader == "" {
					firstHeader = string(record.Name)
				}
				name := strings.Split(string(record.Name), " ")[0]
				gotOrder = append(gotOrder, name)
			}

			if !reflect.DeepEqual(gotOrder, tt.wantOrder) {
				t.Errorf("runNoSort() got order = %v, want %v", gotOrder, tt.wantOrder)
			}

			// If we expect specific annotations, verify them on the first header
			if len(tt.wantFirstHead) > 0 {
				if firstHeader == "" {
					t.Fatalf("expected at least one record in output to inspect header")
				}
				for _, substr := range tt.wantFirstHead {
					if !strings.Contains(firstHeader, substr) {
						t.Errorf("header %q missing expected substring %q", firstHeader, substr)
					}
				}
			}
		})
	}
}

// TestRunDefaultCommand_Stdin verifies that runDefaultCommand calls sortStdin
// when the input is set to "-" (stdin) and that sorting succeeds
func TestRunDefaultCommand_Stdin(t *testing.T) {
	// Prepare temporary stdin FASTQ file
	tmpInFile, err := os.CreateTemp("", "default_cmd_stdin_*.fastq")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpInFile.Name())

	// Two simple records with different qualities so we can verify ordering
	fmt.Fprintf(tmpInFile, "@seq1\nACGT\n+\nIIII\n") // high quality
	fmt.Fprintf(tmpInFile, "@seq2\nACGT\n+\n$$$$\n") // low quality
	tmpInFile.Close()

	// Redirect stdin to the temp file
	oldStdin := os.Stdin
	newStdin, err := os.Open(tmpInFile.Name())
	if err != nil {
		t.Fatal(err)
	}
	os.Stdin = newStdin
	defer func() {
		os.Stdin = oldStdin
		newStdin.Close()
	}()

	// Prepare temporary output file
	tmpOutFile, err := os.CreateTemp("", "default_cmd_stdout_*.fastq")
	if err != nil {
		t.Fatal(err)
	}
	tmpOutFile.Close()
	defer os.Remove(tmpOutFile.Name())

	// Set global flags used by runDefaultCommand
	inFile = "-"
	outFile = tmpOutFile.Name()
	metric = "avgphred"
	minPhred = DEFAULT_MIN_PHRED
	minQualFilter = -math.MaxFloat64
	maxQualFilter = math.MaxFloat64
	headerMetrics = ""
	ascending = false
	compLevel = 0
	version = false

	// Call runDefaultCommand; on success it should not call exitFunc.
	// exitFunc is already mocked to panic in TestMain, but since we
	// expect no error, we can call directly without a recover wrapper
	runDefaultCommand(nil, nil)

	// Verify output order (seq1 should come before seq2)
	reader, err := fastx.NewReader(seq.DNAredundant, tmpOutFile.Name(), fastx.DefaultIDRegexp)
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	var gotOrder []string
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		name := strings.Split(string(record.Name), " ")[0]
		gotOrder = append(gotOrder, name)
	}

	wantOrder := []string{"seq1", "seq2"}
	if !reflect.DeepEqual(gotOrder, wantOrder) {
		t.Errorf("runDefaultCommand(stdin) produced order %v, want %v", gotOrder, wantOrder)
	}
}

// TestSortStdin_ReadError covers the error path when reading from stdin fails
// inside sortStdin (i.e. the "Error reading record" branch)
func TestSortStdin_ReadError(t *testing.T) {
	// Capture stderr and recover from the expected panic triggered via exitFunc
	oldStderr := os.Stderr
	rErr, wErr, _ := os.Pipe()
	os.Stderr = wErr

	defer func() {
		wErr.Close()
		os.Stderr = oldStderr
	}()

	// Prepare a valid output path (we are testing read errors, not write errors)
	tmpDir, err := os.MkdirTemp("", "sortstdin_readerr_*")
	if err != nil {
		t.Fatal(err)
	}
	defer os.RemoveAll(tmpDir)
	outPath := filepath.Join(tmpDir, "out.fastq")

	// Redirect stdin to a malformed FASTQ file to trigger a read error
	tmpInFile, err := os.CreateTemp("", "sortstdin_err_in_*.fastq")
	if err != nil {
		t.Fatal(err)
	}
	// Write clearly invalid FASTQ content (not following 4-line record structure)
	fmt.Fprintln(tmpInFile, "this_is_not_a_valid_fastq_record")
	tmpInFile.Close()
	defer os.Remove(tmpInFile.Name())

	oldStdin := os.Stdin
	newStdin, err := os.Open(tmpInFile.Name())
	if err != nil {
		t.Fatal(err)
	}
	os.Stdin = newStdin
	defer func() {
		os.Stdin = oldStdin
		newStdin.Close()
	}()

	// Call sortStdin and expect it to call exitFunc(1), which panics
	didPanic := false
	func() {
		defer func() {
			if r := recover(); r != nil {
				// exitFunc is configured to panic with "exit <code>"
				if s, ok := r.(string); ok && strings.HasPrefix(s, "exit ") {
					didPanic = true
				} else {
					t.Fatalf("Unexpected panic: %v", r)
				}
			}
		}()

		sortStdin(outPath, false, AvgPhred, 0, nil, DEFAULT_MIN_PHRED, -math.MaxFloat64, math.MaxFloat64)
	}()

	// Read captured stderr
	wErr.Close()
	errOutput, _ := io.ReadAll(rErr)

	if !didPanic {
		t.Fatalf("Expected sortStdin to call exitFunc and panic, but it did not")
	}
	if !strings.Contains(string(errOutput), "Error reading record") {
		t.Errorf("Expected stderr to contain 'Error reading record', got: %s", string(errOutput))
	}
}
// TestMainCommand tests the main command functionality
func TestMainCommand(t *testing.T) {
	// Create temporary directory for test files
	tmpDir, err := os.MkdirTemp("", "phredsort_test_*")
	if err != nil {
		t.Fatal(err)
	}
	defer os.RemoveAll(tmpDir)

	// Helper function to create test FASTQ file
	createTestFastq := func(name string) string {
		path := filepath.Join(tmpDir, name)
		f, err := os.Create(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()

		// Write some test FASTQ data
		fmt.Fprintf(f, "@seq1\nACGT\n+\nIIII\n")
		fmt.Fprintf(f, "@seq2\nACGT\n+\n$$$$\n")
		return path
	}

	// Helper function to capture stdout/stderr
	captureOutput := func(f func()) (string, string) {
		oldStdout := os.Stdout
		oldStderr := os.Stderr
		rOut, wOut, _ := os.Pipe()
		rErr, wErr, _ := os.Pipe()
		os.Stdout = wOut
		os.Stderr = wErr

		f()

		wOut.Close()
		wErr.Close()
		os.Stdout = oldStdout
		os.Stderr = oldStderr

		stdout, _ := io.ReadAll(rOut)
		stderr, _ := io.ReadAll(rErr)
		return string(stdout), string(stderr)
	}

	tests := []struct {
		name          string
		args          []string
		expectedCode  int
		checkStdout   bool
		checkStderr   bool
		wantStdout    string
		wantStderr    string
		setupFiles    bool
		validateFiles bool
	}{
		{
			name:         "Version flag",
			args:         []string{"--version"},
			expectedCode: 0,
			checkStdout:  true,
			wantStdout:   fmt.Sprintf("phredsort %s\n", VERSION),
		},
		{
			name:         "Missing required flags",
			args:         []string{},
			expectedCode: 1,
			checkStderr:  true,
			wantStderr: red("Error: input and output files are required") + "\n" +
				red("Try 'phredsort --help' for more information") + "\n",
		},
		{
			name:         "Invalid metric",
			args:         []string{"--in", "input.fq", "--out", "output.fq", "--metric", "invalid"},
			expectedCode: 1,
			checkStderr:  true,
			wantStderr:   red("Error: invalid metric 'invalid'. Must be one of: avgphred, maxee, meep, lqcount, lqpercent"),
		},
		{
			name:         "Invalid compression level",
			args:         []string{"--in", "input.fq", "--out", "output.fq", "--compress", "23"},
			expectedCode: 1,
			checkStderr:  true,
			wantStderr:   red("Error: compression level must be between 0 and 22") + "\n",
		},
		{
			name:         "Invalid header metrics",
			args:         []string{"--in", "input.fq", "--out", "output.fq", "--header", "invalid,metrics"},
			expectedCode: 1,
			checkStderr:  true,
			wantStderr:   red("Error: invalid header metric: invalid") + "\n",
		},
		{
			name:          "Basic file processing",
			args:          []string{"--in", "input.fq", "--out", "output.fq"},
			expectedCode:  0,
			setupFiles:    true,
			validateFiles: true,
		},
		{
			name: "Complex command with multiple options",
			args: []string{
				"--in", "input.fq",
				"--out", "output.fq",
				"--metric", "avgphred",
				"--minphred", "20",
				"--minqual", "15",
				"--maxqual", "40",
				"--header", "avgphred,maxee,length",
				"--ascending",
				"--compress", "1",
			},
			expectedCode:  0,
			setupFiles:    true,
			validateFiles: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if tt.setupFiles {
				inFile := createTestFastq("input.fq")
				outFile := filepath.Join(tmpDir, "output.fq")

				// Update args with actual file paths
				for i, arg := range tt.args {
					switch arg {
					case "input.fq":
						tt.args[i] = inFile
					case "output.fq":
						tt.args[i] = outFile
					}
				}
			}

			// Reset os.Args and set test arguments
			oldArgs := os.Args
			os.Args = append([]string{"phredsort"}, tt.args...)

			// Capture exit code
			var exitCode int
			oldExit := exitFunc
			exitFunc = func(code int) {
				exitCode = code
				panic(fmt.Sprintf("exit %d", code))
			}
			defer func() {
				exitFunc = oldExit
				os.Args = oldArgs

				if r := recover(); r != nil {
					if exitStr, ok := r.(string); !ok || !strings.HasPrefix(exitStr, "exit ") {
						t.Errorf("Unexpected panic: %v", r)
					}
				}
			}()

			// Run main and capture output
			stdout, stderr := captureOutput(func() {
				defer func() {
					if r := recover(); r != nil {
						if exitStr, ok := r.(string); !ok || !strings.HasPrefix(exitStr, "exit ") {
							panic(r)
						}
					}
				}()
				main()
			})

			// Check exit code
			if exitCode != tt.expectedCode {
				t.Errorf("Expected exit code %d, got %d", tt.expectedCode, exitCode)
			}

			// Check stdout if required
			if tt.checkStdout && stdout != tt.wantStdout {
				t.Errorf("Expected stdout:\n%s\nGot:\n%s", tt.wantStdout, stdout)
			}

			// Check stderr if required (including non-printable characters)
			if tt.checkStderr && stderr != tt.wantStderr {
				t.Errorf("Expected stderr:\n%q\nGot:\n%q", tt.wantStderr, stderr)
			}

			// Validate output file if required
			if tt.validateFiles {
				outPath := filepath.Join(tmpDir, "output.fq")
				if _, err := os.Stat(outPath); os.IsNotExist(err) {
					t.Error("Expected output file was not created")
				}

				// Optional: Add more specific file content validation here
				// For example, check if the file is valid FASTQ format
				reader, err := fastx.NewReader(seq.DNAredundant, outPath, fastx.DefaultIDRegexp)
				if err != nil {
					t.Errorf("Failed to read output file: %v", err)
				}
				defer reader.Close()

				// Try to read at least one record
				_, err = reader.Read()
				if err != nil && err != io.EOF {
					t.Errorf("Failed to read record from output file: %v", err)
				}
			}
		})
	}
}

func TestMain(m *testing.M) {
	// Store original exit function
	originalExit := exitFunc
	// Replace exit function with mock during tests
	exitFunc = func(code int) {
		panic(fmt.Sprintf("exit %d", code))
	}
	defer func() { exitFunc = originalExit }()

	m.Run()
}

// TestParsePreSortRecord verifies parsing of pre-computed header metrics for headersort mode
func TestParsePreSortRecord(t *testing.T) {
	tests := []struct {
		name        string
		header      string
		metric      QualityMetric
		wantID      string
		wantQual    float64
		wantSize    int
		wantHasQual bool
		wantHasSize bool
	}{
		{
			name:        "Space-separated maxee with size",
			header:      ">seq1 maxee=2.5 size=100",
			metric:      MaxEE,
			wantID:      "seq1",
			wantQual:    2.5,
			wantSize:    100,
			wantHasQual: true,
			wantHasSize: true,
		},
		{
			name:        "Semicolon-separated maxee with size",
			header:      ">seq2;maxee=0.5;size=200",
			metric:      MaxEE,
			// For semicolon-separated headers without spaces, the entire header
			// (minus the leading '>') is treated as the ID
			wantID:      "seq2;maxee=0.5;size=200",
			wantQual:    0.5,
			wantSize:    200,
			wantHasQual: true,
			wantHasSize: true,
		},
		{
			name:        "Missing metric but with size",
			header:      ">seq3 size=300",
			metric:      MaxEE,
			wantID:      "seq3",
			wantQual:    0.0,
			wantSize:    300,
			wantHasQual: false,
			wantHasSize: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			record := &fastx.Record{
				Name: []byte(tt.header),
			}

			got, err := parsePreSortRecord(record, tt.metric)
			if err != nil {
				t.Fatalf("parsePreSortRecord() unexpected error: %v", err)
			}

			if got.ID != tt.wantID {
				t.Errorf("ID = %q, want %q", got.ID, tt.wantID)
			}
			if got.HasQual != tt.wantHasQual {
				t.Errorf("HasQual = %v, want %v", got.HasQual, tt.wantHasQual)
			}
			if got.HasSize != tt.wantHasSize {
				t.Errorf("HasSize = %v, want %v", got.HasSize, tt.wantHasSize)
			}
			if tt.wantHasQual && math.Abs(got.Quality-tt.wantQual) > 1e-6 {
				t.Errorf("Quality = %v, want %v", got.Quality, tt.wantQual)
			}
			if tt.wantHasSize && got.Size != tt.wantSize {
				t.Errorf("Size = %d, want %d", got.Size, tt.wantSize)
			}
		})
	}
}

// TestRunPresort verifies headersort mode (runPresort) using pre-computed metrics in headers
func TestRunPresort(t *testing.T) {
	tmpDir, err := os.MkdirTemp("", "headersort_test_*")
	if err != nil {
		t.Fatal(err)
	}
	defer os.RemoveAll(tmpDir)

	createInput := func(name, content string) string {
		path := filepath.Join(tmpDir, name)
		if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
			t.Fatalf("failed to write test input: %v", err)
		}
		return path
	}

	tests := []struct {
		name        string
		content     string
		metric      QualityMetric
		ascending   bool
		minQual     float64
		maxQual     float64
		wantOrder   []string
		wantErr     bool
		errContains string
	}{
		{
			name: "MaxEE default ordering (best to worst)",
			content: "" +
				">seq1 maxee=2.0 size=100\nACGT\n" +
				">seq2 maxee=0.5 size=200\nACGT\n",
			metric:    MaxEE,
			ascending: false,
			minQual:   -math.MaxFloat64,
			maxQual:   math.MaxFloat64,
			// For MaxEE, lower is better; default ordering should be best-to-worst
			wantOrder: []string{"seq2", "seq1"},
		},
		{
			name: "MaxEE ascending (worst to best)",
			content: "" +
				">seq1 maxee=2.0\nACGT\n" +
				">seq2 maxee=0.5\nACGT\n",
			metric:    MaxEE,
			ascending: true,
			minQual:   -math.MaxFloat64,
			maxQual:   math.MaxFloat64,
			// Ascending flips the default ordering.
			wantOrder: []string{"seq1", "seq2"},
		},
		{
			name: "Quality filtering",
			content: "" +
				">seq1 maxee=0.5\nACGT\n" +
				">seq2 maxee=2.0\nACGT\n" +
				">seq3 maxee=5.0\nACGT\n",
			metric:    MaxEE,
			ascending: false,
			minQual:   1.0,
			maxQual:   4.0,
			// Only seq2 (maxee=2.0) falls within [1.0, 4.0]
			wantOrder: []string{"seq2"},
		},
		{
			name: "Missing required metric produces error",
			content: "" +
				">seq1 size=100\nACGT\n",
			metric:      MaxEE,
			ascending:   false,
			minQual:     -math.MaxFloat64,
			maxQual:     math.MaxFloat64,
			wantErr:     true,
			errContains: "record missing required quality metric",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			inPath := createInput("input.fasta", tt.content)
			outPath := filepath.Join(tmpDir, "output.fasta")

			err := runPresort(inPath, outPath, tt.metric, tt.ascending, tt.minQual, tt.maxQual)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("runPresort() expected error, got nil")
				}
				if tt.errContains != "" && !strings.Contains(err.Error(), tt.errContains) {
					t.Fatalf("runPresort() error = %q, want to contain %q", err.Error(), tt.errContains)
				}
				return
			}
			if err != nil {
				t.Fatalf("runPresort() unexpected error: %v", err)
			}

			// Read and verify output order
			reader, err := fastx.NewDefaultReader(outPath)
			if err != nil {
				t.Fatalf("failed to create reader: %v", err)
			}
			defer reader.Close()

			var gotOrder []string
			for {
				record, err := reader.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("failed to read record: %v", err)
				}
				// ID is the first token of the header
				name := strings.Split(string(record.Name), " ")[0]
				gotOrder = append(gotOrder, name)
			}

			if !reflect.DeepEqual(gotOrder, tt.wantOrder) {
				t.Errorf("runPresort() got order = %v, want %v", gotOrder, tt.wantOrder)
			}
		})
	}
}

// captureStdout is a helper that captures standard output produced by f
// (used for testing the help function)
func captureStdout(f func()) string {
	oldStdout := os.Stdout
	rOut, wOut, _ := os.Pipe()
	os.Stdout = wOut

	f()

	wOut.Close()
	os.Stdout = oldStdout

	outBytes, _ := io.ReadAll(rOut)
	return string(outBytes)
}

// TestHelpFuncRoot verifies that the custom help function prints the root help
// information, including version and subcommand details
func TestHelpFuncRoot(t *testing.T) {
	cmd := &cobra.Command{
		Use: "phredsort",
	}

	output := captureStdout(func() {
		helpFunc(cmd, nil)
	})

	if !strings.Contains(output, "phredsort v."+VERSION) {
		t.Errorf("root help output missing version string, got:\n%s", output)
	}
	if !strings.Contains(output, "Subcommands:") {
		t.Errorf("root help output missing 'Subcommands:' section, got:\n%s", output)
	}
	if !strings.Contains(output, "Quality metrics:") {
		t.Errorf("root help output missing 'Quality metrics:' section, got:\n%s", output)
	}
}

// TestHelpFuncSort verifies that the custom help function prints detailed help for the `sort` subcommand
func TestHelpFuncSort(t *testing.T) {
	cmd := &cobra.Command{
		Use: "sort",
	}

	output := captureStdout(func() {
		helpFunc(cmd, nil)
	})

	if !strings.Contains(output, "phredsort sort - Sorts FASTQ based on computed quality metrics") {
		t.Errorf("sort help output missing sort description, got:\n%s", output)
	}
	if !strings.Contains(output, "Input FASTQ file (required, use '-' for stdin)") {
		t.Errorf("sort help output missing input flag description, got:\n%s", output)
	}
	if !strings.Contains(output, "Quality metric (avgphred, maxee, meep, lqcount, lqpercent)") {
		t.Errorf("sort help output missing metric flag description, got:\n%s", output)
	}
}

// TestHelpFuncHeaderSort verifies that the custom help function prints detailed
// help for the `headersort` subcommand, including supported header formats
func TestHelpFuncHeaderSort(t *testing.T) {
	cmd := &cobra.Command{
		Use: "headersort",
	}

	output := captureStdout(func() {
		helpFunc(cmd, nil)
	})

	if !strings.Contains(output, "phredsort headersort - Sorts sequences using header quality metrics") {
		t.Errorf("headersort help output missing headersort description, got:\n%s", output)
	}
	if !strings.Contains(output, "Header metric to use (avgphred, maxee, meep, lqcount, lqpercent)") {
		t.Errorf("headersort help output missing metric flag description, got:\n%s", output)
	}
	if !strings.Contains(output, `">seq1 maxee=2.5 size=100"`) {
		t.Errorf("headersort help output missing space-separated header format example, got:\n%s", output)
	}
	if !strings.Contains(output, `">seq1;maxee=2.5;size=100"`) {
		t.Errorf("headersort help output missing semicolon-separated header format example, got:\n%s", output)
	}
}

