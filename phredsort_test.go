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

// Test quality metric calculations
func TestQualityMetricCalculations(t *testing.T) {
	tests := []struct {
		name     string
		qual     []byte
		metric   QualityMetric
		minPhred int
		want     float64
	}{
		{
			name:     "AvgPhred - All high quality",
			qual:     []byte("IIIII"), // ASCII 73 = Phred 40
			metric:   AvgPhred,
			minPhred: DEFAULT_MIN_PHRED,
			want:     40.0,
		},
		{
			name:     "AvgPhred - Mixed quality",
			qual:     []byte("I$$I$"), // Mix of Phred 40 and 3
			metric:   AvgPhred,
			minPhred: DEFAULT_MIN_PHRED,
			want:     5.21791,
		},
		{
			name:     "AvgPhred - Empty quality",
			qual:     []byte{},
			metric:   AvgPhred,
			minPhred: DEFAULT_MIN_PHRED,
			want:     0.0,
		},
		{
			name:     "MaxEE - All high quality",
			qual:     []byte("IIIII"),
			metric:   MaxEE,
			minPhred: DEFAULT_MIN_PHRED,
			want:     0.0005,
		},
		{
			name:     "MaxEE - Empty quality",
			qual:     []byte{},
			metric:   MaxEE,
			minPhred: DEFAULT_MIN_PHRED,
			want:     math.Inf(1),
		},
		{
			name:     "LQCount - No low quality bases",
			qual:     []byte("IIIII"),
			metric:   LQCount,
			minPhred: 30,
			want:     0,
		},
		{
			name:     "LQCount - All low quality bases",
			qual:     []byte("$$$$$"),
			metric:   LQCount,
			minPhred: 30,
			want:     5,
		},
		{
			name:     "LQPercent - Half low quality bases",
			qual:     []byte("II$$$"),
			metric:   LQPercent,
			minPhred: 30,
			want:     60.0,
		},
		{
			name:     "MaxEE - Single very low quality base",
			qual:     []byte("$"), // ASCII 36 = Phred 3
			metric:   MaxEE,
			minPhred: DEFAULT_MIN_PHRED,
			want:     0.5011872336272722, // Error prob for Phred 3
		},
		{
			name:     "Meep - Mixed quality",
			qual:     []byte("I$$I$"), // Mix of Phred 40 and 3
			metric:   Meep,
			minPhred: DEFAULT_MIN_PHRED,
			want:     30.07523, // (sum of error probs * 100) / len
		},
		{
			name:     "LQCount - Custom minPhred threshold",
			qual:     []byte("BBBBB"), // ASCII 66 = Phred 33
			metric:   LQCount,
			minPhred: 35,
			want:     5, // All bases below Phred 35
		},
		{
			name:     "LQPercent - Single base at threshold",
			qual:     []byte("0"), // ASCII 48 = Phred 15
			metric:   LQPercent,
			minPhred: 15,
			want:     0.0, // Base exactly at threshold
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			record := createTestRecord("test", "ACGT", string(tt.qual))
			got := calculateQuality(record, tt.metric, tt.minPhred)

			// Use approximate comparison for floating point values
			if math.Abs(got-tt.want) > 0.00001 {
				t.Errorf("calculateQuality() = %v, want %v", got, tt.want)
			}
		})
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

// Add test for error probability initialization
func TestErrorProbabilitiesInit(t *testing.T) {
	// Test a few key values from the pre-computed errorProbs array
	tests := []struct {
		phred byte
		want  float64
	}{
		{33, 1},        // Phred 0  (33 - 33 = 0)
		{43, 0.1},      // Phred 10 (43 - 33 = 10)
		{53, 0.01},     // Phred 20 (53 - 33 = 20)
		{63, 0.001},    // Phred 30 (63 - 33 = 30)
		{73, 0.0001},   // Phred 40 (73 - 33 = 40)
		{83, 0.00001},  // Phred 50 (83 - 33 = 50)
		{93, 0.000001}, // Phred 60 (93 - 33 = 60)

	}

	for _, tt := range tests {
		t.Run(fmt.Sprintf("Phred%d", tt.phred-PHRED_OFFSET), func(t *testing.T) {
			if got := errorProbs[tt.phred]; math.Abs(got-tt.want) > 1e-10 {
				t.Errorf("errorProbs[%d] = %v, want %v", tt.phred, got, tt.want)
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
