package main

import (
	"fmt"
	"math"
	"os"
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
