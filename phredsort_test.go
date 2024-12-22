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

