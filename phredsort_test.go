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
			want:     0.0001, // Very low error probability
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
			if math.Abs(got-tt.want) > 0.1 {
				t.Errorf("calculateQuality() = %v, want %v", got, tt.want)
			}
		})
	}
}

