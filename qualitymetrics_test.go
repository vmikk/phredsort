package main

import (
	"fmt"
	"math"
	"testing"
)

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

// TestErrorProbabilitiesInit tests a few key values from the pre-computed errorProbs array
func TestErrorProbabilitiesInit(t *testing.T) {
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

