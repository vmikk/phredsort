package main

import (
	"math"

	"github.com/shenwei356/bio/seqio/fastx"
)

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

// Count the number of low quality bases
func countLowQualityBases(qual []byte, minPhred int) float64 {
	if len(qual) == 0 {
		return math.Inf(1)
	}

	count := 0
	for _, q := range qual {
		if int(q)-PHRED_OFFSET < minPhred {
			count++
		}
	}
	return float64(count)
}

// Calculate the percentage of low quality bases
func calculateLQPercent(qual []byte, minPhred int) float64 {
	if len(qual) == 0 {
		return math.Inf(1)
	}

	count := 0
	for _, q := range qual {
		if int(q)-PHRED_OFFSET < minPhred {
			count++
		}
	}
	return (float64(count) * 100) / float64(len(qual))
}

// A common type for quality calculator functions
type QualityCalculator func([]byte, int) float64

// Wrapper functions to standardize the interface
// (to have the same signature `func([]byte, int) float64`)
func avgPhredWrapper(qual []byte, _ int) float64 {
	return calculateAvgPhred(qual)
}

func maxEEWrapper(qual []byte, _ int) float64 {
	return calculateMaxEE(qual)
}

func meepWrapper(qual []byte, _ int) float64 {
	return calculateMeep(qual)
}

func lqCountWrapper(qual []byte, minPhred int) float64 {
	return countLowQualityBases(qual, minPhred)
}

func lqPercentWrapper(qual []byte, minPhred int) float64 {
	return calculateLQPercent(qual, minPhred)
}

// Map of metric types to their calculator functions
var qualityCalculators = map[QualityMetric]QualityCalculator{
	AvgPhred:  avgPhredWrapper,
	MaxEE:     maxEEWrapper,
	Meep:      meepWrapper,
	LQCount:   lqCountWrapper,
	LQPercent: lqPercentWrapper,
}

// Calculate sequence quality
func calculateQuality(record *fastx.Record, metric QualityMetric, minPhred int) float64 {
	if calcFunc, exists := qualityCalculators[metric]; exists {
		return calcFunc(record.Seq.Qual, minPhred)
	}
	return 0
}
