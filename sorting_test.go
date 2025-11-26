package main

import (
	"reflect"
	"sort"
	"testing"
)

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

