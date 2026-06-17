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

func TestQualityIndexListUsesFloat64Precision(t *testing.T) {
	names := []string{"seq1_worse", "seq2_better"}
	items := []QualityIndex{
		{Index: 0, Value: 1000000.0001},
		{Index: 1, Value: 1000000.0},
	}

	list := NewQualityIndexList(items, names, false, MaxEE)
	sort.Sort(list)

	got := []string{
		names[list.items[0].Index],
		names[list.items[1].Index],
	}
	want := []string{"seq2_better", "seq1_worse"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("QualityIndexList order = %v, want %v", got, want)
	}
}

func TestHeaderSortIndexListUsesFloat64Precision(t *testing.T) {
	ids := []string{"seq1_worse", "seq2_better"}
	items := []HeaderSortIndex{
		{Index: 0, Quality: 1000000.0001},
		{Index: 1, Quality: 1000000.0},
	}

	list := NewHeaderSortIndexList(items, ids, false, MaxEE)
	sort.Sort(list)

	got := []string{
		ids[list.items[0].Index],
		ids[list.items[1].Index],
	}
	want := []string{"seq2_better", "seq1_worse"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("HeaderSortIndexList order = %v, want %v", got, want)
	}
}

func TestAscendingTieBreakKeepsNaturalOrder(t *testing.T) {
	names := []string{"seq2", "seq10", "seq1"}
	items := []QualityIndex{
		{Index: 0, Value: 40},
		{Index: 1, Value: 40},
		{Index: 2, Value: 40},
	}

	list := NewQualityIndexList(items, names, true, AvgPhred)
	sort.Sort(list)

	got := make([]string, len(list.items))
	for i, item := range list.items {
		got[i] = names[item.Index]
	}
	want := []string{"seq1", "seq2", "seq10"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("ascending equal-quality order = %v, want %v", got, want)
	}
}

func TestDuplicateNamesAreNotLessThanEachOther(t *testing.T) {
	names := []string{"dup", "dup"}
	list := NewQualityIndexList([]QualityIndex{
		{Index: 0, Value: 40},
		{Index: 1, Value: 40},
	}, names, false, AvgPhred)

	if list.Less(0, 1) || list.Less(1, 0) {
		t.Fatalf("duplicate equal-quality names must compare equal")
	}

	floatList := NewQualityFloatList([]QualityFloat{
		{Name: "dup", Value: 40, Metric: AvgPhred},
		{Name: "dup", Value: 40, Metric: AvgPhred},
	}, false)
	if floatList.Less(0, 1) || floatList.Less(1, 0) {
		t.Fatalf("duplicate equal-quality names must compare equal in QualityFloatList")
	}
}

func TestQualityFloatListMatchesQualityIndexList(t *testing.T) {
	floatItems := []QualityFloat{
		{Name: "seq2", Value: 0.5, Metric: MaxEE},
		{Name: "seq10", Value: 0.5, Metric: MaxEE},
		{Name: "seq1", Value: 2.0, Metric: MaxEE},
	}
	indexNames := []string{"seq2", "seq10", "seq1"}
	indexItems := []QualityIndex{
		{Index: 0, Value: 0.5},
		{Index: 1, Value: 0.5},
		{Index: 2, Value: 2.0},
	}

	floatList := NewQualityFloatList(floatItems, true)
	indexList := NewQualityIndexList(indexItems, indexNames, true, MaxEE)
	sort.Sort(floatList)
	sort.Sort(indexList)

	gotFloat := make([]string, len(floatList.items))
	for i, item := range floatList.items {
		gotFloat[i] = item.Name
	}
	gotIndex := make([]string, len(indexList.items))
	for i, item := range indexList.items {
		gotIndex[i] = indexNames[item.Index]
	}

	if !reflect.DeepEqual(gotFloat, gotIndex) {
		t.Fatalf("QualityFloatList order = %v, QualityIndexList order = %v", gotFloat, gotIndex)
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
