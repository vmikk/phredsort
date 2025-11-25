package main

import (
	"os"
	"reflect"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

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

