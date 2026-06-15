package main

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

func writeFastqRecords(t *testing.T, path string, records []*fastx.Record) {
	t.Helper()

	fh, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()

	for _, record := range records {
		fmt.Fprintf(fh, "@%s\n%s\n+\n%s\n", record.Name, record.Seq.Seq, record.Seq.Qual)
	}
}

func readFastxIDs(t *testing.T, path string) []string {
	t.Helper()

	reader, err := fastx.NewReader(seq.DNAredundant, path, fastx.DefaultIDRegexp)
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	var ids []string
	for {
		record, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		ids = append(ids, strings.Split(string(record.Name), " ")[0])
	}
	return ids
}

func expectExitWithFastqError(t *testing.T, fn func()) {
	t.Helper()

	oldStderr := os.Stderr
	rErr, wErr, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stderr = wErr
	defer func() {
		os.Stderr = oldStderr
		rErr.Close()
	}()

	didPanic := false
	func() {
		defer func() {
			if r := recover(); r != nil {
				if s, ok := r.(string); ok && strings.HasPrefix(s, "exit ") {
					didPanic = true
					return
				}
				t.Fatalf("unexpected panic: %v", r)
			}
		}()
		fn()
	}()

	wErr.Close()
	errOutput, _ := io.ReadAll(rErr)
	if !didPanic {
		t.Fatalf("expected exitFunc panic")
	}
	if !strings.Contains(string(errOutput), computedQualityFastqError) {
		t.Fatalf("stderr = %q, want FASTQ error %q", string(errOutput), computedQualityFastqError)
	}
}

func TestChunkedStorageBoundaries(t *testing.T) {
	storage := NewChunkedStorage(0, 4)
	payloads := [][]byte{
		[]byte("abcd"),  // exact chunk fill
		[]byte("ef"),    // chunk rollover
		[]byte("ghijk"), // larger than chunk size
		{},              // explicit zero-length record support
		[]byte("z"),
	}

	indices := make([]int, 0, len(payloads))
	for _, payload := range payloads {
		indices = append(indices, storage.Append(payload))
	}

	for i, idx := range indices {
		got := storage.Get(idx)
		if !bytes.Equal(got, payloads[i]) {
			t.Fatalf("storage.Get(%d) = %q, want %q", idx, got, payloads[i])
		}
	}
	if storage.Len() != len(payloads) {
		t.Fatalf("storage.Len() = %d, want %d", storage.Len(), len(payloads))
	}
}

func TestSortRecordsCompressedMatchesUncompressed(t *testing.T) {
	tmpDir := t.TempDir()
	inputPath := filepath.Join(tmpDir, "input.fastq")
	outPlain := filepath.Join(tmpDir, "plain.fastq")
	outCompressed := filepath.Join(tmpDir, "compressed.fastq")

	records := []*fastx.Record{
		createTestRecord("read10", "GATTACA", "5555555"),
		createTestRecord("read2", "TT", "!!"),
		createTestRecord("read3", "GGG", "@@@"),
		createTestRecord("read1", "ACGTAC", "IIIIII"),
	}
	writeFastqRecords(t, inputPath, records)

	sortRecords(inputPath, outPlain, false, AvgPhred, 0, nil, DEFAULT_MIN_PHRED, -math.MaxFloat64, math.MaxFloat64)
	sortRecords(inputPath, outCompressed, false, AvgPhred, 1, nil, DEFAULT_MIN_PHRED, -math.MaxFloat64, math.MaxFloat64)

	plainBytes, err := os.ReadFile(outPlain)
	if err != nil {
		t.Fatal(err)
	}
	compressedBytes, err := os.ReadFile(outCompressed)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(plainBytes, compressedBytes) {
		t.Fatalf("compressed output differs from uncompressed\nplain:\n%s\ncompressed:\n%s", plainBytes, compressedBytes)
	}

	gotIDs := readFastxIDs(t, outCompressed)
	wantIDs := []string{"read1", "read3", "read10", "read2"}
	if !reflect.DeepEqual(gotIDs, wantIDs) {
		t.Fatalf("sorted IDs = %v, want %v", gotIDs, wantIDs)
	}
}

func TestSortRecordsRejectsFasta(t *testing.T) {
	for _, compLevel := range []int{0, 1} {
		t.Run(fmt.Sprintf("compress_%d", compLevel), func(t *testing.T) {
			tmpDir := t.TempDir()
			inputPath := filepath.Join(tmpDir, "input.fasta")
			outputPath := filepath.Join(tmpDir, "output.fastq")
			if err := os.WriteFile(inputPath, []byte(">seq1\nACGT\n>seq2\nTT\n"), 0o644); err != nil {
				t.Fatal(err)
			}

			expectExitWithFastqError(t, func() {
				sortRecords(inputPath, outputPath, false, AvgPhred, compLevel, nil, DEFAULT_MIN_PHRED, -math.MaxFloat64, math.MaxFloat64)
			})
		})
	}
}

func TestRunNoSortRejectsFasta(t *testing.T) {
	tmpDir := t.TempDir()
	inputPath := filepath.Join(tmpDir, "input.fasta")
	outputPath := filepath.Join(tmpDir, "output.fastq")
	if err := os.WriteFile(inputPath, []byte(">seq1\nACGT\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	err := runNoSort(inputPath, outputPath, AvgPhred, nil, DEFAULT_MIN_PHRED, -math.MaxFloat64, math.MaxFloat64)
	if err == nil {
		t.Fatalf("expected FASTQ-only error")
	}
	if !strings.Contains(err.Error(), computedQualityFastqError) {
		t.Fatalf("runNoSort() error = %q, want %q", err.Error(), computedQualityFastqError)
	}
}

func TestRunPresortAcceptsFasta(t *testing.T) {
	tmpDir := t.TempDir()
	inputPath := filepath.Join(tmpDir, "input.fasta")
	outputPath := filepath.Join(tmpDir, "output.fasta")
	content := ">seq2 maxee=2.0\nACGT\n>seq1 maxee=0.5\nTT\n"
	if err := os.WriteFile(inputPath, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}

	if err := runPresort(inputPath, outputPath, MaxEE, false, -math.MaxFloat64, math.MaxFloat64); err != nil {
		t.Fatalf("runPresort() error = %v", err)
	}

	outBytes, err := os.ReadFile(outputPath)
	if err != nil {
		t.Fatal(err)
	}
	var gotIDs []string
	for _, line := range strings.Split(string(outBytes), "\n") {
		if strings.HasPrefix(line, ">") {
			gotIDs = append(gotIDs, strings.Fields(strings.TrimPrefix(line, ">"))[0])
		}
	}
	wantIDs := []string{"seq1", "seq2"}
	if !reflect.DeepEqual(gotIDs, wantIDs) {
		t.Fatalf("headersort FASTA IDs = %v, want %v", gotIDs, wantIDs)
	}
}

func TestSortRecordsFilteringBothModes(t *testing.T) {
	for _, compLevel := range []int{0, 1} {
		t.Run(fmt.Sprintf("compress_%d", compLevel), func(t *testing.T) {
			tmpDir := t.TempDir()
			inputPath := filepath.Join(tmpDir, "input.fastq")
			outputPath := filepath.Join(tmpDir, "output.fastq")
			records := []*fastx.Record{
				createTestRecord("high", "ACGT", "IIII"),
				createTestRecord("medium", "ACGT", "@@@@"),
				createTestRecord("low", "ACGT", "$$$$"),
				createTestRecord("edge", "ACGT", "5555"),
			}
			writeFastqRecords(t, inputPath, records)

			sortRecords(inputPath, outputPath, false, AvgPhred, compLevel, nil, DEFAULT_MIN_PHRED, 20, 35)

			gotIDs := readFastxIDs(t, outputPath)
			wantIDs := []string{"medium", "edge"}
			if !reflect.DeepEqual(gotIDs, wantIDs) {
				t.Fatalf("filtered sorted IDs = %v, want %v", gotIDs, wantIDs)
			}
		})
	}
}
