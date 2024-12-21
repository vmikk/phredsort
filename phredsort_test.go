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

