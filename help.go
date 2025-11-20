package main

import (
	"fmt"

	"github.com/spf13/cobra"
)

// Custom help function used
// It provides nicely formatted help messages for the root command and other subcommands
func helpFunc(cmd *cobra.Command, args []string) {

	// Specialized help for subcommands
	switch cmd.Name() {
	case "headersort":
		fmt.Printf(`
%s

%s
  Sort FASTA/FASTQ sequences using pre-computed quality scores stored
  directly in sequence headers. Supports both space-separated
  (">seq1 maxee=2") and semicolon-separated (">seq1;maxee=2") formats.

%s
  %s
  %s
  %s
  %s
  %s
  %s

%s
  %s

%s
  %s
  %s
  %s

`,
			bold(getColorizedLogo()+" phredsort headersort - Sorts sequences using header quality metrics"),
			bold(yellow("Description:")),
			bold(yellow("Flags:")),
			cyan("-i, --in")+" <string>      : Input FASTA/FASTQ file (required)",
			cyan("-o, --out")+" <string>     : Output FASTA/FASTQ file (required)",
			cyan("-s, --metric")+" <string>  : Header metric to use (avgphred, maxee, meep, lqcount, lqpercent) (default, 'avgphred')",
			cyan("-a, --ascending")+" <bool> : Sort in ascending order of the header metric (default, false)",
			cyan("-m, --minqual")+" <float>  : Minimum header metric value for filtering (optional)",
			cyan("-M, --maxqual")+" <float>  : Maximum header metric value for filtering (optional)",
			bold(yellow("Examples:")),
			cyan("phredsort headersort -i input.fasta -o output.fasta --metric maxee"),
			bold(yellow("Supported header formats:")),
			`  ">seq1 maxee=2.5 size=100"`,
			`  ">seq1;maxee=2.5;size=100"`,
			`  (other metrics use the same "name=value" syntax)`,
		)
		return
	case "sort":
		fmt.Printf(`
%s

%s
  Sort FASTQ sequences by calculating quality metrics from base quality
  scores. Can read from a file or stdin and write to a file or stdout.

%s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s

%s
  %s
  %s

`,
			bold(getColorizedLogo()+" phredsort sort - Sorts FASTQ based on computed quality metrics"),
			bold(yellow("Description:")),
			bold(yellow("Flags:")),
			cyan("-i, --in")+" <string>      : Input FASTQ file (required, use '-' for stdin)",
			cyan("-o, --out")+" <string>     : Output FASTQ file (required, use '-' for stdout)",
			cyan("-s, --metric")+" <string>  : Quality metric (avgphred, maxee, meep, lqcount, lqpercent) (default, 'avgphred')",
			cyan("-m, --minqual")+" <float>  : Minimum quality threshold for filtering (optional)",
			cyan("-M, --maxqual")+" <float>  : Maximum quality threshold for filtering (optional)",
			cyan("-p, --minphred")+" <int>   : Quality threshold for 'lqcount' and 'lqpercent' metrics (default, 15)",
			cyan("-H, --header")+" <string>  : Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')",
			cyan("-a, --ascending")+" <bool> : Sort sequences in ascending order of quality (default, false)",
			cyan("-c, --compress")+" <int>   : Memory compression level for stdin-based mode (0=disabled, 1-22; default, 1)",
			cyan("-v, --version")+"          : Show version information",
			bold(yellow("Examples:")),
			cyan("phredsort sort --metric avgphred --in input.fq.gz --out output.fq.gz"),
			cyan("cat input.fq | phredsort sort --compress 0 -i - -o - > sorted.fq"),
		)
		return
	case "nosort":
		fmt.Printf(`
%s

%s
  Estimate FASTQ sequence quality metrics without sorting. Records are streamed
  from input to output in their original order, with optional filtering and
  header annotation.

%s
  %s
  %s
  %s
  %s
  %s
  %s

%s
  %s
  %s

`,
			bold(getColorizedLogo()+" phredsort nosort - Estimates FASTQ quality without sorting"),
			bold(yellow("Description:")),
			bold(yellow("Flags:")),
			cyan("-i, --in")+" <string>      : Input FASTQ file (required, use '-' for stdin)",
			cyan("-o, --out")+" <string>     : Output FASTQ file (required, use '-' for stdout)",
			cyan("-s, --metric")+" <string>  : Quality metric (avgphred, maxee, meep, lqcount, lqpercent) (default, 'avgphred')",
			cyan("-m, --minqual")+" <float>  : Minimum quality threshold for filtering (optional)",
			cyan("-M, --maxqual")+" <float>  : Maximum quality threshold for filtering (optional)",
			cyan("-p, --minphred")+" <int>   : Quality threshold for 'lqcount' and 'lqpercent' metrics (default, 15)",
			bold(yellow("Examples:")),
			cyan("phredsort nosort --metric avgphred --in input.fq.gz --out output.fq.gz"),
			cyan("cat input.fq | phredsort nosort --metric maxee --maxqual 1 -i - -o - > output.fq"),
		)
		return
	}

	// Default: root command help
	fmt.Printf(`
%s

%s
  %s
  %s
  %s
  %s
  %s

%s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s
  %s

%s
  %s
  %s
  %s

%s
  # File-based mode (reads from a file, lower memory usage)
  %s

  # Stdin-based mode (reads from stdin, higher memory usage)
  %s

  # Remove sequences with average Phred quality score below 20, 
  # add several quality metrics to sequence headers
  %s

  # Sort sequences using pre-computed quality scores (stored in headers)
  # (e.g., ">seq1 maxee=1 avgphred=15.5" or ">seq1;maxee=1;avgphred=15.5")
  %s

  # Estimate quality metrics without sorting, add quality metrics to headers, and filter sequences
  %s

%s
  https://github.com/vmikk/phredsort

`,
		bold(getColorizedLogo()+" phredsort v."+VERSION+" - Sorts FASTQ based on different sequence quality metrics"),
		bold(yellow("Quality metrics:")),
		cyan("avgphred")+"  : average Phred quality score",
		cyan("maxee")+"     : maximum expected error (absolute number)",
		cyan("meep")+"      : maximum expected error (percentage per sequence length)",
		cyan("lqcount")+"   : number of bases below quality threshold (default, 15)",
		cyan("lqpercent")+" : percentage of bases below quality threshold",
		bold(yellow("Flags:")),
		cyan("-i, --in")+" <string>      : Input FASTQ file (required, use '-' for stdin)",
		cyan("-o, --out")+" <string>     : Output FASTQ file (required, use '-' for stdout)",
		cyan("-s, --metric")+" <string>  : Quality metric (avgphred, maxee, meep, lqcount, lqpercent) (default, 'avgphred')",
		cyan("-m, --minqual")+" <float>  : Minimum quality threshold for filtering (optional)",
		cyan("-M, --maxqual")+" <float>  : Maximum quality threshold for filtering (optional)",
		cyan("-p, --minphred")+" <int>   : Quality threshold for 'lqcount' and 'lqpercent' metrics (default, 15)",
		cyan("-H, --header")+" <string>  : Comma-separated list of metrics to add to headers (e.g., 'avgphred,maxee,length')",
		cyan("-a, --ascending")+" <bool> : Sort sequences in ascending order of quality (default, false)",
		cyan("-c, --compress")+" <int>   : Memory compression level for stdin-based mode (0=disabled, 1-22; default, 1)",
		cyan("-h, --help")+"             : Show help message",
		cyan("-v, --version")+"          : Show version information",
		bold(yellow("Subcommands:")),
		cyan("sort")+"       : Sort sequences by computing quality metrics from base qualities",
		cyan("nosort")+"     : Estimate quality and optionally filter/annotate without sorting",
		cyan("headersort")+" : Sort sequences using pre-computed quality scores in headers",
		bold(yellow("Usage examples:")),
		cyan("phredsort --metric avgphred --in input.fq.gz --out output.fq.gz"),
		cyan("cat input.fq | phredsort --compress 0 -i - -o - > sorted.fq"),
		cyan("phredsort -i inp.fq.gz -o out.fq.gz --metric avgphred --minqual 20 --header avgphred,maxee,lqpercent,length"),
		cyan("phredsort headersort -i inp.fq.gz -o out.fq.gz --metric maxee"),
		cyan("phredsort nosort --metric maxee --maxqual 1 --header maxee -i - -o - > output.fq"),
		bold(yellow("More information:")),
	)
}
