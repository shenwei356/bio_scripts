package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/codegangsta/cli"
	"github.com/shenwei356/pmap"
)

var key_index_pattern *regexp.Regexp = regexp.MustCompile(`^[\d,]+$`)

func main() {
	app := cli.NewApp()
	app.Name = "csv_grep"
	app.Version = "1.0"
	app.Author = "Wei Shen <https://github.com/shenwei356/bio_scripts>"
	app.Usage = "grep for csv format"

	app.Flags = []cli.Flag{
		// csv file
		cli.StringFlag{
			Name:  "k, key",
			Value: "1",
			Usage: "column number of key in csvfile. Multiple values shoud be separated by comma [1]",
		},
		cli.BoolFlag{
			Name:  "H, ignoretitle",
			Usage: "ignore title",
		},
		cli.StringFlag{
			Name:  "fs",
			Value: ",",
			Usage: "field separator [,]",
		},
		cli.StringFlag{
			Name:  "fs-out",
			Usage: "field separator of ouput [same as --fs]",
		},
		cli.BoolFlag{
			Name:  "t, tab",
			Usage: `field separator is "\t". Quote char is "\t"`,
		},

		// pattern file
		cli.StringFlag{
			Name:  "p, pattern",
			Usage: "query pattern",
		},
		cli.StringFlag{
			Name:  "pf, patternfile",
			Usage: "pattern file",
		},
		cli.StringFlag{
			Name:  "pk",
			Value: "1",
			Usage: "column number of key in pattern file. Multiple values shoud be separated by comma [1]",
		},
		cli.StringFlag{
			Name:  "pfs",
			Value: ",",
			Usage: "field separator of pattern file [,]",
		},

		// matching option
		cli.BoolFlag{
			Name:  "r, use-regexp",
			Usage: "use regular expression",
		},
		cli.BoolFlag{
			Name:  "d, speedup",
			Usage: "delete matched pattern when matching one record",
		},
		cli.BoolFlag{
			Name:  "i, invert",
			Usage: "invert match (do not match)",
		},

		// performance
		cli.IntFlag{
			Name:  "j, ncpus",
			Value: runtime.NumCPU(),
			Usage: fmt.Sprintf("CPU number [%d]", runtime.NumCPU()),
		},
		cli.IntFlag{
			Name:  "c, chunksize",
			Value: 1000,
			Usage: "chunk size [1000]",
		},

		// outpute
		cli.StringFlag{
			Name:  "o, outfile",
			Usage: "output file [stdout]",
		},
		cli.BoolFlag{
			Name:  "vv, verbose",
			Usage: "verbosely print information",
		},
	}

	app.Action = func(c *cli.Context) {
		if c.String("pattern") == "" && c.String("patternfile") == "" {
			fmt.Fprintln(os.Stderr, "one or both of option -p and -pf needed")
			os.Exit(1)
		}

		verbose := c.Bool("verbose")

		ncpus := c.Int("ncpus")
		if ncpus < 1 {
			ncpus = 1
		}
		runtime.GOMAXPROCS(ncpus)

		chunksize := c.Int("chunksize")
		if chunksize < 1 {
			chunksize = 1
		}

		ignoretitle := c.Bool("ignoretitle")

		key_index := get_key_index(c.String("key"))
		pk_index := get_key_index(c.String("pk"))

		fs := string_of_one_char_to_rune(c.String("fs"))
		pfs := string_of_one_char_to_rune(c.String("pfs"))
		if c.Bool("tab") {
			fs = []rune("\t")[0]
			pfs = []rune("\t")[0]
		}
		var fsout rune
		if c.String("fs-out") == "" {
			fsout = fs
		} else {
			fsout = string_of_one_char_to_rune(c.String("fs-out"))
		}

		// patterns
		patterns := pmap.NewParallelMap()
		if c.String("pattern") != "" {
			patterns.Set(c.String("pattern"), nil)
		}
		patternfile := c.String("patternfile")
		if patternfile != "" {
			for _, pattern := range read_patterns(patternfile, pfs, pk_index) {
				patterns.Set(pattern, nil)
			}
		}

		useregexp := c.Bool("use-regexp")
		invert := c.Bool("invert")
		speedup := c.Bool("speedup")
		// compile pattern
		if useregexp {
			var wg sync.WaitGroup
			tokens := make(chan int, ncpus)
			for k, _ := range patterns.Map {
				tokens <- 1
				wg.Add(1)
				key, _ := k.(string)
				go func(m *pmap.ParallelMap, key string) {
					defer func() {
						wg.Done()
						<-tokens
					}()
					r, err := regexp.Compile(key)
					if err != nil {
						fmt.Fprintf(os.Stderr, "fail to compile regexp: %s\n", k)
						os.Exit(1)
					}
					m.Set(key, r)
				}(patterns, key)
			}
			wg.Wait()
		}
		if verbose {
			fmt.Fprintf(os.Stderr, "   patterns: %d\n", len(patterns.Map))
			fmt.Fprintf(os.Stderr, "   chunsize: %d\n", chunksize)
			fmt.Fprintf(os.Stderr, "concurrency: %d\n", ncpus)
		}

		// outfile
		outfile := c.String("outfile")
		var fhout *os.File
		var err error
		if outfile != "" {
			fhout, err = os.Create(outfile)
			if err != nil {
				fmt.Fprintf(os.Stderr, "fail to write file: %s\n", outfile)
			}
		} else {
			fhout = os.Stdout
		}
		// writter
		writter := csv.NewWriter(fhout)
		writter.Comma = fsout
		// writter.LazyQuotes = true

		// reducer
		out := make(chan Chunk, ncpus)
		done := make(chan int)
		go func() {
			buffer := make(map[int][]Record, ncpus)
			var id int = 0
			for chunk := range out {
				buffer[chunk.Id] = chunk.Data
				if _, ok := buffer[id]; ok {
					for _, record := range buffer[id] {
						writter.Write(record.Content)
					}
					delete(buffer, id)
					id++
				}
			}
			// sort by id
			ids := make([]int, len(buffer))
			i := 0
			for id, _ := range buffer {
				ids[i] = id
				i++
			}
			sort.Ints(ids)
			for _, id := range ids {
				for _, record := range buffer[id] {
					writter.Write(record.Content)
				}
			}
			writter.Flush()
			if outfile != "" {
				fhout.Close()
			}
			done <- 1
		}()

		// mapper
		var wg sync.WaitGroup
		tokens := make(chan int, ncpus)
		csvfiles := c.Args()
		if len(csvfiles) == 0 {
			csvfiles = []string{"stdin"}
		}
		for _, file := range csvfiles {
			var fh *os.File
			var err error
			if file == "stdin" {
				fh = os.Stdin
			} else {
				fh, err = os.Open(file)
				if err != nil {
					fmt.Fprintln(os.Stderr, err)
					os.Exit(1)
				}
				defer fh.Close()
			}

			ch := csv_reader(fh, fs, key_index, chunksize, ncpus, ignoretitle)

			for chunk := range ch {
				tokens <- 1
				wg.Add(1)
				go func(chunk Chunk, patterns *pmap.ParallelMap, useregep bool, invert bool, speedup bool) {
					defer func() {
						wg.Done()
						<-tokens
					}()
					chunk2 := check_chunk(chunk, patterns, useregep, invert, speedup)
					out <- chunk2
				}(chunk, patterns, useregexp, invert, speedup)
			}
		}
		wg.Wait()
		close(out)
		<-done
		patterns.Stop()
	}
	app.Run(os.Args)
}

func check_chunk(chunk Chunk, patterns *pmap.ParallelMap, useregep bool, invert bool, speedup bool) Chunk {
	chunk2 := Chunk{chunk.Id, make([]Record, 0)}
	for _, record := range chunk.Data {
		hit := false
		if useregep {
			for k, r := range patterns.Map {
				re, _ := r.(*regexp.Regexp)
				if re.MatchString(record.Key) {
					hit = true
					if speedup {
						patterns.Delete(k)
					}
					break
				}
			}
		} else {
			if _, ok := patterns.Get(record.Key); ok {
				hit = true
				if speedup {
					patterns.Delete(record.Key)
				}
			}
		}
		if invert {
			if hit {
				continue
			}
		} else {
			if !hit {
				continue
			}
		}
		chunk2.Data = append(chunk2.Data, record)
	}
	return chunk2
}

func read_patterns(file string, fs rune, key_index []int) []string {
	patterns := make([]string, 0)

	fh, err := os.Open(file)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	defer fh.Close()

	ch := csv_reader(fh, fs, key_index, 1000, 1, false)
	for {
		select {
		case chunk := <-ch:
			if len(chunk.Data) == 0 {
				return patterns
			}
			for _, record := range chunk.Data {
				patterns = append(patterns, record.Key)
			}
		}
	}
	return patterns
}

type Record struct {
	Key     string
	Content []string
}

type Chunk struct {
	Id   int
	Data []Record
}

func csv_reader(file *os.File, fs rune, key_index []int, chunksize int, buffersize int, ignoretitle bool) chan Chunk {
	ch := make(chan Chunk, buffersize)
	go func() {
		defer func() {
			close(ch)
		}()
		var id int = 0
		chunkdata := make([]Record, 0)

		reader := csv.NewReader(file)
		reader.Comma = fs
		reader.LazyQuotes = true
		once := true
		for {
			record, err := reader.Read()
			if err == io.EOF {
				break
			} else if err != nil {
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}
			if once && ignoretitle {
				once = false
				continue
			}
			n := len(record)
			keys := make([]string, 0)
			for _, i := range key_index {
				if i > n {
					fmt.Fprintf(os.Stderr, "key index (%d) is beyond number of column (%d)\n", i, n)
					os.Exit(1)
				}
				keys = append(keys, record[i-1])
			}
			chunkdata = append(chunkdata, Record{strings.Join(keys, "_"), record})
			if len(chunkdata) == chunksize {
				ch <- Chunk{id, chunkdata}
				chunkdata = make([]Record, 0)
				id += 1
			}
		}
		if len(chunkdata) > 0 {
			ch <- Chunk{id, chunkdata}
		}
	}()
	return ch
}

func string_of_one_char_to_rune(str string) rune {
	if len(str) != 1 {
		fmt.Fprintf(os.Stderr, "string (%s) should be single character\n", str)
		os.Exit(1)
	}
	return []rune(str)[0]
}

func check_key_index(index string) {
	if !key_index_pattern.MatchString(index) {
		fmt.Fprintln(os.Stderr, "invalid key index")
		os.Exit(1)
	}
}

func get_key_index(index string) (indice []int) {
	if !key_index_pattern.MatchString(index) {
		fmt.Fprintln(os.Stderr, "invalid key index")
		os.Exit(1)
	}

	var items []string
	if strings.Contains(index, ",") {
		items = strings.Split(index, ",")
	} else {
		items = []string{index}
	}
	for _, item := range items {
		if i, err := strconv.Atoi(item); err == nil {
			indice = append(indice, i)
		} else {
			fmt.Fprintln(os.Stderr, "invalid key index")
			os.Exit(1)
		}
	}
	return indice
}
