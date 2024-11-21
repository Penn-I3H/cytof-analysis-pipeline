package main

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"os/exec"
	"runtime"
	"strings"
)

type Analysis struct {
	Filename string
}

func main() {
	programLevel := new(slog.LevelVar)
	logger := slog.New(slog.NewJSONHandler(os.Stderr, &slog.HandlerOptions{Level: programLevel}))
	slog.SetDefault(logger)

	inputDir := os.Getenv("INPUT_DIR")

	totalCPU := runtime.NumCPU()
	fmt.Println("Number of Logical CPUs on machine ", totalCPU)
	defaultCPU := runtime.GOMAXPROCS(0)
	fmt.Println("DefaultCPU(s) ", defaultCPU)

	files, err := os.ReadDir(inputDir)
	if err != nil {
		fmt.Println(err)
		log.Fatalf("error reading input files")
	}

	numberOfAnalyses := len(files)
	analyses := make(chan *Analysis, numberOfAnalyses)
	results := make(chan string, numberOfAnalyses)

	log.Println("Starting pipeline")

	NumConcurrentWorkers := 2 // concurrent workers set to 2
	for i := 1; i <= NumConcurrentWorkers; i++ {
		go worker(i, analyses, results)
	}

	// create work
	for _, j := range files {
		analyses <- &Analysis{j.Name()}
	}
	close(analyses)

	// wait for the done signal
	for j := 1; j <= numberOfAnalyses; j++ {
		log.Println(<-results)
	}

	log.Println("Analysis complete")
}

// processes analyses
func worker(w int, analyses <-chan *Analysis, results chan<- string) {
	for r := range analyses {
		log.Printf("processing %s on worker: %v", r.Filename, w)
		err := Process(r.Filename)
		if err != nil {
			results <- err.Error()
		} else {
			results <- fmt.Sprintf("%v done", w)
		}
	}
}

func Process(filename string) error {
	// run pipeline
	cmd := exec.Command("Rscript", "/service/main_parallel.R", filename)
	cmd.Dir = "/service"
	var stdout strings.Builder
	var stderr strings.Builder
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		fmt.Println(stderr.String())
		return err
	}

	fmt.Println(stdout.String())
	return nil
}
