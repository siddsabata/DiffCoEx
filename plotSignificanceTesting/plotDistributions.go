// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"encoding/csv"
	"fmt"
	"image/color"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"
	"time"

	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

func generateNullCorrelations(moduleGenes, c1Genes, c2Genes []string,
	condition1Data, condition2Data map[string][]float64, numPermutations int) ([]float64, []float64) {

	moduleSize := len(moduleGenes)

	// Create channels for parallel processing
	c1Results := make(chan []float64, numPermutations)
	c2Results := make(chan []float64, numPermutations)

	// Create worker pool
	numWorkers := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup

	// Launch workers for condition 1
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			r := rand.New(rand.NewSource(time.Now().UnixNano()))

			permutationsPerWorker := numPermutations / numWorkers
			for i := 0; i < permutationsPerWorker; i++ {
				// Randomly sample genes for null module
				nullGenes := sampleGenes(c1Genes, moduleSize, r)
				nullCorrs := getModuleCorrelations(nullGenes, condition1Data)
				c1Results <- nullCorrs
			}
		}()
	}

	// Launch workers for condition 2
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			r := rand.New(rand.NewSource(time.Now().UnixNano()))

			permutationsPerWorker := numPermutations / numWorkers
			for i := 0; i < permutationsPerWorker; i++ {
				nullGenes := sampleGenes(c2Genes, moduleSize, r)
				nullCorrs := getModuleCorrelations(nullGenes, condition2Data)
				c2Results <- nullCorrs
			}
		}()
	}

	// Close results channels when all workers are done
	go func() {
		wg.Wait()
		close(c1Results)
		close(c2Results)
	}()

	// Collect results
	var c1NullCorrs, c2NullCorrs []float64
	for corrs := range c1Results {
		c1NullCorrs = append(c1NullCorrs, corrs...)
	}
	for corrs := range c2Results {
		c2NullCorrs = append(c2NullCorrs, corrs...)
	}

	return c1NullCorrs, c2NullCorrs
}

func plotConditionDistribution(moduleName, conditionName string, actualCorrs, nullCorrs []float64) {
	p := plot.New()

	// Calculate t-statistic and p-value
	tstat, pval := calculateTStatistic(actualCorrs, nullCorrs)

	// Set plot title and labels with t-statistic and p-value
	p.Title.Text = fmt.Sprintf("Module %s - %s\nt-statistic: %.2f, p-value: %.4f",
		moduleName, conditionName, tstat, pval)
	p.X.Label.Text = "Correlation"
	p.Y.Label.Text = "Density"

	// Create histograms
	nullHist, _ := plotter.NewHist(plotter.Values(nullCorrs), 50)
	actualHist, _ := plotter.NewHist(plotter.Values(actualCorrs), 50)

	// Style the histograms
	nullHist.FillColor = color.RGBA{R: 200, B: 200, A: 255}
	nullHist.LineStyle.Width = vg.Points(1)
	nullHist.Normalize(1)

	actualHist.FillColor = color.RGBA{R: 255, A: 255}
	actualHist.LineStyle.Width = vg.Points(1)
	actualHist.Normalize(1)

	// Add to plot
	p.Add(nullHist, actualHist)

	// Add legend
	p.Legend.Add("Null Distribution", nullHist)
	p.Legend.Add("Actual Correlations", actualHist)
	p.Legend.Top = true

	// Save the plot to the correct directory
	outputPath := filepath.Join("output", "plotting",
		fmt.Sprintf("%s_%s_distribution.png", moduleName, conditionName))

	if err := p.Save(6*vg.Inch, 4*vg.Inch, outputPath); err != nil {
		fmt.Printf("Error saving plot: %v\n", err)
	}
}

func calculateTStatistic(actual, null []float64) (float64, float64) {
	// Check if we have enough data
	if len(actual) < 2 || len(null) < 2 {
		return 0, 1
	}

	// Calculate means of absolute correlations since we care about strength
	// regardless of direction (positive or negative)
	actualAbs := make([]float64, len(actual))
	nullAbs := make([]float64, len(null))
	for i, v := range actual {
		actualAbs[i] = math.Abs(v)
	}
	for i, v := range null {
		nullAbs[i] = math.Abs(v)
	}

	actualMean := stat.Mean(actualAbs, nil)
	nullMean := stat.Mean(nullAbs, nil)

	actualVar := stat.Variance(actualAbs, nil)
	nullVar := stat.Variance(nullAbs, nil)

	if actualVar == 0 || nullVar == 0 {
		return 0, 1
	}

	actualN := float64(len(actual))
	nullN := float64(len(null))

	se := math.Sqrt((actualVar / actualN) + (nullVar / nullN))
	if se == 0 {
		return 0, 1
	}

	t := (actualMean - nullMean) / se
	z := math.Abs(t)
	p := 2 * (1 - normalCDF(z))

	if math.IsNaN(t) || math.IsNaN(p) {
		return 0, 1
	}

	return t, p
}

// Add the normalCDF function if it's not already in your file
func normalCDF(x float64) float64 {
	return 0.5 * (1 + math.Erf(x/math.Sqrt2))
}

func sampleGenes(genes []string, size int, r *rand.Rand) []string {
	sampled := make([]string, size)
	indices := r.Perm(len(genes))
	for i := 0; i < size; i++ {
		sampled[i] = genes[indices[i]]
	}
	return sampled
}

func loadModules(filename string) (map[string]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	// Skip header
	_, err = reader.Read()
	if err != nil {
		return nil, err
	}

	moduleMap := make(map[string]string) // gene -> module
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		moduleMap[record[0]] = record[1]
	}
	return moduleMap, nil
}

func loadExpressionData(filename string) (map[string][]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	data := make(map[string][]float64)

	// Skip header if it exists
	_, err = reader.Read()
	if err != nil {
		return nil, err
	}

	for {
		record, err := reader.Read()
		if err != nil {
			break
		}

		geneName := record[0]
		values := make([]float64, 0)

		// Convert string values to float64, starting from column 1
		for _, val := range record[1:] {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				continue
			}
			values = append(values, f)
		}

		if len(values) > 0 {
			data[geneName] = values
		}
	}
	return data, nil
}

func getModuleCorrelations(genes []string, expressionData map[string][]float64) []float64 {
	// Find minimum sample size across all genes in the dataset
	minSamples := -1
	for _, expr := range expressionData {
		if minSamples == -1 || len(expr) < minSamples {
			minSamples = len(expr)
		}
	}

	// Calculate total number of correlations needed
	numGenes := len(genes)
	numCorrelations := (numGenes * (numGenes - 1)) / 2

	// Create channels for parallel processing
	results := make(chan float64, numCorrelations)
	numWorkers := runtime.GOMAXPROCS(0)
	var wg sync.WaitGroup

	// Calculate how many gene pairs each worker should handle
	pairsPerWorker := numCorrelations / numWorkers
	if pairsPerWorker < 1 {
		pairsPerWorker = 1
	}

	// Launch workers
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		startIdx := w * pairsPerWorker
		endIdx := (w + 1) * pairsPerWorker
		if w == numWorkers-1 {
			endIdx = numCorrelations
		}

		go func(start, end int) {
			defer wg.Done()

			for idx := start; idx < end; idx++ {
				i, j := getGeneIndices(idx)
				if i >= len(genes) || j >= len(genes) {
					continue
				}

				gene1 := genes[i]
				gene2 := genes[j]

				expr1, ok1 := expressionData[gene1]
				expr2, ok2 := expressionData[gene2]

				if ok1 && ok2 {
					expr1 = expr1[:minSamples]
					expr2 = expr2[:minSamples]

					corr := stat.Correlation(expr1, expr2, nil)

					if !math.IsNaN(corr) {
						results <- corr
					}
				}
			}
		}(startIdx, endIdx)
	}

	go func() {
		wg.Wait()
		close(results)
	}()

	var correlations []float64
	for corr := range results {
		correlations = append(correlations, corr)
	}

	return correlations
}

func getGeneIndices(idx int) (i, j int) {
	i = int(math.Floor((-1 + math.Sqrt(1+8*float64(idx))) / 2))
	j = idx - (i * (i + 1) / 2)
	return i, j + i + 1
}
