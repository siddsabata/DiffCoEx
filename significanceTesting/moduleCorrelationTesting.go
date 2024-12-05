// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"encoding/csv"
	"math"
	"os"
	"runtime"
	"strconv"
	"sync"

	"gonum.org/v1/gonum/stat"
)

/*
	Comparing the difference in module expression between two conditions.
	If the p-value is below 0.05, we can say that the module is
	differentially expressed between the two conditions.
*/

type ModuleStats struct {
	Name       string
	TStatistic float64
	PValue     float64
	Size       int
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

func analyzeModule(moduleName string, moduleMap map[string]string, condition1Data, condition2Data map[string][]float64) ModuleStats {
	// Get genes in this module
	var moduleGenes []string
	for gene, module := range moduleMap {
		if module == moduleName {
			moduleGenes = append(moduleGenes, gene)
		}
	}

	// Get correlation values for both conditions
	condition1Corrs := getModuleCorrelations(moduleGenes, condition1Data)
	condition2Corrs := getModuleCorrelations(moduleGenes, condition2Data)

	// Calculate t-statistic and p-value manually
	tstat, pval := calculateTTest(condition1Corrs, condition2Corrs)

	return ModuleStats{
		Name:       moduleName,
		TStatistic: tstat,
		PValue:     pval,
		Size:       len(moduleGenes),
	}
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

			// Calculate which gene pairs this worker should process
			for idx := start; idx < end; idx++ {
				// Convert linear index to i,j coordinates
				i, j := getGeneIndices(idx)
				if i >= len(genes) || j >= len(genes) {
					continue
				}

				gene1 := genes[i]
				gene2 := genes[j]

				// Get expression values
				expr1, ok1 := expressionData[gene1]
				expr2, ok2 := expressionData[gene2]

				if ok1 && ok2 {
					// Use samples up to the minimum sample size
					expr1 = expr1[:minSamples]
					expr2 = expr2[:minSamples]

					// Calculate correlation
					corr := stat.Correlation(expr1, expr2, nil)

					// Only send if correlation is valid
					if !math.IsNaN(corr) {
						results <- corr
					}
				}
			}
		}(startIdx, endIdx)
	}

	// Close results channel when all workers are done
	go func() {
		wg.Wait()
		close(results)
	}()

	// Collect results
	var correlations []float64
	for corr := range results {
		correlations = append(correlations, corr)
	}

	return correlations
}

// Helper function to convert linear index to i,j coordinates
func getGeneIndices(idx int) (i, j int) {
	// Convert linear index back to i,j coordinates
	// This is the inverse of the triangular number formula
	i = int(math.Floor((-1 + math.Sqrt(1+8*float64(idx))) / 2))
	j = idx - (i * (i + 1) / 2)
	return i, j + i + 1
}

func calculateTTest(x, y []float64) (tstat, pval float64) {
	// Check if we have enough data
	if len(x) < 2 || len(y) < 2 {
		return 0, 1 // Return no significance if not enough data
	}

	meanX := stat.Mean(x, nil)
	meanY := stat.Mean(y, nil)
	varX := stat.Variance(x, nil)
	varY := stat.Variance(y, nil)

	// Check for zero variance
	if varX == 0 || varY == 0 {
		return 0, 1 // Return no significance if no variance
	}

	nx := float64(len(x))
	ny := float64(len(y))

	// Calculate pooled standard error
	se := math.Sqrt((varX / nx) + (varY / ny))

	// Check for zero standard error
	if se == 0 {
		return 0, 1
	}

	// Calculate t-statistic
	tstat = (meanX - meanY) / se

	// Calculate approximate p-value using normal distribution
	z := math.Abs(tstat)
	pval = 2 * (1 - normalCDF(z))

	// Check for NaN
	if math.IsNaN(tstat) || math.IsNaN(pval) {
		return 0, 1
	}

	return tstat, pval
}

// normalCDF returns the cumulative distribution function of the standard normal distribution
func normalCDF(x float64) float64 {
	return 0.5 * (1 + math.Erf(x/math.Sqrt(2)))
}

func getUniqueModules(moduleMap map[string]string) map[string]bool {
	modules := make(map[string]bool)
	for _, module := range moduleMap {
		modules[module] = true
	}
	return modules
}
