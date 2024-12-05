// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"math"
	"math/rand"
	"runtime"
	"sync"
	"time"

	"gonum.org/v1/gonum/stat"
)

/*
   Comparing each condition's module correlations against a null distribution.
   For each condition, we generate 1000 random modules of the same size
   and compare the actual module's correlation pattern against these random modules.
   If the p-value is below 0.05, we can say that the module's correlation pattern
   is significantly different from random expectation in that condition.
*/

type NullDistributionStats struct {
	Name             string
	Size             int
	C1NullTStatistic float64
	C1NullPValue     float64
	C2NullTStatistic float64
	C2NullPValue     float64
}

func createNullDistributions(moduleGenes []string, condition1Data, condition2Data map[string][]float64) (float64, float64, float64, float64) {
	const numPermutations = 1000

	// Get all gene names from each condition
	var c1Genes, c2Genes []string
	for gene := range condition1Data {
		c1Genes = append(c1Genes, gene)
	}
	for gene := range condition2Data {
		c2Genes = append(c2Genes, gene)
	}

	moduleSize := len(moduleGenes)

	// Calculate actual correlations for both conditions
	actualC1Corrs := getModuleCorrelations(moduleGenes, condition1Data)
	actualC2Corrs := getModuleCorrelations(moduleGenes, condition2Data)

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

	// Calculate t-statistics and p-values comparing actual vs null distributions
	c1Tstat, c1Pval := compareWithNull(actualC1Corrs, c1NullCorrs)
	c2Tstat, c2Pval := compareWithNull(actualC2Corrs, c2NullCorrs)

	return c1Tstat, c1Pval, c2Tstat, c2Pval
}

func sampleGenes(genes []string, size int, r *rand.Rand) []string {
	sampled := make([]string, size)
	indices := r.Perm(len(genes))
	for i := 0; i < size; i++ {
		sampled[i] = genes[indices[i]]
	}
	return sampled
}

func compareWithNull(actualCorrs, nullCorrs []float64) (float64, float64) {
	// Check if we have enough data
	if len(actualCorrs) < 2 || len(nullCorrs) < 2 {
		return 0, 1
	}

	meanActual := stat.Mean(actualCorrs, nil)
	meanNull := stat.Mean(nullCorrs, nil)
	varActual := stat.Variance(actualCorrs, nil)
	varNull := stat.Variance(nullCorrs, nil)

	if varActual == 0 || varNull == 0 {
		return 0, 1
	}

	nActual := float64(len(actualCorrs))
	nNull := float64(len(nullCorrs))

	se := math.Sqrt((varActual / nActual) + (varNull / nNull))

	if se == 0 {
		return 0, 1
	}

	tstat := (meanActual - meanNull) / se

	// Calculate p-value
	z := math.Abs(tstat)
	pval := 2 * (1 - normalCDF(z))

	return tstat, pval
}

func analyzeModuleNullDistribution(moduleName string, moduleMap map[string]string, condition1Data, condition2Data map[string][]float64) NullDistributionStats {
	// Get genes in this module
	var moduleGenes []string
	for gene, module := range moduleMap {
		if module == moduleName {
			moduleGenes = append(moduleGenes, gene)
		}
	}

	// Calculate t-statistics and p-values for each condition vs its null distribution
	c1NullTstat, c1Pval, c2NullTstat, c2Pval := createNullDistributions(moduleGenes, condition1Data, condition2Data)

	return NullDistributionStats{
		Name:             moduleName,
		Size:             len(moduleGenes),
		C1NullTStatistic: c1NullTstat,
		C1NullPValue:     c1Pval,
		C2NullTStatistic: c2NullTstat,
		C2NullPValue:     c2Pval,
	}
}
