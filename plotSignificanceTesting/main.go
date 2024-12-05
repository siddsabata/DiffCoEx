// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"fmt"
	"log"
	"os"
)

func main() {
	// ./plotSignificanceTesting moduleMap condition1Data condition2Data module
	// Check if correct number of arguments are provided
	if len(os.Args) != 5 {
		fmt.Println("Usage: ./plotSignificanceTesting moduleMap condition1Data condition2Data module")
		fmt.Println("Example: ./plotSignificanceTesting data/golub/golub_diffcoex.csv data/golub/aml_samples.csv data/golub/all_samples.csv M1")
		os.Exit(1)
	}

	// Get file paths from command line arguments
	moduleMapPath := os.Args[1]
	condition1Path := os.Args[2]
	condition2Path := os.Args[3]
	targetModule := os.Args[4]

	// Create output/plotting directory if it doesn't exist
	outputDir := "output/plotting"
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		log.Fatal("Error creating output directory:", err)
	}

	// Load module assignments
	moduleMap, err := loadModules(moduleMapPath)
	if err != nil {
		log.Fatal("Error loading modules:", err)
	}

	// Load expression data for both conditions
	condition1Data, err := loadExpressionData(condition1Path)
	if err != nil {
		log.Fatal("Error loading condition 1 data:", err)
	}

	condition2Data, err := loadExpressionData(condition2Path)
	if err != nil {
		log.Fatal("Error loading condition 2 data:", err)
	}

	// Check if the specified module exists
	if !moduleExists(targetModule, moduleMap) {
		log.Fatalf("Module %s not found in the module map", targetModule)
	}

	fmt.Printf("Plotting distributions for module %s...\n", targetModule)
	plotModuleDistributions(targetModule, moduleMap, condition1Data, condition2Data)
	fmt.Println("Done!")
}

func moduleExists(targetModule string, moduleMap map[string]string) bool {
	for _, module := range moduleMap {
		if module == targetModule {
			return true
		}
	}
	return false
}

func plotModuleDistributions(moduleName string, moduleMap map[string]string, condition1Data, condition2Data map[string][]float64) {
	// Get genes in this module
	var moduleGenes []string
	for gene, module := range moduleMap {
		if module == moduleName {
			moduleGenes = append(moduleGenes, gene)
		}
	}

	// Get actual correlations
	actualC1Corrs := getModuleCorrelations(moduleGenes, condition1Data)
	actualC2Corrs := getModuleCorrelations(moduleGenes, condition2Data)

	// Generate null distributions
	const numPermutations = 1000
	var c1Genes, c2Genes []string
	for gene := range condition1Data {
		c1Genes = append(c1Genes, gene)
	}
	for gene := range condition2Data {
		c2Genes = append(c2Genes, gene)
	}

	// Get null distributions
	c1NullCorrs, c2NullCorrs := generateNullCorrelations(moduleGenes, c1Genes, c2Genes,
		condition1Data, condition2Data, numPermutations)

	// Create plots for each condition using the function from plotDistributions.go
	plotConditionDistribution(moduleName, "condition1", actualC1Corrs, c1NullCorrs)
	plotConditionDistribution(moduleName, "condition2", actualC2Corrs, c2NullCorrs)
}
