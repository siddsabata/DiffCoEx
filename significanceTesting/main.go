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
	"log"
	"os"
	"path/filepath"
	"strconv"
)

func main() {
	// ./significanceTesting moduleMap condition1Data condition2Data
	// Check if correct number of arguments are provided
	if len(os.Args) != 4 {
		fmt.Println("Usage: ./significanceTesting moduleMap condition1Data condition2Data")
		fmt.Println("Example: ./significanceTesting data/golub/golub_diffcoex.csv data/golub/aml_samples.csv data/golub/all_samples.csv")
		os.Exit(1)
	}

	// Get file paths from command line arguments
	moduleMapPath := os.Args[1]
	condition1Path := os.Args[2]
	condition2Path := os.Args[3]

	// Create output/sigTesting directory if it doesn't exist
	outputDir := "output/sigTesting"
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

	fmt.Println("Writing null distribution results...")
	writeNullDistributionResults(moduleMap, condition1Data, condition2Data)

	fmt.Println("Writing module correlation results...")
	writeModuleCorrelationResults(moduleMap, condition1Data, condition2Data)

	fmt.Println("Done!")
}

func writeNullDistributionResults(moduleMap map[string]string, condition1Data, condition2Data map[string][]float64) {
	// Use path/filepath.Join for proper path construction
	outputPath := filepath.Join("output", "sigTesting", "null_distribution_results.csv")
	outputFile, err := os.Create(outputPath)
	if err != nil {
		log.Fatal("Cannot create null distribution output file:", err)
	}
	defer outputFile.Close()

	writer := csv.NewWriter(outputFile)
	defer writer.Flush()

	// Write header
	header := []string{"Module", "Size", "C1_T-Stat", "C1_P-Value", "C2_T-Stat", "C2_P-Value"}
	if err := writer.Write(header); err != nil {
		log.Fatal("Error writing header:", err)
	}

	// Analyze each module and write results
	for module := range getUniqueModules(moduleMap) {
		stats := analyzeModuleNullDistribution(module, moduleMap, condition1Data, condition2Data)

		row := []string{
			stats.Name,
			strconv.Itoa(stats.Size),
			strconv.FormatFloat(stats.C1NullTStatistic, 'f', 6, 64),
			strconv.FormatFloat(stats.C1NullPValue, 'f', 6, 64),
			strconv.FormatFloat(stats.C2NullTStatistic, 'f', 6, 64),
			strconv.FormatFloat(stats.C2NullPValue, 'f', 6, 64),
		}

		if err := writer.Write(row); err != nil {
			log.Fatal("Error writing result:", err)
		}
	}
}

func writeModuleCorrelationResults(moduleMap map[string]string, condition1Data, condition2Data map[string][]float64) {
	// Use path/filepath.Join for proper path construction
	outputPath := filepath.Join("output", "sigTesting", "module_correlation_results.csv")
	outputFile, err := os.Create(outputPath)
	if err != nil {
		log.Fatal("Cannot create module correlation output file:", err)
	}
	defer outputFile.Close()

	writer := csv.NewWriter(outputFile)
	defer writer.Flush()

	// Write header
	header := []string{"Module", "Size", "T-Statistic", "P-Value"}
	if err := writer.Write(header); err != nil {
		log.Fatal("Error writing header:", err)
	}

	// Analyze each module and write results
	for module := range getUniqueModules(moduleMap) {
		stats := analyzeModule(module, moduleMap, condition1Data, condition2Data)

		row := []string{
			stats.Name,
			strconv.Itoa(stats.Size),
			strconv.FormatFloat(stats.TStatistic, 'f', 6, 64),
			strconv.FormatFloat(stats.PValue, 'f', 6, 64),
		}

		if err := writer.Write(row); err != nil {
			log.Fatal("Error writing result:", err)
		}
	}
}
