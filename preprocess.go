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
	"math"
	"os"
	"sort"
	"strconv"

	"gonum.org/v1/gonum/mat"
)

// DataWithGenes holds both the expression data matrix and gene IDs
type DataWithGenes struct {
	Data    *mat.Dense
	GeneIDs []string
}

// ReadData reads and parses the GDS2901.soft file
func ReadData(filePath string) (*DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	reader.FieldsPerRecord = -1 // Allow variable number of fields
	reader.LazyQuotes = true    // Be more permissive with quotes

	var geneIDs []string
	var dataRows [][]float64
	var dataStarted bool
	var numCols int

	// Read the file line by line
	for {
		record, err := reader.Read()
		if err != nil {
			break // End of file or error
		}

		// Skip empty lines
		if len(record) == 0 {
			continue
		}

		// Look for the start of data
		if record[0] == "!dataset_table_begin" {
			dataStarted = true
			// Read the header in the next iteration
			continue
		}

		// If we found the end of data, break
		if record[0] == "!dataset_table_end" {
			break
		}

		// Skip until we find the start of data
		if !dataStarted {
			continue
		}

		// If this is the first line after data start, it's the header
		if numCols == 0 {
			numCols = len(record) - 2 // Subtract ID and description columns
			continue
		}

		// Process data rows
		if len(record) >= 3 { // Ensure we have at least ID, description, and one data point
			geneIDs = append(geneIDs, record[0])
			rowData := make([]float64, numCols)

			for j := 2; j < len(record) && j-2 < numCols; j++ {
				val, err := strconv.ParseFloat(record[j], 64)
				if err != nil {
					val = 0 // Use 0 for invalid values
				}
				rowData[j-2] = val
			}
			dataRows = append(dataRows, rowData)
		}
	}

	// Verify we have data
	if len(dataRows) == 0 {
		return nil, fmt.Errorf("no valid data found in file")
	}

	// Create matrix from data
	matrix := mat.NewDense(len(dataRows), numCols, nil)
	for i, row := range dataRows {
		matrix.SetRow(i, row)
	}

	return &DataWithGenes{
		Data:    matrix,
		GeneIDs: geneIDs,
	}, nil
}

// ReadGolubData reads and parses the Golub data file
func ReadGolubData(filePath string) (*DataWithGenes, *DataWithGenes, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'

	// Read all records
	records, err := reader.ReadAll()
	if err != nil {
		return nil, nil, fmt.Errorf("error reading file: %v", err)
	}

	// Skip header row
	if len(records) < 2 {
		return nil, nil, fmt.Errorf("file does not contain enough data")
	}
	records = records[1:] // Skip header row

	// Process data rows
	var geneIDs []string
	var allRows [][]float64
	var amlRows [][]float64

	for _, record := range records {
		if len(record) < 39 { // Need at least 38 columns plus gene ID
			continue
		}

		geneIDs = append(geneIDs, record[0])

		// Process ALL samples (columns 1-27)
		allData := make([]float64, 27)
		for j := 1; j < 28; j++ {
			val, err := strconv.ParseFloat(record[j], 64)
			if err != nil {
				val = 0 // Use 0 for invalid values
			}
			allData[j-1] = val
		}
		allRows = append(allRows, allData)

		// Process AML samples (columns 28-38)
		amlData := make([]float64, 11)
		for j := 28; j < 39; j++ {
			val, err := strconv.ParseFloat(record[j], 64)
			if err != nil {
				val = 0 // Use 0 for invalid values
			}
			amlData[j-28] = val
		}
		amlRows = append(amlRows, amlData)
	}

	// Create matrix from data
	allMatrix := mat.NewDense(len(allRows), 27, nil)
	for i, row := range allRows {
		allMatrix.SetRow(i, row)
	}

	amlMatrix := mat.NewDense(len(amlRows), 11, nil)
	for i, row := range amlRows {
		amlMatrix.SetRow(i, row)
	}

	return &DataWithGenes{
			Data:    allMatrix,
			GeneIDs: geneIDs,
		}, &DataWithGenes{
			Data:    amlMatrix,
			GeneIDs: geneIDs,
		}, nil
}

// Apply log2 transformation
func applyLog2(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	result := mat.NewDense(rows, cols, nil)
	result.Apply(func(i, j int, v float64) float64 {
		return math.Log2(v + 1)
	}, data)
	return result
}

// NormalizeQuantiles performs quantile normalization
func NormalizeQuantiles(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	result := mat.NewDense(rows, cols, nil)

	// For each column
	for j := 0; j < cols; j++ {
		// Get column data
		col := mat.Col(nil, j, data)

		// Sort values
		sorted := make([]float64, len(col))
		copy(sorted, col)
		sort.Float64s(sorted)

		// Create rank mapping
		rankMap := make(map[float64]float64)
		for i, v := range sorted {
			rankMap[v] = float64(i)
		}

		// Apply normalization
		for i := 0; i < rows; i++ {
			val := data.At(i, j)
			rank := rankMap[val]
			normalizedVal := sorted[int(rank)]
			result.Set(i, j, normalizedVal)
		}
	}

	return result
}

// Extract Eker rat samples for different conditions
func ExtractEkerSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	ekerCols := makeRange(0, 36)
	result := mat.NewDense(rows, len(ekerCols), nil)
	for j, col := range ekerCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

// Extract wild type samples for different conditions
func ExtractWildSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	wildCols := makeRange(36, 72)
	result := mat.NewDense(rows, len(wildCols), nil)
	for j, col := range wildCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

func ExtractALLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	allCols := makeRange(0, 27)
	result := mat.NewDense(rows, len(allCols), nil)
	for j, col := range allCols {
		colData := mat.Col(nil, col, data)
		for i := 0; i < rows; i++ {
			result.Set(i, j, colData[i])
		}
	}
	return result
}

func ExtractAMLSamples(data *mat.Dense) *mat.Dense {
	rows, _ := data.Dims()
	amlCols := 11 // Number of AML samples (columns 28-38 inclusive)
	result := mat.NewDense(rows, amlCols, nil)

	// Copy AML samples (columns 28-38) to the result matrix
	for j := 0; j < amlCols; j++ {
		col := mat.Col(nil, j+27, data) // Start from index 27 (column 28)
		result.SetCol(j, col)
	}

	return result
}

// Helper functions
func removeRow(data *mat.Dense, rowIndex int) *mat.Dense {
	rows, cols := data.Dims()
	if rowIndex < 0 || rowIndex >= rows {
		return data
	}

	newData := mat.NewDense(rows-1, cols, nil)
	for i := 0; i < rowIndex; i++ {
		row := mat.Row(nil, i, data)
		newData.SetRow(i, row)
	}
	for i := rowIndex + 1; i < rows; i++ {
		row := mat.Row(nil, i, data)
		newData.SetRow(i-1, row)
	}
	return newData
}

func makeRange(min, max int) []int {
	a := make([]int, max-min)
	for i := range a {
		a[i] = min + i
	}
	return a
}

func saveToCSV(data *mat.Dense, geneIDs []string, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error creating file: %v", err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	rows, cols := data.Dims()
	for i := 0; i < rows; i++ {
		row := make([]string, cols+1)
		row[0] = geneIDs[i]
		for j := 0; j < cols; j++ {
			row[j+1] = strconv.FormatFloat(data.At(i, j), 'f', -1, 64)
		}
		if err := writer.Write(row); err != nil {
			return fmt.Errorf("error writing row: %v", err)
		}
	}

	return nil
}

// CleanData removes or replaces invalid values in the matrix
func CleanData(data *mat.Dense) *mat.Dense {
	rows, cols := data.Dims()
	cleaned := mat.NewDense(rows, cols, nil)

	for i := 0; i < rows; i++ {
		rowData := mat.Row(nil, i, data)
		validValues := make([]float64, 0)

		// Collect valid values from the row
		for _, val := range rowData {
			if !math.IsNaN(val) && !math.IsInf(val, 0) {
				validValues = append(validValues, val)
			}
		}

		// Calculate row mean from valid values
		var rowMean float64
		if len(validValues) > 0 {
			sum := 0.0
			for _, val := range validValues {
				sum += val
			}
			rowMean = sum / float64(len(validValues))
		}

		// Replace invalid values with row mean
		for j := 0; j < cols; j++ {
			val := data.At(i, j)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				cleaned.Set(i, j, rowMean)
			} else {
				cleaned.Set(i, j, val)
			}
		}
	}

	return cleaned
}

// writeOutput saves the processed data to both diffcoex and coxpress directories
func writeOutput(d *DataWithGenes, filename string) error {
	// Write to diffcoex directory
	if err := saveToCSV(d.Data, d.GeneIDs, "output/diffcoex/"+filename); err != nil {
		return fmt.Errorf("error saving to diffcoex: %v", err)
	}

	// Write to coxpress directory
	if err := saveToCSV(d.Data, d.GeneIDs, "output/coxpress/"+filename); err != nil {
		return fmt.Errorf("error saving to coxpress: %v", err)
	}

	return nil
}

func main() {
	if len(os.Args) != 3 {
		fmt.Println("Usage: ./preprocess <dataset_type> <file_path>")
		fmt.Println("dataset_type: 'rat' or 'golub'")
		os.Exit(1)
	}

	datasetType := os.Args[1]
	filePath := os.Args[2]

	// Create output directories if they don't exist
	for _, dir := range []string{"output/diffcoex", "output/coxpress"} {
		if err := os.MkdirAll(dir, 0755); err != nil {
			log.Fatalf("Error creating output directory %s: %v", dir, err)
		}
	}

	// Process data using both methods
	switch datasetType {
	case "rat":
		if err := processRatData(filePath); err != nil {
			log.Fatalf("Error processing rat data: %v", err)
		}
		fmt.Println("Rat data processing complete! Files saved:")
		fmt.Println("- output/diffcoex/rat_eker_mutants.csv")
		fmt.Println("- output/diffcoex/rat_wild_types.csv")
		fmt.Println("- output/coxpress/rat_eker_mutants.csv")
		fmt.Println("- output/coxpress/rat_wild_types.csv")

	case "golub":
		if err := processGolubData(filePath); err != nil {
			log.Fatalf("Error processing Golub data: %v", err)
		}
		fmt.Println("Golub data processing complete! Files saved:")
		fmt.Println("- output/diffcoex/golub_ALL_samples.csv")
		fmt.Println("- output/diffcoex/golub_AML_samples.csv")
		fmt.Println("- output/coxpress/golub_ALL_samples.csv")
		fmt.Println("- output/coxpress/golub_AML_samples.csv")

	default:
		log.Fatalf("Unknown dataset type: %s. Use 'rat' or 'golub'", datasetType)
	}
}

func processRatData(filePath string) error {

	// Read data
	dataWithGenes, err := ReadData(filePath)
	if err != nil {
		return fmt.Errorf("error reading rat data: %v", err)
	}

	// DiffCoEx preprocessing
	{
		// Create a deep copy for DiffCoEx processing
		diffCoExData := copyDataWithGenes(dataWithGenes)

		// Remove last row and probeset 2475
		diffCoExData.Data = removeRow(diffCoExData.Data, len(diffCoExData.GeneIDs)-1)
		diffCoExData.GeneIDs = append(diffCoExData.GeneIDs[:2474], diffCoExData.GeneIDs[2475:]...)
		diffCoExData.Data = removeRow(diffCoExData.Data, 2474)

		// Process data
		logData := applyLog2(diffCoExData.Data)
		normData := NormalizeQuantiles(logData)

		// Extract conditions
		ekerMutants := ExtractEkerSamples(normData)
		wildTypes := ExtractWildSamples(normData)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(ekerMutants, diffCoExData.GeneIDs, "output/diffcoex/rat_eker_mutants.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx Eker mutants: %v", err)
		}
		if err := saveToCSV(wildTypes, diffCoExData.GeneIDs, "output/diffcoex/rat_wild_types.csv"); err != nil {
			return fmt.Errorf("error saving DiffCoEx wild types: %v", err)
		}
	}

	// coXpress preprocessing
	{
		// Extract conditions directly from raw data
		ekerMutants := ExtractEkerSamples(dataWithGenes.Data)
		wildTypes := ExtractWildSamples(dataWithGenes.Data)

		// Save with gene IDs using descriptive filenames
		if err := saveToCSV(ekerMutants, dataWithGenes.GeneIDs, "output/coxpress/rat_eker_mutants.csv"); err != nil {
			return fmt.Errorf("error saving coXpress Eker mutants: %v", err)
		}
		if err := saveToCSV(wildTypes, dataWithGenes.GeneIDs, "output/coxpress/rat_wild_types.csv"); err != nil {
			return fmt.Errorf("error saving coXpress wild types: %v", err)
		}
	}

	return nil
}

func processGolubData(filePath string) error {

	// Read and split the data
	allData, amlData, err := ReadGolubData(filePath)
	if err != nil {
		return fmt.Errorf("error reading Golub data: %v", err)
	}

	// Save ALL samples
	if err := writeOutput(allData, "golub_ALL_samples.csv"); err != nil {
		return fmt.Errorf("error saving ALL samples: %v", err)
	}

	// Save AML samples
	if err := writeOutput(amlData, "golub_AML_samples.csv"); err != nil {
		return fmt.Errorf("error saving AML samples: %v", err)
	}

	return nil
}

func copyDataWithGenes(d *DataWithGenes) *DataWithGenes {
	newData := mat.NewDense(d.Data.RawMatrix().Rows, d.Data.RawMatrix().Cols, nil)
	newData.Copy(d.Data)

	newGeneIDs := make([]string, len(d.GeneIDs))
	copy(newGeneIDs, d.GeneIDs)

	return &DataWithGenes{
		Data:    newData,
		GeneIDs: newGeneIDs,
	}
}
