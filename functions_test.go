// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"bufio"
	"fmt"
	"io/fs"
	"io/ioutil"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"encoding/csv"

	"math"

	"gonum.org/v1/gonum/mat"
)

// CopyDataWithGenesTest holds the test input and expected output
type CopyDataWithGenesTest struct {
	input  *DataWithGenes
	result *DataWithGenes
}

// GolubDataTest holds the test input and expected outputs
type GolubDataTest struct {
	inputFile   string
	expectedALL *DataWithGenes
	expectedAML *DataWithGenes
}

// RatDataTest holds the test input and expected outputs
type RatDataTest struct {
	inputFile            string
	expectedEkerDiffcoex *DataWithGenes
	expectedWildDiffcoex *DataWithGenes
	expectedEkerCoxpress *DataWithGenes
	expectedWildCoxpress *DataWithGenes
}

// Add this to your existing test_types.go
type CleanDataTest struct {
	name     string
	input    *DataWithGenes
	expected *DataWithGenes
}

// Add this to your existing test types at the top of the file
type ExtractAMLTest struct {
	input    *mat.Dense
	expected *mat.Dense
}

// Add this to your existing test types at the top of the file
type ExtractALLTest struct {
	input    *mat.Dense
	expected *mat.Dense
}

// Add this to your existing test types at the top of the file
type ExtractWildTest struct {
	input    *mat.Dense
	expected *mat.Dense
}

// Add this to your existing test types at the top of the file
type ExtractEkerTest struct {
	input    *mat.Dense
	expected *mat.Dense
}

// Add this struct with the other test types at the top of the file
type NormalizeQuantilesTest struct {
	input    *mat.Dense
	expected *mat.Dense
}

// Add this struct with the other test types at the top of the file
type ApplyLog2Test struct {
	input    *mat.Dense
	expected *mat.Dense
}

// TestCopyDataWithGenes tests the copyDataWithGenes function
func TestCopyDataWithGenes(t *testing.T) {
	tests := ReadCopyDataWithGenesTests("tests/CopyDataWithGenes/")
	for _, test := range tests {
		// Run the test
		result := copyDataWithGenes(test.input)

		// Check matrix dimensions
		r1, c1 := result.Data.Dims()
		r2, c2 := test.result.Data.Dims()
		if r1 != r2 || c1 != c2 {
			t.Errorf("Matrix dimensions mismatch: got (%d,%d), want (%d,%d)", r1, c1, r2, c2)
			continue
		}

		// Check matrix data
		if !mat.Equal(result.Data, test.result.Data) {
			t.Errorf("Matrix data mismatch")
		}

		// Check gene IDs
		if len(result.GeneIDs) != len(test.result.GeneIDs) {
			t.Errorf("GeneIDs length mismatch: got %v, want %v", len(result.GeneIDs), len(test.result.GeneIDs))
			continue
		}
		for i := range result.GeneIDs {
			if result.GeneIDs[i] != test.result.GeneIDs[i] {
				t.Errorf("GeneIDs mismatch at index %d: got %v, want %v", i, result.GeneIDs[i], test.result.GeneIDs[i])
			}
		}
	}
}

// ReadCopyDataWithGenesTests reads test cases from files
func ReadCopyDataWithGenesTests(directory string) []CopyDataWithGenesTest {
	files, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	var tests []CopyDataWithGenesTest
	for _, file := range files {
		inputPath := directory + "/input/" + file.Name()
		outputPath := directory + "/output/output_" + strings.TrimPrefix(file.Name(), "input_")

		test := CopyDataWithGenesTest{
			input:  ReadDataWithGenesFromFile(inputPath),
			result: ReadDataWithGenesFromFile(outputPath),
		}
		tests = append(tests, test)
	}
	return tests
}

// ReadDataWithGenesFromFile reads a DataWithGenes structure from a file (original format)
func ReadDataWithGenesFromFile(file string) *DataWithGenes {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	// Read dimensions
	scanner.Scan()
	dims := strings.Split(strings.TrimSpace(scanner.Text()), " ")
	rows, _ := strconv.Atoi(dims[0])
	cols, _ := strconv.Atoi(dims[1])

	// Read gene IDs
	scanner.Scan()
	geneIDs := strings.Split(strings.TrimSpace(scanner.Text()), " ")

	// Read data
	data := mat.NewDense(rows, cols, nil)
	i := 0
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		values := strings.Split(line, " ")
		row := make([]float64, len(values))
		for j, val := range values {
			row[j], _ = strconv.ParseFloat(val, 64)
		}
		data.SetRow(i, row)
		i++
	}

	return &DataWithGenes{
		Data:    data,
		GeneIDs: geneIDs,
	}
}

// ReadProcessedDataFromFile reads a DataWithGenes structure from a processed file (CSV format)
func ReadProcessedDataFromFile(file string) *DataWithGenes {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	// Skip header line
	scanner.Scan()

	// Read data rows and gene IDs
	var geneIDs []string
	var rows [][]float64

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		values := strings.Split(line, ",")

		// First column is gene ID
		geneIDs = append(geneIDs, values[0])

		// Rest are data values
		row := make([]float64, len(values)-1)
		for j := 1; j < len(values); j++ {
			val, err := strconv.ParseFloat(values[j], 64)
			if err != nil {
				val = 0
			}
			row[j-1] = val
		}
		rows = append(rows, row)
	}

	// Create matrix
	if len(rows) == 0 || len(rows[0]) == 0 {
		panic("Empty data")
	}

	numGenes := len(rows)
	numSamples := len(rows[0])
	data := mat.NewDense(numGenes, numSamples, nil)

	// Fill the matrix
	for i := 0; i < numGenes; i++ {
		data.SetRow(i, rows[i])
	}

	return &DataWithGenes{
		Data:    data,
		GeneIDs: geneIDs,
	}
}

// TestProcessGolubData tests the processGolubData function
func TestProcessGolubData(t *testing.T) {
	tests := []struct {
		inputFile   string
		expectedALL *DataWithGenes
		expectedAML *DataWithGenes
	}{
		{
			inputFile:   "tests/ProcessGolubData/input/input_0.txt",
			expectedALL: ReadProcessedDataFromFile("tests/ProcessGolubData/output/all_0.txt"),
			expectedAML: ReadProcessedDataFromFile("tests/ProcessGolubData/output/aml_0.txt"),
		},
		// Add more test cases as needed
	}

	for i, test := range tests {
		// Create temporary output directories for the test
		tmpDir := fmt.Sprintf("test_output_%d", i)
		for _, dir := range []string{tmpDir + "/diffcoex", tmpDir + "/coxpress"} {
			if err := os.MkdirAll(dir, 0755); err != nil {
				t.Fatalf("Failed to create test directory: %v", err)
			}
		}

		// Temporarily redirect output directory
		os.Rename("output", "output_backup")
		os.Symlink(tmpDir, "output")

		// Run the test
		err := processGolubData(test.inputFile)
		if err != nil {
			t.Errorf("Test %d failed with error: %v", i, err)
			continue
		}

		// Verify ALL samples
		allData := ReadDataWithGenesFromCSV("output/diffcoex/golub_ALL_samples.csv")

		// Debug print
		t.Logf("ALL Data Comparison for test %d:", i)
		r1, c1 := allData.Data.Dims()
		r2, c2 := test.expectedALL.Data.Dims()
		t.Logf("Dimensions - Got: (%d,%d), Expected: (%d,%d)", r1, c1, r2, c2)

		// Print first few values of both matrices
		t.Log("Actual ALL data (first row):")
		if r1 > 0 {
			row := mat.Row(nil, 0, allData.Data)
			t.Logf("%v", row)
		}
		t.Log("Expected ALL data (first row):")
		if r2 > 0 {
			row := mat.Row(nil, 0, test.expectedALL.Data)
			t.Logf("%v", row)
		}

		if !compareDataWithGenes(allData, test.expectedALL, t) {
			t.Errorf("Test %d: ALL samples mismatch", i)
		}

		// Verify AML samples
		amlData := ReadDataWithGenesFromCSV("output/diffcoex/golub_AML_samples.csv")

		// Debug print
		t.Logf("AML Data Comparison for test %d:", i)
		r1, c1 = amlData.Data.Dims()
		r2, c2 = test.expectedAML.Data.Dims()
		t.Logf("Dimensions - Got: (%d,%d), Expected: (%d,%d)", r1, c1, r2, c2)

		// Print first few values of both matrices
		t.Log("Actual AML data (first row):")
		if r1 > 0 {
			row := mat.Row(nil, 0, amlData.Data)
			t.Logf("%v", row)
		}
		t.Log("Expected AML data (first row):")
		if r2 > 0 {
			row := mat.Row(nil, 0, test.expectedAML.Data)
			t.Logf("%v", row)
		}

		if !compareDataWithGenes(amlData, test.expectedAML, t) {
			t.Errorf("Test %d: AML samples mismatch", i)
		}

		// Clean up
		os.Remove("output")
		os.Rename("output_backup", "output")
		os.RemoveAll(tmpDir)
	}
}

// ReadGolubDataTests reads test cases from files
func ReadGolubDataTests(directory string) []GolubDataTest {
	files, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	var tests []GolubDataTest
	for _, file := range files {
		inputPath := directory + "/input/" + file.Name()
		allOutputPath := directory + "/output/all_" + strings.TrimPrefix(file.Name(), "input_")
		amlOutputPath := directory + "/output/aml_" + strings.TrimPrefix(file.Name(), "input_")

		test := GolubDataTest{
			inputFile:   inputPath,
			expectedALL: ReadDataWithGenesFromCSV(allOutputPath),
			expectedAML: ReadDataWithGenesFromCSV(amlOutputPath),
		}
		tests = append(tests, test)
	}
	return tests
}

// ReadDataWithGenesFromCSV reads a DataWithGenes structure from a CSV file
func ReadDataWithGenesFromCSV(file string) *DataWithGenes {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	reader := csv.NewReader(f)
	records, err := reader.ReadAll()
	if err != nil {
		panic(err)
	}

	if len(records) == 0 {
		return &DataWithGenes{
			Data:    mat.NewDense(0, 0, nil),
			GeneIDs: []string{},
		}
	}

	rows := len(records)
	cols := len(records[0]) - 1 // First column is gene ID
	data := mat.NewDense(rows, cols, nil)
	geneIDs := make([]string, rows)

	for i, record := range records {
		geneIDs[i] = record[0]
		for j := 1; j < len(record); j++ {
			val, _ := strconv.ParseFloat(record[j], 64)
			data.Set(i, j-1, val)
		}
	}

	return &DataWithGenes{
		Data:    data,
		GeneIDs: geneIDs,
	}
}

// compareDataWithGenes compares two DataWithGenes structures
func compareDataWithGenes(actual, expected *DataWithGenes, t *testing.T) bool {
	if actual == nil || expected == nil {
		t.Error("One of the DataWithGenes structs is nil")
		return false
	}

	// Compare dimensions
	actualRows, actualCols := actual.Data.Dims()
	expectedRows, expectedCols := expected.Data.Dims()
	if actualRows != expectedRows || actualCols != expectedCols {
		t.Errorf("Matrix dimensions mismatch: got (%d,%d), want (%d,%d)",
			actualRows, actualCols, expectedRows, expectedCols)
		return false
	}

	// Compare gene IDs
	if len(actual.GeneIDs) != len(expected.GeneIDs) {
		t.Errorf("GeneIDs length mismatch: got %d, want %d",
			len(actual.GeneIDs), len(expected.GeneIDs))
		return false
	}
	for i := range actual.GeneIDs {
		if actual.GeneIDs[i] != expected.GeneIDs[i] {
			t.Errorf("GeneID mismatch at position %d: got %s, want %s",
				i, actual.GeneIDs[i], expected.GeneIDs[i])
			return false
		}
	}

	// Compare matrix data
	for i := 0; i < actualRows; i++ {
		for j := 0; j < actualCols; j++ {
			if math.Abs(actual.Data.At(i, j)-expected.Data.At(i, j)) > 1e-10 {
				t.Errorf("Matrix data mismatch at (%d,%d): got %v, want %v",
					i, j, actual.Data.At(i, j), expected.Data.At(i, j))
				return false
			}
		}
	}

	return true
}

// generateTestInputFile creates a test input file with the required number of rows
func generateTestInputFile(filename string) error {
	dir := filepath.Dir(filename)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return fmt.Errorf("failed to create directory: %v", err)
	}

	f, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("failed to create file: %v", err)
	}
	defer f.Close()

	// Write begin marker
	f.WriteString("!dataset_table_begin\n")

	// Write header with 72 samples (36 Eker + 36 Wild)
	f.WriteString("ID\tDescription")
	// Add 36 Eker sample columns
	for i := 1; i <= 36; i++ {
		f.WriteString(fmt.Sprintf("\tEker%d", i))
	}
	// Add 36 Wild sample columns
	for i := 1; i <= 36; i++ {
		f.WriteString(fmt.Sprintf("\tWild%d", i))
	}
	f.WriteString("\n")

	// Write exactly 2476 rows
	for i := 1; i <= 2476; i++ {
		// ID and Description
		row := fmt.Sprintf("GENE%d\tdesc%d", i, i)

		// Eker samples (first 36 columns)
		for j := 0; j < 36; j++ {
			row += fmt.Sprintf("\t%f", float64(i+j))
		}

		// Wild samples (next 36 columns)
		for j := 36; j < 72; j++ {
			row += fmt.Sprintf("\t%f", float64(i+j))
		}

		row += "\n"
		f.WriteString(row)
	}

	// Write end marker without trailing newline
	f.WriteString("!dataset_table_end")
	return nil
}

// TestProcessRatData tests the processRatData function
func TestProcessRatData(t *testing.T) {
	// Create test directory if it doesn't exist
	testDir := "tests/ProcessRatData/input"
	if err := os.MkdirAll(testDir, 0755); err != nil {
		t.Fatalf("Failed to create test directory: %v", err)
	}

	// Generate test input file
	inputFile := filepath.Join(testDir, "input_0.txt")
	if err := generateTestInputFile(inputFile); err != nil {
		t.Fatalf("Failed to generate test input file: %v", err)
	}

	// Verify file was created and has correct number of rows
	data, err := os.ReadFile(inputFile)
	if err != nil {
		t.Fatalf("Failed to read generated file: %v", err)
	}
	lines := strings.Split(string(data), "\n")
	expectedLines := 2476 + 2 + 1 // data rows + begin marker + header
	if len(lines) != expectedLines {
		t.Fatalf("Generated file has %d lines, expected %d", len(lines), expectedLines)
	}

	// Create temporary output directories
	tmpDir := "test_output_0"
	for _, dir := range []string{tmpDir + "/diffcoex", tmpDir + "/coxpress"} {
		if err := os.MkdirAll(dir, 0755); err != nil {
			t.Fatalf("Failed to create test directory: %v", err)
		}
	}

	// Backup and redirect output directory
	if err := os.Rename("output", "output_backup"); err != nil && !os.IsNotExist(err) {
		t.Fatalf("Failed to backup output directory: %v", err)
	}
	if err := os.Symlink(tmpDir, "output"); err != nil {
		t.Fatalf("Failed to create output symlink: %v", err)
	}

	// Run the test
	if err := processRatData(inputFile); err != nil {
		t.Errorf("processRatData failed: %v", err)
	}

	// Clean up
	os.Remove("output")
	if _, err := os.Stat("output_backup"); err == nil {
		os.Rename("output_backup", "output")
	}
	os.RemoveAll(tmpDir)
}

// ReadRatDataTests reads test cases from the specified directory
func ReadRatDataTests(directory string) []RatDataTest {
	files, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	var tests []RatDataTest
	for _, file := range files {
		inputPath := directory + "/input/" + file.Name()
		baseNum := strings.TrimPrefix(file.Name(), "input_")
		baseNum = strings.TrimSuffix(baseNum, ".txt")

		test := RatDataTest{
			inputFile:            inputPath,
			expectedEkerDiffcoex: ReadDataWithGenesFromCSV(directory + "/output/eker_diffcoex_" + baseNum + ".txt"),
			expectedWildDiffcoex: ReadDataWithGenesFromCSV(directory + "/output/wild_diffcoex_" + baseNum + ".txt"),
			expectedEkerCoxpress: ReadDataWithGenesFromCSV(directory + "/output/eker_coxpress_" + baseNum + ".txt"),
			expectedWildCoxpress: ReadDataWithGenesFromCSV(directory + "/output/wild_coxpress_" + baseNum + ".txt"),
		}
		tests = append(tests, test)
	}
	return tests
}

// TestCleanData tests the CleanData function
func TestCleanData(t *testing.T) {
	tests := ReadCleanDataTests("tests/CleanData")

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			// Run the test
			result := CleanData(test.input.Data)

			// Compare dimensions
			r1, c1 := result.Dims()
			r2, c2 := test.expected.Data.Dims()
			if r1 != r2 || c1 != c2 {
				t.Errorf("Matrix dimensions mismatch: got (%d,%d), want (%d,%d)", r1, c1, r2, c2)
				return
			}

			// Compare values with rounding
			for i := 0; i < r1; i++ {
				for j := 0; j < c1; j++ {
					got := result.At(i, j)
					want := test.expected.Data.At(i, j)

					if math.IsNaN(got) && math.IsNaN(want) {
						continue // Both NaN is okay
					}

					// Round both values to 3 decimal places
					gotRounded := roundFloat(got, 3)
					wantRounded := roundFloat(want, 3)

					if gotRounded != wantRounded {
						t.Errorf("Value mismatch at (%d,%d): got %v, want %v", i, j, gotRounded, wantRounded)
					}
				}
			}
		})
	}
}

// ReadCleanDataTests reads test cases from the specified directory
func ReadCleanDataTests(directory string) []CleanDataTest {
	// Read input files
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]CleanDataTest, numFiles)
	for i, inputFile := range inputFiles {
		// Read test name from filename
		tests[i].name = strings.TrimSuffix(inputFile.Name(), ".txt")

		// Read input matrix and gene IDs
		tests[i].input = ReadMatrixFromFile(directory + "/input/" + inputFile.Name())
	}

	// Read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		// Read expected output matrix and gene IDs
		tests[i].expected = ReadMatrixFromFile(directory + "/output/" + outputFile.Name())
	}

	return tests
}

// ReadMatrixFromFile reads a matrix and gene IDs from a file
func ReadMatrixFromFile(file string) *DataWithGenes {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)

	// Read gene IDs from first line
	scanner.Scan()
	geneIDs := strings.Split(scanner.Text(), "\t")

	// Read matrix dimensions from second line
	scanner.Scan()
	dims := strings.Split(scanner.Text(), "\t")
	rows, _ := strconv.Atoi(dims[0])
	cols, _ := strconv.Atoi(dims[1])

	// Read matrix data
	data := make([]float64, rows*cols)
	idx := 0
	for scanner.Scan() {
		values := strings.Split(scanner.Text(), "\t")
		for _, val := range values {
			if val == "NaN" {
				data[idx] = math.NaN()
			} else {
				f, _ := strconv.ParseFloat(val, 64)
				data[idx] = f
			}
			idx++
		}
	}

	return &DataWithGenes{
		Data:    mat.NewDense(rows, cols, data),
		GeneIDs: geneIDs,
	}
}

// ReadDirectory reads in a directory and returns a slice of fs.DirEntry objects containing file info for the directory
func ReadDirectory(dir string) []fs.DirEntry {
	// read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// roundFloat rounds a float64 to a specified number of decimal places
func roundFloat(val float64, precision uint) float64 {
	ratio := math.Pow(10, float64(precision))
	return math.Round(val*ratio) / ratio
}

// TestExtractAMLSamples tests the ExtractAMLSamples function
func TestExtractAMLSamples(t *testing.T) {
	// Read in all tests from the tests/ExtractAML directory and run them
	tests := ReadExtractAMLTests("tests/ExtractAML/")
	for _, test := range tests {
		// Run the test
		result := ExtractAMLSamples(test.input)

		// Compare dimensions
		r1, c1 := result.Dims()
		r2, c2 := test.expected.Dims()
		if r1 != r2 || c1 != c2 {
			t.Errorf("Matrix dimensions mismatch: got (%d,%d), want (%d,%d)", r1, c1, r2, c2)
			continue
		}

		// Compare values
		for i := 0; i < r1; i++ {
			for j := 0; j < c1; j++ {
				got := roundFloat(result.At(i, j), 4)
				want := roundFloat(test.expected.At(i, j), 4)
				if got != want {
					t.Errorf("Value mismatch at (%d,%d): got %v, want %v", i, j, got, want)
				}
			}
		}
	}
}

// ReadExtractAMLTests reads test cases from the specified directory
func ReadExtractAMLTests(directory string) []ExtractAMLTest {
	// Read input files
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]ExtractAMLTest, numFiles)

	// Read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadAMLFromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadAMLFromFile(outputPath)
	}

	return tests
}

// ReadMatrixFromFile reads a matrix from a file
func ReadAMLFromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	dims := strings.Split(scanner.Text(), "\t")
	rows, _ := strconv.Atoi(dims[0])
	cols, _ := strconv.Atoi(dims[1])

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Split(scanner.Text(), "\t")
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}

func TestMakeRange(t *testing.T) {
	// Read test cases from tests/MakeRange directory
	files, err := os.ReadDir("tests/makeRange/input")
	if err != nil {
		t.Fatalf("Failed to read test directory: %v", err)
	}

	for _, file := range files {
		inputPath := filepath.Join("tests/makeRange/input", file.Name())
		outputPath := filepath.Join("tests/makeRange/output",
			strings.Replace(file.Name(), "input_", "output_", 1))

		// Read input numbers (min and max)
		min, max, err := readMinMax(inputPath)
		if err != nil {
			t.Errorf("Error reading input file %s: %v", file.Name(), err)
			continue
		}

		// Read expected output
		expected, err := readIntSlice(outputPath)
		if err != nil {
			t.Errorf("Error reading output file %s: %v", file.Name(), err)
			continue
		}

		// Calculate result
		result := makeRange(min, max)

		// Compare with expected output
		if !sliceEqual(result, expected) {
			t.Errorf("Test case %s failed: got %v, want %v", file.Name(), result, expected)
		}
	}
}

func readMinMax(path string) (min, max int, err error) {
	file, err := os.Open(path)
	if err != nil {
		return 0, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if scanner.Scan() {
		min, err = strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return 0, 0, err
		}
	}
	if scanner.Scan() {
		max, err = strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return 0, 0, err
		}
	}
	return min, max, scanner.Err()
}

func readIntSlice(path string) ([]int, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var numbers []int
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		num, err := strconv.Atoi(strings.TrimSpace(scanner.Text()))
		if err != nil {
			return nil, err
		}
		numbers = append(numbers, num)
	}
	return numbers, scanner.Err()
}

func sliceEqual(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// ReadMatrix reads a matrix from a given file path.
// The first line should contain the row index to remove, followed by the matrix data.
func ReadMatrix(filePath string) (*mat.Dense, int, error) {
	inputData, err := ioutil.ReadFile(filePath)
	if err != nil {
		return nil, 0, err
	}

	lines := strings.Split(string(inputData), "\n")
	if len(lines) < 2 {
		return nil, 0, fmt.Errorf("Input file must contain at least one row of data and one row index")
	}

	// Read the row to remove
	rowToRemove, err := strconv.Atoi(strings.TrimSpace(lines[0]))
	if err != nil {
		return nil, 0, fmt.Errorf("Invalid row index: %v", err)
	}

	// Parse the matrix data
	var dataRows [][]float64
	for _, line := range lines[1:] {
		if line == "" {
			continue
		}
		var row []float64
		for _, val := range strings.Fields(line) {
			num, _ := strconv.ParseFloat(val, 64)
			row = append(row, num)
		}
		dataRows = append(dataRows, row)
	}

	// Create a matrix from the parsed data
	data := mat.NewDense(len(dataRows), len(dataRows[0]), nil)
	for i, row := range dataRows {
		data.SetRow(i, row)
	}

	return data, rowToRemove, nil
}

// MatrixEqual compares two matrices for equality.
func MatrixEqual(a, b *mat.Dense) bool {
	if a == nil || b == nil {
		return a == b
	}
	aRows, aCols := a.Dims()
	bRows, bCols := b.Dims()
	if aRows != bRows || aCols != bCols {
		return false
	}
	for i := 0; i < aRows; i++ {
		for j := 0; j < aCols; j++ {
			if math.Abs(a.At(i, j)-b.At(i, j)) > 1e-9 { // Use a tolerance for floating-point comparison
				return false
			}
		}
	}
	return true
}

func TestRemoveRow(t *testing.T) {
	// Read test cases from tests/RemoveRow directory
	files, err := os.ReadDir("tests/RemoveRow/input")
	if err != nil {
		t.Fatalf("Failed to read test directory: %v", err)
	}

	for _, file := range files {
		inputPath := filepath.Join("tests/RemoveRow/input", file.Name())
		outputPath := filepath.Join("tests/RemoveRow/output",
			strings.Replace(file.Name(), "input_", "output_", 1))

		// Read input matrix and row to remove
		matrix, rowToRemove, err := ReadMatrix(inputPath)
		if err != nil {
			t.Errorf("Error reading input file %s: %v", file.Name(), err)
			continue
		}

		// Read expected output matrix
		expectedOutput, _, err := ReadMatrix(outputPath)
		if err != nil {
			t.Errorf("Error reading output file %s: %v", file.Name(), err)
			continue
		}

		// Remove the specified row
		result := removeRow(matrix, rowToRemove)

		// Compare with expected output
		if !MatrixEqual(result, expectedOutput) {
			t.Errorf("Test case %s failed: matrices are not equal", file.Name())
		}
	}
}

func TestExtractALLSamples(t *testing.T) {
	// Read test cases from tests/ExtractALL directory
	tests := ReadALLTests("tests/ExtractALL")

	for i, test := range tests {
		result := ExtractALLSamples(test.input)
		if !mat.Equal(result, test.expected) {
			t.Errorf("Test case %d failed: matrices are not equal", i)
		}
	}
}

func ReadALLTests(directory string) []ExtractALLTest {
	// Read input files
	inputFiles, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	// Read output files
	outputFiles, err := os.ReadDir(directory + "/output")
	if err != nil {
		panic(err)
	}

	// Sort files to ensure matching pairs
	sort.Sort(byName(inputFiles))
	sort.Sort(byName(outputFiles))

	tests := make([]ExtractALLTest, len(inputFiles))
	if len(inputFiles) != len(outputFiles) {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadALLFromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadALLFromFile(outputPath)
	}

	return tests
}

// ReadALLFromFile reads a matrix from a file
func ReadALLFromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	dims := strings.Split(scanner.Text(), "\t")
	rows, _ := strconv.Atoi(dims[0])
	cols, _ := strconv.Atoi(dims[1])

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Split(scanner.Text(), "\t")
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}

// byName implements sort.Interface for []os.DirEntry based on Name()
type byName []os.DirEntry

func (f byName) Len() int           { return len(f) }
func (f byName) Less(i, j int) bool { return f[i].Name() < f[j].Name() }
func (f byName) Swap(i, j int)      { f[i], f[j] = f[j], f[i] }

func TestExtractWildSamples(t *testing.T) {
	tests := ReadWildTests("tests/ExtractWild")

	for i, test := range tests {
		result := ExtractWildSamples(test.input)
		if !mat.Equal(result, test.expected) {
			// Add debug printing
			fmt.Printf("\nTest case %d:\n", i)
			fmt.Println("Input matrix:")
			matPrint(test.input)
			fmt.Println("\nExpected matrix:")
			matPrint(test.expected)
			fmt.Println("\nResult matrix:")
			matPrint(result)
			t.Errorf("Test case %d failed: matrices are not equal", i)
		}
	}
}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

func ReadWildTests(directory string) []ExtractWildTest {
	// Read input files
	inputFiles, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	// Read output files
	outputFiles, err := os.ReadDir(directory + "/output")
	if err != nil {
		panic(err)
	}

	// Sort files to ensure matching pairs
	sort.Sort(byName(inputFiles))
	sort.Sort(byName(outputFiles))

	tests := make([]ExtractWildTest, len(inputFiles))
	if len(inputFiles) != len(outputFiles) {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadWildFromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadWildFromFile(outputPath)
	}

	return tests
}

// ReadWildFromFile reads a matrix from a file
func ReadWildFromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	firstLine := scanner.Text()
	dims := strings.Fields(firstLine) // Changed from Split("\t") to Fields()
	if len(dims) != 2 {
		panic(fmt.Sprintf("Invalid dimensions line in file %s: %s", filename, firstLine))
	}

	rows, err := strconv.Atoi(dims[0])
	if err != nil {
		panic(fmt.Sprintf("Invalid row dimension in file %s: %v", filename, err))
	}

	cols, err := strconv.Atoi(dims[1])
	if err != nil {
		panic(fmt.Sprintf("Invalid column dimension in file %s: %v", filename, err))
	}

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Fields(scanner.Text()) // Changed from Split("\t") to Fields()
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}

func TestExtractEkerSamples(t *testing.T) {
	tests := ReadEkerTests("tests/ExtractEker")

	for i, test := range tests {
		result := ExtractEkerSamples(test.input)
		if !mat.Equal(result, test.expected) {
			// Add debug printing
			fmt.Printf("\nTest case %d:\n", i)
			fmt.Println("Input matrix:")
			matPrint(test.input)
			fmt.Println("\nExpected matrix:")
			matPrint(test.expected)
			fmt.Println("\nResult matrix:")
			matPrint(result)
			t.Errorf("Test case %d failed: matrices are not equal", i)
		}
	}
}

func ReadEkerTests(directory string) []ExtractEkerTest {
	// Read input files
	inputFiles, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	// Read output files
	outputFiles, err := os.ReadDir(directory + "/output")
	if err != nil {
		panic(err)
	}

	// Sort files to ensure matching pairs
	sort.Sort(byName(inputFiles))
	sort.Sort(byName(outputFiles))

	tests := make([]ExtractEkerTest, len(inputFiles))
	if len(inputFiles) != len(outputFiles) {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadEkerFromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadEkerFromFile(outputPath)
	}

	return tests
}

// ReadEkerFromFile reads a matrix from a file
func ReadEkerFromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	firstLine := scanner.Text()
	dims := strings.Fields(firstLine) // Changed from Split("\t") to Fields()
	if len(dims) != 2 {
		panic(fmt.Sprintf("Invalid dimensions line in file %s: %s", filename, firstLine))
	}

	rows, err := strconv.Atoi(dims[0])
	if err != nil {
		panic(fmt.Sprintf("Invalid row dimension in file %s: %v", filename, err))
	}

	cols, err := strconv.Atoi(dims[1])
	if err != nil {
		panic(fmt.Sprintf("Invalid column dimension in file %s: %v", filename, err))
	}

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Fields(scanner.Text()) // Changed from Split("\t") to Fields()
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}

func TestNormalizeQuantiles(t *testing.T) {
	tests := ReadNormalizeQuantilesTests("tests/NormalizeQuantiles")

	for i, test := range tests {
		result := NormalizeQuantiles(test.input)

		// Compare dimensions
		r1, c1 := result.Dims()
		r2, c2 := test.expected.Dims()
		if r1 != r2 || c1 != c2 {
			t.Errorf("Test case %d: Matrix dimensions mismatch: got (%d,%d), want (%d,%d)",
				i, r1, c1, r2, c2)
			continue
		}

		// Compare values with tolerance for floating point comparison
		for i := 0; i < r1; i++ {
			for j := 0; j < c1; j++ {
				got := result.At(i, j)
				want := test.expected.At(i, j)
				if math.Abs(got-want) > 1e-10 {
					t.Errorf("Value mismatch at (%d,%d): got %v, want %v", i, j, got, want)

					// Add debug printing
					fmt.Printf("\nTest case %d failed:\n", i)
					fmt.Println("Input matrix:")
					matPrint(test.input)
					fmt.Println("\nExpected matrix:")
					matPrint(test.expected)
					fmt.Println("\nResult matrix:")
					matPrint(result)
					return
				}
			}
		}
	}
}

func ReadNormalizeQuantilesTests(directory string) []NormalizeQuantilesTest {
	// Read input files
	inputFiles, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	// Read output files
	outputFiles, err := os.ReadDir(directory + "/output")
	if err != nil {
		panic(err)
	}

	// Sort files to ensure matching pairs
	sort.Sort(byName(inputFiles))
	sort.Sort(byName(outputFiles))

	tests := make([]NormalizeQuantilesTest, len(inputFiles))
	if len(inputFiles) != len(outputFiles) {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadNormalizeQuantilesFromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadNormalizeQuantilesFromFile(outputPath)
	}

	return tests
}

// ReadNormalizeQuantilesFromFile reads a matrix from a file
func ReadNormalizeQuantilesFromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	firstLine := scanner.Text()
	dims := strings.Fields(firstLine)
	if len(dims) != 2 {
		panic(fmt.Sprintf("Invalid dimensions line in file %s: %s", filename, firstLine))
	}

	rows, err := strconv.Atoi(dims[0])
	if err != nil {
		panic(fmt.Sprintf("Invalid row dimension in file %s: %v", filename, err))
	}

	cols, err := strconv.Atoi(dims[1])
	if err != nil {
		panic(fmt.Sprintf("Invalid column dimension in file %s: %v", filename, err))
	}

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Fields(scanner.Text())
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}

func TestApplyLog2(t *testing.T) {
	tests := ReadApplyLog2Tests("tests/ApplyLog2")

	for i, test := range tests {
		result := applyLog2(test.input)

		// Compare dimensions
		r1, c1 := result.Dims()
		r2, c2 := test.expected.Dims()
		if r1 != r2 || c1 != c2 {
			t.Errorf("Test case %d: Matrix dimensions mismatch: got (%d,%d), want (%d,%d)",
				i, r1, c1, r2, c2)
			continue
		}

		// Compare values with tolerance for floating point comparison
		for i := 0; i < r1; i++ {
			for j := 0; j < c1; j++ {
				got := result.At(i, j)
				want := test.expected.At(i, j)
				if math.Abs(got-want) > 1e-10 {
					t.Errorf("Value mismatch at (%d,%d): got %v, want %v", i, j, got, want)

					// Add debug printing
					fmt.Printf("\nTest case %d failed:\n", i)
					fmt.Println("Input matrix:")
					matPrint(test.input)
					fmt.Println("\nExpected matrix:")
					matPrint(test.expected)
					fmt.Println("\nResult matrix:")
					matPrint(result)
					return
				}
			}
		}
	}
}

func ReadApplyLog2Tests(directory string) []ApplyLog2Test {
	// Read input files
	inputFiles, err := os.ReadDir(directory + "/input")
	if err != nil {
		panic(err)
	}

	// Read output files
	outputFiles, err := os.ReadDir(directory + "/output")
	if err != nil {
		panic(err)
	}

	// Sort files to ensure matching pairs
	sort.Sort(byName(inputFiles))
	sort.Sort(byName(outputFiles))

	tests := make([]ApplyLog2Test, len(inputFiles))
	if len(inputFiles) != len(outputFiles) {
		panic("Error: number of input and output files do not match")
	}

	// Process each test case
	for i := range tests {
		// Read input matrix
		inputPath := directory + "/input/" + inputFiles[i].Name()
		tests[i].input = ReadApplyLog2FromFile(inputPath)

		// Read expected output matrix
		outputPath := directory + "/output/" + outputFiles[i].Name()
		tests[i].expected = ReadApplyLog2FromFile(outputPath)
	}

	return tests
}

// ReadApplyLog2FromFile reads a matrix from a file
func ReadApplyLog2FromFile(filename string) *mat.Dense {
	file, err := os.Open(filename)
	if err != nil {
		panic(fmt.Sprintf("Error opening file %s: %v", filename, err))
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read dimensions from first line
	scanner.Scan()
	firstLine := scanner.Text()
	dims := strings.Fields(firstLine)
	if len(dims) != 2 {
		panic(fmt.Sprintf("Invalid dimensions line in file %s: %s", filename, firstLine))
	}

	rows, err := strconv.Atoi(dims[0])
	if err != nil {
		panic(fmt.Sprintf("Invalid row dimension in file %s: %v", filename, err))
	}

	cols, err := strconv.Atoi(dims[1])
	if err != nil {
		panic(fmt.Sprintf("Invalid column dimension in file %s: %v", filename, err))
	}

	// Read matrix data
	data := make([]float64, 0, rows*cols)
	for scanner.Scan() {
		values := strings.Fields(scanner.Text())
		for _, val := range values {
			f, err := strconv.ParseFloat(val, 64)
			if err != nil {
				panic(fmt.Sprintf("Error parsing float in file %s: %v", filename, err))
			}
			data = append(data, f)
		}
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf("Data length mismatch in file %s: got %d elements, expected %d",
			filename, len(data), rows*cols))
	}

	return mat.NewDense(rows, cols, data)
}
