// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"bufio"
	"io/fs"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"testing"
)

type NormalCDFTest struct {
	x        float64
	expected float64
}

type TStatisticTest struct {
	actual    []float64
	null      []float64
	expectedT float64
	expectedP float64
}

type ModuleCorrelationTest struct {
	genes          []string
	expressionData map[string][]float64
	expected       []float64
}

type GeneIndicesTest struct {
	idx       int
	expectedI int
	expectedJ int
}

// TestNormalCDF tests the normalCDF function
func TestNormalCDF(t *testing.T) {
	tests := ReadNormalCDFTests("Tests/NormalCDF")

	for _, test := range tests {
		result := normalCDF(test.x)
		if math.Abs(result-test.expected) > 1e-3 {
			t.Errorf("normalCDF(%v) = %v, want %v", test.x, result, test.expected)
		}
	}
}

// TestCalculateTStatistic tests the calculateTStatistic function
func TestCalculateTStatistic(t *testing.T) {
	tests := ReadTStatisticTests("Tests/TStatistic")

	for i, test := range tests {
		resultT, resultP := calculateTStatistic(test.actual, test.null)
		if math.Abs(resultT-test.expectedT) > 1e-2 || math.Abs(resultP-test.expectedP) > 1e-2 {
			t.Errorf("Test %d: calculateTStatistic() = (%v, %v), want (%v, %v)",
				i, resultT, resultP, test.expectedT, test.expectedP)
		}
	}
}

// TestGetModuleCorrelations tests the getModuleCorrelations function
func TestGetModuleCorrelations(t *testing.T) {
	tests := ReadModuleCorrelationTests("Tests/ModuleCorrelation")

	for i, test := range tests {
		result := getModuleCorrelations(test.genes, test.expressionData)

		// Sort both slices to ensure consistent comparison
		sort.Float64s(result)
		sort.Float64s(test.expected)

		if len(result) != len(test.expected) {
			t.Errorf("Test %d: got %d correlations, want %d", i, len(result), len(test.expected))
			continue
		}

		for j := range result {
			if math.Abs(result[j]-test.expected[j]) > 1e-6 {
				t.Errorf("Test %d: correlation %d = %v, want %v", i, j, result[j], test.expected[j])
			}
		}
	}
}

// TestGetGeneIndices tests the getGeneIndices function
func TestGetGeneIndices(t *testing.T) {
	tests := ReadGeneIndicesTests("Tests/GeneIndices")

	for i, test := range tests {
		resultI, resultJ := getGeneIndices(test.idx)
		if resultI != test.expectedI || resultJ != test.expectedJ {
			t.Errorf("Test %d: getGeneIndices(%v) = (%v, %v), want (%v, %v)",
				i, test.idx, resultI, resultJ, test.expectedI, test.expectedJ)
		}
	}
}

// ReadNormalCDFTests takes in an input directory and returns a slice
// of NormalCDFTest objects
func ReadNormalCDFTests(directory string) []NormalCDFTest {
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]NormalCDFTest, numFiles)
	for i, inputFile := range inputFiles {
		tests[i] = ReadNormalCDFTest(directory + "/input/" + inputFile.Name())
	}

	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		tests[i].expected = ReadFloatFromFile(directory + "/output/" + outputFile.Name())
	}

	return tests
}

// ReadTStatisticTests takes in an input directory and returns a slice
// of TStatisticTest objects
func ReadTStatisticTests(directory string) []TStatisticTest {
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]TStatisticTest, numFiles)
	for i, inputFile := range inputFiles {
		tests[i] = ReadTStatisticTest(directory + "/input/" + inputFile.Name())
	}

	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		tests[i].expectedT, tests[i].expectedP = ReadTwoFloatsFromFile(directory + "/output/" + outputFile.Name())
	}

	return tests
}

// ReadModuleCorrelationTests reads test cases from files
func ReadModuleCorrelationTests(directory string) []ModuleCorrelationTest {
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]ModuleCorrelationTest, numFiles)
	for i, inputFile := range inputFiles {
		tests[i] = ReadModuleCorrelationTest(directory + "/input/" + inputFile.Name())
	}

	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		tests[i].expected = ReadCorrelationsFromFile(directory + "/output/" + outputFile.Name())
	}

	return tests
}

// ReadGeneIndicesTests reads test cases from files
func ReadGeneIndicesTests(directory string) []GeneIndicesTest {
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]GeneIndicesTest, numFiles)
	for i, inputFile := range inputFiles {
		tests[i] = ReadGeneIndicesTest(directory + "/input/" + inputFile.Name())
	}

	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match!")
	}

	for i, outputFile := range outputFiles {
		tests[i].expectedI, tests[i].expectedJ = ReadTwoIntsFromFile(directory + "/output/" + outputFile.Name())
	}

	return tests
}

// ReadNormalCDFTest reads a single NormalCDFTest from a file
func ReadNormalCDFTest(file string) NormalCDFTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var test NormalCDFTest

	scanner.Scan()
	test.x, _ = strconv.ParseFloat(scanner.Text(), 64)

	return test
}

// ReadTStatisticTest reads a single TStatisticTest from a file
func ReadTStatisticTest(file string) TStatisticTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var test TStatisticTest

	// Read actual values
	scanner.Scan()
	actualStr := strings.Split(scanner.Text(), " ")
	test.actual = make([]float64, len(actualStr))
	for i, str := range actualStr {
		test.actual[i], _ = strconv.ParseFloat(str, 64)
	}

	// Read null values
	scanner.Scan()
	nullStr := strings.Split(scanner.Text(), " ")
	test.null = make([]float64, len(nullStr))
	for i, str := range nullStr {
		test.null[i], _ = strconv.ParseFloat(str, 64)
	}

	return test
}

// ReadModuleCorrelationTest reads a single test case from a file
func ReadModuleCorrelationTest(file string) ModuleCorrelationTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var test ModuleCorrelationTest

	// Read genes
	scanner.Scan()
	test.genes = strings.Split(scanner.Text(), " ")

	// Read expression data
	test.expressionData = make(map[string][]float64)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		parts := strings.Split(line, ":")
		if len(parts) != 2 {
			continue
		}

		gene := strings.TrimSpace(parts[0])
		valuesStr := strings.Fields(parts[1])
		values := make([]float64, len(valuesStr))

		for i, v := range valuesStr {
			values[i], _ = strconv.ParseFloat(v, 64)
		}

		test.expressionData[gene] = values
	}

	return test
}

// ReadGeneIndicesTest reads a single test case from a file
func ReadGeneIndicesTest(file string) GeneIndicesTest {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var test GeneIndicesTest

	scanner.Scan()
	test.idx, _ = strconv.Atoi(scanner.Text())

	return test
}

// ReadDirectory reads in a directory and returns a slice of fs.DirEntry objects containing file info for the directory
func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// ReadFloatFromFile reads in a single float from a file
func ReadFloatFromFile(file string) float64 {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	//read in the line
	scanner.Scan()
	line := scanner.Text()

	//convert the line to a float using strconv
	value, err := strconv.ParseFloat(line, 64)
	if err != nil {
		panic(err)
	}

	return value
}

// ReadTwoFloatsFromFile reads two floats from a file
func ReadTwoFloatsFromFile(file string) (float64, float64) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Scan()
	values := strings.Split(scanner.Text(), " ")

	val1, _ := strconv.ParseFloat(values[0], 64)
	val2, _ := strconv.ParseFloat(values[1], 64)

	return val1, val2
}

// ReadCorrelationsFromFile reads expected correlations from a file
func ReadCorrelationsFromFile(file string) []float64 {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Scan()

	valuesStr := strings.Split(scanner.Text(), " ")
	values := make([]float64, len(valuesStr))

	for i, v := range valuesStr {
		values[i], _ = strconv.ParseFloat(v, 64)
	}

	return values
}

// ReadTwoIntsFromFile reads two integers from a file
func ReadTwoIntsFromFile(file string) (int, int) {
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Scan()
	values := strings.Split(scanner.Text(), " ")

	val1, _ := strconv.Atoi(values[0])
	val2, _ := strconv.Atoi(values[1])

	return val1, val2
}

//-----------------------------------------------------------------------------------
