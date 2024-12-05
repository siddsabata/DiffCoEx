// Jason Hyun (jasonhyu)
// Siddharth Sabata (ssabata)
// Darrick Lo (ddlo)
// Katie Wang (kcw2)

// Dec 1, 2024

// NOTE: Generative AI used to produce following code:

package main

import (
	"bufio"
	"math"
	"os"
	"strconv"
	"strings"
	"testing"
)

func readInputFile(filename string) ([]float64, []float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()

	var actualCorrs, nullCorrs []float64
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Split(line, ": ")
		if len(parts) != 2 {
			continue
		}
		values := strings.Split(parts[1], ", ")
		for _, v := range values {
			f, err := strconv.ParseFloat(v, 64)
			if err != nil {
				return nil, nil, err
			}
			if strings.HasPrefix(parts[0], "actualCorrs") {
				actualCorrs = append(actualCorrs, f)
			} else if strings.HasPrefix(parts[0], "nullCorrs") {
				nullCorrs = append(nullCorrs, f)
			}
		}
	}
	return actualCorrs, nullCorrs, scanner.Err()
}

func readOutputFile(filename string) (float64, float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, 0, err
	}
	defer file.Close()

	var tstat, pval float64
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Split(line, ": ")
		if len(parts) != 2 {
			continue
		}
		f, err := strconv.ParseFloat(parts[1], 64)
		if err != nil {
			return 0, 0, err
		}
		if strings.HasPrefix(parts[0], "t-statistic") {
			tstat = f
		} else if strings.HasPrefix(parts[0], "p-value") {
			pval = f
		}
	}
	return tstat, pval, scanner.Err()
}

func roundToFourDecimalPlaces(value float64) float64 {
	return math.Round(value*10000) / 10000
}

func TestCompareWithNullFromFile(t *testing.T) {
	tests := []struct {
		inputFile  string
		outputFile string
	}{
		{"testing/CompareWithNull/Input/input1.txt", "testing/CompareWithNull/Output/output1.txt"},
		{"testing/CompareWithNull/Input/input2.txt", "testing/CompareWithNull/Output/output2.txt"},
		{"testing/CompareWithNull/Input/input3.txt", "testing/CompareWithNull/Output/output3.txt"},
		{"testing/CompareWithNull/Input/input4.txt", "testing/CompareWithNull/Output/output4.txt"},
	}

	for _, tt := range tests {
		t.Run(tt.inputFile, func(t *testing.T) {
			actualCorrs, nullCorrs, err := readInputFile(tt.inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file: %v", err)
			}

			expectedT, expectedP, err := readOutputFile(tt.outputFile)
			if err != nil {
				t.Fatalf("Failed to read output file: %v", err)
			}

			tstat, pval := compareWithNull(actualCorrs, nullCorrs)
			tstat = roundToFourDecimalPlaces(tstat)
			pval = roundToFourDecimalPlaces(pval)

			if tstat != expectedT || pval != expectedP {
				t.Errorf("compareWithNull() = (%v, %v), want (%v, %v)", tstat, pval, expectedT, expectedP)
			}
		})
	}
}
