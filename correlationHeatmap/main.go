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
	"sort"

	"log"
	"os"
	"path/filepath"
	"strconv"

	"math/rand"
	"time"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette/moreland"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

// MatrixGrid wraps a *mat.Dense and implements the plotter.GridXYZ interface
type MatrixGrid struct {
	*mat.Dense
}

func (m MatrixGrid) Dims() (c, r int) {
	return m.Dense.Dims()
}

func (m MatrixGrid) Z(c, r int) float64 {
	return m.At(r, c)
}

func (m MatrixGrid) X(c int) float64 {
	return float64(c)
}

func (m MatrixGrid) Y(r int) float64 {
	return float64(r)
}

// ReadCSV reads a CSV file, excluding the first row and first column, and returns a 2D array of float64 and the first column as an array of strings.
func ReadCSV(fileName string) ([][]float64, []string, error) {
	file, err := os.Open(fileName)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, nil, err
	}

	var data [][]float64
	var genes []string

	// Process all rows including header
	for _, row := range records {
		var rowData []float64
		for j, value := range row {
			if j == 0 {
				// Save the first column as genes (including the header)
				genes = append(genes, value)
				continue
			}

			floatValue, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, nil, err
			}
			rowData = append(rowData, floatValue)
		}
		data = append(data, rowData)
	}

	return data, genes, nil // Return all genes including the header
}

// MergeMatrices merges two 2D arrays such that the upper triangular part (including the diagonal) comes from the first matrix and the lower triangular part comes from the second matrix.
func MergeMatrices(matrix1, matrix2 [][]float64) [][]float64 {
	rows := len(matrix1)
	cols := len(matrix1[0])
	merged := make([][]float64, rows)
	for i := range merged {
		merged[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			if i <= j {
				merged[i][j] = matrix1[i][j]
			} else {
				merged[i][j] = matrix2[i][j]
			}
		}
	}

	return merged
}

func main() {
	// ./correlationHeatmap condition1Data condition2Data
	// Check if correct number of arguments are provided
	if len(os.Args) != 3 {
		log.Fatalf("Usage: %s condition1Data condition2Data\n", os.Args[0])
	}

	// Get file names from command line arguments
	condition1File := os.Args[1]
	condition2File := os.Args[2]

	// Read the CSV files
	matrix1, genes, err := ReadCSV(condition1File)
	if err != nil {
		log.Fatalf("Error reading %s: %v", condition1File, err)
	}

	matrix2, _, err := ReadCSV(condition2File)
	if err != nil {
		log.Fatalf("Error reading %s: %v", condition2File, err)
	}

	// Find minimum number of columns between the two matrices
	minCols := len(matrix1[0])
	if len(matrix2[0]) < minCols {
		minCols = len(matrix2[0])
	}

	// Create output/plotting directory if it doesn't exist
	outputDir := "output/plotting"
	if err := os.MkdirAll(outputDir, 0755); err != nil {
		log.Fatalf("Error creating output directory: %v", err)
	}

	// Create a color palette
	palette := moreland.Kindlmann().Palette(256)

	var mergedMatrix [][]float64
	n := len(matrix1)

	if n > 50 {
		// Create a list of indices and shuffle it
		rand.Seed(time.Now().UnixNano())
		indices := make([]int, n)
		for i := range indices {
			indices[i] = i
		}
		rand.Shuffle(len(indices), func(i, j int) {
			indices[i], indices[j] = indices[j], indices[i]
		})

		// Take first 50 indices and sort them
		indices = indices[:50]
		sort.Ints(indices)

		// Create matrices with selected genes
		matrix1Dense := mat.NewDense(50, minCols, nil)
		matrix2Dense := mat.NewDense(50, minCols, nil)

		// Fill the dense matrices using selected genes
		for i := 0; i < 50; i++ {
			idx := indices[i]
			for j := 0; j < minCols; j++ {
				matrix1Dense.Set(i, j, matrix1[idx][j])
				matrix2Dense.Set(i, j, matrix2[idx][j])
			}
		}

		// Update genes list to match selected indices
		newGenes := make([]string, 50)
		for i, idx := range indices {
			newGenes[i] = genes[idx]
		}
		genes = newGenes

		// Calculate correlations for the 50 selected genes
		matrix1Corr := mat.NewDense(50, 50, nil)
		matrix2Corr := mat.NewDense(50, 50, nil)

		// Calculate correlations for condition 1
		for i := 0; i < 50; i++ {
			row1i := mat.Row(nil, i, matrix1Dense)
			for j := 0; j < 50; j++ {
				row1j := mat.Row(nil, j, matrix1Dense)
				corr := stat.Correlation(row1i, row1j, nil)
				matrix1Corr.Set(i, j, corr)
			}
		}

		// Calculate correlations for condition 2
		for i := 0; i < 50; i++ {
			row2i := mat.Row(nil, i, matrix2Dense)
			for j := 0; j < 50; j++ {
				row2j := mat.Row(nil, j, matrix2Dense)
				corr := stat.Correlation(row2i, row2j, nil)
				matrix2Corr.Set(i, j, corr)
			}
		}

		// Convert mat.Dense to [][]float64 for merging
		matrix1CorrSlice := make([][]float64, 50)
		matrix2CorrSlice := make([][]float64, 50)
		for i := 0; i < 50; i++ {
			matrix1CorrSlice[i] = make([]float64, 50)
			matrix2CorrSlice[i] = make([]float64, 50)
			for j := 0; j < 50; j++ {
				matrix1CorrSlice[i][j] = matrix1Corr.At(i, j)
				matrix2CorrSlice[i][j] = matrix2Corr.At(i, j)
			}
		}

		mergedMatrix = MergeMatrices(matrix1CorrSlice, matrix2CorrSlice)
		n = 50 // Update n for plotting
	} else {
		// For small matrices, trim to minimum columns and merge directly
		matrix1Trimmed := make([][]float64, n)
		matrix2Trimmed := make([][]float64, n)
		for i := 0; i < n; i++ {
			matrix1Trimmed[i] = matrix1[i][:minCols]
			matrix2Trimmed[i] = matrix2[i][:minCols]
		}
		mergedMatrix = MergeMatrices(matrix1Trimmed, matrix2Trimmed)
	}

	// Flatten the merged matrix for mat.Dense
	flattenedMatrix := make([]float64, n*n)
	for i := range mergedMatrix {
		for j := range mergedMatrix[i] {
			flattenedMatrix[i*n+j] = mergedMatrix[i][j]
		}
	}

	// Create merged heatmap
	pMerged := plot.New()
	pMerged.Title.Text = "Merged Heat Map"
	pMerged.X.Label.Text = "Genes"
	pMerged.Y.Label.Text = "Genes"

	gridMerged := MatrixGrid{mat.NewDense(n, n, flattenedMatrix)}
	hMerged := plotter.NewHeatMap(gridMerged, palette)
	pMerged.Add(hMerged)

	// Add gene labels
	pMerged.NominalX(genes...)
	pMerged.NominalY(genes...)
	pMerged.X.Tick.Label.Rotation = math.Pi * -0.5
	pMerged.X.Tick.Label.YAlign = draw.YCenter
	pMerged.X.Tick.Label.XAlign = draw.XRight

	// Save the merged heatmap to output/plotting directory
	if err := pMerged.Save(10*vg.Inch, 10*vg.Inch, filepath.Join(outputDir, "heatmap_merged.png")); err != nil {
		panic(err)
	}
}
