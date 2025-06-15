# UPGMA Phylogenetic Tree Constructor

A Python implementation of the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm for constructing phylogenetic trees from either DNA sequences or distance matrices.

## Overview

This tool constructs evolutionary trees that show relationships between different species or sequences. UPGMA is a simple clustering method that assumes a constant rate of evolution (molecular clock hypothesis).

## Features

- **Two input methods**: DNA sequences (MSA) or pre-computed distance matrices
- **Automatic distance calculation**: Hamming distance for sequence data
- **Multiple output formats**: Newick format, tabular format, and visual plots
- **Interactive interface**: Easy-to-use command-line interface
- **Built-in examples**: Ready-to-run examples for learning

## Requirements

```bash
pip install numpy matplotlib
```

## Usage

### Running the Program

```bash
python main.py
```

The program will present you with three options:

```
1. Enter aligned sequences (MSA)
2. Enter distance matrix  
3. Run examples
```

### Option 1: Using DNA Sequences

When you choose option 1, you'll enter aligned DNA sequences in this format:

```
Enter sequence (or 'done'): Human ATGCGTACGTACGT
Enter sequence (or 'done'): Chimp ATGCGTACGTACGA
Enter sequence (or 'done'): Mouse ATGCGTACGTACCG
Enter sequence (or 'done'): done
```

**Important**: All sequences must be the same length (Multiple Sequence Alignment - MSA).

### Option 2: Using Distance Matrix

For option 2, you'll input a symmetric distance matrix:

```
Number of species/taxa: 3
Species 1 name: Human
Species 2 name: Chimp  
Species 3 name: Mouse

Row 1 (Human): 0.0 0.2 0.6
Row 2 (Chimp): 0.2 0.0 0.5
Row 3 (Mouse): 0.6 0.5 0.0
```

### Option 3: Examples

Choose option 3 to run built-in examples that demonstrate both input methods.

## Output Files

The program generates three types of output:

1. **`.newick`** - Standard phylogenetic tree format
2. **`.txt`** - Tabular representation with node details
3. **`.png`** - Visual tree diagram (dendrogram)

## How UPGMA Works

1. **Start**: Each sequence/species begins as its own cluster
2. **Find closest pair**: Identify the two clusters with minimum distance
3. **Merge clusters**: Combine them into a new cluster
4. **Update distances**: Calculate distances to the new cluster using weighted averages
5. **Repeat**: Continue until only one cluster remains (the root)

### Distance Calculation

- **For sequences**: Uses Hamming distance (proportion of differing positions)
- **For matrices**: Uses provided distances directly

## Example Output

### Newick Format
```
((Human:0.050000,Chimp:0.050000):0.150000,(Mouse:0.100000,Rat:0.100000):0.100000):0.000000;
```

### Visual Tree
The program creates a horizontal dendrogram showing evolutionary relationships with branch lengths proportional to genetic distances.

## File Structure

```
main.py           # Main program with UPGMA implementation
README.md         # This documentation
example1_tree.*   # Example outputs from sequence data
example2_tree.*   # Example outputs from distance matrix
```

## Educational Notes

### When to Use UPGMA
- When evolutionary rates are relatively constant
- For molecular data with clock-like behavior
- For educational purposes and simple analyses

### Limitations
- Assumes molecular clock (constant evolutionary rate)
- Not suitable for data with varying evolutionary rates
- Can produce incorrect topologies for non-clock data

### Alternative Methods
- **Neighbor-Joining**: Better for non-clock data
- **Maximum Likelihood**: More sophisticated statistical approach
- **Maximum Parsimony**: Minimizes evolutionary changes

## Code Structure

### Main Classes

- `TreeNode`: Represents nodes in the phylogenetic tree
- `UPGMATreeBuilder`: Main class implementing the UPGMA algorithm

### Key Methods

- `load_sequences_from_fasta()`: Load sequences from FASTA format
- `load_sequences_from_list()`: Load sequences from tuples
- `load_distance_matrix()`: Load pre-computed distances
- `calculate_distance_matrix()`: Compute Hamming distances
- `build_upgma_tree()`: Execute UPGMA algorithm
- `plot_tree()`: Generate visual representation

## Tips for Students

1. **Start with examples**: Run option 3 first to see how it works
2. **Use short sequences**: For manual verification, use sequences of 10-20 nucleotides
3. **Check alignment**: Ensure all sequences have the same length
4. **Interpret results**: Closer species have shorter branch lengths
5. **Compare methods**: Try the same data with different phylogenetic methods

## Common Errors

- **"All sequences must have the same length"**: Ensure proper alignment
- **"At least 2 sequences required"**: You need minimum 2 sequences/species
- **"Distance matrix must be symmetric"**: Check that matrix[i][j] = matrix[j][i]

## License

Educational use only. Created for bioinformatics coursework.

---

**Author**: Bioinformatics Tool  
**Version**: 1.0  
**Course**: Medical Bioinformatics
