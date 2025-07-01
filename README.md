# panMetal - Multiple Sequence Alignment Distance Calculator

A Python implementation of the four MSA distance metrics from Blackburne & Whelan (2012).

## Citation

**This software implements algorithms from:**
```
Blackburne, B.P. and Whelan, S. (2012). Measuring the distance between 
multiple sequence alignments. Bioinformatics, 28(4), 495-502.
https://doi.org/10.1093/bioinformatics/btr701
```

**If you use panMetal in your research, please cite:**
```
McInerney, J.O. (2025). panMetal: a python program to measure the distance 
between multiple sequence alignments. https://github.com/mol-evol/panMetal/
```

## Overview

panMetal provides tools to compare Multiple Sequence Alignments (MSAs) using four different distance metrics that capture various aspects of alignment similarity:

1. **dSSP** - Symmetrized Sum-of-Pairs metric (ignores gaps)
2. **dseq** - Simple gap information metric
3. **dpos** - Positional gap information metric
4. **devol** - Evolutionary gap information metric (uses phylogenetic trees)

## Features

- Four complementary distance metrics for comprehensive alignment comparison
- Support for phylogenetic tree input with proper Dollo parsimony (via DendroPy)
- Flexible input formats (edge dictionaries or tree objects)
- Pure Python implementation with minimal dependencies
- Comprehensive test suite with 40+ tests
- Example files included (example_tree.nwk)

## Installation

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MetAL.git
cd MetAL

# Basic usage (no additional dependencies required)
python metal.py
```

### Full Installation (with tree support)

```bash
# Install DendroPy for phylogenetic tree handling
pip install dendropy

# Or install all dependencies
pip install -r requirements.txt
```

## Quick Start

### Command Line Usage

```bash
# Compare two alignments
python metal.py -a1 examples/alignment1.fasta -a2 examples/alignment2.fasta

# Include phylogenetic tree for devol metric
python metal.py -a1 examples/alignment1.fasta -a2 examples/alignment2.fasta -t examples/example_tree.nwk

# Calculate specific metrics only
python metal.py -a1 examples/alignment1.fasta -a2 examples/alignment2.fasta -m dSSP dseq

# Get help
python metal.py --help
```

### Output Format

```
Reading alignment 1: examples/alignment1.fasta
Reading alignment 2: examples/alignment2.fasta

Calculating distances...
==================================================
dSSP: 0.083333
dseq: 0.333333
dpos: 0.333333
devol: 0.333333
==================================================
```

### Python API Usage

```python
from metal import MSAAlignment, MSADistanceCalculator, read_alignment

# Read alignments from files
alignment1 = read_alignment("alignment1.fasta")
alignment2 = read_alignment("alignment2.fasta")

# Calculate distances
calculator = MSADistanceCalculator()
distances = calculator.calculate_all_distances(alignment1, alignment2)

# Access individual metrics
print(f"dSSP: {distances['dSSP']:.4f}")
print(f"dseq: {distances['dseq']:.4f}")
print(f"dpos: {distances['dpos']:.4f}")
```

## Distance Metrics Explained

### dSSP (Symmetrized Sum-of-Pairs)
- Compares homology relationships between non-gap characters
- Ignores gap positions entirely
- Useful for comparing core structural alignments

### dseq (Simple Gap Information)
- Treats gaps as distinct characters labeled by sequence
- Each gap is labeled as G₀, G₁, etc. for sequences 0, 1, etc.
- Captures basic gap distribution patterns

### dpos (Positional Gap Information)
- Labels gaps with both sequence and position information
- Gaps labeled as G₁₅ (gap in sequence 1, following position 5)
- Captures more detailed gap placement information

### devol (Evolutionary Gap Information)
- Uses phylogenetic information to infer gap evolution
- Applies Dollo parsimony (gaps gained once, lost multiple times)
- Labels gaps with evolutionary edge information
- Most biologically informative when tree is available

## Command-Line Options

```
usage: metal.py [-h] -a1 ALIGN1 -a2 ALIGN2 [-t TREE] [-f FORMAT]
                [-m {dSSP,dseq,dpos,devol,all} [{dSSP,dseq,dpos,devol,all} ...]]

panMetal - MSA Distance Calculator

optional arguments:
  -h, --help            show this help message and exit
  -a1 ALIGN1, --align1 ALIGN1
                        First alignment file
  -a2 ALIGN2, --align2 ALIGN2
                        Second alignment file
  -t TREE, --tree TREE  Phylogenetic tree file (Newick format) for devol metric
  -f FORMAT, --format FORMAT
                        Alignment format (default: fasta)
  -m METRICS, --metrics METRICS
                        Metrics to calculate (default: all)
  --cite                Show citation information
```

## Examples

### Example 1: Comparing Similar Alignments

```python
# Nearly identical alignments
alignment1 = MSAAlignment(["ACGT", "ACGT", "ACGT"])
alignment2 = MSAAlignment(["ACGT", "ACGT", "AC-T"])

distances = calculator.calculate_all_distances(alignment1, alignment2)
# All distances will be close to 0
```

### Example 2: Comparing Different Gap Patterns

```python
# Different gap placements
alignment1 = MSAAlignment(["A--CGT", "AC--GT", "ACGT--"])
alignment2 = MSAAlignment(["--ACGT", "-AC-GT", "ACG--T"])

distances = calculator.calculate_all_distances(alignment1, alignment2)
# Gap-aware metrics (dseq, dpos) will show larger distances
```

### Example 3: Using Tree Files

```python
# Load tree from file
tree = dendropy.Tree.get(path="my_phylogeny.nwk", schema="newick")

# Ensure sequence names match tree taxa
alignment1 = MSAAlignment(sequences1, tree_taxa_names)
alignment2 = MSAAlignment(sequences2, tree_taxa_names)

# Calculate with evolutionary information
distances = calculator.calculate_all_distances(alignment1, alignment2, tree)
```

## File Formats

### Input Alignments
- Simple list of strings format
- Each string represents one sequence
- Gaps represented as '-'
- All sequences must have equal length

### Tree Formats (via DendroPy)
- Newick format (.nwk, .tree)
- NEXUS format (.nex, .nexus)
- PhyloXML format (.xml)

Example Newick tree:
```
((Species1:0.1,Species2:0.1):0.05,Species3:0.15);
```

## Limitations

1. **Equal Length Requirement**: Both alignments must have the same number of sequences and alignment length
2. **Sequence Order**: Sequences must be in the same order in both alignments
3. **Gap Character**: Only '-' is recognized as a gap character
4. **Tree Taxa Matching**: For devol metric, sequence names must match tree taxa labels


## Testing

The project includes a comprehensive test suite:

```bash
# Run all tests
python -m unittest discover

# Run specific test modules
python -m unittest test_metal.py         # Core functionality tests
python -m unittest test_tree_handler.py  # Tree handling tests (requires DendroPy)
python -m unittest test_integration.py   # Integration tests with examples

# Run doctests
python -m doctest metal.py
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Troubleshooting

### ImportError: No module named 'dendropy'
- Install DendroPy: `pip install dendropy`
- Or use the tool without tree support (dSSP, dseq, dpos metrics only)

### ValueError: Sequences not in tree
- Ensure sequence names in alignment match taxa labels in tree
- Check tree file format and loading

### Different alignment lengths
- This tool requires alignments to have the same dimensions
- Consider trimming or padding alignments to match

## Contact

For questions or issues, please open an issue on GitHub or contact [your email].