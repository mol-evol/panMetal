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
- Command-line interface with FASTA alignment input
- Pure Python implementation with minimal dependencies
- Example alignments and trees included

## Installation

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/mol-evol/panMetal.git
cd panMetal

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
======================================================================
panMetal v1.0 - MSA Distance Calculator
https://github.com/mol-evol/panMetal/
----------------------------------------------------------------------
Implementation of algorithms from:
Blackburne & Whelan (2012) Bioinformatics 28(4):495-502
======================================================================

Reading alignment 1: examples/alignment1.fasta
Reading alignment 2: examples/alignment2.fasta

Calculating distances...
==================================================
dSSP: 0.083333
dseq: 0.333333
dpos: 0.333333
==================================================

Please cite:
- Blackburne & Whelan (2012) Bioinformatics 28(4):495-502
- McInerney (2025) panMetal: https://github.com/mol-evol/panMetal/
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
                        Alignment format (currently only 'fasta' supported)
  -m METRICS, --metrics METRICS
                        Metrics to calculate (default: all)
  --cite                Show citation information
```

## Examples

The `examples/` directory contains test alignments and phylogenetic trees:
- `alignment1.fasta`, `alignment2.fasta` - Small test alignments (3 sequences)
- `example_tree.nwk` - Phylogenetic tree for the small alignments
- `primate_align1.fasta`, `primate_align2.fasta` - Larger test alignments (12 primate sequences)
- `primate_tree.nwk` - Phylogenetic tree for the primate alignments

### Example 1: Basic Usage

```bash
# Compare two small alignments
python metal.py -a1 examples/alignment1.fasta -a2 examples/alignment2.fasta

# Output:
# dSSP: 0.083333
# dseq: 0.333333
# dpos: 0.333333
```

### Example 2: Including Phylogenetic Tree

```bash
# Include tree for evolutionary distance (devol)
python metal.py -a1 examples/alignment1.fasta -a2 examples/alignment2.fasta -t examples/example_tree.nwk

# Output includes all four metrics:
# dSSP: 0.083333
# dseq: 0.333333
# dpos: 0.333333
# devol: 0.333333
```

### Example 3: Larger Alignments

```bash
# Compare primate alignments
python metal.py -a1 examples/primate_align1.fasta -a2 examples/primate_align2.fasta -t examples/primate_tree.nwk

# Calculate specific metrics only
python metal.py -a1 examples/primate_align1.fasta -a2 examples/primate_align2.fasta -m dSSP dseq
```

## File Formats

### Input Alignments
- FASTA format only
- Gaps represented as '-'
- All sequences must have equal length
- Both alignments must have the same sequences

### Tree Format
- Newick format (.nwk files)
- Taxa names in tree must match sequence names in alignments

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

Basic tests are included:

```bash
# Run tests
python test.py
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