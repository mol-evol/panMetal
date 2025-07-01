"""
panMetal: Python implementation of MSA distance metrics

This software implements the four MSA distance metrics described in:
Blackburne, B.P. and Whelan, S. (2012). Measuring the distance between 
multiple sequence alignments. Bioinformatics, 28(4), 495-502.
https://doi.org/10.1093/bioinformatics/btr701

Citation for this software:
McInerney, J.O. (2025). panMetal: a python program to measure the distance 
between multiple sequence alignments. https://github.com/mol-evol/panMetal/

This module implements four metrics to compare Multiple Sequence Alignments (MSAs):
1. dSSP: Symmetrized Sum-of-Pairs metric
2. dseq: Simple gap information metric
3. dpos: Positional gap information metric
4. devol: Evolutionary gap information metric (requires phylogenetic tree)
"""

import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Set, Optional, Union
import itertools
import argparse
import sys
import os

# Try to import tree handler
try:
    from tree_handler import PhylogeneticTreeHandler, DENDROPY_AVAILABLE
except ImportError:
    DENDROPY_AVAILABLE = False
    PhylogeneticTreeHandler = None


class MSAAlignment:
    """Represents a Multiple Sequence Alignment
    
    Examples:
        >>> sequences = ["ACGT", "AC-T", "A-GT"]
        >>> alignment = MSAAlignment(sequences, ["Human", "Mouse", "Rat"])
        >>> alignment.num_sequences
        3
        >>> alignment.get_column(2)
        ['G', '-', 'G']
    """

    def __init__(self, sequences: List[str], sequence_names: Optional[List[str]] = None):
        """
        Initialize MSA alignment

        Args:
            sequences: List of aligned sequences (with gaps as '-')
            sequence_names: Optional list of sequence names
            
        Raises:
            ValueError: If sequences have different lengths
            
        Examples:
            >>> align = MSAAlignment(["ACGT", "AC-T"])
            >>> align.sequences
            ['ACGT', 'AC-T']
        """
        self.sequences = [seq.upper() for seq in sequences]
        self.num_sequences = len(sequences)
        self.alignment_length = len(sequences[0]) if sequences else 0
        self.sequence_names = sequence_names or [f"Seq{i+1}" for i in range(self.num_sequences)]

        # Validate alignment
        if not all(len(seq) == self.alignment_length for seq in self.sequences):
            raise ValueError("All sequences must have the same length")

    def get_column(self, col_idx: int) -> List[str]:
        """Get a specific column from the alignment"""
        return [seq[col_idx] for seq in self.sequences]

    def get_character(self, seq_idx: int, pos_idx: int) -> str:
        """Get character at specific sequence and position"""
        return self.sequences[seq_idx][pos_idx]


class MSADistanceCalculator:
    """Calculator for MSA distance metrics"""

    def __init__(self):
        self.gap_char = '-'

    def _create_homology_sets_ssp(self, alignment: MSAAlignment) -> Dict[Tuple[int, int], Set[Tuple[int, int]]]:
        """
        Create homology sets for SSP metric (ignores gaps)

        Returns:
            Dictionary mapping (seq_idx, pos_idx) to set of homologous positions
        """
        homology_sets = {}

        for col_idx in range(alignment.alignment_length):
            column = alignment.get_column(col_idx)
            # Get positions of non-gap characters in this column
            non_gap_positions = [(seq_idx, col_idx) for seq_idx, char in enumerate(column)
                               if char != self.gap_char]

            # Each non-gap character is homologous to all others in the same column
            for seq_idx, pos_idx in non_gap_positions:
                homology_sets[(seq_idx, pos_idx)] = set(non_gap_positions)

        return homology_sets

    def _create_homology_sets_seq(self, alignment: MSAAlignment) -> Dict[Tuple[int, int], Set[str]]:
        """
        Create homology sets for seq metric (simple gap labeling)

        Returns:
            Dictionary mapping (seq_idx, pos_idx) to set of homology labels
        """
        homology_sets = {}

        for col_idx in range(alignment.alignment_length):
            column = alignment.get_column(col_idx)
            homology_set = set()

            for seq_idx, char in enumerate(column):
                if char == self.gap_char:
                    homology_set.add(f"G{seq_idx}")
                else:
                    homology_set.add(char)

            # Assign the same homology set to all positions in this column
            for seq_idx in range(alignment.num_sequences):
                homology_sets[(seq_idx, col_idx)] = homology_set.copy()

        return homology_sets

    def _create_homology_sets_pos(self, alignment: MSAAlignment) -> Dict[Tuple[int, int], Set[str]]:
        """
        Create homology sets for pos metric (positional gap information)

        Returns:
            Dictionary mapping (seq_idx, pos_idx) to set of homology labels
        """
        homology_sets = {}

        for col_idx in range(alignment.alignment_length):
            column = alignment.get_column(col_idx)
            homology_set = set()

            for seq_idx, char in enumerate(column):
                if char == self.gap_char:
                    # Find position of real character to the left
                    left_pos = self._find_left_character_position(alignment, seq_idx, col_idx)
                    homology_set.add(f"G{seq_idx}{left_pos}")
                else:
                    homology_set.add(char)

            # Assign the same homology set to all positions in this column
            for seq_idx in range(alignment.num_sequences):
                homology_sets[(seq_idx, col_idx)] = homology_set.copy()

        return homology_sets

    def _find_left_character_position(self, alignment: MSAAlignment, seq_idx: int, col_idx: int) -> int:
        """Find the position of the nearest real character to the left"""
        for pos in range(col_idx - 1, -1, -1):
            if alignment.get_character(seq_idx, pos) != self.gap_char:
                return pos
        return 0  # If no character found to the left, use position 0

    def _create_homology_sets_evol(self, alignment: MSAAlignment, tree_edges: Union[Dict[int, str], 'dendropy.Tree']) -> Dict[Tuple[int, int], Set[str]]:
        """
        Create homology sets for evol metric (evolutionary gap information)

        Args:
            alignment: MSA alignment
            tree_edges: Either a dictionary mapping edge indices to edge labels,
                       or a dendropy.Tree object for proper Dollo parsimony

        Returns:
            Dictionary mapping (seq_idx, pos_idx) to set of homology labels
        """
        homology_sets = {}
        
        # Check if we have a tree object for proper Dollo parsimony
        if DENDROPY_AVAILABLE and hasattr(tree_edges, 'taxon_namespace'):
            # Use proper Dollo parsimony with the tree
            handler = PhylogeneticTreeHandler()
            
            # Create column-specific indel edge mapping
            indel_edges = {}
            for col_idx in range(alignment.alignment_length):
                column = alignment.get_column(col_idx)
                gap_patterns = {}
                
                # Create gap pattern for this column
                for seq_idx, seq_name in enumerate(alignment.sequence_names):
                    gap_patterns[seq_name] = (column[seq_idx] == self.gap_char)
                
                # Infer indel edge using Dollo parsimony
                try:
                    indel_edges[col_idx] = handler.infer_dollo_parsimony(tree_edges, gap_patterns)
                except (ValueError, AttributeError):
                    # Fallback if sequence names don't match tree taxa
                    indel_edges[col_idx] = "unknown"
        else:
            # Use simplified implementation with edge dictionary
            indel_edges = None

        for col_idx in range(alignment.alignment_length):
            column = alignment.get_column(col_idx)

            # Get indel edge for this column
            if indel_edges is not None:
                indel_edge = indel_edges[col_idx]
            else:
                # Use simplified implementation
                gap_pattern = [char == self.gap_char for char in column]
                indel_edge = self._infer_indel_edge(gap_pattern, tree_edges)

            homology_set = set()
            for seq_idx, char in enumerate(column):
                if char == self.gap_char:
                    # Find position of real character to the left
                    left_pos = self._find_left_character_position(alignment, seq_idx, col_idx)
                    homology_set.add(f"G{seq_idx}{left_pos}({indel_edge})")
                else:
                    homology_set.add(char)

            # Assign the same homology set to all positions in this column
            for seq_idx in range(alignment.num_sequences):
                homology_sets[(seq_idx, col_idx)] = homology_set.copy()

        return homology_sets

    def _infer_indel_edge(self, gap_pattern: List[bool], tree_edges: Dict[int, str]) -> str:
        """
        Simplified indel edge inference using Dollo parsimony

        Args:
            gap_pattern: Boolean list indicating gap positions
            tree_edges: Dictionary of tree edges

        Returns:
            Edge label where indel likely occurred
        """
        # Simplified implementation - in practice, this would require
        # full phylogenetic analysis with Dollo parsimony
        # For now, return a default edge based on gap pattern
        num_gaps = sum(gap_pattern)
        edge_idx = num_gaps % len(tree_edges) if tree_edges else 0
        return tree_edges.get(edge_idx, "0")

    def _jaccard_distance(self, set1: Set, set2: Set) -> float:
        """Calculate Jaccard distance between two sets"""
        if not set1 and not set2:
            return 0.0

        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))

        if union == 0:
            return 0.0

        return 1.0 - (intersection / union)

    def _hamming_distance(self, set1: Set, set2: Set) -> float:
        """Calculate normalized Hamming distance between two sets"""
        if not set1 and not set2:
            return 0.0
        
        # Use the maximum set size for normalization to handle unequal sizes
        max_size = max(len(set1), len(set2))
        if max_size == 0:
            return 0.0

        intersection = len(set1.intersection(set2))
        return 1.0 - (intersection / max_size)

    def calculate_dSSP(self, alignment1: MSAAlignment, alignment2: MSAAlignment) -> float:
        """
        Calculate dSSP (Symmetrized Sum-of-Pairs) distance between two alignments
        
        The dSSP metric compares homology relationships between non-gap characters,
        ignoring gap positions entirely. It uses Jaccard distance on homology sets.

        Args:
            alignment1: First MSA alignment
            alignment2: Second MSA alignment

        Returns:
            dSSP distance (0.0 to 1.0)
            
        Examples:
            >>> calc = MSADistanceCalculator()
            >>> align1 = MSAAlignment(["AC-GT", "A--GT"])
            >>> align2 = MSAAlignment(["ACG-T", "A-G-T"])
            >>> dist = calc.calculate_dSSP(align1, align2)
            >>> 0.0 <= dist <= 1.0
            True
        """
        homology_sets1 = self._create_homology_sets_ssp(alignment1)
        homology_sets2 = self._create_homology_sets_ssp(alignment2)

        total_distance = 0.0
        num_characters = 0

        # Calculate distance for each character position
        for seq_idx in range(alignment1.num_sequences):
            for pos_idx in range(alignment1.alignment_length):
                if alignment1.get_character(seq_idx, pos_idx) != self.gap_char:
                    key = (seq_idx, pos_idx)
                    if key in homology_sets1 and key in homology_sets2:
                        distance = self._jaccard_distance(homology_sets1[key], homology_sets2[key])
                        total_distance += distance
                        num_characters += 1

        return total_distance / num_characters if num_characters > 0 else 0.0

    def calculate_dseq(self, alignment1: MSAAlignment, alignment2: MSAAlignment) -> float:
        """
        Calculate dseq distance between two alignments

        Args:
            alignment1: First MSA alignment
            alignment2: Second MSA alignment

        Returns:
            dseq distance (0.0 to 1.0)
        """
        homology_sets1 = self._create_homology_sets_seq(alignment1)
        homology_sets2 = self._create_homology_sets_seq(alignment2)

        total_distance = 0.0
        num_characters = 0

        # Calculate distance for each character position
        for seq_idx in range(alignment1.num_sequences):
            for pos_idx in range(alignment1.alignment_length):
                key = (seq_idx, pos_idx)
                if key in homology_sets1 and key in homology_sets2:
                    distance = self._hamming_distance(homology_sets1[key], homology_sets2[key])
                    total_distance += distance
                    num_characters += 1

        return total_distance / num_characters if num_characters > 0 else 0.0

    def calculate_dpos(self, alignment1: MSAAlignment, alignment2: MSAAlignment) -> float:
        """
        Calculate dpos distance between two alignments

        Args:
            alignment1: First MSA alignment
            alignment2: Second MSA alignment

        Returns:
            dpos distance (0.0 to 1.0)
        """
        homology_sets1 = self._create_homology_sets_pos(alignment1)
        homology_sets2 = self._create_homology_sets_pos(alignment2)

        total_distance = 0.0
        num_characters = 0

        # Calculate distance for each character position
        for seq_idx in range(alignment1.num_sequences):
            for pos_idx in range(alignment1.alignment_length):
                key = (seq_idx, pos_idx)
                if key in homology_sets1 and key in homology_sets2:
                    distance = self._hamming_distance(homology_sets1[key], homology_sets2[key])
                    total_distance += distance
                    num_characters += 1

        return total_distance / num_characters if num_characters > 0 else 0.0

    def calculate_devol(self, alignment1: MSAAlignment, alignment2: MSAAlignment,
                       tree_edges: Union[Dict[int, str], 'dendropy.Tree']) -> float:
        """
        Calculate devol distance between two alignments

        Args:
            alignment1: First MSA alignment
            alignment2: Second MSA alignment
            tree_edges: Either dictionary mapping edge indices to edge labels,
                       or a dendropy.Tree object for proper Dollo parsimony

        Returns:
            devol distance (0.0 to 1.0)
        """
        homology_sets1 = self._create_homology_sets_evol(alignment1, tree_edges)
        homology_sets2 = self._create_homology_sets_evol(alignment2, tree_edges)

        total_distance = 0.0
        num_characters = 0

        # Calculate distance for each character position
        for seq_idx in range(alignment1.num_sequences):
            for pos_idx in range(alignment1.alignment_length):
                key = (seq_idx, pos_idx)
                if key in homology_sets1 and key in homology_sets2:
                    distance = self._hamming_distance(homology_sets1[key], homology_sets2[key])
                    total_distance += distance
                    num_characters += 1

        return total_distance / num_characters if num_characters > 0 else 0.0

    def calculate_all_distances(self, alignment1: MSAAlignment, alignment2: MSAAlignment,
                              tree_edges: Optional[Union[Dict[int, str], 'dendropy.Tree']] = None) -> Dict[str, float]:
        """
        Calculate all four distance metrics between two alignments

        Args:
            alignment1: First MSA alignment
            alignment2: Second MSA alignment
            tree_edges: Optional tree edges for devol calculation

        Returns:
            Dictionary with all four distance metrics
        """
        distances = {
            'dSSP': self.calculate_dSSP(alignment1, alignment2),
            'dseq': self.calculate_dseq(alignment1, alignment2),
            'dpos': self.calculate_dpos(alignment1, alignment2)
        }

        if tree_edges is not None:
            distances['devol'] = self.calculate_devol(alignment1, alignment2, tree_edges)

        return distances


def read_fasta(filename: str) -> Tuple[List[str], List[str]]:
    """
    Read sequences from a FASTA file.
    
    Args:
        filename: Path to FASTA file
        
    Returns:
        Tuple of (sequences, names)
    """
    sequences = []
    names = []
    current_seq = []
    current_name = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    sequences.append(''.join(current_seq))
                    names.append(current_name)
                current_name = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_name is not None:
            sequences.append(''.join(current_seq))
            names.append(current_name)
    
    return sequences, names


def read_alignment(filename: str, format: str = 'fasta') -> MSAAlignment:
    """
    Read alignment from file.
    
    Args:
        filename: Path to alignment file
        format: File format (currently only 'fasta' supported)
        
    Returns:
        MSAAlignment object
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Alignment file not found: {filename}")
    
    if format.lower() == 'fasta':
        sequences, names = read_fasta(filename)
    else:
        raise ValueError(f"Unsupported format: {format}")
    
    return MSAAlignment(sequences, names)


def print_header():
    """Print program header with citations"""
    print("\n" + "="*70)
    print("panMetal v1.0 - MSA Distance Calculator")
    print("https://github.com/mol-evol/panMetal/")
    print("-"*70)
    print("Implementation of algorithms from:")
    print("Blackburne & Whelan (2012) Bioinformatics 28(4):495-502")
    print("="*70 + "\n")


def print_citation_info():
    """Print full citation information"""
    print("\nCitation Information:")
    print("="*70)
    print("\nOriginal Algorithm:")
    print("Blackburne, B.P. and Whelan, S. (2012). Measuring the distance between")
    print("multiple sequence alignments. Bioinformatics, 28(4), 495-502.")
    print("https://doi.org/10.1093/bioinformatics/btr701")
    print("\nThis Software:")
    print("McInerney, J.O. (2025). panMetal: a python program to measure the distance")
    print("between multiple sequence alignments. https://github.com/mol-evol/panMetal/")
    print("="*70 + "\n")


def main():
    """Main command-line interface"""
    parser = argparse.ArgumentParser(
        description='panMetal - MSA Distance Calculator\nImplementation of Blackburne & Whelan (2012) algorithms',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare two alignments
  python metal.py -a1 align1.fasta -a2 align2.fasta
  
  # Include phylogenetic tree for devol metric
  python metal.py -a1 align1.fasta -a2 align2.fasta -t tree.nwk
  
  # Show citation information
  python metal.py --cite

GitHub: https://github.com/mol-evol/panMetal/
        """
    )
    
    parser.add_argument('-a1', '--align1', required=False,
                        help='First alignment file')
    parser.add_argument('-a2', '--align2', required=False,
                        help='Second alignment file')
    parser.add_argument('-t', '--tree',
                        help='Phylogenetic tree file (Newick format) for devol metric')
    parser.add_argument('-f', '--format', default='fasta',
                        help='Alignment format (default: fasta)')
    parser.add_argument('-m', '--metrics', nargs='+',
                        choices=['dSSP', 'dseq', 'dpos', 'devol', 'all'],
                        default=['all'],
                        help='Metrics to calculate (default: all)')
    parser.add_argument('--cite', action='store_true',
                        help='Show citation information')
    
    args = parser.parse_args()
    
    # Handle --cite flag
    if args.cite:
        print_citation_info()
        sys.exit(0)
    
    # Check required arguments if not --cite
    if not args.align1 or not args.align2:
        parser.error("Arguments -a1 and -a2 are required unless using --cite")
    
    try:
        # Print header
        print_header()
        
        # Read alignments
        print(f"Reading alignment 1: {args.align1}")
        alignment1 = read_alignment(args.align1, args.format)
        
        print(f"Reading alignment 2: {args.align2}")
        alignment2 = read_alignment(args.align2, args.format)
        
        # Validate alignments have same dimensions
        if alignment1.num_sequences != alignment2.num_sequences:
            raise ValueError(f"Alignments have different number of sequences: "
                           f"{alignment1.num_sequences} vs {alignment2.num_sequences}")
        
        if alignment1.alignment_length != alignment2.alignment_length:
            raise ValueError(f"Alignments have different lengths: "
                           f"{alignment1.alignment_length} vs {alignment2.alignment_length}")
        
        # Check sequence names match
        if set(alignment1.sequence_names) != set(alignment2.sequence_names):
            print("Warning: Sequence names differ between alignments")
        
        # Load tree if provided
        tree = None
        if args.tree:
            if not DENDROPY_AVAILABLE:
                print("Warning: DendroPy not installed. Cannot use tree file.")
                print("Install with: pip install dendropy")
            else:
                print(f"Reading tree: {args.tree}")
                handler = PhylogeneticTreeHandler()
                tree = handler.load_tree(args.tree, schema="newick")
        
        # Calculate distances
        calculator = MSADistanceCalculator()
        
        print("\nCalculating distances...")
        print("=" * 50)
        
        if 'all' in args.metrics:
            distances = calculator.calculate_all_distances(alignment1, alignment2, tree)
            for metric, value in distances.items():
                print(f"{metric}: {value:.6f}")
        else:
            for metric in args.metrics:
                if metric == 'dSSP':
                    value = calculator.calculate_dSSP(alignment1, alignment2)
                elif metric == 'dseq':
                    value = calculator.calculate_dseq(alignment1, alignment2)
                elif metric == 'dpos':
                    value = calculator.calculate_dpos(alignment1, alignment2)
                elif metric == 'devol':
                    if tree is None:
                        print("devol: Requires tree file (use -t option)")
                        continue
                    value = calculator.calculate_devol(alignment1, alignment2, tree)
                print(f"{metric}: {value:.6f}")
        
        print("=" * 50)
        
        # Print citation reminder
        print("\nPlease cite:")
        print("- Blackburne & Whelan (2012) Bioinformatics 28(4):495-502")
        print("- McInerney (2025) panMetal: https://github.com/mol-evol/panMetal/")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def example_usage():
    """Example usage of the MSA distance calculator"""

    print("=== Example 1: Original alignments from paper ===")
    # Example alignments from the paper
    sequences1 = [
        "AC-GT",
        "A--GT",
        "-CTGT"
    ]

    sequences2 = [
        "ACG-T",
        "A-G-T",
        "-C-GT"
    ]

    # Create MSA alignment objects
    alignment1 = MSAAlignment(sequences1, ["Seq1", "Seq2", "Seq3"])
    alignment2 = MSAAlignment(sequences2, ["Seq1", "Seq2", "Seq3"])

    # Create distance calculator
    calculator = MSADistanceCalculator()

    # Example tree edges (simplified)
    tree_edges = {0: "edge1", 1: "edge2", 2: "edge3"}

    # Calculate all distances
    distances = calculator.calculate_all_distances(alignment1, alignment2, tree_edges)

    print("MSA Distance Metrics:")
    print("-" * 30)
    for metric, distance in distances.items():
        print(f"{metric}: {distance:.4f}")

    print("\n=== Example 2: Identical alignments (should be 0.0) ===")
    # Test with identical alignments
    identical_distances = calculator.calculate_all_distances(alignment1, alignment1, tree_edges)
    for metric, distance in identical_distances.items():
        print(f"{metric}: {distance:.4f}")

    print("\n=== Example 3: More different alignments ===")
    # Test with more divergent alignments
    sequences3 = [
        "ACGT-",
        "AC-GT", 
        "A--GT"
    ]
    
    sequences4 = [
        "--ACGT",
        "-AC-GT",
        "A---GT"
    ]
    
    alignment3 = MSAAlignment(sequences3, ["Seq1", "Seq2", "Seq3"])
    alignment4 = MSAAlignment(sequences4, ["Seq1", "Seq2", "Seq3"])
    
    divergent_distances = calculator.calculate_all_distances(alignment3, alignment4, tree_edges)
    for metric, distance in divergent_distances.items():
        print(f"{metric}: {distance:.4f}")

    return distances


def example_with_tree():
    """Example using a phylogenetic tree object"""
    
    if not DENDROPY_AVAILABLE:
        print("DendroPy not installed. Install with: pip install dendropy")
        print("Falling back to simple edge dictionary example.")
        return example_usage()
    
    print("=== Example with Phylogenetic Tree ===")
    
    # Create sequences and alignments
    sequences1 = [
        "AC-GT",
        "A--GT",
        "-CTGT"
    ]
    
    sequences2 = [
        "ACG-T",
        "A-G-T",
        "-C-GT"
    ]
    
    # Create MSA alignment objects with matching sequence names
    seq_names = ["Seq1", "Seq2", "Seq3"]
    alignment1 = MSAAlignment(sequences1, seq_names)
    alignment2 = MSAAlignment(sequences2, seq_names)
    
    # Create or load a phylogenetic tree
    try:
        import dendropy
        
        # Example 1: Create tree from Newick string
        tree_string = "((Seq1:0.1,Seq2:0.1):0.05,Seq3:0.15);"
        tree = dendropy.Tree.get(data=tree_string, schema="newick")
        
        # Example 2: Load tree from file (uncomment to use)
        # tree = dendropy.Tree.get(path="my_tree.nwk", schema="newick")
        
        # Create distance calculator
        calculator = MSADistanceCalculator()
        
        # Calculate all distances using the tree
        distances = calculator.calculate_all_distances(alignment1, alignment2, tree)
        
        print("MSA Distance Metrics with Proper Dollo Parsimony:")
        print("-" * 50)
        for metric, distance in distances.items():
            print(f"{metric}: {distance:.4f}")
        
        # Show edge extraction
        handler = PhylogeneticTreeHandler()
        edges = handler.extract_edges(tree)
        print("\nTree edges:")
        for idx, edge in edges.items():
            print(f"  {idx}: {edge}")
            
    except Exception as e:
        print(f"Error: {e}")
        print("Falling back to simple example")
        return example_usage()
    
    return distances


if __name__ == "__main__":
    # If command line arguments provided, use CLI
    if len(sys.argv) > 1:
        main()
    else:
        # Otherwise run examples
        print("No arguments provided. Running example...")
        print("For command-line usage, run: python metal.py --help\n")
        try:
            example_with_tree()
        except:
            example_usage()
