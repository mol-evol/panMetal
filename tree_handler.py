"""
Phylogenetic tree handling for MSA distance calculations.
Provides tree loading, edge extraction, and Dollo parsimony implementation.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict, deque

try:
    import dendropy
    DENDROPY_AVAILABLE = True
except ImportError:
    DENDROPY_AVAILABLE = False
    print("Warning: DendroPy not installed. Tree handling functionality will be limited.")
    print("Install with: pip install dendropy")


class PhylogeneticTreeHandler:
    """Handles phylogenetic trees for MSA distance calculations"""
    
    def __init__(self):
        if not DENDROPY_AVAILABLE:
            raise ImportError("DendroPy is required for tree handling. Install with: pip install dendropy")
    
    def load_tree(self, tree_input, schema="newick"):
        """
        Load a phylogenetic tree from file or string.
        
        Args:
            tree_input: Path to tree file or tree string
            schema: Tree format ('newick', 'nexus', 'phyloxml')
            
        Returns:
            dendropy.Tree object
        """
        if isinstance(tree_input, str) and "\n" not in tree_input and "(" not in tree_input:
            # Likely a file path
            return dendropy.Tree.get(path=tree_input, schema=schema)
        else:
            # Tree string
            return dendropy.Tree.get(data=tree_input, schema=schema)
    
    def extract_edges(self, tree) -> Dict[int, str]:
        """
        Extract edges from a phylogenetic tree.
        
        Args:
            tree: dendropy.Tree object
            
        Returns:
            Dictionary mapping edge indices to edge labels
        """
        edges = {}
        edge_index = 0
        
        # Ensure all nodes have labels
        for i, node in enumerate(tree.preorder_node_iter()):
            if not node.label:
                if node.is_leaf():
                    # Try to use taxon label for leaves
                    node.label = node.taxon.label if node.taxon else f"leaf_{i}"
                else:
                    node.label = f"node_{i}"
        
        # Extract edges
        for edge in tree.preorder_edge_iter():
            if edge.tail_node is not None:  # Skip root edge
                parent_label = edge.tail_node.label
                child_label = edge.head_node.label
                edges[edge_index] = f"{parent_label}->{child_label}"
                edge_index += 1
        
        return edges
    
    def infer_dollo_parsimony(self, tree, gap_patterns: Dict[str, List[bool]]) -> str:
        """
        Infer indel placement using Dollo parsimony.
        Under Dollo parsimony, a character (gap) can be gained only once but lost multiple times.
        
        Args:
            tree: dendropy.Tree object
            gap_patterns: Dictionary mapping taxon names to gap presence (True/False)
            
        Returns:
            Edge label where the indel most likely occurred
        """
        # Map taxa to their gap states
        for leaf in tree.leaf_node_iter():
            if leaf.taxon and leaf.taxon.label in gap_patterns:
                leaf.gap_state = gap_patterns[leaf.taxon.label]
            else:
                leaf.gap_state = False
        
        # Reconstruct ancestral states using Dollo parsimony
        # Bottom-up pass: determine possible states at internal nodes
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                # Leaf states are already set
                node.can_be_absent = not node.gap_state
                node.can_be_present = node.gap_state
            else:
                # Internal node: check children
                children = node.child_nodes()
                
                # Can be absent if all children can be absent
                node.can_be_absent = all(child.can_be_absent for child in children)
                
                # Can be present if any child requires presence
                node.can_be_present = any(child.can_be_present for child in children)
        
        # Top-down pass: assign final states and find gain edge
        root = tree.seed_node
        root.final_state = root.can_be_present  # Root state
        gain_edge = None
        
        for node in tree.preorder_node_iter():
            if node.parent_node is None:
                continue
                
            parent = node.parent_node
            
            # Determine node state based on parent and constraints
            if parent.final_state:
                # Parent has gap - child inherits unless it must be absent
                node.final_state = node.can_be_present
            else:
                # Parent lacks gap - child gains only if required
                if node.can_be_present and not node.can_be_absent:
                    node.final_state = True
                    # This is the gain edge!
                    if gain_edge is None:
                        edge = node.edge
                        parent_label = edge.tail_node.label if edge.tail_node else "root"
                        child_label = edge.head_node.label
                        gain_edge = f"{parent_label}->{child_label}"
                else:
                    node.final_state = False
        
        return gain_edge or "root"
    
    def create_dollo_parsimony_calculator(self, tree, alignment_sequences: List[str], 
                                        sequence_names: List[str]) -> Dict[int, str]:
        """
        Create a function that can calculate Dollo parsimony for each column.
        
        Args:
            tree: dendropy.Tree object
            alignment_sequences: List of aligned sequences
            sequence_names: Names corresponding to tree taxa
            
        Returns:
            Function that takes a column index and returns the inferred indel edge
        """
        # Validate that sequence names match tree taxa
        tree_taxa = set(taxon.label for taxon in tree.taxon_namespace)
        seq_taxa = set(sequence_names)
        
        if not seq_taxa.issubset(tree_taxa):
            missing = seq_taxa - tree_taxa
            raise ValueError(f"Sequences not in tree: {missing}")
        
        # Create mapping
        self.tree = tree
        self.seq_name_to_index = {name: i for i, name in enumerate(sequence_names)}
        self.alignment = alignment_sequences
        
        def get_indel_edge_for_column(col_idx: int) -> str:
            """Get the inferred indel edge for a specific column"""
            gap_patterns = {}
            
            for taxon in tree.taxon_namespace:
                if taxon.label in self.seq_name_to_index:
                    seq_idx = self.seq_name_to_index[taxon.label]
                    is_gap = self.alignment[seq_idx][col_idx] == '-'
                    gap_patterns[taxon.label] = is_gap
            
            return self.infer_dollo_parsimony(tree, gap_patterns)
        
        # Create a mapping for all columns
        edge_map = {}
        for col_idx in range(len(alignment_sequences[0])):
            edge_map[col_idx] = get_indel_edge_for_column(col_idx)
        
        return edge_map


def create_example_tree():
    """Create an example tree for testing"""
    if not DENDROPY_AVAILABLE:
        return None
    
    # Example Newick tree string
    tree_string = "((Seq1:0.1,Seq2:0.1):0.05,Seq3:0.15);"
    tree = dendropy.Tree.get(data=tree_string, schema="newick")
    
    return tree


def demo_tree_handling():
    """Demonstrate tree handling capabilities"""
    if not DENDROPY_AVAILABLE:
        print("DendroPy not installed. Cannot run demo.")
        return
    
    handler = PhylogeneticTreeHandler()
    
    # Create example tree
    tree = create_example_tree()
    print("Example tree created")
    
    # Extract edges
    edges = handler.extract_edges(tree)
    print("\nExtracted edges:")
    for idx, edge in edges.items():
        print(f"  {idx}: {edge}")
    
    # Example gap patterns for Dollo parsimony
    gap_patterns = {
        "Seq1": True,
        "Seq2": True,
        "Seq3": False
    }
    
    # Infer indel edge
    indel_edge = handler.infer_dollo_parsimony(tree, gap_patterns)
    print(f"\nInferred indel edge: {indel_edge}")
    
    return edges


if __name__ == "__main__":
    demo_tree_handling()