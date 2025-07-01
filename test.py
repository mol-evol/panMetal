"""
Basic tests for MetAL MSA distance calculator
Run with: python test.py
"""

import unittest
from metal import MSAAlignment, MSADistanceCalculator, read_alignment
import os


class TestMetAL(unittest.TestCase):
    """Basic tests for MetAL functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.calculator = MSADistanceCalculator()
        
        # Example alignments from paper
        self.sequences1 = ["AC-GT", "A--GT", "-CTGT"]
        self.sequences2 = ["ACG-T", "A-G-T", "-C-GT"]
        
        self.alignment1 = MSAAlignment(self.sequences1, ["Seq1", "Seq2", "Seq3"])
        self.alignment2 = MSAAlignment(self.sequences2, ["Seq1", "Seq2", "Seq3"])
    
    def test_alignment_creation(self):
        """Test creating alignment objects"""
        self.assertEqual(self.alignment1.num_sequences, 3)
        self.assertEqual(self.alignment1.alignment_length, 5)
        self.assertEqual(self.alignment1.sequences[0], "AC-GT")
    
    def test_identical_alignments(self):
        """Test that identical alignments have distance 0"""
        distances = self.calculator.calculate_all_distances(
            self.alignment1, self.alignment1
        )
        
        self.assertEqual(distances['dSSP'], 0.0)
        self.assertEqual(distances['dseq'], 0.0)
        self.assertEqual(distances['dpos'], 0.0)
    
    def test_distance_calculations(self):
        """Test distance calculations match expected values"""
        distances = self.calculator.calculate_all_distances(
            self.alignment1, self.alignment2
        )
        
        # Check all metrics are present
        self.assertIn('dSSP', distances)
        self.assertIn('dseq', distances)
        self.assertIn('dpos', distances)
        
        # Check values are in valid range
        for metric, value in distances.items():
            self.assertGreaterEqual(value, 0.0)
            self.assertLessEqual(value, 1.0)
        
        # Test known value from paper
        self.assertAlmostEqual(distances['dSSP'], 0.0833, places=4)
    
    def test_file_reading(self):
        """Test reading alignments from files"""
        # Check if example files exist
        if os.path.exists("examples/alignment1.fasta"):
            align1 = read_alignment("examples/alignment1.fasta")
            align2 = read_alignment("examples/alignment2.fasta")
            
            self.assertEqual(align1.num_sequences, 3)
            self.assertEqual(align2.num_sequences, 3)
            
            # Calculate distances
            distances = self.calculator.calculate_all_distances(align1, align2)
            self.assertIn('dSSP', distances)
    
    def test_hamming_distance_fix(self):
        """Test that Hamming distance handles unequal set sizes"""
        # This was the original bug - ensure it's fixed
        set1 = {'A', 'C', 'G'}
        set2 = {'A', 'C'}
        
        # Should not raise ValueError
        distance = self.calculator._hamming_distance(set1, set2)
        self.assertIsInstance(distance, float)


if __name__ == '__main__':
    print("Running MetAL tests...")
    unittest.main(verbosity=2)