import unittest
from Main import DNA, Person  # Import the main module


class TestMutationAnalysis(unittest.TestCase):

    def test_frameshift_detection_1(self):
        wild_type_1 = "ATGCGTA"
        mutant_1 = "ATGC"
        person_1 = Person(wild_type_1, mutant_1)
        self.assertEqual(person_1.check_frameshift(), "Mutation detected but is in-frame")

    def test_frameshift_detection_2(self):
        wild_type_2 = "ATGCGTA"
        mutant_2 = "ATGCGTAA"
        person_2 = Person(wild_type_2, mutant_2)
        self.assertEqual(person_2.check_frameshift(), "Frameshift mutation detected!")

    def test_insertion_detection_1(self):
        wild_type_3 = "ATGCGTA"
        mutant_3 = "ATGCGTA"
        person_3 = Person(wild_type_3, mutant_3)
        self.assertEqual(person_3.compare_dna_sequence_insertion(), "No insertion detected")

    def test_insertion_detection_2(self):
        wild_type_4 = "ATGCGTA"
        mutant_4 = "ATGCGTAA"
        person_4 = Person(wild_type_4, mutant_4)
        self.assertEqual(person_4.compare_dna_sequence_insertion(), "Insertion detected at the end of the sequence at nt 8")

    def test_deletion_detection(self):
        wild_type_5 = "ATGCGTA"
        mutant_5 = "ATCGTA"
        person_5 = Person(wild_type_5, mutant_5)
        self.assertEqual(person_5.compare_dna_sequence_deletion(), "Deleted nucleotides: G at nt 3")

    def test_deletion_detection_2(self):
        wild_type_6 = "ATGCGTA"
        mutant_6 = "ATGCGTA"
        person_6 = Person(wild_type_6, mutant_6)
        self.assertEqual(person_6.compare_dna_sequence_deletion(), "No deletion detected")

    def test_mutation_type_1(self):
        wild_type_7 = "TACCTA"
        mutant_7 = "TACCAA"
        person_7 = Person(wild_type_7, mutant_7)
        self.assertEqual(person_7.detect_mutation_type(), "Missense mutation (amino acid change detected).")

    def test_mutation_type_2(self):
        wild_type_8 = "TACCTG"
        mutant_8 = "TACCTA"
        person_8 = Person(wild_type_8, mutant_8)
        self.assertEqual(person_8.detect_mutation_type(), "No change in protein sequence")

    def test_no_mutation_detected(self):
        """Test case with no mutation"""
        wild_type_4 = "ATGCGTAC"
        mutant_4 = "ATGCGTAC"
        person_4 = Person(wild_type_4, mutant_4)
        self.assertEqual(person_4.check_frameshift(), "No mutation detected.")
        self.assertEqual(person_4.detect_mutation_type(), "No change in protein sequence")

    def test_no_mutation_detected_2(self):
        """Test case with no mutation"""
        wild_type_4 = "TAGTGAAAAAAA"
        mutant_4 = "TAGTGAAAAAAA"
        person_4 = Person(wild_type_4, mutant_4)
        self.assertEqual(person_4.check_frameshift(), "No mutation detected.")
        self.assertEqual(person_4.detect_mutation_type(), "No change in protein sequence")

    def test_transcription_and_translation(self):
        wild_type_dna = DNA("TACTGCCCCCCC")
        wild_type_mrna = wild_type_dna.transcription()

        if wild_type_mrna:
            mutant_protein = wild_type_mrna.translation()
            self.assertEqual(wild_type_mrna.seq, "AUGACGGGGGGG")
            self.assertEqual(mutant_protein, ['Met-Thr-Gly-Gly'])  # Check this protein sequence carefully

    def test_transcription_and_translation_2(self):
        mutant_dna = DNA("TACTGAAAAAAA")
        mutant_mrna = mutant_dna.transcription()

        if mutant_mrna:
            mutant_protein = mutant_mrna.translation()
            self.assertEqual(mutant_mrna.seq, "AUGACUUUUUUU")
            self.assertEqual(mutant_protein, ['Met-Thr-Phe-Phe'])  # Check this protein sequence carefully

    def test_protein_comparison(self):
        """Test protein comparison between wild-type and mutant"""
        wild_type_7 = "TACCTA"
        mutant_7 = "TACCAA"
        person_7 = Person(wild_type_7, mutant_7)

        # Assuming that Person class has a method to get proteins directly
        wt_proteins = person_7.get_all_proteins()[0]  # Getting wild-type protein
        mt_proteins = person_7.get_all_proteins()[1]  # Getting mutant protein

        self.assertEqual(wt_proteins, ['Met-Asp'])  # Expected wild-type protein
        self.assertEqual(mt_proteins, ['Met-Val'])

    def test_protein_comparison_2(self):
        wild_type_7 = "TACACG"
        mutant_7 = "TACAGT"
        person_7 = Person(wild_type_7, mutant_7)

        # Assuming that Person class has a method to get proteins directly
        wt_proteins = person_7.get_all_proteins()[0]  # Getting wild-type protein
        mt_proteins = person_7.get_all_proteins()[1]  # Getting mutant protein

        self.assertEqual(wt_proteins, ['Met-Cys'])  # Expected wild-type protein
        self.assertEqual(mt_proteins, ['Met-Ser'])

if __name__ == "__main__":
    unittest.main()
