# Final-Project
Main.py is the driver file
- It includes classes and functions to handle DNA sequence processing, mRNA transcription, and protein translation
- The user can input wild-type and mutant DNA sequences to compare and the program outputs the differences and protein translation result

class DNA
- creates DNA objects where DNA strands undergo transcription to make a complementary version of DNA which is an mRNA
- contains dictionary of key-value pairs of complementary bases

class mRNA
- mRNA created from the DNA clas sis then used to undergo translation
- characters of the mRNA is read in thee nucleotides (characters) which are matched using the dictionary codon table
- mRNA strand is then converted to sequence of Amino Acids
- dictionary codon_table{...}

class Protein
- Initializes the protein object with the amino acid sequence
- The protein sequence is of type str and in the format "aminoacid1-aminoacid2"
- protein = [] which stores amino acids as the protein sequence is built

class Person
- Creates the Person object with wild-type and mutant DNA sequences
- The wild type DNA sequence is compared to the mutant DNA sequence of a person
- returns the index of mutation and what the mutation is

Text Files
- test_sequences.txt: This file contains sample DNA sequences for testing the program, including both wild-type and mutant sequences. 
- codon_table.txt: Contains the codon-to-amino acid mapping used for protein translation.

Functions:
- def transcription(self) -> str:
  - purpose: simulate DNA transcription inside a biological cell
  - input: string of nucleotide bases
  - output: string of nucleotide bases complementary to input sequence
- def translation(self) -> str:
  - simulate DNA translation inside a biological cell
  - input: complementary string of nucleotide bases
  - output: corresponding amino acids the sequence codes for (string) -> str:
- compare_dna_sequence_deletion(self, wild_type, mutant) -> str:
  - compare dna sequences of a normal individual with the dna sequence of a diseased/affected individual
  - input: wild-type DNA stand (string) that acts as a negative control and a mutated DNA strand
  - output: index where mutations are located and the deleted/inserted DNA strand and what the protein the strand was supposed to code for










Who did what:
Noa- 
File Handling
Class MRNA
Class Person-  def check_frameshift(self),  def compare_dna_sequence_insertion(self):,  def get_all_proteins(self):

Ryan-
Testing
Class DNA
Class Person- def compare_dna_sequence_deletion(self):, def detect_mutation_type(self):

Both of us worked on the output print statements
