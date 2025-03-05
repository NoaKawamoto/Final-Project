#DNA class
#initializes the creation of the DNA sequence object
#ensures that the DNA sequence is in uppercase

#def transcription(self):
#converts DNA sequence to complementary mRNA strand (transcription
#try and except block will be used to catch possible errors in inputting sequences that are invalid
# it will contain complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        

#class MRNA
#will contain the dictionary of key-value pairs of codon sequence and the amino acid the codon corresponds to
#codon_table = {'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu',
 #                  'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu', 'AUU': 'Ile', 'AUC': 'Ile',
   #                'AUA': 'Ile', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
  #                 'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AUG': 'Met',
   #                'AGU': 'Ser', 'AGC': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro',
   #                'CCG': 'Pro', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
 #                  'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'UAU': 'Tyr',
#                   'UAC': 'Tyr', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
#                   'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp',
#                   'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'UGU': 'Cys', 'UGC': 'Cys',
#                   'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
 #                  'AGA': 'Arg', 'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly',
 #                  'GGG': 'Gly', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'}

#initializes the creation of MRNA object
#Protein variable is declared outside the for loop which will be appended with Amino Acids of corresponding codons
#protein = []
       
              
 #there will be these two classes 
#class Protein:
 #   def __init__(self, sequence):
  #      self.sequence = sequence
  #  def __str__(self):
  #      return self.sequence

#class Person:
#    def __init__(self, wild_type, mutant):
 #       self.wild_type = wild_type
 #       self.mutant = mutant
        
  #def compare_dna_sequence_deletion(self, wild_type, mutant):
# it will return "No deletion in DNA sequence"    or     return f"Deleted sequence: {deleted_seq}, at index: {nt}"

# there will be the following inputs 
#wild_type_seq = input("Enter wild-type DNA sequence (5' to 3'): ")
#mutant_seq = input("Enter mutated DNA sequence (5' to 3'): ")

#the follwoing will be difined and printed 
#print("\n=== DNA Sequence Comparison ===")
#print("Wild-Type DNA: ", wild_type_seq)
#print("Mutant DNA:    ", mutant_seq)
#print(deletion_result)
#mutant_dna = DNA(mutant_seq)
#mutant_mrna = mutant_dna.transcription()
#mutant_protein = mutant_mrna.translation()
#print("\n=== Transcription and Translation ===")
#print("mRNA:          ", mutant_mrna.seq)
#print("Protein:       ", mutant_protein if mutant_protein else "No protein translated")
