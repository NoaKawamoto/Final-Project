class DNA:
#initializes the creation of the DNA sequence object
#ensures that the DNA sequence is in uppercase

    def __init__(self, seq):
        self.seq = seq.upper()

#converts DNA sequence to complementary mRNA strand (transcription
#try and except block to catch possible errors in inputting sequences that are invalid

    def transcription(self):
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            mRNA_seq = ''.join(complementary[base] for base in self.seq)
            return MRNA(mRNA_seq)
        except KeyError as e:
            print('DNA sequence is invalid. Incorrect bases. Try again.', e)

class MRNA:

#dictionary of key-value pairs of codon sequence and the amino acid the codon corresponds to

    codon_table = {'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu',
                   'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu', 'AUU': 'Ile', 'AUC': 'Ile',
                   'AUA': 'Ile', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
                   'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AUG': 'Met',
                   'AGU': 'Ser', 'AGC': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro',
                   'CCG': 'Pro', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
                   'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 'UAU': 'Tyr',
                   'UAC': 'Tyr', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                   'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp',
                   'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu', 'UGU': 'Cys', 'UGC': 'Cys',
                   'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
                   'AGA': 'Arg', 'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly',
                   'GGG': 'Gly', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'}

#initializes the creation of MRNA object

    def __init__(self, seq):
        self.seq = seq

#Protein variable is declared outside the for loop which will be appended with Amino Acids of corresponding codons
#seq[r:r+3] goes over sequence of codons in increments of 3 (2 values after r)
    def translation(self):
        protein = []
        start_index = -1  # To track the start codon position

        # Step 1: Find the first "AUG" (start codon)
        for i in range(len(self.seq) - 2):
            if self.seq[i:i + 3] == "AUG":
                start_index = i
                break

                # If no start codon is found, return an invalid sequence message
        if start_index == -1:
            return "Invalid Sequence (No Start Codon)"

        # Step 2: Process codons from the start codon in steps of 3
        for base in range(start_index, len(self.seq) - 2, 3):
            codon = self.seq[base:base + 3]
            amino_acid = self.codon_table.get(codon, '')

            if amino_acid == "Stop":
                break  # Stop translation if a stop codon is reached

            if amino_acid:
                protein.append(amino_acid)  # Add amino acid to protein chain

        return Protein("-".join(protein)) if protein else "Invalid Sequence"

class Protein:
    def __init__(self, sequence):
        self.sequence = sequence
    def __str__(self):
        return self.sequence

class Person:
    def __init__(self, wild_type, mutant):
        self.wild_type = wild_type
        self.mutant = mutant

    def compare_dna_sequence_deletion(self, wild_type, mutant):
        self.wild_type = wild_type
        self.mutant = mutant
        len_wt = len(wild_type)
        len_mut = len(mutant)
        if len_wt == len_mut:
            return "No deletion in DNA sequence"

        for nt in range(len(mutant)):
            if wild_type[nt] != mutant[nt]:
                break
        else:
            nt = len(mutant)  # If loop completes, deletion is at the end

        # Extract deleted sequence
        deleted_seq = wild_type[nt: len_wt - (len_mut - nt)]
        return f"Deleted sequence: {deleted_seq}, at index: {nt}"

wild_type_seq = input("Enter wild-type DNA sequence (5' to 3'): ").strip().upper()
mutant_seq = input("Enter mutated DNA sequence (5' to 3'): ").strip().upper()

person = Person(wild_type_seq, mutant_seq)
deletion_result = person.compare_dna_sequence_deletion(wild_type_seq, mutant_seq)

print("\n=== DNA Sequence Comparison ===")
print("Wild-Type DNA: ", wild_type_seq)
print("Mutant DNA:    ", mutant_seq)
print(deletion_result)

mutant_dna = DNA(mutant_seq)
mutant_mrna = mutant_dna.transcription()
mutant_protein = mutant_mrna.translation()

print("\n=== Transcription and Translation ===")
print("mRNA:          ", mutant_mrna.seq)
print("Protein:       ", mutant_protein if mutant_protein else "No protein translated")


