# Class representing a DNA sequence and its transcription process
class DNA:
    #Initializes the DNA sequence and converts it to uppercase.

    def __init__(self, seq):
        self.seq = seq.upper()  # Convert sequence to uppercase

    def transcription(self):
        # Transcribes the DNA to mRNA by replacing bases: A->U, T->A, C->G, G->C
        # Purpose: Transcribes the DNA sequence into an mRNA sequence.
        # Input: None (operates on self.seq)
        # Output: MRNA object with the transcribed sequence.
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            mRNA_seq = ''.join(complementary[base] for base in self.seq)
            return MRNA(mRNA_seq)  # Return MRNA object with transcribed sequence
        except KeyError:
            print("Invalid DNA sequence: contains non-standard base(s).")
            return None


# Class representing an mRNA sequence and its translation to protein
class MRNA:
    codon_table = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro',
        'CCA': 'Pro', 'CCG': 'Pro', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'GCU': 'Ala', 'GCC': 'Ala',
        'GCA': 'Ala', 'GCG': 'Ala', 'UAU': 'Tyr', 'UAC': 'Tyr', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg',
        'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    def __init__(self, seq):
        self.seq = seq
        # Purpose: Initializes an mRNA sequence.
        # Input: seq (str) - mRNA sequence
    def translation(self):
        # Translates mRNA sequence into a protein sequence
        # Purpose: Translates the mRNA sequence into a protein sequence.
        # Input: uses self.seq
        # Output: List of translated proteins as strings
        proteins = []
        for start_index in range(len(self.seq) - 2):
            if self.seq[start_index:start_index + 3] == "AUG":  # Start translation at AUG
                protein = []
                for base in range(start_index, len(self.seq) - 2, 3):
                    codon = self.seq[base:base + 3]
                    amino_acid = self.codon_table.get(codon, '')

                    if amino_acid == "Stop":
                        break  # Stop translation if stop codon is reached
                    if amino_acid:
                        protein.append(amino_acid)

                if protein:
                    proteins.append("-".join(protein))  # Join amino acids with a hyphen

        return proteins if proteins else ["No valid proteins translated"]


# Class to handle mutation analysis between wild-type and mutant DNA sequences
class Person:
    def __init__(self, wild_type, mutant):
        self.wild_type = wild_type
        self.mutant = mutant
        # Purpose: Initializes a person object with wild-type and mutant DNA sequences.
        # Input: wild_type (str) = Wild-type DNA sequence
        #        mutant (str) = Mutant DNA sequence

    def check_frameshift(self):
        # Purpose: Checks if the mutation is a frameshift mutation.
        # Input: uses self.wild_type and self.mutant
        # Output: A message telling whether a frameshift mutation occurred
        difference = abs(len(self.wild_type) - len(self.mutant))
        if difference == 0:
            return "No mutation detected."
        if difference % 3 == 0:
            return "No frameshift detected (mutation is in-frame)."
        return "Frameshift mutation detected!"

    def compare_dna_sequence_insertion(self):
        # Purpose: Determines the type of mutation (frameshift, missense, or nonsense).
        # Input: uses self.wild_type and self.mutant
        # Output: A string describing the type of mutation
        if len(self.mutant) <= len(self.wild_type):
            return "No insertion detected"

        for i in range(len(self.wild_type)):
            if i >= len(self.mutant) or self.wild_type[i] != self.mutant[i]:
                inserted_seq = self.mutant[i:len(self.mutant) - (len(self.wild_type) - i)]
                return f"Inserted nucleotides: {inserted_seq} at nucleotide {i+1}"

        return f"Insertion detected at the end of the sequence at nucleotide {len(self.wild_type)+1}"

    def compare_dna_sequence_deletion(self):
        if len(self.wild_type) <= len(self.mutant):
            return "No deletion detected"

        for i in range(len(self.mutant)):
            if i >= len(self.wild_type) or self.wild_type[i] != self.mutant[i]:
                deleted_seq = self.wild_type[i:len(self.wild_type) - (len(self.mutant) - i)]
                return f"Deleted nucleotides: {deleted_seq} at nucleotide {i+1}"

        return f"Deletion detected at the end of the sequence at nucleotide {len(self.mutant)+1}"

    def detect_mutation_type(self):
        frameshift_result = self.check_frameshift()
        if "Frameshift mutation detected" in frameshift_result:
            return "Frameshift mutation detected. Point mutation is not possible."

        wt_mrna = DNA(self.wild_type).transcription()
        mut_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mut_mrna:
            return "Invalid DNA sequence."

        wt_protein = wt_mrna.translation()
        mut_protein = mut_mrna.translation()

        if wt_protein == mut_protein:
            return "No change in protein sequence"

        if any("Stop" in p for p in mut_protein) and not any("Stop" in p for p in wt_protein):
            return "Nonsense mutation (premature stop codon detected)."

        return "Missense mutation (amino acid change detected)."

    def get_all_proteins(self):
        # Purpose: Retrieves the translated protein sequences for both wild-type and mutant DNA.
        # Input: uses self.wild_type and self.mutant
        # Output: A tuple of two lists containing the wild-type and mutant protein sequences.
        wt_mrna = DNA(self.wild_type).transcription()
        mt_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mt_mrna:
            return "Invalid DNA sequence.", "Invalid DNA sequence."

        return wt_mrna.translation(), mt_mrna.translation()


# User input
wild_type_seq = input("Enter wild-type DNA sequence (5' to 3'): ").strip().upper()
mutant_seq = input("Enter mutated DNA sequence (5' to 3'): ").strip().upper()

# Create a Person object with the input sequences
person = Person(wild_type_seq, mutant_seq)

# Display mutation analysis results
print("\n--- Mutation Analysis ---")
print("Frameshift:", person.check_frameshift())
print("Mutation Type:", person.detect_mutation_type())

# DNA sequence comparison
print("\n--- DNA Sequence Comparison ---")
print(person.compare_dna_sequence_deletion())
print(person.compare_dna_sequence_insertion())

# Transcription and translation
mutant_dna = DNA(mutant_seq)
mutant_mrna = mutant_dna.transcription()

if mutant_mrna:
    mutant_protein = mutant_mrna.translation()
    print("\n--- Transcription and Translation ---")
    print(f"mRNA: {mutant_mrna.seq}")
    print(f"Protein: {mutant_protein}")
else:
    print("\nInvalid DNA sequence detected. Transcription failed.")

# Protein comparison between wild-type and mutant sequences
wt_proteins, mt_proteins = person.get_all_proteins()

print("\n--- Protein Comparison ---")
print("Wild-Type Proteins:", wt_proteins)
print("Mutant Proteins:   ", mt_proteins)
