class DNA:
    def __init__(self, seq):
        self.seq = seq.upper()
    def transcription(self):
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        mRNA_seq = ''.join(complementary[base] for base in self.seq)
        return MRNA(mRNA_seq)

class MRNA:
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
    def __init__(self, seq):
        self.seq = seq.upper()
    def translation(self):
        protein = []
        found_start_codon = False
        for r in range(0, len(self.seq), 3):
            codon = self.seq[r:r+3]
            amino_acid = self.codon_table.get(codon, '')

            if amino_acid == "Met":
                found_start_codon = True

            if found_start_codon:
                if amino_acid == "Stop":
                        break
                if amino_acid:
                        protein.append(amino_acid)
        return Protein("-".join(protein)) if protein else "Invalid Sequence"


class Protein:
    def __init__(self, sequence):
        self.sequence = sequence
    def __str__(self):
        return self.sequence

class Person:
    def __init__(self, wild_type, mutation):
        self.wild_type = wild_type
        self.mutation = mutation
    def __str__(self):
        return self.wild_type + self.mutation


dna_sequence = input("Enter DNA template strand (5' to 3'): ")
dna = DNA(dna_sequence)
mrna = dna.transcription()
protein = mrna.translation()

print("DNA:", dna.seq)
print("mRNA:", mrna.seq)
print("Protein:", protein if protein else "No valid protein translated")
