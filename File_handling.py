# Class representing a DNA sequence and its transcription process
class DNA:
    def __init__(self, seq: str):
        self.seq = seq.upper()

    def transcription(self):
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            mRNA_seq = ''.join(complementary[base] for base in self.seq)
            return MRNA(mRNA_seq)
        except KeyError:
            raise ValueError("Invalid DNA sequence: contains unrecognized bases (should be A, T, C, G).")


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

    def __init__(self, seq: str):
        self.seq = seq

    def translation(self):
        proteins = []
        for i in range(len(self.seq) - 2):
            if self.seq[i:i+3] == "AUG":
                protein = []
                for base in range(i, len(self.seq) - 2, 3):
                    codon = self.seq[base:base + 3]
                    amino_acid = self.codon_table.get(codon, '')
                    if amino_acid == "Stop":
                        break
                    if amino_acid:
                        protein.append(amino_acid)
                if protein:
                    proteins.append("-".join(protein))
        return proteins if proteins else ["No valid proteins translated"]


class Person:
    # Purpose: Initializes a person object with wild-type and mutant DNA sequences.
    # Input: wild_type (type string) - Wild-type DNA sequence
    #        mutant (type string) - Mutant DNA sequence
    def __init__(self, wild_type: str, mutant: str):
        self.wild_type = wild_type
        self.mutant = mutant

    def get_all_proteins(self):
        # Purpose: Translates both wild-type and mutant DNA to proteins.
        # Input: uses self.wild_type and self.mutant
        # Output: Tuple of two lists (wild-type proteins, mutant proteins)
        wt_mrna = DNA(self.wild_type).transcription()
        mt_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mt_mrna:
            return ["Invalid DNA sequence"], ["Invalid DNA sequence"]
        return wt_mrna.translation(), mt_mrna.translation()

    def compare_proteins(self):
        # Purpose: Compares wild-type and mutant proteins and marks mutations.
        # Input: uses self.wild_type and self.mutant
        # Output: List of strings describing protein differences#
        wt_proteins, mt_proteins = self.get_all_proteins()
        max_len = max(len(wt_proteins), len(mt_proteins))

        protein_differences = []
        for i in range(max_len):
            wt_protein = wt_proteins[i] if i < len(wt_proteins) else "None"
            mt_protein = mt_proteins[i] if i < len(mt_proteins) else "None"

            if wt_protein == "None":
                protein_differences.append(f"+ New protein in mutant: {mt_protein}")
            elif mt_protein == "None":
                protein_differences.append(f"- Protein lost in mutant: {wt_protein}")
            elif wt_protein != mt_protein:
                # Highlight mutated proteins with `*`
                highlighted_mutation = [
                    f"*{a}*" if a != b else a
                    for a, b in zip(wt_protein.split("-"), mt_protein.split("-"))
                ]
                protein_differences.append(f"Wild-Type: {wt_protein} -> Mutant: {'-'.join(highlighted_mutation)}")
            else:
                protein_differences.append(f"  {wt_protein}")

        return protein_differences if protein_differences else ["No protein changes detected."]


def read_sequences_from_file(filename: str):
    # Purpose: Reads wild-type and mutant DNA sequences from a given file.
    # Input: filename (type string) in read mode
    # Output: Tuple (wild-type DNA sequence, mutant DNA sequence) or (None, None) if file not found
    try:
        with open(filename, "r") as file:
            lines = file.readlines()
            if len(lines) < 2:
                raise ValueError("Error: File must contain at least two lines (wild-type and mutant DNA sequences).")
            return lines[0].strip().upper(), lines[1].strip().upper()
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None, None


def write_results_to_file(filename: str, results: list[str]):
    # Purpose: Writes the mutation analysis results to a specified file.
    # Input: filename (str) - Name of the file to write results
    # Output: None (writes results to file)
    try:
        with open(filename, "w") as file:
            file.writelines("\n".join(results))
        print(f"Results successfully saved to '{filename}'.")
    except Exception as e:
        print(f"Error writing to file: {e}")


if __name__ == "__main__":
    input_file = "DNA_sample_insertion.txt"
    output_file = "mutation_results.txt"

    wild_type_seq, mutant_seq = read_sequences_from_file(input_file)

    if wild_type_seq and mutant_seq:
        person = Person(wild_type_seq, mutant_seq)
        protein_changes = person.compare_proteins()

        results = ["\n------ DNA Sequence Comparison ------\n",
                   f"Wild-Type DNA: {wild_type_seq}\n",
                   f"Mutant DNA: {mutant_seq}\n",
                   "\n------ Protein Comparison ------\n",
                   "Wild-Type Proteins:\n" + "\n".join(person.get_all_proteins()[0]),
                   "\nMutant Proteins:\n" + "\n".join(person.get_all_proteins()[1]),
                   "\n------ Mutated Proteins ------\n",
                   "\n".join(protein_changes),]

        write_results_to_file(output_file, results)

        print("\n".join(results))
