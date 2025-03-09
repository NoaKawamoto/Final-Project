class DNA:
    def __init__(self, seq):
        self.seq = seq.upper()

    def transcription(self):
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            mRNA_seq = ''.join(complementary[base] for base in self.seq)
            return MRNA(mRNA_seq)
        except KeyError:
            raise ValueError("Invalid DNA sequence: contains unrecognized bases (should be A, T, C, G).")


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

    def translation(self):
        start_index = self.seq.find("AUG")
        if start_index == -1:
            return ["No valid proteins translated"]

        protein = []
        for base in range(start_index, len(self.seq) - 2, 3):
            codon = self.seq[base:base + 3]
            amino_acid = self.codon_table.get(codon, '')

            if amino_acid == "Stop":
                break
            if amino_acid:
                protein.append(amino_acid)

        return ["-".join(protein)] if protein else ["No valid proteins translated"]


class Person:
    def __init__(self, wild_type, mutant):
        self.wild_type = wild_type
        self.mutant = mutant

    def check_frameshift(self):
        difference = abs(len(self.wild_type) - len(self.mutant))
        if difference == 0:
            return "No mutation detected."
        if difference % 3 == 0:
            return "Mutation detected but is in-frame."
        return "Frameshift mutation detected!"

    def compare_dna_sequence_insertion(self):
        if len(self.mutant) <= len(self.wild_type):
            return "No insertion detected"

        for i in range(len(self.wild_type)):
            if i >= len(self.mutant) or self.wild_type[i] != self.mutant[i]:
                inserted_seq = self.mutant[i:len(self.mutant) - (len(self.wild_type) - i)]
                return f"Inserted nucleotides: {inserted_seq} at nt {i+1}"

        return f"Insertion detected at the end of the sequence at nt {len(self.wild_type)+1}"

    def compare_dna_sequence_deletion(self):
        if len(self.wild_type) <= len(self.mutant):
            return "No deletion detected"

        for i in range(len(self.mutant)):
            if i >= len(self.wild_type) or self.wild_type[i] != self.mutant[i]:
                deleted_seq = self.wild_type[i:len(self.wild_type) - (len(self.mutant) - i)]
                return f"Deleted nucleotides: {deleted_seq} at nt {i+1}"

        return f"Deletion detected at the end of the sequence at nt {len(self.mutant)+1}"

    def detect_mutation_type(self):
        frameshift_result = self.check_frameshift()
        if "Frameshift" in frameshift_result:
            return "Frameshift mutation detected. Point mutation is not possible."

        wt_mrna = DNA(self.wild_type).transcription()
        mut_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mut_mrna:
            return "Invalid DNA sequence."

        wt_protein = wt_mrna.translation()
        mut_protein = mut_mrna.translation()

        if wt_protein == mut_protein:
            return "No change in protein sequence."

        return "Missense mutation (amino acid change detected)."

    def get_all_proteins(self):
        wt_mrna = DNA(self.wild_type).transcription()
        mt_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mt_mrna:
            return "Invalid DNA sequence."

        return wt_mrna.translation(), mt_mrna.translation()


# ** Ensure File Handling is Separate from Unit Testing **
if __name__ == "__main__":
    input_file = "DNA_sample.txt"
    output_file = "mutation_results.txt"

    def read_sequences_from_file(filename):
        try:
            with open(filename, "r") as file:
                lines = file.readlines()
                if len(lines) < 2:
                    raise ValueError("File must contain two lines: wild-type and mutant sequences.")

                return lines[0].strip().upper(), lines[1].strip().upper()
        except FileNotFoundError:
            print(f"Error: The file '{filename}' was not found.")
            return None, None

    def write_results_to_file(filename, results):
        try:
            with open(filename, "w") as file:
                file.write(results)
            print(f"Results saved to '{filename}'.")
        except Exception as e:
            print(f"Error writing to file: {e}")

    wild_type_seq, mutant_seq = read_sequences_from_file(input_file)

    if wild_type_seq and mutant_seq:
        person = Person(wild_type_seq, mutant_seq)

        results = [
            "\n------ DNA Sequence Comparison ---",
            f"Wild-Type DNA: {wild_type_seq}",
            f"Mutant DNA: {mutant_seq}",
            person.compare_dna_sequence_deletion(),
            person.compare_dna_sequence_insertion(),
            "\n--- Mutation Analysis ---",
            person.check_frameshift(),
            person.detect_mutation_type(),
            "\n--- Protein Comparison ---",
            f"Wild-Type Proteins: {', '.join(person.get_all_proteins()[0])}",
            f"Mutant Proteins: {', '.join(person.get_all_proteins()[1])}"
        ]

        write_results_to_file(output_file, "\n".join(results))
