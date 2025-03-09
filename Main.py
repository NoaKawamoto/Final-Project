class DNA:
    def __init__(self, seq):
        self.seq = seq.upper()

    def transcription(self):
        complementary = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        try:
            mRNA_seq = ''.join(complementary[base] for base in self.seq)
            return MRNA(mRNA_seq)
        except KeyError:
            print("Invalid DNA sequence: does not code for AA")
        except ValueError:
            print("Invalid DNA sequence: contains invalid base.")  # Invalid base in the sequence
        except Exception as e:
            raise Exception(f"An error occurred during transcription: {str(e)}")


class MRNA:
    codon_table = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser', 'CCU': 'Pro', 'CCC': 'Pro',
        'CCA': 'Pro', 'CCG': 'Pro', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr', 'GCU': 'Ala', 'GCC': 'Ala',
        'GCA': 'Ala', 'GCG': 'Ala', 'UAU': 'Tyr', 'UAC': 'Tyr', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGG': 'Trp', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg',
        'AGG': 'Arg', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGA': 'Stop'
    }

    def __init__(self, seq):
        self.seq = seq

    def translation(self):
        proteins = []
        for start_index in range(len(self.seq) - 2):
            if self.seq[start_index:start_index + 3] == "AUG":
                protein = []
                for base in range(start_index, len(self.seq) - 2, 3):
                    codon = self.seq[base:base + 3]
                    amino_acid = self.codon_table.get(codon, '')

                    if amino_acid == "Stop":
                        break  # Stop translation if a stop codon is reached
                    if amino_acid:
                        protein.append(amino_acid)

                if protein:
                    proteins.append("-".join(protein))

        return proteins if proteins else ["No valid proteins translated"]


class Person:
    def __init__(self, wild_type, mutant):
        self.wild_type = wild_type
        self.mutant = mutant

    def check_frameshift(self):
        difference = abs(len(self.wild_type) - len(self.mutant))
        if difference == 0:
            return "No mutation detected."
        if difference % 3 == 0:
            return "No frameshift detected (mutation is in-frame)."
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
        # First, check for a frameshift mutation
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

        return "Missense mutation (amino acid change detected)."

    def get_all_proteins(self):
        wt_mrna = DNA(self.wild_type).transcription()
        mt_mrna = DNA(self.mutant).transcription()
        if not wt_mrna or not mt_mrna:
            return "Invalid DNA sequence."

        return wt_mrna.translation(), mt_mrna.translation()


def read_input_file(file_path):
    try:
        with open(file_path, 'r') as file:
            wild_type_seq = file.readline().strip()
            mutant_seq = file.readline().strip()
        return wild_type_seq, mutant_seq
    except FileNotFoundError:
        raise Exception(f"File {file_path} not found.")
    except Exception as e:
        raise Exception(f"Error reading file {file_path}: {str(e)}")


def write_output_file(file_path, content):
    try:
        with open(file_path, 'w') as file:
            file.write(content)
    except Exception as e:
        raise Exception(f"Error writing to file {file_path}: {str(e)}")


# Main program
def main():
    try:
        # Read from input file
        wild_type_seq, mutant_seq = read_input_file('DNA_sample.txt')

        # Create Person object
        person = Person(wild_type_seq, mutant_seq)

        # Mutation analysis
        mutation_analysis = "\n--- Mutation Analysis ---\n"
        mutation_analysis += person.check_frameshift() + "\n"
        mutation_analysis += person.detect_mutation_type() + "\n"

        # DNA Sequence Comparison
        dna_comparison = "\n--- DNA Sequence Comparison ---\n"
        dna_comparison += f"Wild-Type DNA: {wild_type_seq}\n"
        dna_comparison += f"Mutant DNA: {mutant_seq}\n"
        dna_comparison += person.compare_dna_sequence_deletion() + "\n"
        dna_comparison += person.compare_dna_sequence_insertion() + "\n"

        # Transcription and Translation
        mutant_dna = DNA(mutant_seq)
        mutant_mrna = mutant_dna.transcription()
        transcription_translation = "\n--- Transcription and Translation ---\n"
        if mutant_mrna:
            mutant_protein = mutant_mrna.translation()
            transcription_translation += f"mRNA: {mutant_mrna.seq}\n"
            transcription_translation += f"Protein: {mutant_protein}\n"
        else:
            transcription_translation += "Invalid DNA sequence detected. Transcription failed.\n"

        # Protein Comparison
        wt_proteins, mt_proteins = person.get_all_proteins()
        protein_comparison = "\n--- Protein Comparison ---\n"
        protein_comparison += f"Wild-Type Proteins: {wt_proteins}\n"
        protein_comparison += f"Mutant Proteins:   {mt_proteins}\n"

        # Combine all content
        full_content = mutation_analysis + dna_comparison + transcription_translation + protein_comparison

        # Write results to output file
        write_output_file('output.txt', full_content)
        print("Results have been written to 'output.txt'.")

    except Exception as e:
        print(f"Error: {e}")


# Run the program
if __name__ == "__main__":
    main()