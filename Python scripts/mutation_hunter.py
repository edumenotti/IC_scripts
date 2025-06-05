import xml.etree.ElementTree as ET
import argparse
from Bio import AlignIO
from Bio import SeqIO
import os
import tempfile
import hashlib
import csv
import subprocess      
import re
import sys
import cProfile
from collections import defaultdict

class Freq_Calculator:
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-f", "--filename", help="Specify the Blast XML results input file.\n(required)", type=str, required=True)
        parser.add_argument("-wf", "--wholefreq", action='store_true', help='Calculate whole frequency of amino acids in the query sequence.')
        parser.add_argument("-e", "--external", action='store_true', help='Indica que o alinhamento é externo ao banco de dados e fornece o nome do arquivo XML gerado pelo novo blastp.')
        parser.add_argument("-p","--position", type=int, help='Posição específica para calcular a frequência', default=None)
        parser.add_argument("-od", "--output_dir", help="Specify the output directory for results", type=str, default=None)
        parser.add_argument("--alignment_dir", help="Specify the directory for alignment files.", default=os.path.join(os.getcwd(), 'alignment'))
        parser.add_argument("--fasta_dir", help="Specify the directory for fasta files.", default=os.path.join(os.getcwd(), 'fasta_files'))
        parser.add_argument("--blast_result_dir", help="Specify the directory for blast result files.", default=os.path.join(os.getcwd(), 'blast_results'))
        parser.add_argument("--blast2_dir", help="Specify the directory for blast2 output.", default=os.path.join(os.getcwd(), 'blast2'))
        parser.add_argument("--blast_db", help="Specify the directory for the BLASTP database.", default=os.path.join(os.getcwd(), 'blast_db'))
        
        self.args = parser.parse_args()

        # Verifying and creating directories if they don't exist
        for dir_path in [self.args.output_dir, self.args.alignment_dir, self.args.fasta_dir, self.args.blast_result_dir, self.args.blast2_dir]:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

        self.output_dir = self.args.output_dir if self.args.output_dir else os.getcwd()
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def extract_sequences(self, filename):
        unknown_gene_counter = 1
        gene_names = []
        sequence_count = 0

        for record in SeqIO.parse(filename, "fasta"):
            sequence_count += 1  
            name = record.description
            seq = str(record.seq)
            match = re.search(r'(?<=protein=)[^]]*', name)
            if match:
                gene_name = match.group(0)
            else:
                gene_name = f'unknown_gene_{unknown_gene_counter}'
                unknown_gene_counter += 1

            gene_names.append(gene_name)
            output_file = os.path.join(self.fasta_dir, f"{gene_name}.fasta")
            with open(output_file, 'w') as f:
                f.write(f">{name}\n{seq}\n")

            print(f"Processando sequência: {gene_name}")

        print(f"Número total de sequências processadas: {sequence_count}")

        return gene_names   

    def check_ms(self, fasta_file):
        sequence_count = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_count += 1
            if sequence_count > 1:
                return True

        return False
    
    def process_sequences(self, fasta_files, strain_db):
        for fasta_file in fasta_files:
            blast_result_path = self.run_blastp(fasta_file, strain_db)
            
    def run_blastp(self, fasta_file, strain_db):
        unique_id = hashlib.md5(fasta_file.encode()).hexdigest()
        new_xml_filename = f"fasta_{unique_id}.xml"
        blast_result_path = os.path.join(self.blast_result_dir, new_xml_filename)
        if not os.path.exists(blast_result_path):
            subprocess.run([
                "blastp",
                "-query", fasta_file,
                "-db", strain_db,
                "-out", blast_result_path,
                "-outfmt", "5",
                "-num_threads", "8"
            ])
        return blast_result_path
        
    def parse_blast_xml(self, xml_filename, seq_ID=None):
        tree = ET.parse(xml_filename)
        root = tree.getroot()
        
        mutation_count = {}
        mutation_subjects = {}
        mutation_positions = []
        
        for iteration in root.findall('./BlastOutput_iterations/Iteration'):
            for hit in iteration.findall('./Iteration_hits/Hit'):
                cur_hit = Hit()
                if seq_ID:
                    cur_hit.def_ = seq_ID
                else:
                    cur_hit.def_ = hit.find('Hit_def').text
                cur_hit.mismatches = 0 
                cur_hit.mismatch_positions = []

                for hsp in hit.findall('./Hit_hsps/Hsp'):
                    cur_hit.query_seq = hsp.find("Hsp_qseq").text
                    cur_hit.subject_seq = hsp.find("Hsp_hseq").text
                    cur_hit.subject_init = int(hsp.find("Hsp_hit-from").text)

                    for i in range(len(cur_hit.query_seq)):
                        if cur_hit.query_seq[i] != cur_hit.subject_seq[i]:
                            cur_hit.mismatches += 1
                            mismatch_query_aa = cur_hit.query_seq[i]
                            mismatch_ref_aa = cur_hit.subject_seq[i]
                            subject_pos = cur_hit.subject_init + i - cur_hit.subject_seq[:i].count("-")
                            mutation_positions.append(subject_pos)
                            cur_hit.mismatch_positions.append((i+1, mismatch_ref_aa, mismatch_query_aa, subject_pos))
                            
                            mutation = f"{mismatch_ref_aa}{subject_pos}{mismatch_query_aa}"
                            mutation_count[mutation] = mutation_count.get(mutation, 0) + 1
                            if mutation not in mutation_subjects:
                                mutation_subjects[mutation] = []
                            mutation_subjects[mutation].append(cur_hit.def_)
                            
                            

        return mutation_count, mutation_subjects, mutation_positions
        
    def calculate_whole_frequency(self,args):
        # Parse the XML file once to get the query sequence
        tree = ET.parse(self.args.filename)
        root = tree.getroot()

        # Extract the first seq_id for alignment file path
        first_hit_def = root.find(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit[1]/Hit_def").text
        first_seq_id = first_hit_def.split()[0]
        first_alignement = os.path.join(self.alignment_dir, first_hit_def.split()[1] + ".fa.mafft")
        alignment = AlignIO.read(first_alignement, 'clustal')
        
         # Find the sequence in the alignment that matches the seq_ID
        ref_seq = None
        for record in alignment:
            if record.id == first_seq_id:
                ref_seq = record.seq
                break

        # Initialize a list to store amino acid counts for each position
        aa_counts = [{} for _ in range(len(ref_seq))]

        # Count amino acids for each position in one pass
        for record in alignment:
            for pos, aa in enumerate(record.seq):
                if pos < len(aa_counts):  # Verificação adicional para evitar IndexError
                    aa_counts[pos][aa] = aa_counts[pos].get(aa, 0) + 1

        # Write the results to 'wholefrequency.txt'
        with open('wholefrequency.txt', 'w') as f:
            original_pos = 0  # Position in the original sequence without gaps
            for pos, aa in enumerate(ref_seq):
                if aa != '-':  # Skip gaps
                    original_pos += 1  # Increment the position in the original sequence
                    counts = aa_counts[pos]
                    total_count = sum(counts.values())
                    freq = counts.get(aa, 0) / total_count * 100 if total_count > 0 else 0
                    count = counts.get(aa, 0)
                    f.write(f"Position: {original_pos} | Reference Amino Acid: {aa} | Frequency: {freq:.2f}% | Count: {count}\n")
                    f.write("\n")

    def process_input(self):
        if self.args.wholefreq:
            self.calculate_whole_frequency(self.args.filename)

        elif self.args.external or self.args.filename.endswith(('.fa', '.fasta', '.fastq')):
            gene_names = None
            gene_data = {}
            if self.args.filename.endswith(('.fa', '.fasta', '.fastq')):
                if self.check_ms(self.args.filename):
                    # Processamento para arquivos multifasta
                    gene_names = self.extract_sequences(self.args.filename)
                    fasta_files = [os.path.join(self.fasta_dir, f) for f in os.listdir(self.fasta_dir) if f.endswith('.fasta')]
                else:
                    # Processamento para um único arquivo fasta
                    gene_names = [None]
                    fasta_files = [self.args.filename]

                for fasta_file in fasta_files:
                    xml_file = self.run_blastp(fasta_file, self.args.blast_db)
                    gene_name = os.path.basename(fasta_file).split('.')[0] if gene_names != [None] else None
                    gene_data[gene_name] = self.process_ext(xml_file)

            # Tratamento para outros formatos de arquivo
            else:
                gene_names = [None]
                xml_file = self.args.filename
                gene_data[None] = self.process_ext(xml_file)

            return gene_names, gene_data
                
        else:
            mutation_count, mutation_subjects, mutation_positions = self.parse_blast_xml(self.args.filename)
            alignment_file = os.path.join(self.alignment_dir, mutation_subjects[list(mutation_subjects.keys())[0]][0].split(" ")[1] + ".fa.mafft")
            seq_ID = mutation_subjects[list(mutation_subjects.keys())[0]][0].split(" ")[0]
            gene_data[None] = (mutation_count, mutation_subjects, mutation_positions, alignment_file, seq_ID)
            return None, gene_data

            
    def process_ext(self, xml_filename):
        unique_id = hashlib.md5(self.args.filename.encode()).hexdigest()
        new_xml_filename = os.path.join(self.blast2dir, f"new_blast_output_{unique_id}.xml")
            
        # Verifique se o arquivo XML já existe
        if not os.path.exists(new_xml_filename):
            
            # Extraia a sequência do query do arquivo XML original
            tree = ET.parse(xml_filename)
            root = tree.getroot()
            query_seq = root.find(".//BlastOutput_iterations/Iteration/Iteration_hits/Hit[1]/Hit_hsps/Hsp/Hsp_qseq").text
            
            # Crie um arquivo temporário para armazenar a sequência do query
            with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp:
                temp.write(query_seq)
                temp.seek(0)
                
                # Execute um novo BLASTp
                os.system(f"blastp -query {temp.name} -db /home/edumenotti/ic/bancoblast/banco_blast -out {new_xml_filename} -outfmt 5 -num_threads 8")
        
        _, mutation_subjects, _ = self.parse_blast_xml(new_xml_filename)
        if not mutation_subjects:
            print(f"No mutations found for {self.args.filename}")
            sys.exit(0)
        seq_ID = mutation_subjects[list(mutation_subjects.keys())[0]][0].split(" ")[0]
        mutation_count, _, mutation_positions = self.parse_blast_xml(xml_filename, seq_ID=seq_ID)
        alignment_file = os.path.join(self.alignment_dir, mutation_subjects[list(mutation_subjects.keys())[0]][0].split(" ")[1] + ".fa.mafft")
        filename = new_xml_filename
            
        return mutation_count, mutation_subjects, mutation_positions, alignment_file, seq_ID

    def calculate_rare_mutations(self, gene_name, mutation_count, mutation_positions, alignment_file, seq_ID):
        aa_counts = {pos: {} for pos in mutation_positions}
        if os.path.getsize(alignment_file) == 0:
            print(f"The CLUSTAL file {alignment_file} is empty. Terminating process.")
            sys.exit(0)
            
        alignment = AlignIO.read(alignment_file, 'clustal')
        sequence_hash = {record.id: str(record.seq) for record in alignment}

        # Encontrar a sequência de referência
        ref_seq = next((record for record in alignment if record.id == seq_ID), None)
        ref_seq_str = str(ref_seq.seq)

        # Contar aminoácidos para as posições das mutações
        mutation_positions_dict = defaultdict(list)
        for seq_id, seq in sequence_hash.items():
            pos_in_ref_seq = 0
            for pos, aa in enumerate(seq):
                if ref_seq_str[pos] != '-':
                    pos_in_ref_seq += 1

                if pos_in_ref_seq in aa_counts:
                    current_count = aa_counts[pos_in_ref_seq]
                    current_count[aa] = current_count.get(aa, 0) + 1
                    mutation_positions_dict[pos_in_ref_seq].append((seq_id, aa))

        # Calcular informações de mutação
        mutation_info = {}
        for mutation, count in sorted(mutation_count.items(), key=lambda x: int(x[0][1:-1])):
            seq_pos = int(mutation[1:-1])
            query_aa = mutation[0] #Aqui está invertido o que é query e subject para B123A, o original seria query_aa = mutation[0]
            subject_aa = mutation[-1]

            if query_aa == '-' or subject_aa == '-':
                continue
            
            real_count = aa_counts[seq_pos].get(subject_aa, 0)
            total_count = sum(aa_counts[seq_pos].values())
            query_aa_frequency = (aa_counts[seq_pos].get(query_aa, 0) / total_count * 100) if total_count > 0 else 0
            subject_aa_frequency = (aa_counts[seq_pos].get(subject_aa, 0) / total_count * 100) if total_count > 0 else 0

            
            if subject_aa_frequency <= 10 and real_count > 2:    #como a cepa boa está no subject agora, alterei o filtro para ele
                mutation_info[mutation] = {
                    'query_aa_frequency': query_aa_frequency,
                    'subject_aa_frequency': subject_aa_frequency,
                    'real_count': real_count,
                    }

        strains_with_mutations = defaultdict(list)
        for mutation in mutation_info:
            query_aa, pos, subject_aa = mutation[0], int(mutation[1:-1]), mutation[-1]

            for seq_id, aa in mutation_positions_dict[pos]:
                if aa == subject_aa: #Como a cepa boa está no subject, eu alterei para guardar os que tiverem a mutação. Original era query_aa
                    strains_with_mutations[mutation].append(seq_id[:3])

        return mutation_info, strains_with_mutations
    
    def main(self):
        input_filename = os.path.splitext(os.path.basename(self.args.filename))[0]

        output_filename_gapslist = f'{input_filename}_mutations_list.txt'
        output_filename_csv = f'{input_filename}_mutations.csv'

        gene_names, gene_data = self.process_input()

        all_mutation_info = defaultdict(dict)
        all_strains_with_mutations = defaultdict(dict)

        if gene_names and len(gene_names) > 1:
            # Processamento para arquivos multifasta
            for gene_name in gene_names:
                # Obter informações específicas para cada gene
                mutation_count, mutation_subjects, mutation_positions, alignment_file, seq_ID = gene_data[gene_name]

                # Calcular mutações raras e cepas com mutações para cada gene
                gene_mutation_info, gene_strains_with_mutations = self.calculate_rare_mutations(gene_name, mutation_count, mutation_positions, alignment_file, seq_ID)
                all_mutation_info[gene_name] = gene_mutation_info
                all_strains_with_mutations[gene_name] = gene_strains_with_mutations
        else:
            # Processamento para arquivos de sequência única ou externos
            if gene_data.get(None):
                mutation_count, mutation_subjects, mutation_positions, alignment_file, seq_ID = gene_data[None]
                single_mutation_info, single_strains_with_mutations = self.calculate_rare_mutations(None, mutation_count, mutation_positions, alignment_file, seq_ID)
                all_mutation_info['single_sequence'] = single_mutation_info
                all_strains_with_mutations['single_sequence'] = single_strains_with_mutations

        if self.args.filename.endswith(('.fa', '.fasta', '.fastq')):
            if self.check_ms(self.args.filename):
                # Verifique se o diretório de saída existe; crie se não existir

                output_filename_path = os.path.join(self.output_dir, output_filename_gapslist)
                with open(output_filename_path, 'w') as f:
                    for gene_name in gene_names:
                        for mutation, info in all_mutation_info[gene_name].items():
                            f.write(f"Gene: {gene_name} | Mutation: {mutation} | {mutation[0]} Frequency: {info['query_aa_frequency']:.2f}% | {mutation[-1]} Frequency: {info['subject_aa_frequency']:.2f}% | Count: {info['real_count']}\n\n")
                        f.write("\n")
                
                output_filename_csv_path = os.path.join(self.output_dir, output_filename_csv)
                with open(output_filename_csv_path, 'w', newline='') as final_csv_file:
                    writer = csv.writer(final_csv_file)
                    writer.writerow(['Gene', 'Mutation', 'Strains']) 
                    for gene_name in gene_names:
                        individual_csv_file = os.path.join(self.output_dir, f"{gene_name}_mutations.csv")
                        with open(individual_csv_file, 'w', newline='') as csvfile:
                            fieldnames = ['Mutation', 'Strains']
                            individual_writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                            individual_writer.writeheader()
                            for mutation, strains in all_strains_with_mutations[gene_name].items():
                                individual_writer.writerow({
                                    'Mutation': mutation,
                                    'Strains': ';'.join(strains)
                                })
                        with open(individual_csv_file, 'r') as csvfile:
                            next(csvfile)  
                            for line in csvfile:
                                writer.writerow([gene_name] + line.strip().split(','))
                        os.remove(individual_csv_file)
            else:  # para arquivos fasta, fa ou fastq de sequencia unica                          
                output_filename_path = os.path.join(self.output_dir, output_filename_gapslist)
                with open(output_filename_path, 'w') as f:
                    for mutation, info in all_mutation_info['single_sequence'].items():
                        f.write(f"Mutation: {mutation} | {mutation[0]} Frequency: {info['query_aa_frequency']:.2f}% | {mutation[-1]} Frequency: {info['subject_aa_frequency']:.2f}% | Count: {info['real_count']}\n\n")
                
                output_filename_csv_path = os.path.join(self.output_dir, output_filename_csv)
                with open(output_filename_csv_path, 'w', newline='') as csvfile:
                    fieldnames = ['Mutation', 'Strains']
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for mutation, strains in all_strains_with_mutations['single_sequence'].items():
                        writer.writerow({
                            'Mutation': mutation,
                            'Strains': ';'.join(strains)
                        })
        else: # para arquivos em outros formatos que não fasta, fa ou fastq
            output_filename_path = os.path.join(self.output_dir, output_filename_gapslist)            
            with open(output_filename_path, 'w') as f:
                for mutation, info in all_mutation_info['single_sequence'].items():
                    f.write(f"Mutation: {mutation} | {mutation[0]} Frequency: {info['query_aa_frequency']:.2f}% | {mutation[-1]} Frequency: {info['subject_aa_frequency']:.2f}% | Count: {info['real_count']}\n\n")

            output_filename_csv_path = os.path.join(self.output_dir, output_filename_csv)
            with open(output_filename_csv_path, 'w', newline='') as csvfile:
                fieldnames = ['Mutation', 'Strains']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for mutation, strains in all_strains_with_mutations['single_sequence'].items():
                    writer.writerow({
                        'Mutation': mutation,
                        'Strains': ';'.join(strains)
                    })
                   
if __name__ == "__main__":
    freq_calculator = Freq_Calculator()
    freq_calculator.main()
