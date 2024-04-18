import argparse
from Bio import SeqIO

def extract_translated_sequences_with_protein_name(genbank_file):
    """
    Extrai as sequências de aminoácidos traduzidas de todos os CDS de um arquivo GenBank e as salva em arquivos .fasta individuais,
    usando o nome do gene ou locus tag como o nome do arquivo. O cabeçalho de cada arquivo .fasta incluirá [protein="nome do gene"].

    Args:
    - genbank_file (str): Caminho para o arquivo GenBank Flatfile.
    """
    # Ler o arquivo GenBank e extrair as sequências de CDS
    for record in SeqIO.parse(genbank_file, "genbank"):
        # Para cada feature do tipo 'CDS' no arquivo GenBank
        for feature in record.features:
            if feature.type == "CDS":
                # Obter o nome do gene a partir das qualifiers
                gene_name = feature.qualifiers.get('gene', feature.qualifiers.get('locus_tag', ['unknown_gene']))[0]

                # Extrair a sequência de aminoácidos traduzida
                protein_sequence = feature.qualifiers.get('translation', [''])[0]

                # Nome do arquivo baseado no nome do gene ou locus tag
                filename = f"{gene_name}.fasta"

                # Salvar a sequência de aminoácidos em um arquivo .fasta
                with open(filename, "w") as fasta_file:
                    # Incluir [protein="nome do gene"] no cabeçalho
                    fasta_header = f">[protein={gene_name}]\n"
                    fasta_file.write(fasta_header + protein_sequence + "\n")

if __name__ == "__main__":
    # Configurar o parser de argumentos
    parser = argparse.ArgumentParser(description="Extract translated protein sequences from a GenBank file and save them as individual .fasta files with protein name in the header.")
    parser.add_argument("genbank_file", type=str, help="The path to the GenBank file to process.")
    
    # Parsear os argumentos da linha de comando
    args = parser.parse_args()
    
    # Executar a função principal com o arquivo GenBank fornecido
    extract_translated_sequences_with_protein_name(args.genbank_file)
