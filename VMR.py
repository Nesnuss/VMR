#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse
from Bio import Entrez
from Bio import SeqIO
from argparse import RawTextHelpFormatter
import time
from io import StringIO
import subprocess
import http.client 


version = "1.01"
help = """
VMR Program version """+version+""" - 4 Oct 2024
Obtem codigos de acessos proteicos apartir de codigos de acessos nucleotidicos da tabela VMR
(c) 2024. Rafael Santos da Silva & Arthur Gruber

Usage: VMR.py -i <tabela VMR> -o <tabela output>

-i input<VMR table>         VMR table
-o output <file name>      	Output table file
-s <file name>              auxiliary table
"""

parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument('-o')
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
#parser.add_argument('-c', type=int, default= 22)
#parser.add_argument('-w', type=int, default= 27)
parser.add_argument('-s')
args = parser.parse_args()



#função responsavel por pegar acession_code_nuccore e retornar acession_code_protein especificos
def busca_entrez(nuc_acc, query_terms, negative_terms, min_len, max_len):
    Entrez.email = "rafass2003@gmail.com"
    Entrez.api_key = "511ef882e71fdff1e01eaaa3177e47c43e09"

    for attempt in range(3):
        try:
            handle = Entrez.esearch(db="nuccore", term=nuc_acc, usehistory='y')
            record = Entrez.read(handle)
            handle.close()

            webenv = record['WebEnv']
            query_key = record['QueryKey']

            pesquisa_elink = Entrez.elink(dbfrom="nuccore", db='protein', query_key=query_key, webenv=webenv, linkname="nuccore_protein")
            Resultado_elink = Entrez.read(pesquisa_elink)
            pesquisa_elink.close()

            protein_ids = []
            for link in Resultado_elink[0].get("LinkSetDb", [{}])[0].get("Link", []):
                protein_ids.append(link.get("Id", ""))

            filtered_ids = []

            separated_terms = []
            for item in negative_terms:
                if pd.isna(item):
                    continue
                separated_terms.extend([term.strip() for term in str(item).split(",")])

            len_search = f'"{min_len}"[SLEN] : "{max_len}"[SLEN]'
            search = f"({'[uid] OR '.join(protein_ids)}[uid]) AND ({' OR '.join(query_terms)}) AND ({len_search}) NOT ({' OR '.join(separated_terms)})"
            print(search)
            filtro = Entrez.esearch(db="protein", term=search)
            resultado_filtro = Entrez.read(filtro)
            filtro.close()

            for id in resultado_filtro.get("IdList", []):
                filtered_ids.append(id)

            if not filtered_ids:
                print("Nenhum ID filtrado encontrado. Pulando...")
                return []

            with Entrez.efetch(db="protein", id=",".join(filtered_ids), rettype="gb", retmode="text") as pesquisa_efetch:
                response_data = pesquisa_efetch.read()

            seqio = SeqIO.parse(StringIO(response_data), "genbank")
            extracted_ids = [record.id for record in seqio]
            print(f"IDs extraídos: {extracted_ids}")
            return extracted_ids[0]

        except (http.client.IncompleteRead, ValueError) as e:
            print(f"Erro ao buscar dados: {e}. Tentando novamente...")
            time.sleep(1)  

    print("Falha ao obter dados após várias tentativas.")
    return []
# recuperar o fasta do genoma do virus
def fasta(nuc_acc, path):
    Entrez.email = "rafass2003@gmail.com"
    Entrez.api_key = "511ef882e71fdff1e01eaaa3177e47c43e09"

    # Formata o nome do arquivo de saída
    output_file = f"{path}/{nuc_acc}_prot.fasta"

    # Verifica se o arquivo já existe
    if os.path.exists(output_file):
        print(f"O arquivo {output_file} já existe. Usando arquivo existente.")
        return output_file

    handle = Entrez.esearch(db="nuccore", term=nuc_acc, retmax=10000, usehistory='y')
    record = Entrez.read(handle)
    handle.close()

    webenv = record['WebEnv']
    query_key = record['QueryKey']

    pesquisa_elink = Entrez.elink(dbfrom="nuccore", db='protein', query_key=query_key, webenv=webenv, linkname="nuccore_protein", usehistory='y')
    Resultado_elink = Entrez.read(pesquisa_elink)
    pesquisa_elink.close()

    idlist = Resultado_elink[0].get("LinkSetDb", [{}])[0].get("Link", [])

    fasta_records = []
    
    for id in idlist:
        with Entrez.efetch(db="protein", id=id.values(), rettype="fasta", retmode="text") as pesquisa_efetch:
            response_data = pesquisa_efetch.read()
            seqio = SeqIO.parse(StringIO(response_data), "fasta")
            time.sleep(1) 
        for seq_record in seqio:
            fasta_records.append(seq_record)  

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")

    return output_file
# baixa os fastas de proteinas de uma familia 
def refdb(name_protein, taxid, familia, path, pn):
    Entrez.email = "rafass2003@gmail.com"
    Entrez.api_key = "511ef882e71fdff1e01eaaa3177e47c43e09"

    # Formata o nome do arquivo de saída
    if len(pn) > 1: 
        names = f'{pn[0]}_{pn[1]}' 
    else:    
        names = f'{pn[0]}'
    
    output_file = f"{path}/{familia}_{names}.fasta"

    # Verifica se o arquivo já existe
    if os.path.exists(output_file):
        print(f"O arquivo {output_file} já existe. Usando arquivo existente.")
        return output_file

    search = f"{name_protein} AND txid{taxid}[Organism]"
    print(search)

    handle = Entrez.esearch(db="ipg", term=search, retmax=10000, usehistory='y')
    record = Entrez.read(handle)
    handle.close()

    idlist = record['IdList']
    print(idlist)
    
    if not idlist:
        print("Nenhum ID encontrado.")
        return []

    fasta_records = []
    
    for id in idlist:
        print(id)
        print(type(id))
        for attempt in range(3):
            try:
                with Entrez.efetch(db="ipg", id=id, rettype="fasta", retmode="text") as pesquisa_efetch:
                    response_data = pesquisa_efetch.read()
                    seqio = SeqIO.parse(StringIO(response_data.decode('utf-8')), "fasta") 
                    #time.sleep(8)        
                for seq_record in seqio:
                    fasta_records.append(seq_record) 
                break 
            except (http.client.IncompleteRead, ValueError) as e:
                print(f"Erro ao ler dados: {e}. Tentando novamente...")
                time.sleep(1)  

    with open(output_file, "w") as output_handle:
        SeqIO.write(fasta_records, output_handle, "fasta")

    return output_file
# Reculpera o taxid da familia
def txid(family):
    
    Entrez.email = "rafass2003@gmail.com"

    handle = Entrez.esearch(db="taxonomy", term= family)
    record = Entrez.read(handle)
    handle.close()

    taxid =str(record["IdList"][0]) if record["IdList"] else None

    return taxid
#Execulta os comando da função blast_plus
def auxiliar(command):
    """Executa um comando e verifica se foi bem-sucedido."""
    try:
        result = subprocess.run(command, check=True, text=True)
        print(f"Comando executado com sucesso: {' '.join(command)}")
        return result
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar o comando: {' '.join(command)}")
        print(e)
        return None
# Crias linhas de comando que seram rodados na shell
def blast_plus(database_fasta,genome):
    makeblastdb_cmd = [
        "makeblastdb",
        "-in", database_fasta,
        "-dbtype", "prot",
        "-out", database_fasta
    ]
    auxiliar(makeblastdb_cmd)


    blastp_cmd = [
        "blastp",
        "-query", genome,
        "-db", database_fasta,
        "-out", "resultado_blast.txt",
        "-evalue", "1e-20",
        "-num_alignments", "1",
        "-outfmt", ""'6 qseqid '""
    ]
    auxiliar(blastp_cmd)

    if os.path.exists("resultado_blast.txt"):
        with open("resultado_blast.txt", "r") as leitura:
            lines = leitura.readlines()
            if lines:  # Verifica se há linhas no arquivo
                protein = lines[0].strip()
                return protein
            else:
                print("Nenhum resultado encontrado no arquivo resultado_blast.txt.")
                return None
    else:
        print("O arquivo resultado_blast.txt não foi encontrado.")
        return None
    os.remove("resultado_blast.txt")

    return protein

inicio_tempo = time.perf_counter()
# Seleciona os arquivos que seram ultilizados
print("Lendo os dados...")
tabelaX = pd.DataFrame(pd.read_csv(args.i, delimiter=';'))
tabelaY = pd.read_csv(args.s, delimiter=';')
# Selecionando as colunas 
termos_positivos = tabelaY['Positive_terms']
termos_negative = tabelaY['Negative_terms']
min_len = tabelaY['min_length']
max_len = tabelaY['max_length']
# Crindo nova tabela
nova_tabelaX = []

#nao_anotado = []
#Loop central; roda tando o protocolo defalt, quanto o protocolo anotação funcional
for i in range(len(tabelaX)):
    for j in range(len(tabelaY)):
        line = dict(tabelaX.iloc[i])
        code = line.get("Virus GENBANK accession", "")
        familia = line.get("Family", "")
        taxiadi = txid(familia)
        print(f"Processando código: {code}")
        time.sleep(1) 
# Protocolo defalt: reculpera as acession_code_protein
        proteina = busca_entrez(code, [termos_positivos[j]], [termos_negative[j]], min_len[j], max_len[j])
        #print(proteina)
        #print(termos_positivos[j])
# Protocolo anotação funcional: identifica proteinas por similariedade e reculpera o acession_code_protein
        if proteina == []:
            print("DEU ERRADOOOOO!!")

            protein_name = termos_positivos[j].split()
            
            if len(protein_name) > 1: 
                diretorio = f'refdb3/{familia}/{protein_name[0]}_{protein_name[1]}' 
            else:    
                diretorio = f'refdb3/{familia}/{protein_name[0]}'

            genome_file = f'cds_virus_fasta3/{familia}'
            #diretorio = f'refdb/{familia}/"{termos_positivos[j]}"'
            os.makedirs(diretorio, exist_ok=True)
            os.makedirs(genome_file, exist_ok=True)
            fasta_file = fasta(code, genome_file)
            database = refdb(termos_positivos[j],taxiadi,familia, diretorio, protein_name)
            #print(fasta_file)
            #print(database)
            proteina = blast_plus(database,fasta_file)
            print(proteina)
            #print(type(proteina))


        line['min_length'] = min_len[j]
        line['max_length'] = max_len[j]
        line['Negative_Terms'] = termos_negative[j]
        line['Positive_Terms'] = termos_positivos[j]
        line['Protein_codes'] = proteina if proteina else ""
        line['family_taxid'] = taxiadi
        nova_tabelaX.append(line)


print("Salvando os resultados...")
nova_tabelaX = pd.DataFrame(nova_tabelaX)
nova_tabelaX.to_csv(args.o, sep=';', index=False)

df = pd.read_csv(args.o, delimiter=';')

nova_ordem = ['Sort', 'Isolate Sort', 'Realm', 'Subrealm', 'Kingdom', 'Subkingdom', 'Phylum', 'Subphylum', 'Class', 'Subclass', 'Order', 'Suborder', 'Family','family_taxid', 'Subfamily', 'Genus', 'Subgenus', 'Species', 'Exemplar or additional isolate', 'Virus name(s)', 'Virus name abbreviation(s)', 'Virus isolate designation', 'Virus GENBANK accession', 'Positive_Terms','Negative_Terms','min_length','max_length', 'Protein_codes', 'Virus REFSEQ accession', 'Genome coverage', 'Genome composition', 'Host source']  
df = df[nova_ordem]
df.to_csv(args.o, index=False, sep=';')

fim_tempo = time.perf_counter()

tempo_total = fim_tempo - inicio_tempo


# Convertendo o tempo total em horas, minutos e segundos

horas, resto = divmod(tempo_total, 3600)

minutos, segundos = divmod(resto, 60)


print(f"Tempo total de execução: {int(horas)} horas, {int(minutos)} minutos e {int(segundos)} segundos")

print("Processo concluído!")