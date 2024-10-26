#!/usr/bin/env python3

import sys
import argparse
from argparse import RawTextHelpFormatter

version = "1.00"

help = """
VMR Program version """+version+""" - 4 Oct 2024
Obtem codigos de acessos proteicos de codigos de acessos nucleotidicos da tabela VMR
(c) 2024. Rafael Santos da Silva & Arthur Gruber

Usage: VMR.py -i <tabela VMR> -o <tabela output>

-i input<tabela VMR>             Tabela VMR
-o output <file name>      	Output table file
-c <int>					Acession code colunm number
-w <int>					Keywords colunm number
-s <file name>              Sequence FASTA of protein IDs
"""

parser = argparse.ArgumentParser(add_help=False, formatter_class=RawTextHelpFormatter)
parser.add_argument('-i')
parser.add_argument('-o')
parser.add_argument('-h', '--help', action='store_true')
parser.add_argument('-v', '--version', action='store_true')
parser.add_argument('-c', type=int, default= 1)
parser.add_argument('-w', type=int, default= 2)
parser.add_argument('-s')
args = parser.parse_args()



def protein_fetch(accession_code, words):
	proteins=dict()
	from Bio import Entrez
	from Bio import SeqIO
	Entrez.email = "rafass2003@gmail.com"
	handle = Entrez.efetch(db="nucleotide", id=accession_code, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	for f in record.features:		
		if f.type == "CDS":
			for word in words:
				if 'product' in f.qualifiers:
					if word in f.qualifiers["product"][0]:
						#print(f.qualifiers["product"][0])
						#print(f.qualifiers["note"][0])
						proteins[f.qualifiers["protein_id"][0]]=f.qualifiers["translation"][0]
						#break
				if 'note' in f.qualifiers:
					if word in f.qualifiers["note"][0]:
						#print(f.qualifiers["product"][0])
						#print(f.qualifiers["note"][0])
						proteins[f.qualifiers["protein_id"][0]]=f.qualifiers["translation"][0]
						#break
	return proteins

if __name__ == '__main__':
	if not len(sys.argv)>1:
		print(help)
	elif args.help == True:
		print(help)
	elif args.version == True:
		print("""
VMR Program version """+version+""" - 04 Oct 2024
(c) 2024. Rafael Santos da Silva & Arthur Gruber
""")
	else:
		with open(args.i, 'r') as fille:
			lines = fille.readlines()
		with open(args.o, 'w') as arquivo:
			allseq = {}
			for line in lines:
				coluns = line.split(';')
				code = coluns[args.c-1]
				terms = coluns[args.w-1].strip().split(',')
				proteins = protein_fetch(code,terms)
				allseq.update(proteins)
				proteins_code = ','.join(list(proteins.keys()))
				arquivo.write(line.strip() + ';' + proteins_code + '\n')
		if args.s:
			with open(args.s,'w') as fasta:
				for seq in allseq:
					fasta.write('>{}\n{}\n'.format(seq,allseq[seq]))
