import os
import re
import sys
import argparse
import pandas as pd
from Bio import SeqIO
import subprocess as sp
from ete3 import NCBITaxa
sys.path.append(r'/davidb/yatirsolan/scripts/python/bio_utilities')
import mmseq_taxonomy
import phylogenetics
import fasta_import
import databases
import HMM

#####################################################################

def contig_clutser(gnrl_anlz_file):
    scr_path = r'/davidb/yatirsolan/yatirsBioUtils/bioutils/scripts/qmmseqs.cluster'
    base_name = gnrl_anlz_file.replace('.fa','')
    sp.run(f'{scr_path} {gnrl_anlz_file} 8 {base_name}_clu 7.5 0.5 --cluster-mode 2 --cov-mode 1', shell=True, capture_output=True, text=True)
    return f'{base_name}_clu_rep_seq.fasta'

#####################################################################

def mmseq_taxonomy_make(cntgs_fasta):
    tax_res = mmseq_taxonomy.main(cntgs_fasta, queue='dudulight')
    return tax_res

#####################################################################

def WGS_Genomes_tax_annotation(contigs_file):
    ncbi = NCBITaxa()
    contigs = [contig.split('\n')[0] for contig in open(contigs_file).readlines()]
    reg = re.compile(r'(?<=\|)GCA_\d+\.*.(?=\|)')
    contig_gcas = {contig:reg.search(contig).group() for contig in contigs}
    gca_taxid = databases.gca_taxids_dic()
    taxid_lin = phylogenetics.taxid_lineage()

    df_dict = {col:list() for col in ['contigs', 'taxid', 'rank', 'tax']}
    for contig in contigs:
        gca = contig_gcas.get(contig)
        taxid = gca_taxid.get(gca)

        try:
            df_dict.get('rank').append(list(ncbi.get_rank([taxid]).values())[0])
            df_dict.get('tax').append(list(ncbi.get_taxid_translator([taxid]).values())[0])
        except IndexError : 
            continue # in case of a un-identified taxid this contig will be excluded, later this contig will be annotated as 'no_mmseq_result' (in the 'review_summary.py' script)
        df_dict.get('contigs').append(contig)
        df_dict.get('taxid').append(taxid)
    tax_res = contigs_file.replace('.headers', '_TaxResult.tsv')
    pd.DataFrame.from_dict(df_dict).to_csv(tax_res, header=None, index=None, sep='\t')
    return tax_res

#####################################################################

def main(gnrl_anlz_file, cluster=False): # The cluster option is OFF and never was used.
    WGS_Gnms = True if 'WGS_Genomes' in os.getcwd() else False

    cntg_set = set(pd.read_table(gnrl_anlz_file).contig)

    gnrl_anlz_file = os.path.basename(gnrl_anlz_file)
    gnrl_anlz_file = os.path.splitext(gnrl_anlz_file)[0]

    fasta_import_output = fasta_import.main(hdrs=cntg_set, 
                                            output_type='fa', 
                                            output_file_name=gnrl_anlz_file, 
                                            create_report=True, 
                                            WGS_Genomes=WGS_Gnms).get('file_path')

    if cluster:
        clustered_fasta = contig_clutser(fasta_import_output) # clustered_fasta is the base name without the '.fa'
        cntgs_fasta = f'{os.path.splitext(fasta_import_output)[0]}_clstr.fa'
        SeqIO.write(SeqIO.parse(clustered_fasta, 'fasta'), cntgs_fasta, 'fasta') # creating new fasta file of the contigs. spcifically this one is formatted correctly (60 chars per line)
        cntg_set = {rec.name for rec in SeqIO.parse(cntgs_fasta, 'fasta')} # creating cntg_set of the clustered output. (set comprehension...)
    else:
        cntgs_fasta = fasta_import_output

    with open(f'{gnrl_anlz_file}.headers', 'w') as contigs_file: # contig headers for the follow up
        contigs_file.writelines(f'{contig}\n' for contig in cntg_set)

    if WGS_Gnms:
        tax_res = WGS_Genomes_tax_annotation(contigs_file.name) # annotating with the already implemented taxonomy inside the contigs name.
    else:
        tax_res = mmseq_taxonomy_make(cntgs_fasta) # using mmseq to annotate the contigs -> against RefSeq database.

    return tax_res, contigs_file.name.replace('.headers', '.report')

#####################################################################

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-P', '--gnrl_anlz_file', type=str, help='pattern filtration file')
    input_details = argparse.parse_args()
    gnrl_anlz_file = input_details.gnrl_anlz_file
    main(gnrl_anlz_file)