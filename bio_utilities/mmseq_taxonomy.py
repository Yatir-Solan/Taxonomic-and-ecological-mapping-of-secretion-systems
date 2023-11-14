# modules import:
import os
import re
import shutil
import sys
import argparse
import datetime
import random
import sys
import pandas as pd 
from Bio import SeqIO
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/bio_utilities')
sys.path.append(r'/davidb/yatirsolan/scripts/python/cluster_utilities')
import job_runner
#

def target_database(database):
    return {'ref_seq':r'/davidb/yatirsolan/databases/mmseq_for_taxonomy/RefProteome_2021_03/all.RefProteome._mmdb',
            'uniref_100':r'/davidb/bio_db/UniProt/UniRef100/uniref100._mmdb'}.get(database)

def default_blacklist():
    # default values of mmseq, that are ignored in first place!
    return [12908, 28384] 

def main(query, 
         target=target_database(database='ref_seq'), 
         threads=48,
         sensitivity=4, 
         queue='duduheavy', # was changed to duduheavy instead of duduhimem due to overload in phoenix
         limit_ram='N', 
         blacklist=None, # blacklist should be list/set or int/strq
         delete_folder=False, 
         cpu=4, 
         wait=60):

    query_base = os.path.basename(query).split('.f')[0]
    target_base = os.path.basename(target).split('._mm')[0]
    working_folder = os.path.dirname(query)
    mmseq_folder = os.path.join(working_folder, f'Mmseq_vs_{target_base}_{str(random.randint(0, 1000000))}')
    query_db = os.path.join(mmseq_folder, f'{query_base}.mmdb')
    tax_result = os.path.join(mmseq_folder, f'{query_base}_vs_{target_base}_TaxResult._mmdb')

    # this directory is used for search - delete folder for rerunning search
    tmp = f'/davidb/scratch/mmseqs2/{datetime.datetime.today().now().isoformat()}'

    params = f'--threads {threads} -s {sensitivity}'
    if limit_ram != 'N':
        # best to use 250G if needed
        params = f'{params} --split-memory-limit {limit_ram}'

    if not os.path.exists(mmseq_folder):
        os.makedirs(mmseq_folder)
    
    if blacklist: # blacklist is an mmsesq feature that allows the exclusion of taxids from the taxonomy survey.
        # can be seen on terminal -> 'mmseqs taxonomy --help'
        blacklist = blacklist if isinstance(blacklist, list) else [blacklist] # the input of blacklist can be int/list
        blacklist += default_blacklist() 
        blacklist = ','.join([str(el) for el in blacklist])
        blacklist = f'--blacklist {blacklist}'
    else:
        blacklist = ''

#################################################################
# New part - 25.07.23

    mmseq_file = f'{query_base}_vs_{target_base}_TaxResult.tsv'
    difference = True
    prev_query = None

#################################################################
# New part - 25.07.23
    if os.path.isfile(mmseq_file):
        pass
        # old_mmseq_df = pd.read_table(mmseq_file, header=None)
        # mmseq_hdrs = set(old_mmseq_df.loc[:, 0])
        # fa_dict = SeqIO.to_dict(SeqIO.parse(query, 'fasta'))
        # fa_hdrs = set(SeqIO.to_dict(SeqIO.parse(query, 'fasta')).keys())
        # difference = fa_hdrs.difference(mmseq_hdrs)
        # if difference: 
        #     fraction_query = os.path.join(working_folder, 'fraction.fa')
        #     SeqIO.write([fa_dict.get(hdr) for hdr in difference], fraction_query, 'fasta')
        #     prev_query, query = query, fraction_query

#################################################################

    if difference:
        bash_commands = [r"module load MMseqs2/dec_2020",
                         f"mmseqs createdb {query} {query_db}",
                         f"mmseqs taxonomy {query_db} {target} {tax_result} {tmp} {params} {blacklist}",
                         f"mmseqs createtsv {query_db} {tax_result} {tax_result.replace('._mmdb','.tsv')}",
                         f"mmseqs taxonomyreport {target} {tax_result} {tax_result.replace('._mmdb', '.report')}"]
        
        if os.path.exists(query_db):
            bash_commands = list(filter(lambda line:'mmseqs createdb' not in line, bash_commands))

        tax_job_id = job_runner.job_submit(script='\n'.join(bash_commands), job_name='mmseq_taxonomy_assay', queue=queue, cpu=cpu)
        job_runner.wait_for_job_to_end(tax_job_id, wait)
        shutil.copy(os.path.join(mmseq_folder, mmseq_file), mmseq_file)

#################################################################
# New part - 25.07.23

        # if prev_query:
        #     new_mmseq_df = pd.read_table(mmseq_file, header=None)
        #     pd.concat([old_mmseq_df, new_mmseq_df]).to_csv(mmseq_file, header=None, index=None, sep='\t')

#################################################################

        os.remove(f'mmseq_taxonomy_assay.o{tax_job_id}')
    
    if delete_folder:
        shutil.rmtree(mmseq_folder) # deleting the folder

    return mmseq_file

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-QRY', '--query', required=True, type=str, help='a path to a query contig given as fasta file')
    argparse.add_argument('-TRG', '--target', default=target_database(database='ref_seq'), type=str, help='a path to a protein db created by mmseqs, default is RefProteome_2021_03, aditional option is uniref100 -> /davidb/bio_db/UniProt/UniRef100/uniref100._mmdb')
    argparse.add_argument('-TDS', '--threads', default=48, type=int, help='number of threads [default 48]')
    argparse.add_argument('-S', '--sensitivity', default=4, type=float, help='mmseq sensitivity. [default 4]')
    argparse.add_argument('-C', '--cpu', default=4, type=int, help='number of CPUs of the mmseq job that will be submitted.[default 4]')
    argparse.add_argument('-Q', '--queue', default="duduhimem", type=str, help='queue to run on. [default duduhimem with all RAM]')
    argparse.add_argument('-L', '--limit_ram', default='N', type=str, help='limit RAM usage, N ignores that flag. example: 10K, 100M, 100G [default N]')
    argparse.add_argument('-B', '--blacklist', default=None, type=str, help='comma separated list of ignored taxa in LCA computation')
    argparse.add_argument('-D', '--delete_folder', type=str, choices=['True','False'], default='False', help="True-deletes the folder of mmseq, False-Don't [default is False]")
    input_details = argparse.parse_args()
    #
    query = input_details.query
    target = input_details.target
    threads = input_details.threads
    sensitivity = input_details.sensitivity
    cpu = input_details.cpu
    queue = input_details.queue
    limit_ram = input_details.limit_ram
    blacklist = re.sub(r'[\[\]\s]', '', input_details.blacklist).split(',') if input_details.blacklist else input_details.blacklist # converting pythonic string that looks like list to be a real list : '[123,456,789]' -> ['123', '456', '789'] or '123' -> ['123']
    delete_folder = {'True':True,'False':False}.get(input_details.delete_folder)
    #
    main(query, target, threads, sensitivity, queue, limit_ram, blacklist, delete_folder, cpu)
