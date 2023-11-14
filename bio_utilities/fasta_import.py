# modules import:
import re
import os
import sys
import math
import time
import json
import glob
import random
import string
import argparse
import itertools
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/cluster_utilities')
import job_runner
import databases

#####################################################################

def arrange_as_lst(hdrs):
    if type(hdrs) is str:
        if os.path.isfile(hdrs):
            cntg_hdrs_lst = [contig.replace('\n', '') for contig in open(hdrs, 'r').readlines()] # it means that it is a file of headers, and in that case we should make a list out of it.
        else:
            cntg_hdrs_lst = re.sub(r'[\[\]\s]', '', hdrs).split(',')
    else:
        cntg_hdrs_lst = list(hdrs)
    return cntg_hdrs_lst

#####################################################################

def output_type_dict():
    return {'fa':'.fa', 
            'faa':'.proteins.faa', 
            'genes':'.genes.ffn'}

#####################################################################

def get_pathes():
    denovo_reg, mgnify_reg, wgs_reg = list(databases.rglr_expr().values())[:-1:]
    denovo_path = r'/davidb/assemblies/all_min10k/denovo.sample_filePrefix.tsv'
    df = pd.read_table(denovo_path, names=['sample', 'path'])
    df['sample'] = df['sample'].apply(lambda x:re.search(denovo_reg, x).group())
    mgnify_path = r'/davidb/assemblies/all_min10k/mgnify.sample_filePrefix.tsv'
    mf = pd.read_table(mgnify_path, names=['sample', 'path'])
    mf['sample'] = mf['sample'].apply(lambda x:re.search(mgnify_reg, x).group() if 'ERZ' in x else x)
    wgs_genomes_path = r'/davidb/assemblies/all_min10k/wgs_genomes.sample_filePrefix.tsv'
    wgf = pd.read_table(wgs_genomes_path, names=['sample', 'path'])
    wgf['sample'] = wgf['sample'].apply(lambda x:re.search(wgs_reg, x).group())
    wgs_uncultured_path = r'/davidb/assemblies/all_min10k/wgs_uncultured.sample_filePrefix.tsv'
    wuf = pd.read_table(wgs_uncultured_path, names=['sample', 'path'])
    wuf['sample'] = wuf['sample'].apply(lambda x:re.search(wgs_reg, x).group())
    wgs_metagenomes_path = r'/davidb/assemblies/all_min10k/wgs_metagenomes.sample_filePrefix.tsv'
    wmf = pd.read_table(wgs_metagenomes_path, names=['sample', 'path'])
    wmf['sample'] = wmf['sample'].apply(lambda x:re.search(wgs_reg, x).group())
    main_df = pd.concat([df, mf, wgf, wuf, wmf])
    return dict(main_df.values)

#####################################################################

def fasta_source(contig, path_dic, output_type):
    def gem_source(contig, output_type):
        gem_output_type = {'.fa':'fna', 
                           '.proteins.faa':'faa', 
                           '.genes.ffn':'ffn'}.get(output_type_dict().get(output_type))
        
        for file in glob.glob(os.path.join(r'/davidb/bio_db/GEM', gem_output_type, f"{contig.split('_')[0]}*")): # let say we have a contig "2004178001_GEM_gws2_d1_ctg_0448" -> searching for directories starting with "2004178001..."
            with open(file, 'r') as fsta_srce:
                for line in fsta_srce:
                    if line.startswith('>'):
                        if contig in line.split()[0][1::]:
                            return fsta_srce.name
                        
    def all_other_source(contig, path_dic, output_type, source):
        match = re.search(databases.rglr_expr().get(source), contig).group()
        fst_source = f'{path_dic.get(match)}{output_type_dict().get(output_type)}'
        if os.path.isfile(fst_source):
            return fst_source
        elif output_type == 'fa': # this condition deals with an option (that was encountered) in which the 'fa' file is abscent. if by chance the gff file is present that a solution.
            return f'{path_dic.get(match)}.gff' # the function get_rcrds can do its job on gff files as well.
        return 'file_not_found'

    source = databases.genomic_source_classifire(contig)
    if source is 'GEM':
        return gem_source(contig, output_type)
    else:
        return all_other_source(contig, path_dic, output_type, source)

#####################################################################

def get_rcrds(contig, fasta_to_extract_seq_from):
    hits = list()
    mtgnmic_srce = databases.genomic_source_classifire(contig)
    cntg_rgx = {'DeNovo':re.compile(r'(?<=\.)[SE]RR[0-9]+_ctg_\d+$'),
                'Mgnify':re.compile(r'(?<=\|)ERZ\d+_ctg\d+$'),
                'WGS':re.compile(r'(?<=\|)GCA_.+$'),
                'GEM':contig}.get(mtgnmic_srce)
    if cntg_rgx is contig: # it means that it is a GEM contig
        line_rgx = contig
    else:
        line_rgx = cntg_rgx.search(contig).group().replace('|', '\|').replace('.', '\.')
        line_rgx = f'.+{line_rgx}(_+|$)'
    for req in SeqIO.parse(fasta_to_extract_seq_from, 'fasta'):
        if re.search(line_rgx, req.id):
            hits.append(req)
    return hits

#####################################################################

def get_interval_and_method(cntgs_num):
    # if rec_num is  28 and lower -> run by multiprocessing (every process will run one record) - small tasks...
    # higher -> spread into jobs -> we will use 1000 jobs. above we shall encounter 'to many open files' error.
    method = 'multiprocessing'
    interval = 1
    if cntgs_num > 28:
        method = 'jobs'
        desired_jobs_num = 500
        if cntgs_num > desired_jobs_num:
            interval = math.ceil(cntgs_num/desired_jobs_num) # the ronuded number prevents the to many open files error.
    return interval, method

#####################################################################

def search_and_import(contigs, fasta_dic, report_dic, method, path_dic, path_json, output_type):
    rcrds = list()
    report = True if type(report_dic) is mp.managers.DictProxy else False
    if method is 'multiprocessing':
        for contig in contigs:
            fasta_to_extract_seq_from = fasta_source(contig, path_dic, output_type)
            cntg_rcrds = get_rcrds(contig, fasta_to_extract_seq_from)
            if report:
                report_dic[contig] = fasta_to_extract_seq_from
            for req in cntg_rcrds:
                rcrds.append(req)
        fasta_dic[mp.current_process().pid] = rcrds # dictionary : key:process_id -> value:list of biopython seq records
    
    if method is 'jobs':
        lst_as_str_contigs = f"[{','.join(contigs)}]"
        output_file_name = ''.join(random.choices(string.ascii_letters + string.digits, k=20))

        job_id = job_runner.job_submit(script = f"python /davidb/yatirsolan/scripts/python/bio_utilities/fasta_import.py -F '{lst_as_str_contigs}' -T {output_type} -O {output_file_name} -J {path_json} -R {report} --contigs_subgroup True", 
                                       cpu = 1,
                                       job_name = 'fasta_import_job')

        fasta_dic[job_id] = f'{output_file_name}{output_type_dict().get(output_type)}'
        if report:
            report_dic[job_id] = f'{output_file_name}.report'

#####################################################################

def subgroup_contig_search(hdrs, path_json, output_type, output_file_name, create_report):
    cntg_hdrs_lst = arrange_as_lst(hdrs) 
    path_dic = json.load(open(path_json)) 
    rcrds = list()
    subgroup_report_dic = dict()
    for contig in cntg_hdrs_lst:
        fasta_to_extract_seq_from = fasta_source(contig, path_dic, output_type)
        subgroup_report_dic[contig] = fasta_to_extract_seq_from
        cntg_rcrds = get_rcrds(contig, fasta_to_extract_seq_from)
        for req in cntg_rcrds:
            rcrds.append(req)
        
    SeqIO.write(rcrds, f'{output_file_name}{output_type_dict().get(output_type)}', 'fasta')
    if create_report:
        with open(f'{output_file_name}.report', 'w') as contig_report:
            json.dump(subgroup_report_dic, contig_report, indent=0)

#####################################################################

def main_paralel(hdrs, output_type, output_file_name=None, create_report=True):
    cntg_hdrs_lst = arrange_as_lst(hdrs)  
    cntgs_num = len(cntg_hdrs_lst)
    rec_interval, method = get_interval_and_method(cntgs_num)
    path_dic = get_pathes()
    if method is 'jobs':
        with open('path.json', 'w') as json_file:
            json.dump(path_dic, json_file, indent=1)
    processes = list()
    cntgs_lst = list()
    with mp.Manager() as man:
        fasta_dic = man.dict() # dictionary the it is shared through all the process memory
        report_dic = man.dict() if create_report else None # dictionary the it is shared through all the process memory
        for count, contig in enumerate(cntg_hdrs_lst):
            cntgs_lst.append(contig)
            if len(cntgs_lst) == rec_interval or count+1 == cntgs_num:
                search_and_import_proc = mp.Process(target=search_and_import,
                                                    kwargs={'contigs':cntgs_lst, 
                                                            'fasta_dic':fasta_dic, 
                                                            'report_dic':report_dic,
                                                            'method':method, 
                                                            'path_dic':path_dic, 
                                                            'path_json':'path.json',
                                                            'output_type':output_type})

                search_and_import_proc.start()
                processes.append(search_and_import_proc)
                cntgs_lst.clear() # cleaning the list for the next round

        for pr in processes:
            pr.join()
        
        if method is 'jobs':  
            job_runner.wait_for_job_to_end(list(fasta_dic.keys()), 30)
            file_to_remove = [f'fasta_import_job.o{job_id}' for job_id in fasta_dic.keys()] # list with all the output files of the jobs.
            file_to_remove.append(json_file.name) # the json file should also be erased.

            for file_name in file_to_remove:
                os.remove(file_name) # deleting all the jobs finish files...

        output_rcrds = list()
        for rcrds_lst in fasta_dic.values():
            if type(rcrds_lst) is str:
                fasta_to_remove = rcrds_lst
                rcrds_lst = SeqIO.parse(rcrds_lst, 'fasta')
                os.remove(fasta_to_remove)
            for req in rcrds_lst:
                output_rcrds.append(req)
        
        if type(hdrs) is str and os.path.isfile(hdrs):
            base_name = hdrs.replace('.headers', '')
        else:
            base_name = output_file_name if output_file_name else f'fasta_import_{str(random.randint(0, 1000))}'
        SeqIO.write(output_rcrds, f'{base_name}{output_type_dict().get(output_type)}', 'fasta')

        if create_report:
            with open(f'{base_name}.report', 'w') as report_main_file:
                report_main_file.write('contig\tpath\n')
                if method is 'multiprocessing':
                    for cntg, path in report_dic.items():
                        report_main_file.write(f'{cntg}\t{path}\n')

                elif method is 'jobs':
                    for rprt_json in report_dic.values():
                        for cntg, path in json.load(open(rprt_json)).items():
                            report_main_file.write(f'{cntg}\t{path}\n')
                        os.remove(rprt_json) 
        return {'file_path':os.path.realpath(f'{base_name}{output_type_dict().get(output_type)}'), 'biopython_records':output_rcrds}

#####################################################################

def main(hdrs, output_type, output_file_name=None, create_report=True, WGS_Genomes=False):
    contigs = arrange_as_lst(hdrs)  
    path_dic = get_pathes()
    report_dic = dict()

    if type(hdrs) is str and os.path.isfile(hdrs):
        base_name = hdrs.replace('.headers', '')
    else:
        base_name = output_file_name if output_file_name else f'fasta_import_{str(random.randint(0, 1000))}'

    report_dic = {contig:fasta_source(contig, path_dic, output_type) for contig in contigs}

    if WGS_Genomes:
        output_rcrds = list()
    else:
        output_rcrds = [get_rcrds(contig, fasta_to_extract_seq_from) for contig, fasta_to_extract_seq_from in report_dic.items()]
        output_rcrds = list(itertools.chain.from_iterable(output_rcrds))                     

    SeqIO.write(output_rcrds, f'{base_name}{output_type_dict().get(output_type)}', 'fasta')

    if create_report:
        pd.DataFrame.from_dict(report_dic, orient='index').reset_index().rename(columns={'index':'contig', 0:'path'}).to_csv(f'{base_name}.report', sep='\t', index=False)

    return {'file_path':os.path.realpath(f'{base_name}{output_type_dict().get(output_type)}'), 'biopython_records':output_rcrds}

#####################################################################

def pullseq(headers, target, output=None):
    def hdrs_arngmt(hdrs):
        if type(hdrs) is str:
            if os.path.isfile(hdrs):
                headers_lst = [hd.replace('\n', '') for hd in open(hdrs, 'r').readlines()]
            else:
                headers_lst = list(hdrs)
        else:
            headers_lst = hdrs  
        headers_lst = [hd.split()[0].replace('>', '') for hd in headers_lst]
        return headers_lst
    def get_seq(lst, trgt_seq_dic):
        seq_lst = list()
        for hd in lst:
            if hd in trgt_seq_dic.keys():
                seq_lst.append(trgt_seq_dic[hd])
        return seq_lst
    headers = hdrs_arngmt(headers)
    if type(target) is not dict:
        target = SeqIO.to_dict(SeqIO.parse(target, 'fasta'))
    seq_lst = get_seq(headers, target)
    if output:
        SeqIO.write(seq_lst, output, 'fasta') # if output is True a fasta file wil be produced.
    return seq_lst # list of biopython seq_records.
        
#####################################################################

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-F', '--hdrs', type=str, required=True, help='hdrs')
    argparse.add_argument('-T', '--output_type', type=str, choices=['fa','faa','genes'], required=True, help='fa-nucleotides, faa-proteins, genes-genes, gff-gff')
    argparse.add_argument('-O', '--output_file_name', type=str, help='name of output file in case a list or string is given')
    argparse.add_argument('-R', '--create_report', type=str, choices=['True','False'], default='True', help="True-create report file with paths, False-Don't [default is True]")
    argparse.add_argument('-J', '--json_path', type=str, help="this feature is oriented only for the specific script")
    argparse.add_argument('-CSG', '--contigs_subgroup', type=str, choices=['True','False'], default='False', help="this feature is oriented only for the specific script")
    
    input_details = argparse.parse_args() 
    hdrs = input_details.hdrs
    output_type = input_details.output_type
    output_file_name = input_details.output_file_name
    create_report = {'True':True, 'False':False}.get(input_details.create_report)
    json_path = input_details.json_path
    contigs_subgroup = {'True':True, 'False':False}.get(input_details.contigs_subgroup)

    if contigs_subgroup:
        subgroup_contig_search(hdrs, json_path, output_type, output_file_name, create_report)
    else:
        main(hdrs, output_type, output_file_name, create_report)