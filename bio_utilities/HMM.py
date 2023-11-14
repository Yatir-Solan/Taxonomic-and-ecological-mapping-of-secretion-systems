import argparse
import re
import os
import shutil
import multiprocessing as mp
import sys
import pandas as pd
# python scripts import:
sys.path.append('/davidb/yatirsolan/scripts/python/cluster_utilities')
import job_runner
import databases

class Hmmsearch_line:
    def __init__(self,line):
        self.Line = line.split()
        self.targetName = self.Line[0]
        self.targetAccession = self.Line[1]
        self.queryName = self.Line[2]
        self.queryAccession = self.Line[3]
        self.Evalue = float(self.Line[4])
        self.score = float(self.Line[5])
        self.bias = float(self.Line[6])
        self.domainEvalue = float(self.Line[7])
        self.domainScore = float(self.Line[8])
        self.domainBias = float(self.Line[9])
        self.exp = float(self.Line[10])
        self.reg = float(self.Line[11])
        self.clu = float(self.Line[12])
        self.ov = float(self.Line[13])
        self.env = float(self.Line[14]) 
        self.dom = float(self.Line[15])
        self.rep = float(self.Line[16])
        self.inc = float(self.Line[17])
        self.targetDescription = ' '.join(self.Line[18::])
    
class Hmm:
    def __init__(self, hmm_file):
        self.main_file_name = hmm_file
        self.parse_dic = self.parse() 

    def parse(self):
        hmms_dic = dict()
        hmm_lines = list()
        with open(self.main_file_name) as f: 
            for line in f: 
                hmm_lines.append(line) 
                if line[:4:] == 'NAME':
                    hmm_name = line.split()[1].replace('\n','') 
                if line[:3:] == 'ACC':
                    acc = line.split()[1].replace('\n','')
                    hmm_name = f'{hmm_name}_{acc}'
                if '//' in line: 
                    hmms_dic[hmm_name] = hmm_lines
                    hmm_lines = list()
        return hmms_dic

def hmmsearch_path():
    return r'/davidb/local/bin/hmmsearch'

def hmmbuild_path():
    return r'/davidb/local/bin/hmmbuild'
    
def hmmsearch(hmm_file, target=databases.genomic_datasets(), evalue=1e-4):
    def prll_hmmsrch(hmm_objct, trgt_name, trgt_path,evalue):
        # target_name = re.sub(r'all_|[\._]proteins\.faa|\.min10k.*\.faa','',os.path.basename(trgt_path))
        hmm_srch_lst = list() # using this list to create pd df to concate all of them.
        job_dic = dict()
        work_dir = f'{qry_hmm_name}_vs_{trgt_name}_hmmsearch'
        os.mkdir(work_dir) # the directory in which all the the search will be inside.
        for hmm_name in hmm_objct.parse_dic:
            with open(os.path.join(work_dir,f'{hmm_name}.Hmm'),'w') as sub_hmm_file:
                sub_hmm_file.writelines(hmm_objct.parse_dic.get(hmm_name))
            res_name = os.path.join(work_dir,f'{hmm_name}_vs_{trgt_name}.hmmsearch') 
            job_id = job_runner.job_submit(f'{hmmsearch_path()} -o /dev/null --tblout {res_name} --notextw -E {evalue} --cpu 4 {sub_hmm_file.name} {trgt_path}',job_name=hmm_name,queue='duduheavy')
            job_dic[job_id] = hmm_name
            hmm_srch_lst.append(res_name)
        job_runner.wait_for_job_to_end(list(job_dic.keys()),180)
        job_runner.report_file_delete(job_dic)
        
        with open(f'{qry_hmm_name}_vs_{trgt_name}.hmmsearch','a') as res_file:
            for hmm_srch_file in hmm_srch_lst:
                with open(hmm_srch_file,'r') as f:
                    for line in f:
                        res_file.write(line)
        shutil.rmtree(work_dir)

    hmm_objct = Hmm(hmm_file)
    qry_hmm_name = re.sub(r'\.[Hh][M|m]{2}','',os.path.basename(hmm_objct.main_file_name))
    processes = list()
    if isinstance(target, str):
        target = {re.sub(r'all_|[\._]proteins\.faa|\.min10k.*\.faa','',os.path.basename(target)):target}
    for trgt_name,trgt_path in target.items():
        search = mp.Process(target=prll_hmmsrch,kwargs={'hmm_objct':hmm_objct,'trgt_name':trgt_name,'trgt_path':trgt_path,'evalue':evalue})
        search.start()
        processes.append(search)
    for pr in processes:
        pr.join()
    return {key:f'{qry_hmm_name}_vs_{key}.hmmsearch' for key in target.keys()}    

def hmmsearch_to_tsv(hmmsearch_file, del_orgl_file=False):
    dic = {key:list() for key in ['targetName','targetAccession','queryName','queryAccession','E-value', 
                                  'score','bias','domainE-value','domainScore','domainBias', 
                                  'exp','reg','clu','ov','env','dom','rep','inc','targetDescription']}

    for line in list(filter(lambda line:line[0]!='#', open(hmmsearch_file,'r').readlines())):
        Line = Hmmsearch_line(line)
        dic.get('targetName').append(Line.targetName)
        dic.get('targetAccession').append(Line.targetAccession)
        dic.get('queryName').append(Line.queryName)
        dic.get('queryAccession').append(Line.queryAccession)
        dic.get('E-value').append(Line.Evalue)
        dic.get('score').append(Line.score)
        dic.get('bias').append(Line.bias)
        dic.get('domainE-value').append(Line.domainEvalue)
        dic.get('domainScore').append(Line.domainScore)
        dic.get('domainBias').append(Line.domainBias)
        dic.get('exp').append(Line.exp)
        dic.get('reg').append(Line.reg)
        dic.get('clu').append(Line.clu)
        dic.get('ov').append(Line.ov)
        dic.get('env').append(Line.env)
        dic.get('dom').append(Line.dom)
        dic.get('rep').append(Line.rep)
        dic.get('inc').append(Line.inc)
        dic.get('targetDescription').append(Line.targetDescription)
    pd.DataFrame.from_dict(dic).to_csv(f'{hmmsearch_file}.tsv', index=None, sep='\t')
    if del_orgl_file:
        os.remove(hmmsearch_file)
    return f'{hmmsearch_file}.tsv'

def hmm_naming(hmm_file, names_dic, replace_file=False):
    # the function recives hmm and a names dictionary (key=unwanted_name : value=wanted_name). the dictionary supposed to look like - /davidb/yatirsolan/secretion_systems/T3SS/knowledge/T3SS_naming.json 
    # the function create another hmm file in which the name section was changed to an informative : 'NAME  K012345.2' -> 'NAME  TssA'
    l = list() 
    r = re.compile(r'(?<=^NAME\s{2}).+(?=\n$)') 
    for line in open(hmm_file).readlines(): 
        if  line.startswith('NAME'): 
            unwanted_name = r.search(line).group() 
            line = line.replace(unwanted_name, names_dic.get(unwanted_name, unwanted_name)) 
            # line = line.replace(unwanted_name,names_dic.get(unwanted_name.split('.')[0],unwanted_name)) 
        l.append(line) 
    new_file_name = hmm_file if replace_file else f'{hmm_file}_named'
    with open(new_file_name,'w') as f: 
        f.writelines(l)

def hmm_fetch(hmm_lst, output_file_name, hmm_rep=databases.kegg_hmms_database()):
    kegg = 'Kegg' in hmm_rep # True/Flase
    hmm_rep_dic = dict()
    hmm_rep_lines_lst = open(hmm_rep).readlines()
    r_name = re.compile(r'(?<=^NAME\s{2}).+(?=\n$)')
    r_ac = re.compile(r'(?<=^ACC\s{3}).+(?=\n$)')
    for i, line in enumerate(hmm_rep_lines_lst):
        if line[:4:] == 'NAME':
            name = r_name.search(hmm_rep_lines_lst[i]).group() 
            acc = r_ac.search(hmm_rep_lines_lst[i+1]).group() 
            hmm_rep_dic[acc] = {'name':name,'start':i-1}
        if '//' in line : 
            hmm_rep_dic.get(acc)['end'] = i

    acc_rep_dic = {hmm_acc.split('.')[0] if kegg else hmm_acc : None for hmm_acc in hmm_rep_dic.keys()}

    for hmm_acc in hmm_rep_dic.keys():
        base_hmm_acc = hmm_acc.split('.')[0] if kegg else hmm_acc
        if acc_rep_dic.get(base_hmm_acc):
            acc_rep_dic.get(base_hmm_acc).append(hmm_acc)
        else:
            acc_rep_dic[base_hmm_acc] = [hmm_acc]
    with open(output_file_name,'w') as fetched:
        for hmm in sorted(hmm_lst):
            key = lambda x:int(x.split('.')[1]) if kegg else None
            for acc in sorted(acc_rep_dic.get(hmm), key=key):
                coor = hmm_rep_dic.get(acc)
                strt = coor.get('start')
                end = coor.get('end')
                fetched.writelines(hmm_rep_lines_lst[strt:end+1:])

def add_accessions(hmm_file, prefix=str()):
    # the function gets HMM file and adds ACC line (with the same name) below the NAME line.
    # if the accession line is present, no alteration will occure to the specific HMM.
    # finally the a new file overides the older one.
    changes = 0
    hmm_file_lines_lst = open(hmm_file).readlines() 
    r_name = re.compile(r'(?<=^NAME\s{2}).+(?=\n$)') 
    accessions = list()
    for i,line in enumerate(hmm_file_lines_lst): 
        if line[:4:] == 'NAME':
            if hmm_file_lines_lst[i+1][:3:] == 'ACC':
                continue
            changes += 1
            name = r_name.search(hmm_file_lines_lst[i]).group() 
            accession = f'ACC   {prefix}{name}\n' 
            accessions.append((accession,i))

    for count,(accession,idx) in enumerate(accessions): 
        hmm_file_lines_lst.insert(idx+count+1,accession)

    print(changes)

    if changes: # it means that there was accession additions into the list. otherwise, the file will not be altered.
        with open(hmm_file, 'w') as f:
            f.writelines(hmm_file_lines_lst)

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-H', '--hmm_file', type=str, help='hmm file to search with (the query)')
    argparse.add_argument('-T', '--target', type=str, help='specific proteins database to search in, the default is all the metagenomic data of the lab')
    input_details = argparse.parse_args()
    hmm_file = input_details.hmm_file
    target = input_details.target
    if target:
        hmmsearch(hmm_file, target)
    else:
        hmmsearch(hmm_file)