# modules import:
import os
import re
import argparse
import sys
import pandas as pd
pd.options.display.max_colwidth = 1000
from Bio import SeqIO
import sys
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/bio_utilities')
import phylogenetics
#

def main(general_analyzer, tax_res, report_file, lineage_dic=phylogenetics.taxid_lineage()):
    msdf = pd.read_table(tax_res, header=None) # creating a dataframe out of the mmseq result.tsv 
    msdf['det'] = msdf.apply(lambda x:f"[{'-'.join(str(x.iloc[i]) for i in list(msdf.columns)[1::])}]",axis=1)
    mmseq_dic = {contig:{'taxid':taxid, 'det':det} for contig, taxid, det in msdf.loc[:,[0,1,'det']].values}
    
    rf = pd.read_table(report_file)
    cntg_len_dic = {rec.name:len(rec.seq) for rec in SeqIO.parse(report_file.replace('.report','.fa'),'fasta')}
    gf = pd.read_table(general_analyzer)
   
    gf_lst = [{'contig':contig,
               'pattern':pattern,
               'mmseq_tax':mmseq_dic.get(contig,{}).get('det','no_mmseq_result'),
               'mmseq_lineage':lineage_dic.get(mmseq_dic.get(contig,{}).get('taxid'),'no_lineage').replace('; ','->'),
               'length':cntg_len_dic.get(contig,'not_found'),
               'ORFs':len(re.sub(r"['\[\]\s]","",orfs).split(','))} 
               for contig, pattern, orfs in gf.loc[:,['contig','pattern','ORFs']].values] # it is a list of dictionaries - each of them represent specific candidtae.

    opt = {col:list() for col in ['contig','pattern','mmseq_tax','contig_length',
                                   'ORFs','unmatched_ORFs','mmseq_lineage']}

    for contig_desc in gf_lst: # contig_desc is a dictionary.
        
        # contig name
        contig = contig_desc.get('contig')
        opt.get('contig').append(contig)

        # pattern
        opt.get('pattern').append(contig_desc.get('pattern'))
        
        # mmseq_tax 
        opt.get('mmseq_tax').append(contig_desc.get('mmseq_tax'))
        
        # contig_length column
        opt.get('contig_length').append(contig_desc.get('length'))

        # ORFs & unmatched_ORFs columns
        with open(re.sub(r'\.fa$', '.proteins.faa', rf[rf.loc[:,'contig']==contig].loc[:,'path'].values.tolist()[0]).replace('fna','faa')) as rec_file:
            sum_ORFs = len(list(filter(lambda x : x[0] == '>' and contig in x ,rec_file.readlines())))
        opt.get('ORFs').append(sum_ORFs)
        opt.get('unmatched_ORFs').append(sum_ORFs - contig_desc.get('ORFs'))

        # mmseq_lineage column
        opt.get('mmseq_lineage').append(contig_desc.get('mmseq_lineage'))
        #

    opt_df = pd.DataFrame.from_dict(opt) 
    opt_df[opt_df['mmseq_lineage']!='no_lineage'] # filteres contigs that have lineages that are not in the database.
    opt_df[opt_df['mmseq_tax'].apply(lambda x:x.split('-',1)[0].replace('[','') not in phylogenetics.uninformative_taxids().keys())] # filters contigs that annotated to be irelevant. like 'bacteria' or so...
    opt_df['tmp_col'] = opt_df['pattern'].apply(lambda x:len(x))
    opt_df.sort_values(by='tmp_col', ascending=False, inplace=True)
    opt_df.drop(columns='tmp_col', inplace=True) 
    opt_df.to_csv(f"{os.path.splitext(report_file)[0]}_review_summary.tsv", index=None, sep='\t') # '...10_min_core_prts.report' -> '...10_min_core_prts_review_summary.tsv'
    return f"{os.path.splitext(report_file)[0]}_review_summary.tsv"

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-G', '--general_analyzer', type=str, help='the general analyzer file that was created by operon_finder.py and was filtered in structure filtration.py')
    argparse.add_argument('-TS', '--tax_res', type=str, help='the mmseq result as tsv file')
    argparse.add_argument('-R', '--report_file', type=str, help='a file that contains the path of the contigs source')
    input_details = argparse.parse_args()

    general_analyzer = input_details.general_analyzer
    tax_res = input_details.tax_res
    report_file = input_details.report_file

    main(general_analyzer, tax_res, report_file)
