# modules import:
import os
import sys
import glob
import json
import pickle
import shutil
import random
import argparse
import datetime
import itertools
import pandas as pd
from Bio import SeqIO
from ete3 import NCBITaxa
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/bio_utilities')
sys.path.append(r'/davidb/yatirsolan/scripts/python/cluster_utilities')
from phylogenetics import taxid_to_name
import mmseq_taxonomy
import fasta_import
import phylogenetics
import job_runner
import databases

##########################################################################################################

def system_paths():
    return {'T3SS':{'T3SS':r'/davidb/yatirsolan/secretion_systems/T3SS/T3SS_vs_Metagenomics/review_summary/T3SS_summary.tsv'},   

            'T4SS':{'T4SSA':r'/davidb/yatirsolan/secretion_systems/T4SS/T4SSA_vs_Metagenomics/review_summary/T4SSA_summary.tsv',
                    'T4SSB':r'/davidb/yatirsolan/secretion_systems/T4SS/T4SSB_vs_Metagenomics/review_summary/T4SSB_summary.tsv'},

            'T6SS':{'T6SSi':r'/davidb/yatirsolan/secretion_systems/T6SS/T6SSi_vs_Metagenomics/review_summary/T6SSi_summary.tsv',
                    'T6SSii':r'/davidb/yatirsolan/secretion_systems/T6SS/T6SSii_vs_Metagenomics/review_summary/T6SSii_summary.tsv',
                    'T6SSiii':r'/davidb/yatirsolan/secretion_systems/T6SS/T6SSiii_vs_Metagenomics/review_summary/T6SSiii_summary.tsv'}}

##########################################################################################################
                
def phylogenetic_distribution(tax_info, rnk, ncbi, tax_map_dic=phylogenetics.taxid_mapping()):
    # the function mainley recieves a tax_info dictionary/or summary path with taxonomy annotation already implemented, and summarize within a dictionary the phylogenetic representation of every taxids in the desired rank (family/order).
    # e.g. - {3321378:[contig 1, contig 2, ...], 13345169:[contig 5, contig 6, contig 7,...]}
    # unmapped contigs are also delivered as an output - some of them due to mmseq taxonomy assay could have not match them to anything, some due to some unupdated taxids.
    # e.g. - {contig 10:'no_mmseq_result', contig 17:'no_mmseq_result', contig 20:{2321108: 'Arcobacteraceae'}} # here for instance, the taxid 2321108 was not found in my taxids databases.

    sngl_sys_distrbtion = dict()
    sngl_sys_unmapped = dict()

    if isinstance(tax_info, str): # means that and summary path was given, and the next line is parsing this file in to a dicitonary.
        tax_info = {cntg:tx.split('-',1)[0].replace('[','') for cntg,tx in pd.read_table(tax_info).loc[:,['contig', 'mmseq_tax']].values}
    elif isinstance(tax_info, pd.core.frame.DataFrame):
        # the other option is that tax_info is already a dataframe of the mmseq results file.
        tax_info = {cntg:str(tx) for cntg, tx in tax_info.iloc[:,[0,1]].values} # added the str(tx) recently, worth while to check if it is always needed.

    for cntg, tx in tax_info.items(): 
        res = tax_map_dic.get(tx, dict()).get(rnk)
        if res:
            if res not in sngl_sys_distrbtion:
                sngl_sys_distrbtion[res] = list()
            sngl_sys_distrbtion.get(res).append(cntg)
        else:
            if tx == 'no_mmseq_result':
                sngl_sys_unmapped.update({cntg:tx}) # {contig 10:'no_mmseq_result'}
            else:
                sngl_sys_unmapped.update({cntg:ncbi.get_taxid_translator([tx])}) # {2321108: 'Arcobacteraceae'}
                
    return sngl_sys_distrbtion, sngl_sys_unmapped

##########################################################################################################

def wgs_only(contigs):
    return all(databases.genomic_source_classifire(c) == 'WGS' for c in contigs)

##########################################################################################################

def bash_commands(query, taxid, target=mmseq_taxonomy.target_database(database='uniref_100')):
    query_base = os.path.basename(query).split('.f')[0]
    target_base = os.path.basename(target).split('._mm')[0]
    working_folder = os.path.dirname(query)
    mmseq_folder = os.path.join(working_folder, f'Mmseq_vs_{target_base}_{str(random.randint(0, 1000000))}')
    query_db = os.path.join(mmseq_folder, f'{query_base}.mmdb')
    tax_result = os.path.join(mmseq_folder, f'{query_base}_vs_{target_base}_TaxResult._mmdb')
    threads=64 # was 48
    sensitivity=4
    blacklist=phylogenetics.customize_blacklist()+[taxid]
    tmp = f'/davidb/scratch/mmseqs2/{datetime.datetime.today().now().isoformat()}'
    params = f'--threads {threads} -s {sensitivity}'

    if blacklist: # blacklist is an mmsesq feature that allows the exclusion of taxids from the taxonomy survey.
        # can be seen on terminal -> 'mmseqs taxonomy --help'
        blacklist = blacklist if isinstance(blacklist, list) else [blacklist] # the input of blacklist can be int/list
        blacklist += mmseq_taxonomy.default_blacklist() 
        blacklist = ','.join([str(el) for el in blacklist])
        blacklist = f'--blacklist {blacklist}'
    else:
        blacklist = ''

    if not os.path.exists(mmseq_folder):
        os.makedirs(mmseq_folder)

    bash_commands = [r"module load MMseqs2/dec_2020",
                     f"mmseqs createdb {query} {query_db}",
                     f"mmseqs taxonomy {query_db} {target} {tax_result} {tmp} {params} {blacklist}",
                     f"mmseqs createtsv {query_db} {tax_result} {tax_result.replace('._mmdb','.tsv')}",
                     f"mmseqs taxonomyreport {target} {tax_result} {tax_result.replace('._mmdb', '.report')}"]
    
    return ('\n'.join(bash_commands), mmseq_folder)
    
##########################################################################################################

def descendancy(rank_taxid, taxids, rnk, tax_map_dic, ncbi):
# this function gets a tuple which is a single dictionary element ({key:value}) - the key is a rank in taxid, and the values are list of taxids
# once one of the taxids within the taxids in the values, is equal or a desendent of the key taxid, the element will be concluded. 
# it means that the possession veryfication was complited regarding this specific taxid of the rank (family/order)
    return any((taxid == rank_taxid or
                taxid in [str(txd) for txd in ncbi.get_descendant_taxa(rank_taxid, collapse_subspecies=False)] or # 'get_descendant_taxa' returnes integers, therefor it is crucial to convert the output to strings.
                tax_map_dic.get(taxid,dict()).get(rnk) == rank_taxid) 
                for taxid in taxids)

##########################################################################################################

def main(mtdata, data_presentation_dir, rnk, system, specified_tree_name=str()):
    assay_dir = os.path.join(data_presentation_dir, rnk, system) # e.g. - /davidb/yatirsolan/data_presentation/family/T6SS
    if not os.path.isdir(assay_dir):
        os.mkdir(assay_dir)
    os.chdir(assay_dir)

    basename = os.path.join(assay_dir, f"{os.path.basename(mtdata).split('_phylogenetic_mtdta')[0]}_{system}") # e.g. - '/davidb/yatirsolan/review_tree/family/T6SS/review_family_T6SS'
    if 'mtdta_' in mtdata and not specified_tree_name:
        specified_tree_name = f"_{mtdata.split('mtdta_')[-1].replace('.tsv','')}"

    system_dic = system_paths().get(system) # here is where the data of the research pipeline is entering this script.
    sys_pos_df = pd.read_table(mtdata, dtype={'rnk_txd':str}) # This the table, of the families in the tree, it is not related to the research pipeline itself.
    ncbi = NCBITaxa()
    false_positives = list()
    systems_distribution = dict()
    gca_taxids_dic = databases.gca_taxids_dic()
    tax_map_dic = phylogenetics.taxid_mapping()
    main_fasta_file = f'{system}_wgs_only_all_cntgs.fa'
    main_fasta_dict = dict()
    wgs_only_mmseq_file = f"{basename.rsplit('_',1)[0]}_wgs_only_TaxResult.tsv"
    wgs_only_mmseq_df = pd.read_table(wgs_only_mmseq_file, header=None) if os.path.exists(wgs_only_mmseq_file) else None
    cmplmntry_mmseq_file = None

    if os.path.exists(wgs_only_mmseq_file) and not os.path.exists(main_fasta_file):
        main_fasta_file = fasta_import.main(hdrs=set(wgs_only_mmseq_df.iloc[:,0]), 
                                            output_type='fa', 
                                            output_file_name=f'{system}_wgs_only_all_cntgs', 
                                            create_report=False).get('file_path')

    for system_sbtype, smry_path in system_dic.items():
        print(system_sbtype)
        systems_distribution[system_sbtype], _ = phylogenetic_distribution(smry_path, rnk, ncbi)
        wgs_only_distribution = dict(filter(lambda elem: wgs_only(elem[1]), systems_distribution.get(system_sbtype).items()))
        wgs_only_cntgs = set(itertools.chain.from_iterable(wgs_only_distribution.values())) # leaving only the contigs themself - later the script uses this set to import their sequences for an mmseq taxonomy assay.

        if os.path.exists(main_fasta_file):
            main_fasta_dict = SeqIO.to_dict(SeqIO.parse(main_fasta_file, 'fasta'))
            cntgs_fasta_cmplmntry_mmseq = fasta_import.main(hdrs=filter(lambda cntg:cntg not in main_fasta_dict, wgs_only_cntgs), 
                                                            output_type='fa', 
                                                            output_file_name='tmp', 
                                                            create_report=False).get('file_path')
            if os.path.exists(cntgs_fasta_cmplmntry_mmseq):
                SeqIO.write(list(SeqIO.parse(cntgs_fasta_cmplmntry_mmseq, 'fasta')) + list(SeqIO.parse(main_fasta_file, 'fasta')), main_fasta_file, 'fasta')
                main_fasta_dict = SeqIO.to_dict(SeqIO.parse(main_fasta_file, 'fasta'))
                os.remove(cntgs_fasta_cmplmntry_mmseq)

        cntgs_for_cmplmntry_mmseq = wgs_only_cntgs

        if os.path.exists(wgs_only_mmseq_file):
            cntgs_for_cmplmntry_mmseq = set(filter(lambda cntg:cntg not in set(wgs_only_mmseq_df.iloc[:,0]), wgs_only_cntgs))

        if cntgs_for_cmplmntry_mmseq: # an empty set/list/dictionary is considered False in python...
            cntgs_fasta_cmplmntry_mmseq = fasta_import.main(hdrs=cntgs_for_cmplmntry_mmseq, 
                                                            output_type='fa', 
                                                            output_file_name=f"{basename.rsplit('_',1)[0]}_complementary_mmseq", 
                                                            create_report=False).get('file_path')
            
            main_fasta_dict = {**main_fasta_dict, 
                               **SeqIO.to_dict(SeqIO.parse(cntgs_fasta_cmplmntry_mmseq, 'fasta'))}
            
            SeqIO.write(main_fasta_dict.values(), main_fasta_file, 'fasta')
            cmplmntry_mmseq_file = mmseq_taxonomy.main(query=cntgs_fasta_cmplmntry_mmseq, 
                                                       blacklist=phylogenetics.customize_blacklist(), 
                                                       delete_folder=False)
            if os.path.exists(cntgs_fasta_cmplmntry_mmseq):
                os.remove(cmplmntry_mmseq_file)

        wgs_only_mmseq_df = pd.concat([wgs_only_mmseq_df]+[pd.read_table(cmplmntry_mmseq_file, header=None) if cmplmntry_mmseq_file else None], ignore_index=True)
        if cntgs_for_cmplmntry_mmseq:
            os.remove(cmplmntry_mmseq_file)
        mmseq_results_dic = dict(wgs_only_mmseq_df.iloc[:, [0, 1]].values)
        own_mpd_cntgs = list(filter(lambda cntg:gca_taxids_dic.get(databases.gca_from_accession(cntg)) == mmseq_results_dic.get(cntg), wgs_only_cntgs))

        valid_taxids, _ = phylogenetic_distribution(wgs_only_mmseq_df[~(wgs_only_mmseq_df.iloc[:,0].isin(own_mpd_cntgs))], 
                                                    rnk, 
                                                    ncbi, 
                                                    tax_map_dic)
        own_mpd_cntgs_not_to_chck = list()

        for own_mapped_cntg in own_mpd_cntgs:
            own_mapped_cntg_taxid = gca_taxids_dic.get(databases.gca_from_accession(own_mapped_cntg))
            rank_taxid = tax_map_dic.get(str(own_mapped_cntg_taxid)).get(rnk)
            if rank_taxid in valid_taxids:
                own_mpd_cntgs_not_to_chck.append(own_mapped_cntg)
        own_mpd_cntgs = list(filter(lambda x:x not in own_mpd_cntgs_not_to_chck, own_mpd_cntgs))
        wgs_only_mmseq_df = wgs_only_mmseq_df[~(wgs_only_mmseq_df.iloc[:,0].isin(own_mpd_cntgs))]

        if own_mpd_cntgs:
            own_mpd_cntgs_fasta_file = f'{system}_own_mapped.fa'
            SeqIO.write([main_fasta_dict.get(cntg) for cntg in own_mpd_cntgs], own_mpd_cntgs_fasta_file, 'fasta')
            GCAs = [databases.gca_from_accession(cntg) for cntg in own_mpd_cntgs]
            taxids = [gca_taxids_dic.get(gca) for gca in GCAs]
            target = mmseq_taxonomy.target_database(database='uniref_100')
            own_mpd_cntgs = mmseq_taxonomy.main(query=own_mpd_cntgs_fasta_file, 
                                                target=target, 
                                                blacklist=taxids, 
                                                delete_folder=False)
            
            own_mpd_mmseq_df = pd.read_table(own_mpd_cntgs, header=None)
            wgs_only_mmseq_df = pd.concat([wgs_only_mmseq_df, own_mpd_mmseq_df])
        wgs_only_mmseq_df.to_csv(wgs_only_mmseq_file, sep='\t', index=None, header=None)
        mmseq_results_dic = {cntg:str(txd) for cntg, txd in wgs_only_mmseq_df.iloc[:,[0,1]].values}

        false_positives = dict(filter(lambda item: not descendancy(rank_taxid=item[0], # taxid of the family/order
                                                                   taxids=[mmseq_results_dic.get(contig) for contig in item[1]], # list of taxids in any rank. hopefully most are species.
                                                                   rnk=rnk,
                                                                   tax_map_dic=tax_map_dic, # phylogenetics.taxid_mapping() 
                                                                   ncbi=ncbi), wgs_only_distribution.items()))

        systems_distribution[system_sbtype] = dict(filter(lambda item:item[0] not in false_positives, systems_distribution.get(system_sbtype).items()))
        sys_pos_df[system_sbtype] = sys_pos_df.rnk_txd.apply(lambda rnk_txd: 1 if rnk_txd in systems_distribution.get(system_sbtype) else 0) # adding system_sbtype possession annotation. 
        systems_distribution[system_sbtype] = {taxid_to_name(txid, ncbi):cntgs for txid, cntgs in systems_distribution.get(system_sbtype).items()} # changing the dictionary keys from taxids into textual taxas, for the pickle file later to be written.
        
    with open(f'{basename}_systems_distribution{specified_tree_name}.pkl', 'wb') as f:
        pickle.dump(systems_distribution, f)
    
    sys_pos_df.drop('rnk_txd', axis=1, inplace=True)
    sys_pos_df = pd.melt(sys_pos_df.loc[:,['rnk_txn']+list(system_dic.keys())], 
                         id_vars=['rnk_txn'], 
                         var_name='system', 
                         value_name='possession') # reducing df columns to only 'rnk_txn' and the systems possession, then melting the df.
    sys_pos_df['possession'] = sys_pos_df.apply(lambda x:x.loc['system'] if x.loc['possession'] == 1 else 'none', axis=1)
    sys_pos_df.to_csv(f'{basename}_system_possession_mtdta{specified_tree_name}.tsv', sep='\t', index=None)

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-S', '--system', type=str, choices=['T3SS', 'T4SS', 'T6SS'], required=True, help='system_name')
    input_details = argparse.parse_args()
    mtdata = r'/davidb/yatirsolan/review_tree/family/review_family_phylogenetic_mtdta.tsv'
    data_presentation_dir = r'/davidb/yatirsolan/data_presentation' # The main insertion directory.
    rnk = 'family'
    system = input_details.system 
    main(mtdata, data_presentation_dir, rnk, system)