# modules import:
import os
import argparse
import sys
import re
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/secretion_systems_mapping_pipeline')
import operon_finder
import structure_filtration
import tax_annotation
import review_summary
import summarizing_files
#

def main(hmmsearch, system):

    ptrns, gnrl_anlyz = operon_finder.main(
        HmmSearchFile = hmmsearch,
        Evalue = 1e-6,
        max_space_betw_ORFs = 15,
        min_operon_length = 3,
        min_unq_prts_in_operon = 2,
        gff_file = None)

    strctr_fltrd_ptrns = structure_filtration.main(
        operon_patterns_file = ptrns,
        general_analyzer_file = gnrl_anlyz,
        system = system)
    
    if strctr_fltrd_ptrns == 'no_systems_passed_filtration':
        return

    tax_res, cntgs_rprt = tax_annotation.main(
        gnrl_anlz_file = gnrl_anlyz,
        cluster = False)

    final_file = review_summary.main(
        general_analyzer = gnrl_anlyz,
        tax_res = tax_res,
        report_file = cntgs_rprt)
        
if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-S', '--system', type=str, choices=['T3SS', 'T4SSA', 'T4SSB', 'T6SSi', 'T6SSii', 'T6SSiii'], required=True, help='system_name')
    input_details = argparse.parse_args()
    system = input_details.system 
    system_dir = system[:4:]
    general_path = f'/davidb/yatirsolan/secretion_systems/{system_dir}/{system}_vs_Metagenomics/'
    assay = 'review'
    paths = [os.path.join(general_path, src, assay) for src in ['Mgnify', 
                                                                'GEM', 
                                                                'WGS_Metagenomes', 
                                                                'WGS_Uncultured', 
                                                                'WGS_Genomes']]
    
    for db in paths:
        os.chdir(db)
        for f in os.listdir(db):
            r = re.search(r'.+\.[Hh]mm[Ss]earch$', f)
            if r:
                hmmsearch = f
        main(hmmsearch, system)

    summarizing_files.main(general_path, assay, system)