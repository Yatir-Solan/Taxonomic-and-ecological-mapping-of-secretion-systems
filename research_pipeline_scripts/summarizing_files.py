import os
import argparse
import pandas as pd

def main(general_path, assay, system):
    # This function units all the research results, from all the databses into a single file.
    os.chdir(general_path)
    dir_name = f'{assay}_summary'
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    for ptrn in ['_summary.tsv', '.report', '_gnrl_anlyz.tsv', '_min_core_prts.tsv']:
        df_lst = list()
        for data in ['GEM', 'Mgnify', 'WGS_Uncultured', 'WGS_Metagenomes', 'WGS_Genomes']:
            os.chdir(os.path.join(general_path, data, assay))
            for file in os.listdir():
                if file.endswith(ptrn):
                    df = pd.read_table(os.path.realpath(file))
                    if ptrn == '_min_core_prts.tsv':
                        df = df[~(df.patterns.str.startswith('parameters:'))] # rejects the last raw in the patterns_files.
                    df_lst.append(df)
        os.chdir(general_path)
        df = pd.concat(df_lst, ignore_index=True)
        if ptrn == '_summary.tsv':
            df['tmp_col'] = df['pattern'].apply(lambda x:len(x))
            df.sort_values(by='tmp_col', ascending=False, inplace=True)
            df.drop(columns='tmp_col', inplace=True) 
        df.to_csv(os.path.join(dir_name, f'{system}{ptrn}'), sep='\t', index=None)

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-P', '--general_path', type=str, help='general_path')
    argparse.add_argument('-A', '--assay', type=str, help='mostly review')
    argparse.add_argument('-S', '--system', type=str, choices=['T3SS', 'T4SSA', 'T4SSB', 'T6SSi', 'T6SSii', 'T6SSiii'], required=True, help='system_name')
    input_details = argparse.parse_args()
    general_path = input_details.general_path
    assay = input_details.assay
    system = input_details.system
    main(general_path, assay, system)