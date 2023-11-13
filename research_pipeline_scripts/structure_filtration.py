# modules import:
import os 
import re
import sys
import json
import argparse
import itertools
import collections
import pandas as pd
# python scripts import:
sys.path.append(r'/davidb/yatirsolan/scripts/python/bio_utilities')
sys.path.append(r'/davidb/yatirsolan/thesis_work/figures/general')
import databases
#

def structure_default_values(system):
    structure_dic = \
            {'T3SS':{'core_proteins': ['sctC', 'sctJ', 'sctN', 'sctQ', 'sctR', 'sctS', 'sctT', 'sctU', 'sctV'], # all the searched proteins can be seen here -> set(json.load(open('/davidb/yatirsolan/secretion_systems/T3SS/knowledge/T3SS_naming.json')).values())
                     'minimal_core_proteins':7,
                     'minimal_general_matches':7,
                     'exclusion_features':{'proteins':['FlgB', 'FlgC', 'FliE'],'fragment':any}, # flagellar specific
                     'inclusion_features':{'proteins':['sctC'], 'fragment':all}}, # T3SS specific

             'T4SSA':{'core_proteins':list(),
                     'minimal_core_proteins':0,
                     'minimal_general_matches':10,
                     'exclusion_features':{'proteins':json.load(open(r'/davidb/yatirsolan/secretion_systems/T4SS/knowledge/T4SS_rocha_names.json')).get('conjugational'), 'fragment':10},
                     'inclusion_features':{'proteins':['virb4'], 'fragment':all}},

             'T4SSB':{'core_proteins':['DotC', 'DotD', 'IcmG_DotF', 'IcmE_DotG', 'IcmK_DotH'], # based on -> https://www.pnas.org/doi/full/10.1073/pnas.1404506111
                     'minimal_core_proteins':4, # based on -> https://www.pnas.org/doi/full/10.1073/pnas.1404506111.
                     'minimal_general_matches':7, # the system is much larger, but in Legionella for instance there is one operon that is exactley 7 proteins.
                     'exclusion_features':None,
                     'inclusion_features':{'proteins':['IcmS', 'IcmW', 'IcmQ', 'IcmF', 'IcmH_DotU'], 'fragment':any}}, # T4SSB specific - proteins that are abscent of homologous both in IncI r64 conjugational plasmid and T4SSA. based on -> https://pubmed.ncbi.nlm.nih.gov/15652976, Following a meeting in 20.07.23 I omitted DotK from the inclusion feature.
            
             'T6SSi':{'core_proteins':json.load(open(r'/davidb/yatirsolan/secretion_systems/T6SS/hmm_repository/rocha/T6SS_naming.json')).get('T6SS_i'), 
                     'minimal_core_proteins':10,
                     'minimal_general_matches':10,
                     'exclusion_features':None,
                     'inclusion_features':None},

             'T6SSii':{'core_proteins':json.load(open(r'/davidb/yatirsolan/secretion_systems/T6SS/hmm_repository/rocha/T6SS_naming.json')).get('T6SS_ii'), 
                        'minimal_core_proteins':10,
                        'minimal_general_matches':10,
                        'exclusion_features':None,
                        'inclusion_features':None},

             'T6SSiii':{'core_proteins':json.load(open(r'/davidb/yatirsolan/secretion_systems/T6SS/hmm_repository/rocha/T6SS_naming.json')).get('T6SS_iii'), 
                        'minimal_core_proteins':10,
                        'minimal_general_matches':10,
                        'exclusion_features':None,
                        'inclusion_features':None}
            }
    #  list(itertools.chain.from_iterable()), # itertools.chain.from_iterable just flats 2d list in 1d list to put all the hmm names of all subtypes together.
    return structure_dic.get(system)

### the function gets list and check out if thers a minimal number of the desired core proteins.
###
def filter_by_core(pattern, core_proteins, minimal_core_proteins):
    return sum([core_protein in pattern for core_protein in core_proteins]) >= minimal_core_proteins # ensure appearence of core protein in the pattern.

### the function gets list of exclusive proteins and ensures their abscence (return False in case of the oposite [of any of them]).
###
def filter_by_exclusion(pattern, exclusion_features):
    if not exclusion_features:
        return True
    if callable(exclusion_features.get('fragment')): # checks if the fragment value is a function (callable...)
        func = exclusion_features.get('fragment') # func should be rather any or all
        return not func(exclusion_feature in pattern for exclusion_feature in exclusion_features.get('proteins')) 
    else:
        valid_exc_freq_dic = dict(collections.Counter([False if p in exclusion_features.get('proteins') else True for p in filter(lambda prot:prot!='#',pattern)])) # creates a dictionery with the frequencies of all the excluded and not-excluded proteins as False and True appearances.
        return (valid_exc_freq_dic.get(False,0)/sum(valid_exc_freq_dic.values()))*100 <= exclusion_features.get('fragment') # now checks that False appearances (excluded proteins) fragment in respect to the total system is lower than the value given.

### the function gets list of inclusive proteins and ensures their presence (return False in case of the oposite).
###
def filter_by_inclusion(pattern,inclusion_features):
    if not inclusion_features:
        return True
    if callable(inclusion_features.get('fragment')): # checks if the fragment value is a function (callable...)
        func = inclusion_features.get('fragment') # func should be rather any or all
        return func(inclusion_feature in pattern for inclusion_feature in inclusion_features.get('proteins')) #  ensure presence of all inclusion features (protein) in the pattern. generator is used to avoid entire pattern check after the first inclusion feature is abscent.

### the fuction gets a list and returns a dictionary that for all the core proteins it says if it is in the list
###
def pattern_info(pattern, core_proteins_lst):
    pattern_set = set(filter(lambda x:x!='#', pattern))
    return {core_protein:core_protein in pattern_set for core_protein in core_proteins_lst}
     
###
def mixed_filter(pattern,systems_prots_dic,uniq_sys=2):
    sec_systems = {sec_sys for sec_sys in systems_prots_dic.values()} # set comprehension for the secretion systems exists in the assay.
    pattern_set = set(filter(lambda x:x!='#', pattern))
    check_dic = {sec_sys:sec_sys in {systems_prots_dic.get(p) for p in pattern_set} for sec_sys in sec_systems} # the dictionary consits True/False in respect to the presence of proteins of the systems check{'T3SS': True, 'T6SS': False}
    return sum(check_dic.values()) >= uniq_sys
    
###
def flattened_contigs(fl):
    mltpl_cntgs_fl = fl[fl.loc[:,'num of patterns'] > 1] # Leaves only raws with duplicates patterns (a pattern that is maintained by more than one contig).
    mltpl_cntgs_fl.drop(columns=['num of patterns'], inplace=True) # Drop the num_of_patterns column.
    fltned_df_lst = list()

    for _, sub_df in mltpl_cntgs_fl.groupby(by='patterns'):
        pat, number_of_matches, pattern_length, contigs = sub_df.values[0]
        contigs = contigs.split(' , ')
        fltn_df_len = len(contigs)
        fltned_df_lst.append(pd.DataFrame(data={'patterns': [pat for _ in range(fltn_df_len)],
                                                'number of matches':[number_of_matches for _ in range(fltn_df_len)],
                                                'pattern length':[pattern_length for _ in range(fltn_df_len)],
                                                'contigs':contigs}))

    mltple_cntgs_df = pd.concat(fltned_df_lst + [fl[~(fl.loc[:,'num of patterns'] > 1)].drop(columns=['num of patterns'])], ignore_index=True)

    return mltple_cntgs_df
    # consize_by_contigs_df_lst = list()

#------------------------

def uniq(s): # The function returns a set of the unique values in the operons.
    return set(filter(lambda x:x!='#', re.split('--///--|--', s)))
    # return set(filter(lambda x:x!='#', s.split('--')))

def share_elements(s1, s2): # The function checks if two operons share any elements
    return len(uniq(s1).intersection(uniq(s2))) > 0

#------------------------

def unit_operons_by_contigs(fl, exclusion_features):

    exclusion_features_proteins = exclusion_features.get('proteins') if exclusion_features else list()
    
    def for_sort_by_contigs(s): # The function needed in order to find the longest operon, that is absent of exclusion features.
        if len(uniq(s).intersection(exclusion_features_proteins)) > 0: 
            return 0 
        else: 
            return len(uniq(s))
        
    def sum_united_elements(lst, index_to_sum): # The function sums the values that are being united to be part of the first value, and then deletes the values that has been united.
        lst[0] = sum([lst[i] for i in [0] + index_to_sum])
        lst = list(filter(lambda x:x, [None if i in index_to_sum else oprn for i, oprn in enumerate(lst)]))
        return lst

    def unit_operons_by_contigs_main(patterns, number_of_matches, pattern_length, exclusion_features_proteins):
        patterns = sorted(patterns, key=for_sort_by_contigs, reverse=True) # The operons are now sorted from big to small    
        index_to_drop = list()
        for i in range(1, len(patterns)):
            if len(uniq(patterns[i]).intersection(exclusion_features_proteins)) > 0:
                continue 
            if not share_elements(patterns[0], patterns[i]):
                patterns[0] = f'{patterns[0]}--///--{patterns[i]}'
                index_to_drop.append(i)

        patterns = list(filter(lambda x:x, [None if i in index_to_drop else oprn for i, oprn in enumerate(patterns)])) # filters the seperated version of the systems that were united. 
        number_of_matches = sum_united_elements(number_of_matches, index_to_drop)
        pattern_length = sum_united_elements(pattern_length, index_to_drop) 

        return patterns, number_of_matches, pattern_length

    consize_by_contigs_df_lst = list()

    for contig, sub_df in fl.groupby(by='contigs'):
        patterns, number_of_matches, pattern_length = [sub_df.loc[:,column_name].to_list() for column_name in ['patterns', 'number of matches', 'pattern length']]
        patterns, number_of_matches, pattern_length = unit_operons_by_contigs_main(patterns, number_of_matches, pattern_length, exclusion_features_proteins)
        consize_by_contigs_df_lst.append(pd.DataFrame(data={'patterns':patterns,
                                                            'number of matches':number_of_matches,
                                                            'pattern length':pattern_length,
                                                            'contigs':[contig for _ in patterns]}))
        
    return pd.concat(consize_by_contigs_df_lst, ignore_index=True) # a df containing all contigs, and there patterns. if the patterns are from different contig part they will be seperated by --///-- . e.g -> 'DotB--DotC--DotD--///--IcmJ_DotN--#--IcmE_DotG--IcmK_DotH--IcmL_DotI'

    #------------------------

def unit_operons_by_genomes(fl, exclusion_features, inclusion_features):

    def unit_operons_by_genomes_main(check_df, 
                                     inclusion_features_elements, inclusion_features_fragment,
                                     exclusion_features_elements, exclusion_features_fragment):

        res_df = check_df.copy()
        
        check_df['pattern_elements'] = check_df.patterns.apply(lambda pattern:uniq(pattern))
        check_df['pattern_size'] = check_df.patterns.apply(lambda pattern:len(uniq(pattern)))

        if callable(exclusion_features_fragment):
            check_df['exclusion_features'] = check_df.pattern_elements.apply(lambda pat_elm:exclusion_features_fragment([elm in exclusion_features_elements for elm in pat_elm]))
        else:
            check_df['exclusion_features'] = check_df.pattern_elements.apply(lambda pat_elm:dict(collections.Counter([False if elm in exclusion_features_elements else True for elm in pat_elm])))
            check_df['exclusion_features'] = check_df.exclusion_features.apply((lambda valid_exc_freq_dic: valid_exc_freq_dic.get(False, 0)/sum(valid_exc_freq_dic.values())))
            check_df['exclusion_features'] = check_df.exclusion_features*100
            check_df['exclusion_features'] = check_df.exclusion_features > exclusion_features_fragment

        check_df['inclusion_features'] = check_df.pattern_elements.apply(lambda pat_elm:inclusion_features_fragment([elm in inclusion_features_elements for elm in pat_elm]))
        
        check_df = check_df[check_df.inclusion_features]
        check_df = check_df[~(check_df.exclusion_features)]

        check_df = check_df.sort_values(['contigs', 'pattern_size'], ascending=[True, False])
        check_df.drop(['pattern_elements', 'pattern_size', 'exclusion_features', 'inclusion_features'], inplace=True, axis=1)

        grouped = check_df.groupby('contigs')

        if len(grouped) != 2: 
            return res_df

        contig_df1, contig_df2 = grouped
        contig1, contig2 = contig_df1[0], contig_df2[0]
        pattern1, pattern2 = contig_df1[1].patterns.to_list()[0], contig_df2[1].patterns.to_list()[0]
        number_of_matches1, number_of_matches2 = contig_df1[1].loc[:,'number of matches'].to_list()[0], contig_df2[1].loc[:,'number of matches'].to_list()[0]
        pattern_length1, pattern_length2 = contig_df1[1].loc[:,'pattern length'].to_list()[0], contig_df2[1].loc[:,'pattern length'].to_list()[0]

        if not share_elements(pattern1, pattern2):
            contig_df1 = pd.DataFrame(data={'patterns':contig_df1[1].patterns.to_list(),
                                            'number of matches':contig_df1[1].loc[:, 'number of matches'].to_list(),
                                            'pattern length':contig_df1[1].loc[:, 'pattern length'].to_list(),
                                            'contigs':[contig1 for _ in range(len(contig_df1[1]))],
                                            'GCA':[gca for _ in range(len(contig_df1[1]))]})
            
            cross_contigs_pattern = f'{pattern1}--///CCP/{contig2}/CCP///--{pattern2}' # CCP - Cross Contigs Patterns 
            contig_df1.loc[0, ['patterns', 'number of matches', 'pattern length']] = cross_contigs_pattern, (number_of_matches1+number_of_matches2), (pattern_length1+pattern_length2)
            
            contig_df2 = contig_df2[1].iloc[1::,:]

            check_df = pd.concat([contig_df1, contig_df2], ignore_index=True)
            res_df = check_df

        return res_df

    inclusion_features_elements, inclusion_features_fragment = inclusion_features.values() if type(inclusion_features) is dict else (list(), lambda _:True)
    exclusion_features_elements, exclusion_features_fragment = exclusion_features.values() if type(exclusion_features) is dict else (list(), lambda _:False)

    consize_by_contigs_df_lst = list()

    fl['GCA'] = fl.contigs.map(databases.gca_from_accession)
    for gca, sub_df in fl.groupby(by='GCA'):
        sub_df = sub_df.sort_values('number of matches', ascending=False).iloc[:2:,:]
        consize_by_contigs_df_lst.append(unit_operons_by_genomes_main(sub_df, 
                                                                      inclusion_features_elements, inclusion_features_fragment,
                                                                      exclusion_features_elements, exclusion_features_fragment))

    fl = pd.concat(consize_by_contigs_df_lst, ignore_index=True)
    fl.drop(columns=['GCA'], inplace=True, axis=1)

    return fl

    #------------------------

###
def main(operon_patterns_file, general_analyzer_file, system): 
    core_proteins_lst, minimal_core_proteins, minimal_general_matches, exclusion_features, inclusion_features = structure_default_values(system).values()

    ###
    fl = pd.read_table(operon_patterns_file)
    fl = flattened_contigs(fl)
    fl = unit_operons_by_contigs(fl, exclusion_features)
    if 'WGS_' in os.path.basename(operon_patterns_file):
        fl = unit_operons_by_genomes(fl, exclusion_features, inclusion_features)
    fl['patterns_as_lst'] = fl['patterns'].apply(lambda pattern : re.split('--///--|--', re.sub(r'[><]', '', (re.sub(r'--///CCP/.*/CCP///--', '--', pattern))))) # this line splits alson patterns from multiple DNA area in the contig and also fron other contigs from the same genome.

    fl = fl[[filter_by_core(pattern, core_proteins_lst, minimal_core_proteins) for pattern in fl.patterns_as_lst]] # filters the data frame with the help of filter_by_core function.
    if exclusion_features:
        fl = fl[[filter_by_exclusion(pattern, exclusion_features) for pattern in fl.patterns_as_lst]] # filters the data frame with the help of filter_by_exclusion function.
    if inclusion_features:
        fl = fl[[filter_by_inclusion(pattern, inclusion_features) for pattern in fl.patterns_as_lst]] # filters the data frame with the help of filter_by_inclusion function.

    fl.dropna(inplace=True)
    fl['core proteins'] = fl['patterns_as_lst'].apply(lambda pattern:pattern_info(pattern, core_proteins_lst)) # core proteins is now a dictionary.
    fl['sum of core proteins'] = fl['core proteins'].apply(lambda d:sum(d.values())) # sum of core proteins is sum of the True 'values' in core proteins dict.
    fl.rename(columns={'number of matches':'matches'}, inplace=True)
    fl = fl[fl['matches'] >= minimal_general_matches]
    fl.sort_values(['sum of core proteins', 'matches'], ascending=False, ignore_index=True, inplace=True)
    fl.drop(columns=['patterns_as_lst', 'core proteins', 'sum of core proteins'], inplace=True)
    fl_output_name = f"{operon_patterns_file.replace('.tsv','')}_{minimal_general_matches}_min_gnrl_prts_{minimal_core_proteins}_min_core_prts.tsv"

    if fl.empty:
        return 'no_systems_passed_filtration'
    
    fl.to_csv(fl_output_name, index=None, sep='\t')
    with open(fl_output_name, 'a') as out_put_file:
        out_put_file.write('\n')
        out_put_file.write(f'parameters: minimal_general_matches-{minimal_general_matches}, core_proteins_lst-{core_proteins_lst}, minimal_core_proteins_matches-{minimal_core_proteins}')

    contigs_to_keep = set(fl.contigs) # the set here is just for me to emphasize that these is a set of contigs. as matter effect the df is already absenct of contig duplications.
    for pat in fl.patterns:
        search_res = re.search(r'(?<=--///CCP/).*(?=/CCP///--)', pat) # regex to extract the additional the contig name from the pattern itself.
        if search_res: 
            contigs_to_keep.add(search_res.group()) 
    patterns_to_keep = set(itertools.chain.from_iterable([re.sub(r'--///CCP/.*/CCP///--', '--///--', sys).split('--///--') for sys in fl.patterns])) # regex to find and replace the CCP marker.
    
    ###
    gf = pd.read_table(general_analyzer_file)

    gf = gf[gf.contig.isin(contigs_to_keep)]
    pattern_filter_func = lambda pat, lst : (pat in lst) or ('--'.join(pat.split('--')[::-1]) in lst)
    gf = gf[[pattern_filter_func(ptrn, patterns_to_keep) for ptrn in gf.pattern]]

    gf['patterns_as_lst'] = gf['pattern'].apply(lambda pattern:pattern.split('--'))
    gf['matches'] = gf['patterns_as_lst'].apply(lambda ptrn_lst:len(list(filter(lambda x:x!='#', ptrn_lst))))
    gf.sort_values(by=['contig', 'matches'], ascending=[True, False], ignore_index=True, inplace=True)
    gf.drop(columns=['matches', 'patterns_as_lst'], inplace=True)
    gf.to_csv(general_analyzer_file, index=None, sep='\t') # filtered overide of the orignial 
    
    return general_analyzer_file

### main
if __name__ == "__main__":
    argparse = argparse.ArgumentParser()
    argparse.add_argument('-P', '--operon_patterns_file', type=str, help='tsv file, one of the outputs of operon_sorting.py hopefully with tax implemented within it already')
    argparse.add_argument('-G', '--general_analyzer_file', type=str, help='tsv file, one of the outputs of operon_sorting.py')
    argparse.add_argument('-S', '--system', type=str, choices=['T3SS','T4SSA','T4SSB','T6SS'], required=True, help='system_name')

    input_details = argparse.parse_args()
    operon_patterns_file = input_details.operon_patterns_file
    general_analyzer_file = input_details.general_analyzer_file
    system = input_details.system
    
    main(operon_patterns_file, general_analyzer_file, system)
