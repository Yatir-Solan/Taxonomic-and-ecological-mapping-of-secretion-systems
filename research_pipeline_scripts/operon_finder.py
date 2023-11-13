# modules import:
import re
import argparse
import operator
import pandas as pd

Gff_line_direc_dict = dict()

# the function is activated in respect to whethere a gff file was given in the input or not.
# if it is given the func will produce dictionary (Gff_line_direc_dict) of ORFs (as keys) and translation direction (as values). 
       
def produce_direc_dict(gff_file):
    direc_df = pd.read_table(gff_file,header=None).iloc[:, [8,6]]
    direc_df.columns = ['ORF', 'direc']
    direc_df.iloc[:,0] = direc_df.iloc[:,0].apply(lambda desc : desc[3:desc.find(';'):])
    return dict(direc_df.values)

class HmmSearchLine:
    def __init__(self, line):
        self.Line = line.split()
        self.target_name = self.Line[0]
        self.accession1 = self.Line[1]
        self.query_name = self.Line[2]
        self.accession2 = self.Line[3]
        self.full_sequence_Evalue = float(self.Line[4])
        self.full_sequence_score = float(self.Line[5])
        self.bias1 = float(self.Line[6])
        self.best_domain_Evalue = float(self.Line[7])
        self.best_domain_score = float(self.Line[8])
        self.bias2 = float(self.Line[9])
        self.exp = float(self.Line[10])
        self.reg = float(self.Line[11])
        self.clu = float(self.Line[12])
        self.ov = float(self.Line[13])
        self.env = float(self.Line[14])
        self.dom = float(self.Line[15])
        self.rep = float(self.Line[16])
        self.inc = float(self.Line[17])
        self.description = ' '.join(self.Line[18::])
        # metagenomic paramters:
        self.sample = None
        self.contig = None
        self.ORF = None
        self.direc = None

    def __hash__(self):
        return hash(self.query_name)

    def __eq__(self,other):
        return self.query_name == other.query_name and self.accession2 == other.accession2 
        # this function was modified (12.10.21) in order for it to work on proteins in which their HMMs include 
        # more then one HMM, in that case the accession2 was considred a second identity to be compare with. 

    def __repr__(self):
        repr_str = f'(ORF:{str(self.ORF)},{self.query_name[self.query_name.find("-")+1::]},{self.accession2},{str(self.full_sequence_Evalue)},{str(self.full_sequence_score)},{str(self.direc)})'
        return repr_str
        
    # the function get HmmSearchLine object, and return it's ORF and Contig as tupple. 
    # (of course it has to be a metgenomic format). it contains some REGEX expressions.
    def get_contig_and_ORF_from_targetName(self):
        target_name = self.target_name
        if 'ctg' in target_name: # meaning that it is a metagenomnic annotation.
            contig_and_ORF = re.search(r'ctg_{0,1}[0-9]+_[0-9]+',target_name).group()
            contig_only = re.search(r'(?<=ctg)(_{0,1})[0-9]+(?=_[0-9]+)',contig_and_ORF).group().replace('_','') # not in use as for now
            ORF_only = contig_and_ORF[len('ctg_')+len(contig_only)+1::] if contig_and_ORF.count('_') > 1 else contig_and_ORF[len('ctg')+len(contig_only)+1::]
            full_contig = target_name[:-len(ORF_only)-1:]
        elif 'contig' in target_name:
            contig_and_ORF = re.search(r'contig_{0,1}[0-9]+_[0-9]+',target_name).group()
            contig_only = re.search(r'(?<=contig)(_{0,1})[0-9]+(?=_[0-9]+)',contig_and_ORF).group().replace('_','') # not in use as for now
            ORF_only = target_name.split("_")[-1]
            full_contig = target_name[:-len(ORF_only)-1:]
        elif 'GCA' in target_name : # meaning that it is an ncbi wgs annotation.
            contig_and_ORF = target_name.split('|')[-1]
            contig_only = contig_and_ORF.split('_')[0]
            ORF_only = str(contig_and_ORF.split('_')[1])
            full_contig = target_name[:-len('_'+ORF_only):]
        else:
            contig_and_ORF = target_name.rsplit('|',1)[-1]
            contig_only = contig_and_ORF.rsplit('_',1)[0]
            ORF_only = str(contig_and_ORF.rsplit('_',1)[-1])
            full_contig = target_name.rsplit('_',1)[0]
        return (str(full_contig),int(ORF_only))

# the function gets HmmSearchFile and return list of HmmSearchLine objects, which there Evalue is significant (lower) than the cutoff. 
def get_lines_from_hmmSearchFile(HmmSearchFile, Evalue_cutoff):
    lines = list()
    with open(HmmSearchFile, 'r') as file:
        for line in file:
            if line.startswith(('#', 'targetName')): # suits both to tsv and regular text
                continue
            Line = HmmSearchLine(line)
            if Line.full_sequence_Evalue <= Evalue_cutoff :
                lines.append(Line)
    return lines

def main_arranger(list_of_HmmSearchLine_objects):
    HmmSearchLines = list_of_HmmSearchLine_objects
    contig_dict = dict()
    for Line in HmmSearchLines :
        Line.contig,Line.ORF = Line.get_contig_and_ORF_from_targetName() # specific function of the HmmSearchLine object - gets the contig and ORF number. the one that uses regex...
        Line.direc = Gff_line_direc_dict[Line.target_name] if Line.target_name in Gff_line_direc_dict.keys() else None # specific function of the HmmSearchLine object - extracts the direction of translation.
        # the next if an else statement sort all the ORF as HmmSearchLine objects as the values for each contig which are the keys. the HmmSearchLines are stored as list.
        if Line.contig in contig_dict:
            contig_dict[Line.contig].append(Line)
        else:
            contig_dict[Line.contig] = [Line]
    for key in contig_dict.keys():
        contig_dict[key] = sorted(contig_dict.get(key), key=operator.attrgetter('ORF')) # so, the lists are now sorted based on their ORFs.
    return contig_dict

class Contig:
    ### static section ###
    def list_to_list_of_lists(lst):
        list_of_lists = []
        for elem in lst:
            list_of_lists.append([elem])
        return list_of_lists

    # recursive function that gets list of list of sorted ORFs, and arrange them as operon candidate. 
    # (each ORF is represented as HmmSearchLine object, thay are sorted in respect to their ORF).
    def contig_arranger(contig_list, max_space_betw_ORFs):
        def stop_cond(lst, sbORFs):
            for i in range(len(lst)-1):
                if lst[i+1][0].ORF - lst[i][-1].ORF <= max_space_betw_ORFs :
                    return False
            return True
        if stop_cond(contig_list, max_space_betw_ORFs):
            return contig_list
        for i in range(len(contig_list)):
            if contig_list[i+1][0].ORF - contig_list[i][-1].ORF <= max_space_betw_ORFs :
                contig_list[i].append(contig_list[i+1][0])
                del contig_list[i+1]
                return Contig.contig_arranger(contig_list, max_space_betw_ORFs)
    # end of static section

    def __init__(self, contig_name, list_of_ORFs, min_operon_length, max_space_betw_ORFs, min_unq_prts_in_operon): # list of ORF's - each ORF is actually an HmmSearchLine object.
        self.contig_name = contig_name
        self.taxonomy = None # was thought by me that sometime when everything will run autimatically I shall put the analyzed taxonomy here (probabely wont happen...)
        self.min_operon_length = min_operon_length
        self.min_unq_prts_in_operon = min_unq_prts_in_operon
        self.list_of_operons = [Operon(self, operon_cand) for operon_cand in Contig.contig_arranger(Contig.list_to_list_of_lists(list_of_ORFs), max_space_betw_ORFs)]
        self.contig_fitler()
        self.operon_num = len(self.list_of_operons)
        
        
    # this function filters the contig out of operons that are smaller then desired or does not have enough unique desired protiens. (desired means - by the script input) 
    def contig_fitler(self):
        copy_of_list_of_operons = self.list_of_operons.copy()
        for operon in self.list_of_operons:
            if operon.operon_length < self.min_operon_length or operon.unique_proteins_num < self.min_unq_prts_in_operon :
                copy_of_list_of_operons.remove(operon)
        self.list_of_operons = copy_of_list_of_operons
                
class Operon:
    # static #    
    # the function adds '#' to represented un-idetified proteins matches.
    def add_number_signs(sublist):
        new_sublist = sublist.copy()
        space_count = 0
        for i in range(len(sublist)-1):
            space = abs(sublist[i+1].ORF-sublist[i].ORF)
            if space > 1 :
                for j in range(space-1):
                    space_count += 1
                    new_sublist.insert(i+space_count, '#')
        return new_sublist
    
    def number_to_strings(sublist):
        new_sublist = sublist.copy()
        for i in range(len(new_sublist)):
            new_sublist[i] = str(new_sublist[i])
        return new_sublist
        
    def ORF_filter(sublist): # when two hits are for the same ORF the function discard the less significant one (E-value). when equal E-value, the ORF score will be considerd. return : (filtered operon,list of deleted ORFs).
        ORF_dict = dict()
        final_sublist = list()
        for i in range(len(sublist)):
            if sublist[i].ORF in ORF_dict.keys():
                ORF_dict[sublist[i].ORF].append(sublist[i])
            else:
                ORF_dict[sublist[i].ORF] = [sublist[i]]
        for ORF in ORF_dict.keys():
            final_sublist.append(min(ORF_dict.get(ORF), key=operator.attrgetter('full_sequence_Evalue')))     
        final_sublist = sorted(final_sublist, key=operator.attrgetter('ORF'))
        deleted_ORFs = list(filter(lambda HmmSearcLine : HmmSearcLine not in final_sublist, sublist))
        return (final_sublist, deleted_ORFs) 
        
    # end of static #
    
    def __init__(self, contig, sublist):
        self.contig = contig
        self.contig_name = self.contig.contig_name
        self.original_sublist = sublist
        self.sublist, self.deleted_ORFs = Operon.ORF_filter(self.original_sublist) # list of HmmSearchLine objects orderd as ORFs in an operon state -> filtered to be without duplications.
        self.target_names_list = [ORF.target_name for ORF in self.sublist]
        self.unique_proteins_num = len(set(self.sublist))
        self.tmp = Operon.add_number_signs(self.sublist)
        self.pattern = Operon_pattern(self.contig_name, self.tmp)
        self.operon_for_repr = Operon.number_to_strings(self.tmp)
        self.operon_length = len(self.sublist)
      
    def __repr__(self):
        return '--'.join(self.operon_for_repr)
    
class Operon_pattern:
    def __init__(self, contig_name, list_of_ORF_and_numberSigns):
        self.contig_name = contig_name
        if len(Gff_line_direc_dict.keys()) == 0:
            self.pattern = [ORF.query_name if type(ORF) is HmmSearchLine else ORF for ORF in list_of_ORF_and_numberSigns]
            self.pattern_list = [elem[elem.find('-')+1::] for elem in self.pattern]
        else:
            self.pattern = [('>>' if ORF.direc=='+' else '<<') + ORF.query_name + ('>>' if ORF.direc=='+' else '<<') if type(ORF) is HmmSearchLine else ORF for ORF in list_of_ORF_and_numberSigns]
            self.pattern_list = [elem[:elem.find('LOG'):]+elem[elem.find('-')+1::] for elem in self.pattern]
        self.pattern_length = len(self.pattern)
        self.num_of_hits_in_pattern = len(list(filter(lambda x:x!='#' and x!='<<' and x!='>>',self.pattern)))
        
    def __eq__(self, other):
        pattern1 = [elem.replace('<<','').replace('>>','') for elem in self.pattern_list]
        pattern2 = [elem.replace('<<','').replace('>>','') for elem in other.pattern_list]
        return pattern1 == pattern2 or pattern1[::-1] == pattern2

    def __hash__(self):
        return hash(''.join(sorted([elem.replace('<<','').replace('>>','') for elem in self.pattern_list])))
        
    def __repr__(self):
        return '--'.join(self.pattern_list)

def main(HmmSearchFile, Evalue, max_space_betw_ORFs, min_operon_length, min_unq_prts_in_operon, gff_file=None):
    files_base_name = f'{HmmSearchFile}_'
    if gff_file != None:
        Gff_line_direc_dict = produce_direc_dict(gff_file) # if a gff file is given in the script run, the function will be activated.
    contigs_dict = main_arranger(get_lines_from_hmmSearchFile(HmmSearchFile, Evalue))
    dict_for_general = {col:list() for col in ['contig', 'operon length', 'operon', 'pattern', 'deleted duplications', 'ORFs']}
    patterns_dict = dict()
    main_operon_lst = list()
    for key in contigs_dict :
        contig = Contig(key, contigs_dict[key], min_operon_length, max_space_betw_ORFs, min_unq_prts_in_operon)
        if contig.operon_num > 0 :
            main_operon_lst += contig.list_of_operons
            
    for operon in sorted(main_operon_lst, key=operator.attrgetter('operon_length'), reverse=True):
        if operon.pattern in patterns_dict.keys():
            patterns_dict[operon.pattern][0] += 1
            patterns_dict[operon.pattern][1].append(operon.pattern.contig_name)
        else:
            patterns_dict[operon.pattern] = [1, [operon.pattern.contig_name]]
        dict_for_general['contig'].append(str(operon.contig_name))
        dict_for_general['operon length'].append(operon.operon_length)
        dict_for_general['operon'].append(str(operon))
        dict_for_general['pattern'].append(str(operon.pattern))
        dict_for_general['deleted duplications'].append(str(operon.deleted_ORFs))
        dict_for_general['ORFs'].append(operon.target_names_list)
        
    dict_for_patterns = {col:list() for col in ['patterns', 'num of patterns', 'number of matches', 'pattern length', 'contigs']}
    
    for elem in sorted(patterns_dict.items(), key=lambda x:x[0].num_of_hits_in_pattern, reverse=True):
        dict_for_patterns['patterns'].append(elem[0])
        dict_for_patterns['num of patterns'].append(elem[1][0])
        dict_for_patterns['number of matches'].append(elem[0].num_of_hits_in_pattern)
        dict_for_patterns['pattern length'].append(elem[0].pattern_length)
        dict_for_patterns['contigs'].append(' , '.join(elem[1][1]))     
    pf = pd.DataFrame.from_dict(dict_for_patterns)
    gf = pd.DataFrame.from_dict(dict_for_general)
    pf.to_csv(f'{files_base_name}oprn_ptrns.tsv', index=None, sep='\t')
    gf.to_csv(f'{files_base_name}oprn_gnrl_anlyz.tsv', index=None, sep='\t')
    return f'{files_base_name}oprn_ptrns.tsv', f'{files_base_name}oprn_gnrl_anlyz.tsv'

if __name__ == '__main__':
    argparse = argparse.ArgumentParser()
    pd.options.display.max_colwidth = 400000
    argparse.add_argument('-F', '--HmmSearchFile', type=str, help='HmmSearchFile (original Hmmer file and not a modified TSV file...)')
    argparse.add_argument('-E', '--Evalue', type=float, help='Evalue cutoff')
    argparse.add_argument('-S', '--max_space_betw_ORFs', type=int, help='max_space_betw_ORFs')
    argparse.add_argument('-M', '--min_operon_length', type=int, help='min_operon_length')
    argparse.add_argument('-G', '--gff_file', type=str, default=None, help='gff_file')
    argparse.add_argument('-U', '--min_unq_prts_in_operon', default='', type=int, help='min_unq_prts_in_operon')
    input_details = argparse.parse_args()

    HmmSearchFile = input_details.HmmSearchFile
    Evalue = input_details.Evalue
    max_space_betw_ORFs = input_details.max_space_betw_ORFs
    min_operon_length = input_details.min_operon_length
    gff_file = input_details.gff_file
    min_unq_prts_in_operon = input_details.min_unq_prts_in_operon

    main(HmmSearchFile, Evalue, max_space_betw_ORFs, min_operon_length, min_unq_prts_in_operon, gff_file)