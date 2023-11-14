import re
from dna_features_viewer import GraphicFeature,GraphicRecord

def main(path):
    if 'GEM' in path:
        res_dic = GEM_parse(path)
    else:
        res_dic = gff_parse(path)
    return res_dic

def gff_parse(gff_path):
    output_dic = {}
    id_reg = r'(?<=ID=).*'
    with open(gff_path) as gff_file:
        for line in gff_file:
            if line[0] == '>':
                break

            if '##gff-version' in line or '##FASTA' in line:
                continue

            if line[:2:] == '##': # the dictionary contigs are imported from the first region of the gff file.
                line = line.split()
                output_dic[line[1]] = list()
                continue

            line = line.split()
            output_dic.get(line[0]).append({'id':re.search(id_reg,line[8].split(';',1)[0]).group(),'source':line[1],'annotation':line[2],'start':int(line[3]),'end':int(line[4]),'strand':line[6]})  
    return output_dic

def GEM_parse(gem_faa_path):
    output_dic = {}
    with open(gem_faa_path) as gem_faa:
        for line in gem_faa:

            if line[0] != '>':
                continue

            line = line.split()
            contig = line[0].replace('>','').rsplit('_',1)[0] 
            orf = line[0].replace('>','')

            if contig not in output_dic.keys():
                output_dic[contig] = list() 

            output_dic.get(contig).append({'id':orf,'source':'GEM','annotation':'GEM_orf_annotation','start':int(line[3]),'end':int(line[5]),'strand':line[7]})
    return output_dic