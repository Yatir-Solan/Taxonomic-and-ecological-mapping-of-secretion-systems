# modules import:
import re
import os
import json
import pandas as pd
#

########################### genomic datasets ########################### genomic datasets ########################### genomic datasets ###########################

def genomic_datasets():
    return {'denovo':r'/davidb/assemblies/analysis/subdatasets/all_denovo.min10k.faa',
            'Mgnify':r'/davidb/assemblies/analysis/subdatasets/all_Mgnify.min10k.faa',
            'GEM':r'/davidb/bio_db/GEM/all_GEM_proteins.faa',
            'WGS_Uncultured':r'/davidb/assemblies/analysis/subdatasets/all_WGS_UnculturedGenomes.min10k.faa',
            'WGS_Metagenomes':r'/davidb/assemblies/analysis/subdatasets/all_WGS_Metagenomes.min10k.faa',
            'WGS_Genomes':r'/davidb/assemblies/analysis/subdatasets/all_WGS_Genomes.min10k.cln2.faa'} # a filtered version of wgs_genomes. within this version extrimly large proteins were filtered.
            # 'WGS_Genomes':r'/davidb/assemblies/analysis/subdatasets/all_WGS_Genomes.min10k.faa'}

def rglr_expr():
    # regex for every dataset type given as dictionary.
    # the regex itself is a specific classifier for every database.
    return {'DeNovo':r'(?<=\.)[SE]RR[0-9]+',
            'Mgnify':r'(?<=\|)ERZ[0-9]+',
            'WGS':r'(?<=\|)GCA_[0-9]+\.*[0-9]*(?=\|)*',
            'GEM':r'_GEM_'}

def genomic_source_classifire(contig):
    for source, reg in rglr_expr().items():
        match = re.search(reg, contig)
        if match:
            return source
    return None

########################### refseq section ########################### refseq section ########################### refseq section ###########################

def ref_seq_file_paths(key):
    return {'bacteria_SeqIO_dic_pickle':r'/davidb/yatirsolan/review_tree/Bacteria.faa_SeqIO_dictionary.pkl',
            'all_fasta':r'/davidb/bio_db/UniProt/RefProteome/2021_03/all.RefProteome.faa',
            'all_fasta_for_blast':r'/davidb/yatirsolan/databases/RefProteome/all_RefProteome_2021_03/all.RefProteome.faa',
            'bacteria_fasta':r'/davidb/bio_db/UniProt/RefProteome/2021_03/Bacteria.faa',
            'details_json':r'/davidb/yatirsolan/databases/RefProteome/all.RefProteome.faa.json'}.get(key)
            
def ref_seq_details(ref_seq_path=ref_seq_file_paths('all_fasta')):
    # the func return a dictionary consits relevant details (tax, taxid, functional annotation) of any header in the refseq database.
    refseq_json_path = ref_seq_file_paths('details_json')
    
    if os.path.isfile(refseq_json_path):
        refseq_details_dic = json.load(open(refseq_json_path)) 
        # creating a dictionary out the refseq data -> using the json file that was maid previously exactley as the 'else' section would have done.
        # 'tr|K8E1J6|K8E1J6_CARML', {'tax': 'Carnobacterium maltaromaticum', 'id': '1234679', 'id': '1234679', 'annot': 'Aldehyde dehydrogenase'}
    else:
        refseq_details_dic = dict()
        with open(ref_seq_path) as ref_seq_file:
            r_dic = {'annot':re.compile(r'(?<=\s).+(?= OS=)'),
                     'tax':re.compile(r'(?<=OS=).+(?= OX=)'),
                     'id':re.compile(r'(?<=OX=)\d+\S')}
            for line in ref_seq_file:
                if line[0] == '>':
                    hd = line.split()[0].replace('>','')
                    refseq_details_dic[hd] = dict()
                    for t,r in r_dic.items():
                        refseq_details_dic.get(hd)[t] = r.search(line).group() 
                        # creating a dictionary out the refseq data -> /
                        # 'tr|K8E1J6|K8E1J6_CARML', {'tax': 'Carnobacterium maltaromaticum', 'id': '1234679', 'id': '1234679', 'annot': 'Aldehyde dehydrogenase'} 
    return refseq_details_dic

########################### GCA ########################### GCA ########################### GCA ###########################

def gca_taxids_dic():
    gca_file = r'/davidb/bio_db/NCBI/WGS/latest.assembly_summary_genbank.wLineage.no_fungi_metazoa_viridiplantae.tab'
    return dict(pd.read_table(gca_file, usecols=['# assembly_accession', 'taxid']).values)

def gca_from_accession(accession):
    res = re.search(r'(?<=\|)GCA_\d+\.*.(?=\|)',accession)
    res = re.search(r'((?<=\|)GCA_\d+\.*.(?=\|))|((?<=\|)GCA_\d+\.*.$)',accession)
    if res:
        return res.group()
    return None

########################### kegg ########################### kegg ########################### kegg ###########################

def kegg_hmms_database():
    return r'/davidb/bio_db/Kegg/2021-05-14/kg.05_21.ren4prok.2.hmm.db'

########################### GEM ############################ GEM ############################ GEM #############################