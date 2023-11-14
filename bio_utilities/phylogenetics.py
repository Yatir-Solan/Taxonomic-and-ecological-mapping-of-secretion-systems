# modules import:
import json
import pandas as pd
from ete3 import NCBITaxa
# python scripts import:

def taxid_lineage():
    return dict(pd.read_table(r'/davidb/bio_db/Taxonomy/2022-05-30/taxid_lineage.tab', header=None).values)

def taxid_mapping():
    return json.load(open(r'/davidb/yatirsolan/ncbi_taxids/ncbi_taxids_phylo_ranks_map.json')) # here, all items are string with the execption of None values.

def taxid_to_name(taxid, ncbi=NCBITaxa()):
    try :
        return list(ncbi.get_taxid_translator([taxid]).values())[0]
    except IndexError:
        return None

def name_to_taxid(name, ncbi=NCBITaxa()):
    try :
        return str(list(ncbi.get_name_translator([name]).values())[0][0])
    except IndexError:
        return None

def uninformative_taxids():
    return {'2759': 'cellular organisms->Eukaryota',
            '2': 'cellular organisms->Bacteria',
            '1427524': 'unclassified entries->unclassified sequences->mixed sample',
            '131567': 'cellular organisms'}

def customize_blacklist():
    # this function returns taxids values that could not have been found within functions as 'ncbi.get_descendant_taxa' or others.
    # for instance the taxid 177416 of one of francisella strain, is constantly in mmseq results and the result of {34064: 'Francisellaceae'} 
    # in the function ncbi.get_descendant_taxa don't include it for some reason. in that case I want to blacklist this taxid from 
    # the mmseq assay in first place so possession_verifier script will work good.
    return list()
    return [177416]

def ribosomal_HMMs(key):
    return {'HMMs':r'/davidb/yatirsolan/taxonomy_check/Bacterial_check/bacterial_ribosomal_proteins.Hmm',
            'hmm_search':r'/davidb/yatirsolan/taxonomy_check/Bacterial_check/bacterial_ribosomal_proteins_vs_Bacteria.faa.hmmsearch',
            'hmm_search_tsv':r'/davidb/yatirsolan/taxonomy_check/Bacterial_check/bacterial_ribosomal_proteins_vs_Bacteria.faa.hmmsearch.tsv'}.get(key)

def family_to_phylum(fmly_nme, tax_map_dic=taxid_mapping(), proteobacteria='1224'):
    if fmly_nme.isdigit():
        fmly_txID = fmly_nme
    else:
        fmly_txID = name_to_taxid(fmly_nme)
    phylum = tax_map_dic.get(fmly_txID).get('phylum')
    if not phylum:
        return None
    phylum = phylum if phylum!=proteobacteria else tax_map_dic.get(fmly_txID).get('class') # '1224' is Proteobacteria 
    phylum = taxid_to_name(phylum)
    return phylum

def txid_to_species(txid, tax_map_dic=taxid_mapping()):
    species = tax_map_dic.get(txid).get('species')
    if not species:
        return None
    return species

def txid_to_choose(txid, rnk, tax_map_dic=taxid_mapping()):
    rank_det = tax_map_dic.get(txid)
    if not rank_det:
        return None
    rank = rank_det.get(rnk)
    if not rank:
        return None
    return rank