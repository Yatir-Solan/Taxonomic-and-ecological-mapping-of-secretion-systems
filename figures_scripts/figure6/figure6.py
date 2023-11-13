import os 
import json
import pickle
import collections
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
plt.style.use('seaborn-white')

# --------------------------------------------------------------------------------------------------------------------------------------------

def GEM_contigs(main_distribution_dic):
    # the function gets the main distribution dictionary and return all GEM contigs out of it.
    GEM_contigs = list()
    for dictionaries in main_distribution_dic.values():
        for contigs in dictionaries.values():
            GEM_contigs.extend(filter(lambda contig:'_GEM_' in contig, contigs))
    return GEM_contigs

# --------------------------------------------------------------------------------------------------------------------------------------------

GEM_mtdta = r'/davidb/bio_db/GEM/genome_metadata.tsv'
GEM = pd.read_table(GEM_mtdta, usecols=['genome_id', 'metagenome_id', 'ecosystem_category', 'ecosystem_type', 'longitude', 'latitude'])

# --------------------------------------------------------------------------------------------------------------------------------------------

T3SS_main_distribution_dic = pickle.load(open(r'/davidb/yatirsolan/data_presentation/family/T3SS/review_family_T3SS_systems_distribution.pkl','rb'))
T3SS_contigs = GEM_contigs(T3SS_main_distribution_dic)
T3SS_metagenome_ids = [int(contig.split('_GEM')[0]) for contig in T3SS_contigs]
T3SS_metagenome_ids_counter = dict(collections.Counter(T3SS_metagenome_ids))
GEM['T3SS_encounters'] = GEM.metagenome_id.apply(lambda metagenome_id:T3SS_metagenome_ids_counter.get(metagenome_id, 0))
T3SS_metagenome_ids = {int(contig.split('_GEM')[0]) for contig in T3SS_contigs} # no duplicates 
GEM['T3SS_appearance'] = GEM.metagenome_id.apply(lambda metagenome_id: True if metagenome_id in T3SS_metagenome_ids else False)

T4SS_main_distribution_dic = pickle.load(open(r'/davidb/yatirsolan/data_presentation/family/T4SS/review_family_T4SS_systems_distribution.pkl','rb'))
T4SS_contigs = GEM_contigs(T4SS_main_distribution_dic)
T4SS_metagenome_ids = [int(contig.split('_GEM')[0]) for contig in T4SS_contigs]
T4SS_metagenome_ids_counter = dict(collections.Counter(T4SS_metagenome_ids))
GEM['T4SS_encounters'] = GEM.metagenome_id.apply(lambda metagenome_id:T4SS_metagenome_ids_counter.get(metagenome_id, 0))
T4SS_metagenome_ids = {int(contig.split('_GEM')[0]) for contig in T4SS_contigs} # no duplicates 
GEM['T4SS_appearance'] = GEM.metagenome_id.apply(lambda metagenome_id: True if metagenome_id in T4SS_metagenome_ids else False)

T6SS_main_distribution_dic = pickle.load(open(r'/davidb/yatirsolan/data_presentation/family/T6SS/review_family_T6SS_systems_distribution.pkl','rb'))
T6SS_contigs = GEM_contigs(T6SS_main_distribution_dic)
T6SS_metagenome_ids = [int(contig.split('_GEM')[0]) for contig in T6SS_contigs]
T6SS_metagenome_ids_counter = dict(collections.Counter(T6SS_metagenome_ids))
GEM['T6SS_encounters'] = GEM.metagenome_id.apply(lambda metagenome_id:T6SS_metagenome_ids_counter.get(metagenome_id, 0))
T6SS_metagenome_ids = {int(contig.split('_GEM')[0]) for contig in T6SS_contigs} # no duplicates 
GEM['T6SS_appearance'] = GEM.metagenome_id.apply(lambda metagenome_id: True if metagenome_id in T6SS_metagenome_ids else False)

# --------------------------------------------------------------------------------------------------------------------------------------------

GEM = GEM[~((GEM.latitude.isna()) | (GEM.longitude.isna()))] # filter bad coordinates
threshold = 100
ecosystems_count = dict(GEM.ecosystem_category.value_counts())
less_than_thresh_ecosystems = list(dict(filter(lambda items:items[1] < threshold, ecosystems_count.items())).keys())
GEM = GEM[~(GEM.ecosystem_category.isin(less_than_thresh_ecosystems))]

# --------------------------------------------------------------------------------------------------------------------------------------------

def ratio(appearances_col, ecos_GEM):
    ecos_ratio_dic = {grp:sliced_GEM.loc[:,appearances_col].sum()/len(sliced_GEM) for grp, sliced_GEM in ecos_GEM.groupby('ecosystem_category')}    
    ecos_ratio_dic = dict(sorted(ecos_ratio_dic.items(), key=lambda items:items[1], reverse=False)) # sorting subtype dict so that the family with higher systems will appear first.
    ecos_ratio_dic = dict(filter(lambda items:items[1] > 0.005, ecos_ratio_dic.items()))
    return ecos_ratio_dic

# --------------------------------------------------------------------------------------------------------------------------------------------

def dual_plot(geo_GEM, title_lbl, encntrs_cnt_clm, ecos_ratio_dic):
    color_dict = json.load(open(r'/davidb/yatirsolan/thesis_work/figures/echosystems_distribution/ecosystems_color_nrmlzd.json'))
    labels = list(ecos_ratio_dic.keys())
    values = list(ecos_ratio_dic.values())

    if title_lbl in ['T3SS', 'T6SS']:
        fig, (geo_ax, ecos_ax) = plt.subplots(nrows=1, ncols=2, figsize=(15,6), edgecolor='w', gridspec_kw=dict(width_ratios=[5,1]))
        ecos_ax.barh(y=labels, 
                    width=values,
                    height=0.7,
                    color=[color_dict.get(ecosystem) for ecosystem in labels])
        ecos_ax.spines['top'].set_visible(False)
        ecos_ax.spines['right'].set_visible(False)
        ecos_ax.set_xlim(0,1)
        ecos_ax.set_xticks(np.arange(0, 1.1, 0.1))
        ecos_ax.tick_params(labelcolor='black', labelsize=9, width=1, length=5, color='black')
        ecos_ax.set_ylabel(ylabel='Ecosystem')
        ecos_ax.set_xlabel(xlabel='Appearance ratio')
        for i in np.arange(len(labels)):
            plt.text(x=values[i]+.0125, 
                    y=labels[i], 
                    s=round(values[i], 2), 
                    size=9,
                    va='center_baseline')

    elif title_lbl == 'T4SS': # a version to get the bar
        fig, (geo_ax, ecos_ax) = plt.subplots(nrows=1, ncols=2, figsize=(15,6), edgecolor='w', gridspec_kw=dict(width_ratios=[5,1]))

        left, bottom, width, height = ecos_ax.get_position().bounds
        ecos_ax.set_position([left, bottom+1.5, width, height-1.5])

        ecos_ax.barh(y=labels, 
                    width=values,
                    height=0.7,
                    color=[color_dict.get(ecosystem) for ecosystem in labels])
        ecos_ax.spines['top'].set_visible(False)
        ecos_ax.spines['right'].set_visible(False)
        ecos_ax.set_xlim(0,1)
        ecos_ax.set_xticks(np.arange(0, 1.1, 0.1))
        ecos_ax.tick_params(labelcolor='black', labelsize=9, width=1, length=5, color='black')
        ecos_ax.set_ylabel(ylabel='Ecosystem')
        ecos_ax.set_xlabel(xlabel='Appearance ratio')
        for i in np.arange(len(labels)):
            plt.text(x=values[i]+.0125, 
                    y=labels[i], 
                    s=round(values[i], 2), 
                    size=9,
                    va='center_baseline')

    # elif title_lbl == 'T4SS': # a version to get the map
    #     fig, (geo_ax, ecos_ax) = plt.subplots(nrows=1, ncols=2, figsize=(15,6.5), edgecolor='w', gridspec_kw=dict(width_ratios=[7,1]))
    #     fig.subplots_adjust(wspace=1)
    #     labels = labels[::-1]
    #     values = values[::-1]
    #     ecos_ax.bar(x=labels, 
    #                 height=values,
    #                 width=.65,
    #                 color=[color_dict.get(ecosystem) for ecosystem in labels])
    #     plt.xticks(rotation=45, rotation_mode='anchor', ha='right')
    #     ecos_ax.spines['top'].set_visible(False)
    #     ecos_ax.spines['right'].set_visible(False)
    #     ecos_ax.set_ylim(0,1)
    #     ecos_ax.set_yticks(np.arange(0, 1.1, 0.1))
    #     ecos_ax.tick_params(labelcolor='black', labelsize=9, width=1, length=5, color='black')
    #     ecos_ax.set_xlabel(xlabel='Ecosystem')
    #     ecos_ax.set_ylabel(ylabel='Appearance ratio')
    #     for i in np.arange(len(labels)):
    #         plt.text(x=labels[i],
    #                 y=values[i]+.0125, 
    #                 s=round(values[i],2), 
    #                 size=9,
    #                 ha='center')

    map_legend_text_properties = {'family':'sans', 'weight':'semibold', 'style':'normal', 'size':10}
    factor=2

    geo_ax.set_title(title_lbl, pad=20, loc='center', fontdict={'family':'sans', 'color':'black', 'weight':'bold', 'style':'normal', 'size':12})
    map = Basemap(llcrnrlon=-180, urcrnrlon=180,
                llcrnrlat=-90, urcrnrlat=90, 
                projection='robin', 
                lon_0=0, 
                resolution='c',
                ax=geo_ax)
        
    map.drawmapboundary(fill_color='white', linewidth=0)
    map.fillcontinents(color='grey', alpha=0.3)
    map.drawcoastlines(linewidth=0, color="white")

    data = geo_GEM[['ecosystem_category', 'longitude', 'latitude', encntrs_cnt_clm]]
    data = data[data[encntrs_cnt_clm] > 0].groupby(['ecosystem_category', 'longitude', 'latitude'], as_index=False)[encntrs_cnt_clm].sum()
    data['ecosys_coordinate_percentage'] = (data[encntrs_cnt_clm] / data[encntrs_cnt_clm].sum()).mul(100)
    data['ecosystem_category'] = data.apply(lambda df: 'Others' if df.ecosystem_category in less_than_thresh_ecosystems else df.ecosystem_category, axis=1) # Later, 'Others ( < 1% )' should simply be changed to 'Others' 
    data['color'] = data.ecosystem_category.apply(lambda ecosystem:color_dict.get(ecosystem))
    data.sort_values(by=['longitude', 'latitude', 'ecosys_coordinate_percentage'], ascending=[True, False, False], inplace=True)

    # legend
    for size, label in zip([100, 25, 5], ['100', '25', '5']):
        geo_ax.scatter([], [], c='grey', alpha=0.5, s=size*factor, label=label)
    geo_ax.legend(scatterpoints=1, 
                  scatteryoffsets=[.5], 
                  frameon=True, 
                  labelspacing=2.2, 
                  loc='upper left', 
                  ncol=1, 
                  prop=map_legend_text_properties, 
                  bbox_to_anchor=(0, 1.02), 
                  shadow=True, 
                  markerfirst=True)
    #
    map.scatter(latlon=True,
                x=data.longitude, 
                y=data.latitude, 
                s=data[encntrs_cnt_clm].mul(factor),
                alpha=.8, 
                c=data.color)

    plt.tight_layout()

    return fig

# --------------------------------------------------------------------------------------------------------------------------------------------

for title_lbl, encntrs_cnt_clm, appearance_cnt_clm in zip(['T3SS',            'T4SS',            'T6SS'], 
                                                          ['T3SS_encounters', 'T4SS_encounters', 'T6SS_encounters'], 
                                                          ['T3SS_appearance', 'T4SS_appearance', 'T6SS_appearance']):

    f = dual_plot(GEM, title_lbl, encntrs_cnt_clm, ratio(appearance_cnt_clm, GEM))
    f.savefig(f'{title_lbl}_figure6.svg', bbox_inches='tight')
