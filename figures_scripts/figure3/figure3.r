library(treeio)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
library(svglite)
library(castor)
#
tree <- read.tree('C:/Users/DuduB21/OneDrive/code/MSC_thesis/general/trees/T4SS.nw')     

tree$tip.label <- sub('._Incertae_Sedis', '', tree$tip.label)
tree$tip.label <- sub('_Family_', '_', tree$tip.label)
tree$tip.label <- sub('_subdivision_3', '_sub._3', tree$tip.label)

end_of_line <- c((1:88), (97:105),(108:111), (112:164), (254:285), (365:373))                                                                             
start_of_line <- setdiff(c(1:391), end_of_line)

old_labels <- tree$tip.label
new_labels <- tree$tip.label

max_label = max(nchar(tree$tip.label))
new_labels[end_of_line] <- unlist(lapply(X=end_of_line, FUN=function(x) str_pad(string=new_labels[x], width=max_label, side='left', pad='-')))            
new_labels[start_of_line] <- unlist(lapply(X=start_of_line, FUN=function(x) str_pad(string=new_labels[x], width=max_label, side='right', pad='-')))       
new_labels[354] <- str_pad(substring(new_labels[354], 1, max_label-8), width=max_label, side='right', pad=' ')                                            
new_labels[189] <- str_pad(substring(new_labels[189], 1, max_label-7 ), width=max_label, side='right', pad=' ')                                           
new_labels[190] <- str_pad(substring(new_labels[190], 1, max_label-12), width=max_label, side='right', pad=' ')                                             
new_labels[end_of_line] <- sub('-(?=[a-zA-Z])', ' ', new_labels[end_of_line], perl=TRUE)                                                                    
new_labels[start_of_line] <- sub('(?<=[A-Za-z0-9])-', ' ', new_labels[start_of_line], perl=TRUE)                                                            

tree$tip.label <- new_labels    
tree$tip.label.empty <- str_replace(new_labels, '[^-]+', strrep(' ', nchar(new_labels)))                                                                                                                            
chang_labels_df <- data.frame(old_labels, new_labels)

rnk_phy_col_df <- read.csv('C:/Users/DuduB21/OneDrive/code/MSC_thesis/general/review_family_phylogenetic_mtdta.tsv', header=TRUE, sep='\t')               
rnk_phy_col_df$rnk_txn <- sub('._Incertae_Sedis', '', rnk_phy_col_df$rnk_txn)
rnk_phy_col_df$rnk_txn <- sub('_Family_', '_', rnk_phy_col_df$rnk_txn)
rnk_phy_col_df$rnk_txn <- sub('_subdivision_3', '_sub._3', rnk_phy_col_df$rnk_txn)                                                                          

for (family in c(rnk_phy_col_df$rnk_txn)) {
    if (family %in% old_labels) {
        rnk_phy_col_df[rnk_phy_col_df$rnk_txn==family, 'rnk_txn'] <- chang_labels_df[chang_labels_df$old_labels == family, 'new_labels']                  
    }
}

T3SS <- read.csv('C:/Users/DuduB21/OneDrive/code/MSC_thesis/all_together/T3SS.tsv', header=TRUE, sep='\t')                                                
T3SS$system <- factor(T3SS$system, levels=c('T3SS'))                                                                                                      
T3SS$possession <- factor(T3SS$possession, levels=c('T3SS','none'))
T3SS$rnk_txn <- sub('._Incertae_Sedis', '', T3SS$rnk_txn)
T3SS$rnk_txn <- sub('_Family_', '_', T3SS$rnk_txn)
T3SS$rnk_txn <- sub('_subdivision_3', '_sub._3', T3SS$rnk_txn)

for (family in c(T3SS$rnk_txn)) {
    if (family %in% old_labels) {
        T3SS[T3SS$rnk_txn==family,'rnk_txn'] <- chang_labels_df[chang_labels_df$old_labels == family, 'new_labels']                                       
    }
}

T4SS <- read.csv('C:/Users/DuduB21/OneDrive/code/MSC_thesis/all_together/T4SS.tsv', header=TRUE, sep='\t')                                                
T4SS$system <- factor(T4SS$system, levels=c('T4SSA','T4SSB'))                                                                                             
T4SS$possession <- factor(T4SS$possession, levels=c('T4SSA','T4SSB','none'))                                                                              
T4SS$rnk_txn <- sub('._Incertae_Sedis', '', T4SS$rnk_txn)                                                                                                 
T4SS$rnk_txn <- sub('_Family_', '_', T4SS$rnk_txn)                                                                                                        
T4SS$rnk_txn <- sub('_subdivision_3', '_sub._3', T4SS$rnk_txn)                                                                                                    

for (family in c(T4SS$rnk_txn)) {
    if (family %in% old_labels) {
        T4SS[T4SS$rnk_txn==family,'rnk_txn'] <- chang_labels_df[chang_labels_df$old_labels == family, 'new_labels']                                       
    }
}

T6SS <- read.csv('C:/Users/DuduB21/OneDrive/code/MSC_thesis/all_together/T6SS.tsv', header=TRUE, sep='\t')                                                
T6SS$system <- factor(T6SS$system, levels=c('T6SSi','T6SSii','T6SSiii'))                                                                                  
T6SS$possession <- factor(T6SS$possession, levels=c('T6SSi','T6SSii','T6SSiii','none'))                                                                           
T6SS$rnk_txn <- sub('._Incertae_Sedis', '', T6SS$rnk_txn)                                                                                                     
T6SS$rnk_txn <- sub('_Family_', '_', T6SS$rnk_txn)                                                                                                            
T6SS$rnk_txn <- sub('_subdivision_3', '_sub._3', T6SS$rnk_txn)                                                                                                

for (family in c(T6SS$rnk_txn)) {
    if (family %in% old_labels) {
        T6SS[T6SS$rnk_txn==family,'rnk_txn'] <- chang_labels_df[chang_labels_df$old_labels == family, 'new_labels']                                       
    }
}                                                                                                                                                         

# figure making :
fig <- ggtree(tree, layout='circular', size=.15)

# twitching the tree.
fig <- rotate(fig, 564)
fig <- rotate(fig, 649)
fig <- rotate(fig, 492)
fig <- rotate(fig, 493)
fig <- rotate(fig, 666)
fig <- rotate(fig, 462)
fig <- rotate(fig, 394)
fig <- rotate(fig, 472)
fig <- rotate(fig, 651)
fig <- rotate(fig, 652)
fig <- rotate(fig, 473)

fig <- fig %<+% rnk_phy_col_df

fig <- fig + geom_tippoint(aes(colour=phyla_class), 
                           alpha=0.75, 
                           size=.5, 
                           show.legend=TRUE) # Extended Data Figure 1 - TRUE, Figure 3 - FALSE.

fig <- fig + geom_tiplab(colour='#999999', # Figure 3 - Comment lines 109-118,  Extended Data Figure 1 - keep them.
                         geom='text',
                         fontface=3, 
                         size=.7,
                         offset=-.12, 
                         align=TRUE, 
                         linesize=.05, 
                         linetype='longdash',
                         show.legend=FALSE, 
                         family='mono')

fig <- fig + geom_fruit(data=T3SS, 
                        geom=geom_tile,
                        mapping=aes(y=rnk_txn, fill=possession),
                        width=0.05, 
                        color='white',
                        pwidth=0.07,
                        offset=0.16)

fig <- fig + geom_fruit(data=T4SS,
                        geom=geom_tile,
                        mapping=aes(x=system, y=rnk_txn, fill=possession),
                        width=0.05,
                        color='white',
                        pwidth=0.035,
                        offset=0.022)

fig <- fig + geom_fruit(data=T6SS, 
                        geom=geom_tile,
                        mapping=aes(x=system, y=rnk_txn, fill=possession),
                        width=0.05,
                        color='white',
                        pwidth=0.07,
                        offset=0.022)

fig <- fig + scale_fill_manual(name='System Possession',
             values=c('#E6E6E6', '#3C3C3C','#1b9e77', '#ff7f00', '#b2182b', '#ba36e3', '#1f78b4'), 
             na.translate=TRUE,
             guide=guide_legend(keywidth=1.5, keyheight=1.1, order=1, label.theme=element_text(face='italic'))) 

fig + theme(legend.position = 'none')
ggsave('G:/Shared drives/Lab.Yatir/paper/figures/figure3_and_extended_data_figure1.pdf', width=7, height=7) 