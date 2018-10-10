# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 12:29:55 2018

@author: Angela
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 20:19:36 2018

@author: Angela
"""

import pandas as pd
import os
#import sys
import itertools
from matplotlib import pyplot as plt

os.chdir('C:\\Users\\Angela\\Google Drive\\healthy_microbiome\\updated_networks')
#sys.path.append(os.path.abspath("C:\\Users\\Angela-Thinkpad-X1\\Google Drive\\Code\\Chicken"))
#from used_fxns import taxonomy_trim, uniquify
#from used_fxns import within_clust

def taxonomy_trim(taxon_list,sep,level,num_level): 
    """fxn takes in taxon list typically provided in 16s analysis, separator between taxonomy levels, required text level to retrieve (e.g. 'g__' or 'D_5'), and the level designated numerically (e.g.genus would be level 6), and processes the taxon list such taht only the required level is shown"""           
    new_tax = []
    for taxon in taxon_list:
        current = taxon.split(sep)
        if (num_level >= len(current)):
            lastElement = current[len(current) - 1]
            new_tax.append("unclassified_" + lastElement[3:])
        else:
            lastElement = current[num_level]
            if ("uncultured" in lastElement):
                new_tax.append("unclassified_" + current[num_level-1][3:])
            elif ("Ambiguous" in lastElement):
                current = [a for a in current if 'Ambiguous' not in a]
                new_lastElement = current[len(current)-1]
                new_tax.append("unclassified_" + new_lastElement[3:])
            elif (lastElement[:3] == level): #if first 3 elements is equal to level (e.g. "g__" or "D_5")
                new_tax.append(lastElement[5:]) 
            elif (lastElement == "Unassigned"):
                new_tax.append("Unassigned")
            else:
                new_tax.append("unclassified_" + lastElement[3:])
    return(new_tax)

taxa_updated = pd.read_csv('all_taxa_map.txt',sep='\t') ## Read in OTU IDs and their corresponding taxa from SILVA database, map file is for network consisting of all 1572 samples 

# using taxonomy trim_function to get columns w/ just genus, family, and order
taxa_updated['genus'] = taxonomy_trim(taxa_updated['taxonomy'],';','D_5',5)
taxa_updated['family'] = taxonomy_trim(taxa_updated['taxonomy'],';','D_4',4)
taxa_updated['order'] = taxonomy_trim(taxa_updated['taxonomy'],';','D_3',3)

taxa_updated.to_csv("updated_taxa_map.txt",sep='\t') #this will be used for cytoscape network node attributes 

taxa_updated.set_index("#OTU ID",inplace = True)

# Node size is based on the number of samples an OTU is found in

updated_official = pd.read_csv('all_studies_updated100_sparcc.txt',header=0, sep='\t') # Read in OTU table that was submitted to sparcc for network generation 
updated_official.set_index("#OTU ID",inplace = True)

updated_official.index = updated_official.index.str.replace('.','')

updated_official_genus = pd.merge(updated_official, pd.DataFrame(taxa_updated['genus']), left_index = True, right_index = True)

node_size = updated_official.astype(bool).sum(axis=1) # compute the number of samples an OTU has been in
node_size.to_csv('updated_node_size.csv', header=None) # used for node size attribute in cytoscape 


def node_size_calc(otu_table_path):
    """reads in an otu table and computes node size attribute for each otu"""
    otu_table = pd.read_csv(otu_table_path,header=0, sep='\t') 
    otu_table.set_index("#OTU ID",inplace = True)
    
    otu_table.index = otu_table.index.str.replace('.','')
    
    node_size = otu_table.astype(bool).sum(axis=1)
    return(node_size)

## Interactions between clusters (to generate Figure 2b and c)

def uniquify(df_columns):
    """Certain OTUs are found in more than one clusters so will appear as duplicates in clusterONE (cytoscape plugin) outputted file, this fxn will take in a column of OTUs and number the duplicates, ie. if there are 3 Xs in column, it will generate new column numbering them X, X_2, X_3"""
    seen = set()

    for item in df_columns:
        fudge = 1
        newitem = item

        while newitem in seen:
            fudge += 1
            newitem = "{}_{}".format(item, fudge)

        yield newitem
        seen.add(newitem)

def node_map(clusterONE_output, sep):
    """Reads in clusterONE (cytoscape plugin) outputted text file and separator used in that file then returns a mapping file with two columns, first is every OTU found in clusters, second is the cluster that they're found in, third is a column with duplicate OTUs (ones that appear in multiple clusters) numbered"""
    node_attr = pd.read_csv(clusterONE_output,header=None)
    node_attr['node'] = node_attr.index
    node_attr.rename(columns={0:'OTUs'}, inplace=True)
    
    node_attr = pd.DataFrame(node_attr.OTUs.str.split(sep).tolist(), index=node_attr.node).stack()
    node_attr = node_attr.reset_index()[[0, 'node']] 
    node_attr.columns = ['OTUs', 'node'] 
    node_attr.node=node_attr.node + 1
    node_attr['unique'] = list(uniquify(node_attr.OTUs))
    return(node_attr)

#cluster_otus06 = node_map('updated_clusters_0.6overlap.txt', ' ') # 06 because the overlapping threshold used in clusterONE was 0.6
#
cluster_otus08 = node_map('updated_clusters_0.8overlap.txt', ' ') # reads in clusterONE outputted file 

cluster_node_sizes = pd.DataFrame(cluster_otus08.node.value_counts())
cluster_node_sizes.reset_index(inplace=True)
cluster_node_sizes.columns = ['cluster', 'size']
cluster_node_sizes.to_csv('cluster_node_sizes.txt', sep = '\t')

# For cytoscape clustered figure 

interact_data=pd.read_csv('updated_network.txt','\t') # end product of sparcc pipeline, contains every single interaction in network
interact_data.drop(['pval'],axis=1,inplace=True)
interact_data.drop(['abs'],axis=1,inplace=True)

def within_clust(node_attr,data_interact):
    """Includes only interactions within a cluster, this will be used to generate the clustered network where all the otus within the same cluster are grouped toegether"""
    clusters = node_attr['node'].unique()
    within_inter = pd.DataFrame([])
    for i in range(0,len(clusters)+1):
        cluster_num = node_attr.loc[node_attr['node'] == i]
        filt=interact_data.loc[(interact_data['OTU 1'].isin(cluster_num['OTUs'])) & (interact_data['OTU 2'].isin(cluster_num['OTUs'])) ]
        filt['OTU1 attr'] = i
        filt['OTU2 attr'] = i
        within_inter = within_inter.append(filt)   
    within_inter['type'] = 'intra'
    return(within_inter)


def btwClusters(words1,words2,data):
    """Computes interactions between two clusters by taking the mean of all the positive interactions and mean of all the negative interactions"""
# isolate all lines in file containing words
# forward interactions going from words1 to words2
    dataA=data.loc[data['OTU 1'].isin(words1)]
    dataAB=dataA.loc[dataA['OTU 2'].isin(words2)]
    
    dataB=data.loc[data['OTU 1'].isin(words2)]
    dataBA=dataB.loc[dataB['OTU 2'].isin(words1)]

    all = [dataAB['corr'], dataBA['corr']]
    result=pd.concat(all)
    pos=result[result>0]
    neg=result[result<0]
    return(list([pos.mean(),neg.mean()]))

def generate_all_files(clusterONEoutput, sep, num_clusters):
    """Generate all the files necessary to create cytoscape network figure in paper (Figure 2a,b,c) to current network, need to enter separator used in cluster file, and nubmer of clusters"""
    
    cluster_otus = node_map(clusterONEoutput, sep)

    cluster_node_sizes = pd.DataFrame(cluster_otus.node.value_counts()) # Get number of otus in each cluster 
    cluster_node_sizes.reset_index(inplace=True)
    cluster_node_sizes.columns = ['cluster', 'size']
    cluster_node_sizes.to_csv('cluster_node_sizes.txt', sep = '\t') # Will be used to determine the size of circles that show the actual cluster numbers (not the individual otus) in figure 2b

# Interactions within clusters
    within_interactions=within_clust(cluster_otus,interact_data)

    within2 = pd.merge(within_interactions, cluster_otus, how='left', left_on=['OTU 1','OTU1 attr'], right_on=['OTUs','node'])
    within2 = pd.merge(within2, cluster_otus,how='left', left_on=['OTU 2','OTU2 attr'], right_on=['OTUs','node'])
    within2.drop(['node_x','node_y'],axis=1,inplace=True)
    within2.drop(['OTUs_x','OTUs_y'],axis=1,inplace=True)
    within2.columns.values[6:8] = ['OTU1u','OTU2u']
    within2.to_csv('Clusters_within_interactions.txt','\t',index=False) #these interactions will be used to group clusters together then hidden to generate figure 2b

# Interactions between clusters 
    t1=open(clusterONEoutput)
    t2=open(clusterONEoutput)
    t1lines=t1.readlines()
    t2lines=t2.readlines()

    cluster_interact = pd.DataFrame([])
    for num1, line1 in enumerate(t1lines):
        for num2, line2 in enumerate(t2lines):
            if num1==num2:
                continue
            else:
                words1 = line1.split() 
                words2 = line2.split()
                output = btwClusters(words1,words2,interact_data)
                app_data=pd.concat([pd.DataFrame([num1+1, num2+1]).T, pd.DataFrame(output).T],axis=1,ignore_index=True)
                cluster_interact=cluster_interact.append(app_data)

    cluster_interact.columns = ['cluster1','cluster2','positive','negative']
    cluster_interact = cluster_interact.dropna(thresh=3) #drop rows where clusters have no interactions, only nan values
    cluster_interact = cluster_interact.fillna(0)
    cluster_interact['interactions'] = cluster_interact[['positive', 'negative']].sum(axis=1)

    btw_interactions = cluster_interact.drop(['positive','negative'],axis=1)
    btw_interactions['type'] = 'inter'
    btw_interactions.columns.values[0:3] = ['OTU 1', 'OTU 2', 'corr']

    btw_interactions.reset_index(inplace=True)
    btw_interactions.drop(btw_interactions.columns[0], axis=1, inplace=True)
    btw_reduced = btw_interactions[['OTU 1', 'OTU 2']].apply(sorted, axis = 1).drop_duplicates() # step removes duplicate interactions where OTUs in A and B are just reversed

    btw_interactions = pd.merge(btw_reduced, btw_interactions, how = 'left', left_on=['OTU 1', 'OTU 2'], right_on=['OTU 1', 'OTU 2'])

    btw_interactions.to_csv("btw_cluster_interactions.txt", "\t") ## Need to put btw_interactions and within2 together to generate clustered networks where clusters are grouped together by "intra" edges and then cluster names and between cluster interactions are manually dragged to those clusters

# For figure 2C, where the clusters are represented as pie graphs according to taxonomy breakdown 
    pie_otus = pd.DataFrame(columns=range(1,num_clusters+1),index=cluster_otus['OTUs'].unique())
    pie_otus=pie_otus.fillna(0)

    for i, row in pie_otus.iterrows():
        slice_node_attr = cluster_otus[['OTUs', 'node']][cluster_otus['OTUs']==i] 
        for node in slice_node_attr.node:
            row[node]=1
        
    pie_otus.to_csv('cluster_pie_attr.txt','\t')

    subset_taxa_map = taxa_updated.loc[taxa_updated.index.isin(pie_otus.index)]

## Generating the colouring legend (by genus) for clusters and original network   

    taxa_pie = pie_otus.join(taxa_updated[['genus']])

    taxa_pie = taxa_pie.groupby(taxa_pie.genus).sum()
    subset_taxa_map = subset_taxa_map.groupby(['genus','order']).size().reset_index()
    subset_taxa_map.drop(0,axis=1,inplace=True)
    subset_taxa_map.set_index('genus',inplace=True)
    subset_taxa_map['order'] = 'order_' + subset_taxa_map['order'].astype(str)
    
    
    taxa_pie = taxa_pie.join(subset_taxa_map)
    
    taxa_pie["taxa_final"] = ""
    subset_taxa_map["final"] = ""
    for ind,row in taxa_pie.iterrows():
        if taxa_pie.loc[ind,list(range(1,num_clusters + 1))].sum() < 15:
            taxa_pie.ix[ind,'taxa_final'] = taxa_pie.ix[ind,'order']
            subset_taxa_map.ix[ind,"final"] = taxa_pie.ix[ind,'order']
        else:
            taxa_pie.ix[ind,'taxa_final'] =ind
            subset_taxa_map.ix[ind,"final"] = ind
            
    taxa_pie = taxa_pie.groupby(taxa_pie.taxa_final).sum()
    
    
    taxa_pie = taxa_pie.T
    
    subset_taxa_map.reset_index(inplace=True)
    c2_subset_taxa_map = taxa_updated.loc[taxa_updated.index.isin(pie_otus.index)]
    c2_subset_taxa_map.drop('order',axis=1,inplace=True)
    c2_subset_taxa_map.reset_index(inplace=True)
    c2_subset_taxa_map = pd.merge(c2_subset_taxa_map, subset_taxa_map, how='left', left_on = 'genus', right_on='genus')
    c2_subset_taxa_map = pd.merge(cluster_otus, c2_subset_taxa_map, how='left', left_on = 'OTUs', right_on='#OTU ID')
    c2_subset_taxa_map.drop(['OTUs'],axis=1,inplace=True)
    c2_subset_taxa_map.to_csv('c2_taxa_map.txt',sep='\t') #nodes attribute file for figure 2b
    
    #compute edges (number of OTUs shared) for each node (for edges attribute for figure 1c)
    cols = pie_otus.columns
    bt = pie_otus.apply(lambda x: x > 0)
    nonzero_indices = bt.apply(lambda x: list(cols[x.values]), axis=1)
    
    edges_shared_otus=pd.DataFrame([])
    for otu in nonzero_indices:
        edges_shared_otus = edges_shared_otus.append(pd.DataFrame(list(itertools.combinations(otu,2))))
        
    edges_shared_otus.columns = ['cluster1','cluster2']
    edges_shared_otus.reset_index(inplace=True)
    edges_shared_otus.drop('index',axis=1,inplace=True)
    edges_shared_otus.columns = ['c1','c2']
    edges_shared_otus=edges_shared_otus.groupby(['c1','c2']).size().reset_index().rename(columns={0:'count'})
    
    taxa_pie.to_csv('taxa_pie.txt', sep='\t') #used to generate pie charts according to proportion of different taxa in each cluster in cytoscape (for figure 2C)
    edges_shared_otus.to_csv('cluster_edges.txt', sep='\t',index=False) #edge attribute for figure 2C

#os.chdir('C:\\Users\\Angela-Thinkpad-X1\\Google Drive\\healthy_microbiome\\updated_networks\\threshold07new')
#generate_all_files('clusterONE_threshold07.txt', ' ', 10)

os.chdir('C:\\Users\\Angela\\Google Drive\\healthy_microbiome\\updated_networks\\clusterONE_08_positive')
generate_all_files('threshold08_positive.txt', ' ', 13)

def classify_prop(df,node_otus,cluster_num):
    '''function creates input df for stacked bar chart where each cluster is evaluated on how prevalent it is in all samples, ie. each bar is a cluster, and the bar is split up into proprotional categories; a popular cluster must have 100% of its otus present in all samples, while a niche cluster might have a large proportion of the bar in the 0-20% category
    df must be a dataframe produced from 16s processing where columns are samples and rows are OTU abundances, node_otus must have an OTU column and a node column where the node dictates which cluster the OTU belogns to, cluster_num is the current cluster that the fucntion is trying to build a stacked bar for'''
    cluster_prop = pd.DataFrame([0,0,0,0,0],index=['0','0-25','25-50','50-75','75-100'])
    for column in df:
        sample = df[column]#.to_frame()
        sample = sample[sample>5].to_frame()
        curr_cluster = node_otus[node_otus['node']==cluster_num]
        sample_otus = pd.merge(sample, curr_cluster,how='inner',left_index=True,right_on='OTUs')
        num=sample_otus.shape[0]/curr_cluster.shape[0]
        if num == 0:
            cluster_prop.ix['0']=cluster_prop.ix['0']+1
        elif 0<num<0.25:
            cluster_prop.ix['0-25']=cluster_prop.ix['0-25']+1
        elif 0.25 <= num < 0.5:
            cluster_prop.ix['25-50']=cluster_prop.ix['25-50']+1
        elif 0.5 <= num < 0.75:
            cluster_prop.ix['50-75']=cluster_prop.ix['50-75']+1
        elif 0.75 <= num < 1:
            cluster_prop.ix['75-100']=cluster_prop.ix['75-100']+1
    return(cluster_prop)


stacked_clstr_pr = pd.DataFrame([],index=['0','0-25','25-50','50-75','75-100'])

node_attr = node_map('threshold08_positive.txt', ' ')

for x in range(1,14):
    stacked_clstr_pr = pd.concat([stacked_clstr_pr,classify_prop(updated_official,node_attr,x)], axis=1) 

stacked_clstr_pr.columns = range(1,14)

def normalize_raw(df):
    """normalize a raw otu table so that OTUs are represented by their proportions in their samples as opposed to raw counts"""
    for column in df:
        df[column] = df[column] / df[column].sum()
    return(df)

stacked_clstr_pr_norm = normalize_raw(stacked_clstr_pr)
stacked_clstr_pr_norm = stacked_clstr_pr*100
stacked_clstr_pr_norm.T.plot.bar(stacked=True).legend(loc='center left', bbox_to_anchor=(1, 0.5)) # generates figure S6


def bar_clusters(cluster_interact, loop_list):
    """loop list is list of nodes whose interactions with other nodes need to be counted (nodes could be individual clusters or otus), will be used to generate scatter plot in figure 2D
    cluster_interact must contain 3 columns: cluster1, cluster2, and interactions """
    clust_ratio = pd.DataFrame([])
    cluster_interact=cluster_interact.loc[cluster_interact['cluster1'] != cluster_interact['cluster2']]
    print(cluster_interact.shape[0])
    for node in loop_list:
            cluster = cluster_interact.loc[(cluster_interact['cluster1'] == node) | (cluster_interact['cluster2'] == node)]
            pos=cluster[cluster.interactions > 0]
            neg=cluster[cluster.interactions < 0]
            int_ratio = pd.DataFrame([node,pos.shape[0],neg.shape[0]]).T
            clust_ratio=clust_ratio.append(int_ratio)
    clust_ratio.columns=['cluster','positive','negative']
    return(clust_ratio)

otu_feed = interact_data.copy()


def plot_scatter(interact_data, taxa_updated):
    interact_data.columns=['cluster1','cluster2','interactions']
    all_OTUs=list(pd.concat([interact_data['cluster1'],interact_data['cluster2']]).unique())
    clust_ratio = bar_clusters(interact_data,all_OTUs)
    
    clust_ratio = clust_ratio[(clust_ratio != 0).all(1)]
    clust_ratio['ratio']=clust_ratio.negative/clust_ratio.positive
    clust_ratio['degree']=clust_ratio[['positive', 'negative']].sum(axis=1)
    clust_ratio = pd.merge(clust_ratio,taxa_updated, how='left', left_on="cluster", right_index=True)
    
    clust_ratio.reset_index(inplace=True, drop=True)
    return(clust_ratio)

updated_plot = plot_scatter(otu_feed,taxa_updated)

def plot_ratio_degree(clust_ratio, font): 
    """plots scatter plot to generate Figure 2D"""
    for i in range(len(clust_ratio.ratio)):
        x = clust_ratio.ratio[i]
        y = clust_ratio.degree[i]
        plt.plot(x, y, 'ro', markersize=5)
        plt.text(x, y, clust_ratio.genus[i],fontsize=font)
            
    plt.grid(True)

plot_ratio_degree(updated_plot,7)



os.chdir('C:\\Users\\Angela\\Google Drive\\healthy_microbiome\\updated_networks')

def read_taxa(taxa_file):
    """generates taxa mapping file from otu table with taxonomy column"""
    taxa_map = pd.read_csv(taxa_file,'\t')
    taxa_map.set_index('#OTU ID',inplace=True)
    taxa_map['genus'] = taxonomy_trim(taxa_map['taxonomy'],'; ','D_5',5)
    taxa_map.index = taxa_map.index.str.replace('.','')
    return(taxa_map)

def process_otu_tab(otu_tab_doc,taxa_map):
    """reads in otu table and processes it so it can be suitable for merging and scatter plot generation"""
    otu_tab = pd.read_csv(otu_tab_doc,'\t', index_col = 0)
    otu_tab.index = otu_tab.index.str.replace('.','')
    otu_tab = pd.merge(otu_tab, pd.DataFrame(taxa_map['genus']), left_index = True, right_index = True)
    return(otu_tab)


V13_FecCec_samples = pd.DataFrame(node_size_calc('mgp12337_otus.txt'))
V13_FecCec_samples = V13_FecCec_samples.rename(columns = {0:'num_samples'})
V13_FecCec_samples.to_csv('V13_FecCec_sample_num.txt','\t')  # Needed to generate cytoscape file in figure S8
FecCec_taxa = read_taxa('FecalCecal_taxa.txt')
FecCec_taxa['family'] = taxonomy_trim(FecCec_taxa['taxonomy'],'; ','D_4',4)
FecCec_taxa.to_csv("FecCec_taxa_map.txt",sep='\t') # Need this file to generate cytoscape file in figure S8
V13_FecCec_samples_num = pd.merge(V13_FecCec_samples, FecCec_taxa, left_index = True, right_index = True)
V13_FecCec_samples_core = V13_FecCec_samples_num[V13_FecCec_samples_num['num_samples'] > 16]
V13_FecCec_samples_core['genus'] = taxonomy_trim(V13_FecCec_samples_core['taxonomy'], '; ', 'D_5', 5)
V13_FecCec = pd.read_csv('V13_FecCec_network_025.txt','\t')
V13_FecCec.drop(['pval'],axis=1, inplace=True)
V13_FecCec_deg_ratio = plot_scatter(V13_FecCec, FecCec_taxa)
plot_ratio_degree(V13_FecCec_deg_ratio, 14) # Figure S9a


FecCec_otu_tab = process_otu_tab('mgp12337_otus.txt',FecCec_taxa)  # Need later for figure S7
FecCec_otu_tab.drop(['taxonomy'],axis=1,inplace=True)

V13_perform = pd.read_csv('V13_perform_network_025.txt','\t')

perform_taxa = read_taxa('perform_taxa_map.txt')

V13_perform.drop(['pval'],axis=1, inplace=True)
V13_perform_deg_ratio = plot_scatter(V13_perform, perform_taxa)
plot_ratio_degree(V13_perform_deg_ratio, 16) # Figure S9b
V13_perform_otu_tab = process_otu_tab('V13_perform_otus.txt', perform_taxa) #need it for figure S7

V13_perform_samples = pd.DataFrame(node_size_calc('V13_perform_otus.txt'))
V13_perform_samples = V13_perform_samples.rename(columns = {0:'num_samples'})
V13_perform_samples.to_csv('perform_sample_num.txt','\t') # Need this file to generate cytoscape file in figure S8

V68_AGP = pd.read_csv('V68_AGP_network_025.txt','\t')
V68_AGP_taxa = pd.read_csv('V68_AGP_taxa.txt','\t')
V68_AGP_taxa['genus'] = taxonomy_trim(V68_AGP_taxa['taxonomy'],'; ','D_5',5)
V68_AGP_taxa.set_index('#OTU ID', inplace=True)
V68_AGP_taxa['family'] = taxonomy_trim(V68_AGP_taxa['taxonomy'],'; ','D_4',4)
V68_AGP_taxa.to_csv("V68_taxa_map.txt",sep='\t')
V68_AGP.drop(['pval'],axis=1, inplace=True)
V68_AGP_deg_ratio = plot_scatter(V68_AGP, V68_AGP_taxa)
plot_ratio_degree(V68_AGP_deg_ratio,10) # Figure S9c

V68_otu_tab = process_otu_tab('AGP_26425940_otus.txt', V68_AGP_taxa) # Need later for figure S7
V68_samples = pd.DataFrame(node_size_calc('AGP_26425940_otus.txt'))
V68_samples = V68_samples.rename(columns = {0:'num_samples'})
V68_samples.to_csv('V68_sample_num.txt','\t') # Need this file to generate cytoscape file in figure S8

Salvax_taxa = pd.read_csv('Salvax_taxa.txt','\t')
Salvax_taxa['genus'] = taxonomy_trim(Salvax_taxa['taxonomy'],'; ','D_5',5)
Salvax_taxa.set_index('#OTUID', inplace=True)
Salvax_taxa['family'] = taxonomy_trim(Salvax_taxa['taxonomy'],'; ','D_4',4)
V4_Salvax_otu_tab = process_otu_tab('Salvax_otu_min100.txt', Salvax_taxa)  # Need later for figure S7

biodiv_taxa = pd.read_csv('biodiv_taxa.txt','\t')
V3_biodiv_otu_tab = process_otu_tab('biodiv_otu_min100.txt', biodiv_taxa)  # Need later for figure S7

def lacto_stacked_rel(otu_tab):
    """given an otu table, convert to relative abundance, and then loop through each sample to see count relative abundance of lactobacillus, returns a stacked bar chart showing percentage of samples that had lactobacillus greater than X %"""
    otu_tab = otu_tab.groupby(otu_tab.genus).sum()
    otu_rel = normalize_raw(otu_tab)
    lacto_rel = pd.DataFrame([0,0,0,0,0],index=['0-5','5-10','10-15','15-20','>20'])
    for column in otu_rel:
        num=float(otu_rel[column][otu_rel[column].index == 'Lactobacillus'])
        if  0<=num<0.05:
            lacto_rel.ix['0-5']=lacto_rel.ix['0-5']+1
        elif 0.05 <= num < 0.1:
            lacto_rel.ix['5-10']=lacto_rel.ix['5-10']+1
        elif 0.1 <= num < 0.15:
            lacto_rel.ix['10-15']=lacto_rel.ix['10-15']+1
        elif 0.15 <= num < 0.2:
            lacto_rel.ix['15-20']=lacto_rel.ix['15-20']+1
        elif 0.2 <= num :
            lacto_rel.ix['>20']=lacto_rel.ix['>20']+1
    return(lacto_rel)

lacto_rel_all = pd.DataFrame([],index=['0-5','5-10','10-15','15-20','>20'])

lacto_rel_all = pd.concat([lacto_rel_all, lacto_stacked_rel(V3_biodiv_otu_tab)], axis=1)
lacto_rel_all = pd.concat([lacto_rel_all,lacto_stacked_rel(V4_Salvax_otu_tab)], axis =1)
lacto_rel_all = pd.concat([lacto_rel_all,lacto_stacked_rel(V68_otu_tab)], axis =1)
lacto_rel_all = pd.concat([lacto_rel_all,lacto_stacked_rel(V13_perform_otu_tab)], axis=1)
lacto_rel_all = pd.concat([lacto_rel_all,lacto_stacked_rel(FecCec_otu_tab)], axis=1)

lacto_rel_all.columns = ['V3 27129897', 'V4 26835461', 'V6-V8 26425940', 'V1-V3 26925052', 'V1-V3 25887695']
lacto_rel_all.T.plot.bar(stacked=True).legend(loc='center left', bbox_to_anchor=(1, 0.5)) #figure S7a

def lacto_stacked_abs(otu_tab, threshold):
    """given an otu table, loop through each sample and count nubmer of lactobacillus OTUs with read counts over 100 (to prove that they are detecting lactobacillus), returns a stacked bar chart showing number of lactobacillus otus detected"""
# This hopefully proves that it is not because illumina miseq is detecting more species which lowers relative abundance of lactobacillus, it proves that other hypervariable regions are less sensitive to lactobacillus otus, thus making their contribution harder to infer
    otu_tab = otu_tab[otu_tab['genus'] == 'Lactobacillus']
    
    lacto_abs = pd.DataFrame([0,0,0,0,0],index=['0-5','5-10','10-15','15-20','>20'])
    for column in otu_tab.drop(['genus'],axis=1):
        detected=otu_tab[column][otu_tab[column]>threshold]
        num=detected.shape[0]
        if  0<=num<5:
            lacto_abs.ix['0-5']=lacto_abs.ix['0-5']+1
        elif 5 <= num < 10:
            lacto_abs.ix['5-10']=lacto_abs.ix['5-10']+1
        elif 10 <= num < 15:
            lacto_abs.ix['10-15']=lacto_abs.ix['10-15']+1
        elif 15 <= num < 20:
            lacto_abs.ix['15-20']=lacto_abs.ix['15-20']+1
        elif 20 <= num :
            lacto_abs.ix['>20']=lacto_abs.ix['>20']+1
    return(lacto_abs)

lacto_abs_all = pd.DataFrame([],index=['0-5','5-10','10-15','15-20','>20'])

lacto_abs_all = pd.concat([lacto_abs_all, lacto_stacked_abs(V3_biodiv_otu_tab,20)], axis=1)
lacto_abs_all = pd.concat([lacto_abs_all,lacto_stacked_abs(V4_Salvax_otu_tab,20)], axis =1)
lacto_abs_all = pd.concat([lacto_abs_all,lacto_stacked_abs(V68_otu_tab,1)], axis =1)
lacto_abs_all = pd.concat([lacto_abs_all,lacto_stacked_abs(V13_perform_otu_tab,1)], axis=1)
lacto_abs_all = pd.concat([lacto_abs_all,lacto_stacked_abs(FecCec_otu_tab,1)], axis=1)

lacto_abs_all.columns = ['V3 27129897', 'V4 26835461', 'V6-V8 26425940', 'V1-V3 26925052', 'V1-V3 25887695']

lacto_abs_all = normalize_raw(lacto_abs_all)

lacto_abs_all.T.plot.bar(stacked=True).legend(loc='center left', bbox_to_anchor=(1, 0.5)) #figure S7b


##For supplemental figure 1

os.chdir('C:\\Users\\Angela\\Google Drive\\healthy_microbiome\\updated_taxa_charts')
otu_tab_L6 = pd.read_csv('ChickenType_otu_table_sorted_L6_edited.txt', '\t', index_col = 0)
breed_corr = otu_tab_L6.corr('spearman')
breed_corr.to_csv('breed_correlations.txt', sep='\t')


# For supplemental table 3

def find_pop_comb(otu_table,node_otus,cluster_num):
    curr_cluster = node_otus[node_otus['node']==cluster_num]
    base_lim = int(curr_cluster.shape[0]/4)
    tup_otus=()
    for column in otu_table:
        sample = otu_table[column]#.to_frame()
        sample = sample[sample>5].to_frame()
        sample_otus = pd.merge(sample, curr_cluster,how='inner',left_index=True,right_on='OTUs')
        sample_otus = tuple(list(sample_otus['OTUs']))
        tup_otus = tup_otus + (sample_otus,)
    tup_otus = [x for x in tup_otus if x != ()]
    return(tup_otus)

os.chdir('C:\\Users\\Angela\\Google Drive\\healthy_microbiome\\cluster_frequencies_new')

clusters = list(range(1,14)) 
for x in clusters:
    tup_cluster = find_pop_comb(updated_official, node_attr, x)
    f = open('tup_cluster' + str(x) + '.txt', 'w')
    for t in tup_cluster:
        f.write(' '.join(str(s) for s in t) + '\n')
    f.close()
    
    
f = open('tup_cluster' + str(13) + '.txt', 'w')
for t in tup_cluster:
    f.write(' '.join(str(s) for s in t) + '\n') #the file will be then imported into an R script so the arules R package can be used to calculate the most popular OTU combinations
f.close()

