#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP_haplotyping
Created on Mon Dec 21 12:21:06 2020

@author: Zac Scurlock
"""
import pandas as pd
import numpy as np
import sys
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from pairsnp import calculate_snp_matrix, calculate_distance_matrix


input_path=sys.argv[1]
output_path=sys.argv[2]
metadata_path=sys.argv[3]
outbreak_only=sys.argv[4]
outbreaks=sys.argv[5]
n_threshold=sys.argv[6]
remove_n=sys.argv[7]

def sample_indexes(fasta_object):
    array = fasta_object['fasta'].str.find('>')
    indexes=[x for x in range(len(array)) if array[x] == 0]
    
    return indexes

def reformat(fasta_data, index, IDs):
    sample_dict = {'Samples':[]}
    for y in range(len(index)):
        if y == (len(index) -1):
            real = fasta_data.iloc[index[y]+1: len(fasta_data)]
        else:
            real = fasta_data.iloc[index[y]+1: index[y+1]]
        real_1 = real['fasta'].str.cat(sep='\n')
        real_1 = real_1.replace('\n',"")
        sample_dict['Samples'].append({str(IDs[y]):real_1})
   
    
    return sample_dict

def sample_n_filter(data, IDs, N_threshold):
    new_dict = {}
    app = []
    for num, x in enumerate(IDs):
        n_count=[data['Samples'][num][str(x)].lower().count('n')]
        if n_count[0] <= int(N_threshold):
            app.append(num)
            #del sample_dict['Samples'][x], IDs[x]
            new_IDs = [IDs[i] for i in app]
            new_dict['Samples'] = [data['Samples'][i] for i in app]
        else:
            print('Excluding ' + x + ' for ambigious sequences')
    return new_dict, new_IDs


def filtering(sample_data, sample_id):
    test_dict={}   
    for length in range(len(sample_data['Samples'])):
        test_dict.setdefault(len(sample_data['Samples'][length][sample_id[length]]),[]).append(sample_id[length])

    lis = list(map(len, test_dict.values()))
    if len(lis) > 1:
        num = [i for i, x in enumerate(lis) if x !=max(lis)]
        for li in range(len(num)):
            for p in range(len(test_dict[list(test_dict.keys())[num[li]]])):
                pos = [i for i, x in enumerate(sample_id) if x == test_dict[list(test_dict.keys())[num[li]]][p]]
                print("Excluding sample ID from analysis: " + sample_id[pos[0]] + " due to differences in sample length")
                del sample_data['Samples'][pos[0]], sample_id[pos[0]]

    return sample_data, sample_id

def new_SNP(sample_dict, IDs, remove_n):
    bases = ['A', 'T', 'C', 'G',
             'a', 't', 'c', 'g']
    first = []
    genome = pd.DataFrame({'genome':IDs})
    
    for x in range(len(sample_dict['Samples'][0][IDs[0]])):
        for y in range(len(sample_dict['Samples'])):
            first.append(sample_dict['Samples'][y][IDs[y]][x])
        a = pd.Series(data=first)
        if not all([first in bases for first in first]) and remove_n:
            a=pd.Series(data=None)
        if len(a.value_counts()) >= 2:
            genome[str(x+1)] = first
            
        first=[]
    return genome


def append_haplotype(snp_df):
    bases = []
    frame = snp_df.drop('genome', axis=1)
    for per in range(len(frame)):
        haplo =""
        for snp in frame.iloc[per]:
            haplo += snp
        bases.append(haplo)
    snp_df.insert(1, 'haplotype', bases) 
    
    return snp_df

def haplotype_number(data):
    haplo = data['haplotype']
    set_hap = set(haplo)
    dict_hap = {}
    for y,  x in enumerate(set_hap):
        dict_hap.setdefault(x,[]).append('h.' + str(y+1))

    order =[]
    for y in range(len(haplo)):
        for x in dict_hap:
            if x in haplo[y]:
                order.append(dict_hap[x][0])
        
    data.insert(2, 'haplotype_number', order)
    
    return data
    
    

#Nexus file
def position(outbreak_data, test_dict, set_li, x):
    for counter, i in enumerate(set_li):
        if i == test_dict[outbreak_data['sample'].iloc[x]][0]:
            return counter

def create_array(outbreak_data, filtered_ids, traitlabels):
    sample_ids = list(outbreak_data['sample'])
    one_array = pd.DataFrame(data=None)
    trait_list=[]
    array_list=[]
    for sample in filtered_ids:
        if sample in sample_ids:
            ind = sample_ids.index(sample)
            test_dict={outbreak_data['sample'].iloc[ind]:
            [str(outbreak_data['number'].iloc[ind]) + 
             str(outbreak_data['p/s'].iloc[ind])]}
            trait_list.append(test_dict)

            one = traitlabels.index(list(test_dict.values())[0][0])
            arr = np.zeros( (1, len(traitlabels)+2), dtype=int)[0]
            arr[one] = 1
            array_list.append(np.array2string(arr,separator=',')[1:-1])
        else:
            trait_list.append({str(sample):['None']})
            
            arr = np.zeros( (1, len(traitlabels)+2), dtype=int)[0]
            arr[-2] = 1
            array_list.append(np.array2string(arr,separator=',')[1:-1])
       
    one_array['genome']=filtered_ids
    one_array['array']=array_list
    
    return one_array, trait_list

def create_nexus(outbreak_data, output_path, filtered_ids, snps_hap, outbreak_sample):
    #Traits labels
    trait_labels = [(str(outbreak_data['number'].iloc[x]) + 
                     str(outbreak_data['p/s'].iloc[x])) 
                    for x in range(len(outbreak_data))]
    
    label_set = sorted(set(trait_labels))
    
    
    ###START####
    textlist = ['#NEXUS',
                '\n', 
                '\n', 
                'BEGIN TAXA;\n', 
                'DIMENSIONS NTAX=', str(len(filtered_ids)), ';\n\n',
                'TAXLABELS\n']
    
    textlist2 = [';\n\nEND;\n\n',
                 'BEGIN CHARACTERS;\n',
                 'DIMENSIONS NCHAR='+str(len(snps_hap['haplotype'].iloc[0])) +';''\n',
                 'FORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n']
 
    textlist3 = [';\nEND;\n\nBegin Traits;\n',
                 'Dimensions NTraits='+ str(len(label_set)+2) +';\n'
                 'Format labels=yes missing=? separator=Comma;\n',
                'TraitLabels ']
    print('Writing haplotypes to NEXUS file..')
    
    if outbreak_sample == False:
        test_file = open(str(output_path) +".nex", 'w')
    else:
        test_file = open(str(output_path) +".outbreak_nex", 'w')
    test_file.writelines(textlist)
    for element in filtered_ids:
        print(element, file = test_file)
            
    test_file.writelines(textlist2)
    #Haplotype creation
    for p in range(len(snps_hap[['genome', 'haplotype']])):
        test_file.writelines(snps_hap[['genome', 'haplotype']].iloc[p][0])
        test_file.write('\t')
        test_file.writelines(snps_hap[['genome', 'haplotype']].iloc[p][1])
        test_file.write('\n')
    
    #Traits creation
    test_file.writelines(textlist3)
    
    #Traits labels - Caveat - Letter has to be a single one
    for x in label_set:
        number=x[:-1]
        letter=x[-1]
        if letter == 'S':
            let = 'staff'
        else:
            let = 'patient'
        test_file.write('Outbreak_' + number + '_' + let +' ')
    test_file.write('Other Reference;\nMatrix\n')
    
    
    array, trait_list=create_array(outbreak_data, filtered_ids, label_set)
    for p in range(len(array)):
        test_file.writelines(array['genome'].iloc[p])
        test_file.write('\t')
        test_file.writelines(array['array'].iloc[p])
        test_file.write('\n')
    
    test_file.write(';\nEND;')
    test_file.close()
    print(str(output_path) +'.nex is complete')
    
    return trait_list

def pseudosequence(output_path,snps_hap, outbreak_sample):
    print('Producing pseudosequences..')
    if outbreak_sample == False:
        pseudo_fna = open(output_path +'.pseudo.fna', 'w')
    else:
        pseudo_fna = open(output_path +'.outbreak_pseudo.fna', 'w')
    
    for p in range(len(snps_hap[['genome', 'haplotype']])):
        pseudo_fna.write('>')
        pseudo_fna.writelines(snps_hap[['genome', 'haplotype']].iloc[p][0])
        pseudo_fna.write('\n')
        pseudo_fna.writelines(snps_hap[['genome', 'haplotype']].iloc[p][1])
        pseudo_fna.write('\n')
        
    pseudo_fna.close()
    print(str(output_path) + '.pseudo.fna is complete')

def outbreak(used_outbreak_data, filtered_ids, filtered_dict):
    if len(used_outbreak_data) >= 0:

        used_ids = list(used_outbreak_data['sample'])
        outbreak_dict={'Samples':[]}
        for pos in used_ids:
            num = filtered_ids.index(pos)
            outbreak_dict['Samples'].append(filtered_dict['Samples'][num])

    return outbreak_dict, used_ids

def make_graph(output_path, trait_list, outbreak_sample):
    if outbreak_sample == False:
        suffix = '.pseudo.fna'
    else:
        suffix = '.outbreak_pseudo.fna'
        
    seq_input = str(output_path + suffix)
    
    #Calculate initial distance matrix
    sparse_matrix, consensus, seq_names = calculate_snp_matrix(seq_input)
    dist_matrix = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)


#Determining any duplicates within the samples - same haplotypes
    lis=[]
    for y in range(len(dist_matrix)):
        a = [i for i, x in enumerate(dist_matrix) if (x == dist_matrix[y]).all()]
        if len(a) > 1:
            lis.append(a)
    identical_snps = list({tuple(i) for i in lis})


#Calculate the possible samples we could have
    possible_nodes = list(range(len(dist_matrix)))
    for y in range(len(identical_snps)):
        for x in range(1,len(identical_snps[y])):
            possible_nodes.remove(identical_snps[y][x])


#Variable node size to the node that has been concatenated into
    app=[]
    node_size = {node:1 for node in possible_nodes}
    for x in identical_snps:
        app.append(x[1:])
        node_size[x[0]] = len(x)
    
    node_size_list = []  
    for node in node_size.values():
        size = node*300
        node_size_list.append(size)
    
    app_list = [x for y in range(len(app)) for x in app[y]]

    all_in_one = np.delete(dist_matrix, app_list, axis=0)
    all_in_two = np.delete(all_in_one, app_list, axis=1)
    ait = nx.from_numpy_matrix(all_in_two)


    #Labelling
    pre_labels={num:seq_names[num] for num in range(len(seq_names))}
    #pre_labels = {0:'0x', 1:'1x', 2:'2x', 3:'3x', 4:'4x', 5:'5x', 6:'6x', 7:'7x', 8:'8x'} # Obtain from sample names
    
    label_dict={num:pre_labels[num] for num in possible_nodes}
    #networkx relabelling node
    relabel_dict={}
    
    for x in range(len(ait.nodes())):
        relabel_dict[list(ait.nodes())[x]]=list(label_dict.values())[x]    
    nx.relabel_nodes(ait, relabel_dict, False)
    
    #Edge labels
    edgelabels = nx.get_edge_attributes(nx.minimum_spanning_tree(ait), 'weight')
    rounded_edgelabels = {list(edgelabels.keys())[weight]:round(list(edgelabels.values())[weight]) for weight in range(len(edgelabels))}


    #Adding colour
    new_dict = {}
    for x in range(len(trait_list)):
        new_dict[list(trait_list[x].keys())[0]] = list(trait_list[x].values())[0][0]
        
    colours_categories = pd.Series(new_dict, index=list(new_dict.keys()))
    colours_cats=pd.Categorical(colours_categories)
    colour_dict = {colours_categories[x]:colours_cats.codes[x] for x in range(len(colours_categories))}
    
    #Node colour dictionary
    #traits = colours_categories.value_counts().index
    new_dict2 = {colours_categories.index[x]:colour_dict[colours_categories[x]] 
    for x in range(len(colours_categories))}
    
    ##
    con = [new_dict2[x] for x in ait.nodes()]
    codes = np.array(con)
    
    if len(trait_list) == 0:
        vmax = 0
    else:
        vmax=max(codes)
    
    #Colour mapping
    palette = plt.cm.Set2
    cNorm=colors.Normalize(vmin=0, vmax=vmax)
    scalarMap=cmx.ScalarMappable(norm=cNorm, cmap=palette)
    
    print('Producing minimum spanning tree...')
    #Adding a legend
    f = plt.figure(figsize=(8,8))
    ax = f.add_subplot(1,1,1)
    for label in colour_dict:
        ax.plot([0],[0],color=scalarMap.to_rgba(colour_dict[label]),label=label)
    #Drawing graph
    pos=nx.spring_layout(nx.minimum_spanning_tree(ait), scale=2)
    nx.draw(nx.minimum_spanning_tree(ait),pos,with_labels=True,vmin=0, vmax=vmax,node_size=node_size_list, node_color=codes, cmap=palette, font_size=8, ax=ax)
    nx.draw_networkx_edge_labels(nx.minimum_spanning_tree(ait),pos, edge_labels=rounded_edgelabels, font_size=9, font_color='red')
    plt.legend()
    f.tight_layout()
    print('Minimum spanning tree is complete.')
    plt.savefig(str(seq_input + '.png'))

def outbreak_loop(outbreak_data, outbreaks):
    if outbreaks == '-a':
        new_meta=outbreak_data
        print('Using all outbreak samples..')
    elif outbreaks == 'None':
        new_meta = pd.DataFrame(data=None, 
                                columns=['sample', 'number', 'p/s'])
        print('Using no outbreak samples')
    else:
        print('Using ' + outbreaks + ' outbreak samples')
        outbreaks=[int(s) for s in outbreaks.split(',')]
        empty=[]
        pos=[] 
        for num in outbreaks:
            empty.append([i for i, x in enumerate(outbreak_data['number']) if x ==num])
    
        for num in range(len(empty)):
            pos +=empty[num]
        new_meta = outbreak_data.iloc[pos]
    print(new_meta)
    return new_meta

def calculate_outbreak(outbreak_data, outbreaks, output_path, filtered_ids, filtered_dict, remove_n):
    outbreak_dict, outbreak_ids = outbreak(outbreak_data, filtered_ids, filtered_dict)
    if len(outbreak_ids) != 0:
        print("Calculating Outbreak SNPs...")
        outbreak_SNP = new_SNP(outbreak_dict, outbreak_ids, remove_n)
        #print(outbreak_SNP)
        outbreak_hap = append_haplotype(outbreak_SNP)
        #print(outbreak_hap)
        outbreak_haplo = haplotype_number(outbreak_hap)
        trait_list = create_nexus(outbreak_data, output_path, outbreak_ids, outbreak_haplo, True)
        pseudosequence(output_path, outbreak_haplo, True)
        make_graph(output_path, trait_list, True)
        
        csv_id = [str(ids) for ids in filtered_ids 
        if ids in outbreak_ids]
        
        csv_id = [str(ids) + '_' + 'Outbreak_' + 
          str(outbreak_data['number'].iloc[outbreak_ids.index(ids)])
          + str(outbreak_data['p/s'].iloc[outbreak_ids.index(ids)]) for ids in filtered_ids 
          if ids in outbreak_ids]
        outbreak_SNP['genome']=csv_id
        outbreak_haplo.to_csv(str(output_path) + 'outbreak_only.csv', index=False)

def pipeline(input_path, output_path, metadata_path, outbreak_only, outbreaks, N_threshold, remove_n):
    data = pd.read_csv(str(input_path), 
    delimiter='\n', header=None, names=['fasta'])
    print('Fasta import is complete')
    indexes = sample_indexes(data)
    IDs = [data.iloc[x][0].split('/')[0] for x in indexes]
    IDs=[x[1:] for x in IDs]
    sample_dict = reformat(data, indexes, IDs)
    sample_dict, IDs = sample_n_filter(sample_dict,IDs, N_threshold)

    filtered_dict, filtered_ids = filtering(sample_dict, IDs)
    print('Sample filtering is complete')
    outbreak_data = pd.read_csv(str(metadata_path), 
                                delimiter='\t', header=None, 
                                names=['sample', 'number', 'p/s'], index_col=None)
    print('Outbreak_data is loaded')
    
    outbreak_ids = list(outbreak_data['sample'])
    red = outbreak_data['sample'].value_counts()[outbreak_data['sample'].value_counts() >1]
    pos = []
    noted_set=[]
    for num, x in enumerate(outbreak_ids):
        if x in list(red.index):
            noted_set.append(x)
            if noted_set.count(x) >= 2:
                pos.append(num)
                
    outbreak_data.drop(pos, inplace=True)
    outbreak_ids = list(outbreak_data['sample'])

    counts = [counter for counter, ids in enumerate(outbreak_ids) if ids in filtered_ids]
    used_outbreak_data = outbreak_data.iloc[counts]

    #Base ^  
    if outbreak_only == str(False):
        print("Calculating all sample SNPS...")
        snps = new_SNP(filtered_dict, filtered_ids, remove_n)
        snps_hap = append_haplotype(snps)
        snps_haplo = haplotype_number(snps_hap)
        #Nexus/Graph
        trait_list = create_nexus(used_outbreak_data, output_path, filtered_ids, snps_haplo, False)
        pseudosequence(output_path, snps_haplo, False)
        make_graph(output_path, trait_list, False)
        
        print('Written SNP file to ' + str(output_path))
        #Adding trait data to sample names in csv
        csv_id = [str(ids) + '_' + 'Outbreak_' + str(outbreak_data['number'].iloc[outbreak_ids.index(ids)])
        + str(outbreak_data['p/s'].iloc[outbreak_ids.index(ids)]) if ids in outbreak_ids else str(ids) for ids in IDs]
        snps['genome']=csv_id
        snps_haplo.to_csv(str(output_path) + '.csv', index=False)
        
        new_outbreak = outbreak_loop(used_outbreak_data, outbreaks)
        calculate_outbreak(new_outbreak,outbreaks,output_path, filtered_ids, filtered_dict, remove_n)
   
    elif outbreak_only == str(True):
        print("Calculating outbreak sample SNPS")
        new_outbreak = outbreak_loop(used_outbreak_data, outbreaks)    
        calculate_outbreak(new_outbreak,outbreaks,output_path, filtered_ids, filtered_dict, remove_n)


pipeline(input_path, output_path, metadata_path, outbreak_only, outbreaks,N_threshold, remove_n)
