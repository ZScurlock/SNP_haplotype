#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:28:59 2020

@author: Z Scurlock
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

class cmd():
    def __init__(self, input_path, output_path, metadata_path, outbreak_only, outbreaks, n_threshold, remove_n):
        self.input_path = input_path
        self.output_path = output_path
        self.metadata_path = metadata_path
        self.outbreak_only = outbreak_only
        self.outbreaks = outbreaks
        self.n_threshold = n_threshold
        self.remove_n = remove_n
        self.fasta = ''
        self.indexes = []
        self.ids = []
        
    def import_data(self, input_path):
        data = pd.read_csv(str(input_path), 
        delimiter='\n', header=None, names=['fasta'])
        self.fasta=data
        return data
    
    def find_indexes(self, fasta_file):
        array = fasta_file['fasta'].str.find('>')
        indexes=[x for x in range(len(array)) if array[x] == 0]
        self.indexes = indexes
        
    def reformat(self, fasta_data, index, IDs):
        sample_dict = {'Samples':[]}
        for y in range(len(index)):
            if y == (len(index) -1):
                real = fasta_data.iloc[index[y]+1: len(fasta_data)]
            else:
                real = fasta_data.iloc[index[y]+1: index[y+1]]
            real_1 = real['fasta'].str.cat(sep='\n')
            real_1 = real_1.replace('\n',"")
            sample_dict['Samples'].append({str(IDs[y]):real_1})
        self.fasta = sample_dict
    
    def sample_n_filter(self, data, ids, n_threshold):
        new_dict = {}
        app = []
        for num, x in enumerate(ids):
            n_count=[data['Samples'][num][str(x)].lower().count('n')]
            if n_count[0] <= int(n_threshold):
                app.append(num)
                new_ids = [ids[i] for i in app]
                new_dict['Samples'] = [data['Samples'][i] for i in app]
            else:
                print('Excluding ' + x + ' for ambigious sequences')
        self.fasta = new_dict
        self.ids = new_ids
        
    def filtering(self, sample_data, sample_id):
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
    
        self.fasta = sample_data
        self.ids = sample_id

    def find_snps(self, samples, ids, remove_n):
            
        genome = pd.DataFrame({'genome':ids})
        for x in range(len(samples[0].fasta)):
            base=[samples[y].fasta[x] for y in range(len(samples))]
            summ = pd.Series(base)
            if not all([base in samples[0].bases for base in base]) and remove_n:
                summ=pd.Series(data=None)
            if len(summ.value_counts()) >= 2:
                genome[str(x+1)] = base
        
        self.snps = genome
            

    def append_haplotype(self, snps):
        bases = []
        frame = snps.drop('genome', axis=1)
        for per in range(len(frame)):
            haplo =""
            for snp in frame.iloc[per]:
                haplo += snp
            bases.append(haplo)
        snps.insert(1, 'haplotype', bases) 
    
        self.snps = snps

    def haplotype_number(self, data):
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
        
        self.snps = data
        
    def write_csv(self, output_path, outbreak_data, ids, outbreak_ids, snps, outbreak_only):
        if outbreak_only == True:
            csv_id = [str(x) + '_' + 'Outbreak_' + 
          str(outbreak_data['number'].iloc[outbreak_ids.index(x)])
          + str(outbreak_data['p/s'].iloc[outbreak_ids.index(x)]) for x in ids
          if x in outbreak_ids]
            snps_out = pd.DataFrame.copy(snps, deep=True)
            snps_out['genome']=csv_id
            snps_out.to_csv(str(output_path) + 'outbreak_only.csv', index=False)
            print('Written SNP file to ' + str(output_path) + 'outbreak_only.csv')
        else:
            snps.to_csv(str(output_path) + '.csv', index=False)
            print('Written SNP file to ' + str(output_path))
            
#Sample class
class sample():
    def __init__(self, ids):
        self.ids = ids
        self.fasta = ""
        self.n_threshold = ""
        self.bases=['A','T','C','G','a','t','c','g']
        

class outbreak(cmd):
    def __init__(self, input_path, output_path, metadata_path, outbreak_only, outbreaks,n_threshold, remove_n):
        super().__init__(input_path, output_path, metadata_path, outbreak_only,outbreaks, n_threshold,remove_n)
        self.ids = ""
        self.data=""
        self.new_meta=""
        self.used_data=""

    def import_data(self, metadata_path):
        outbreak_data = pd.read_csv(str(metadata_path), 
                                delimiter='\t', header=None, 
                                names=['sample', 'number', 'p/s'], index_col=None)
        print('Outbreak_data is loaded')
        
        self.data = outbreak_data
        self.ids = outbreak_data['sample']
        
    def remove_duplicates(self, data, ids):
        dup = ids.value_counts()[ids.value_counts()>1]
        pos = []
        noted_set=[]
        for num, x in enumerate(ids):
            if x in list(dup.index):
                noted_set.append(x)
                if noted_set.count(x) >= 2:
                    pos.append(num)
        print('Outbreak sample duplicates removed')
        data.drop(pos, inplace=True)
        self.data = data
        self.ids = list(data['sample'])
        
    #Only using metadata there are samples for
    def meta_outbreaks_used(self, data, ids, filtered_ids):
        counts = [counter for counter, x in enumerate(ids) if x in filtered_ids]
        self.used_data = data.iloc[counts]
        
    #Selecting metadata based on outbreak argument
    def outbreak_process(self, outbreak_data, outbreaks):
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
        self.used_data = new_meta
        self.used_ids = list(new_meta['sample'])

        

class results():
    def __init__(self, output_path):
        self.output_path = output_path
        self.trait_list = []
    
    def traits(self, outbreak_data):
        t_list = [(str(outbreak_data['number'].iloc[x]) + 
                     str(outbreak_data['p/s'].iloc[x])) 
                    for x in range(len(outbreak_data))]
        self.trait_list = sorted(set(t_list))
    
    def create_array(self, outbreak_data, filtered_ids, traitlabels):
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

        self.array = one_array
        self.labels = trait_list
    
    def create_nexus(self, outbreak_data, output_path, filtered_ids, snps_hap, 
                     outbreak_only, trait_labels, array):
        #Traits labels
        
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
        
        if outbreak_only == False:
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
        
        for p in range(len(array)):
            test_file.writelines(array['genome'].iloc[p])
            test_file.write('\t')
            test_file.writelines(array['array'].iloc[p])
            test_file.write('\n')
        
        test_file.write(';\nEND;')
        test_file.close()
        print(str(output_path) +'.nex is complete')
        
    def pseudosequence(self, output_path, snps_hap, outbreak_only):
        print('Producing pseudosequences..')
        if outbreak_only == False:
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
    
    
    def make_graph(self, output_path, trait_list, outbreak_only):
        if outbreak_only == False:
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
    
#Start#

test = cmd(input_path, output_path, metadata_path,
           outbreak_only, outbreaks, n_threshold, remove_n)

test.import_data(test.input_path)
test.find_indexes(test.fasta)
test.ids = [test.fasta.iloc[x][0].split('/')[1] for x in test.indexes]
test.reformat(test.fasta, test.indexes, test.ids)
test.sample_n_filter(test.fasta, test.ids, test.n_threshold)
test.filtering(test.fasta, test.ids)


out = outbreak(input_path, output_path, metadata_path, outbreak_only, outbreaks,n_threshold, remove_n)
out.import_data(out.metadata_path)
out.remove_duplicates(out.data, out.ids)
out.meta_outbreaks_used(out.data, out.ids, test.ids)


# End of pre-processing data #
if outbreak_only == str(False):
    samples = [sample(test.ids[i]) for i in range(len(test.ids))]
    
    if len(samples) == len(test.fasta['Samples']):
        for x in range(len(test.ids)):
            samples[x].fasta = test.fasta['Samples'][x][samples[x].ids]
    
    test.find_snps(samples, test.ids, test.remove_n)
    test.append_haplotype(test.snps)
    test.haplotype_number(test.snps)
    
    test.write_csv(test.output_path, out.data, test.ids, out.ids, test.snps, test.outbreak_only)

    ne = results(output_path)
    ne.traits(out.data)
    ne.create_array(out.data, test.ids, ne.trait_list)
    ne.create_nexus(out.data, ne.output_path, test.ids, test.snps,False,
                    ne.trait_list, ne.array)
    ne.pseudosequence(ne.output_path, test.snps, False)
    ne.make_graph(ne.output_path, ne.labels, False)
    
    out.outbreak_process(out.used_data, out.outbreaks)
    outbreak_samp = [sample(out.used_ids[i]) for i in range(len(out.used_ids))]
    
    for x in range(len(out.used_ids)):
        pos = test.ids.index(outbreak_samp[x].ids)
        outbreak_samp[x].fasta = test.fasta['Samples'][pos][outbreak_samp[x].ids]
        
    if len(outbreak_samp)> 0:
        out.find_snps(outbreak_samp, out.used_ids, out.remove_n)
        out.append_haplotype(out.snps)
        out.haplotype_number(out.snps)
        out.write_csv(out.output_path, out.used_data, test.ids, out.used_ids, out.snps, True)

        me = results(output_path)
        me.traits(out.used_data)
        me.create_array(out.used_data, out.used_ids, me.trait_list)
        me.create_nexus(out.used_data, me.output_path, out.used_ids,out.snps, True,
                        me.trait_list, me.array)
        me.pseudosequence(me.output_path, out.snps, True)
        me.make_graph(me.output_path, me.labels, True)
    else:
        print(str(len(outbreak_samp)) + ' outbreak samples have been selected to run. Try changing the outbreak_only or outbreaks arguments')

elif outbreak_only == str(True):
    out.outbreak_process(out.used_data, out.outbreaks)
    outbreak_samp = [sample(out.used_ids[i]) for i in range(len(out.used_ids))]
    
    for x in range(len(out.used_ids)):
        pos = test.ids.index(outbreak_samp[x].ids)
        outbreak_samp[x].fasta = test.fasta['Samples'][pos][outbreak_samp[x].ids]
        
    if len(outbreak_samp) > 0:
        out.find_snps(outbreak_samp, out.used_ids, out.remove_n)
        out.append_haplotype(out.snps)
        out.haplotype_number(out.snps)
        out.write_csv(out.output_path, out.used_data, test.ids, out.used_ids, out.snps, True)

        me = results(output_path)
        me.traits(out.used_data)
        me.create_array(out.used_data, out.used_ids, me.trait_list)
        me.create_nexus(out.used_data, me.output_path, out.used_ids,out.snps, True,
                        me.trait_list, me.array)
        me.pseudosequence(me.output_path, out.snps, True)
        me.make_graph(me.output_path, me.labels, True)
    else:
        print(str(len(outbreak_samp)) + ' outbreak samples have been selected to run. Try changing the outbreak_only or outbreaks arguments')
#End#

