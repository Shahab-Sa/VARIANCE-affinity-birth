# Python script for analyses of "Antibody affinity birth through somatic hypermutation" publication.
# This pipeline is devided to 8 sections. At the beginning of each section there is a comment which indicates which figures of the publication are generated based on that section.


# input sequences for these analyses are uploaded in data folder. The result of each section will be saved in output folder.
print('Running...')
data_folder='../data'

import os
#import sys
import pandas as pd
import numpy as np

import time
import itertools
import matplotlib.pyplot as plt
import glob
import logomaker #https://logomaker.readthedocs.io

# Functions
def display_big():

    # df = pd.DataFrame()
    # pd.options.display.max_colwidth = 2000
    pd.set_option('display.max_rows', 20)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)

modes_dic={ 1 : {'mice' : 'B18-383', 'datasets': ['OVA', 'APC', 'CGG', 'Passenger']},
            2 : {'mice' : 'HA', 'datasets': ['0_1_OVA', '0_1_APC', '0_1_CGG', '1_1_mix']},

            3 : {'mice' : 'B18-383', 'datasets': ['OVA-CTLA4', 'OVA-Isotype']},
            4 : {'mice' : 'HA', 'datasets': ['0_1_OVA-CTLA4', '0_1_OVA-Isotype']},

            5 : {'mice' : 'HA', 'datasets': ['1_1000_CGG-CTLA4', '1_1000_CGG-Isotype']},
          }

del_sign='-'
aas_dic={'AAA':'K','AAC':'N','AAT':'N','AAG':'K','ACA':'T','ACC':'T','ACT':'T','ACG':'T','ATA':'I','ATC':'I',\
        'ATT':'I','ATG':'M','AGA':'R','AGC':'S','AGT':'S','AGG':'R','CAA':'Q','CAC':'H','CAT':'H','CAG':'Q',\
        'CCA':'P','CCC':'P','CCT':'P','CCG':'P','CTA':'L','CTC':'L','CTT':'L','CTG':'L','CGA':'R','CGC':'R',\
        'CGT':'R','CGG':'R','TAA':'*','TAC':'Y','TAT':'Y','TAG':'*','TCA':'S','TCC':'S','TCT':'S','TCG':'S',\
        'TTA':'L','TTC':'F','TTT':'F','TTG':'L','TGA':'*','TGC':'C','TGT':'C','TGG':'W','GAA':'E','GAC':'D',\
        'GAT':'D','GAG':'E','GCA':'A','GCC':'A','GCT':'A','GCG':'A','GTA':'V','GTC':'V','GTT':'V','GTG':'V',\
        'GGA':'G','GGC':'G','GGT':'G','GGG':'G','---':del_sign}
aas_list=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', del_sign]
aas_chemistry_list=['I', 'V', 'L', 'F', 'C', 'M', 'A', 'W', 'G', 'T', 'S', 'Y', 'P', 'H', 'N', 'D', 'Q', 'E', 'K', 'R']
nts_list=['A', 'C', 'G', 'T', '-']

def set_output_folder(section_output):
    output_folder=data_folder+'/output/'+section_output

    if not os.path.isdir(output_folder): # make output folder if it doesn't exist
        os.makedirs(output_folder)
    return(output_folder)

def translate(seq_nt):
    seq_aa=[]
    frameshift, stopcodon=False, False
    codons=[seq_nt[i:i+3] for i in range(0, len(seq_nt), 3)]
    for codon in codons:
        if codon in aas_dic:
            seq_aa.append(aas_dic[codon])
            if aas_dic[codon]=='*': stopcodon=True
        else:
            seq_aa.append('X')
            frameshift=True
    return([seq_aa, frameshift, stopcodon])

def trans_fra_stop(df):
    for i in df.index:
        translation=translate(df.loc[i, 'seq_nt'])
        df.loc[i, 'seq_aa']=''.join(translation[0])
        df.loc[i, 'frameshift']=translation[1]
        df.loc[i, 'stopcodon']=translation[2]
        df.loc[i, 'len_aa']=len(translation[0])
    return(df)

def import_fasta(files_List):
    dfs=pd.DataFrame()
    for file in files_List:
        # print('file:', file)
        file_name=file.split('/')[-1].rstrip('.fasta') # '0-1_APC-VH'
        # print('file_name:', file_name)
        mice=file.split('/')[-2]
        # print('mice:', mice)
        chain=file_name.split('-')[-1]
        # print('chain:', chain)
        # dataset=file_name.split('-')[0]

        dataset=file_name.rstrip('-{}'.format(chain))

        # print('dataset:', dataset)
        # antigen=file_name.split('-')[0].split('_')[-1]
        antigen=dataset.split('_')[-1]


        # print('antigen:', antigen)
        df=pd.DataFrame(pd.read_csv(file, sep='\t',header=None).values.reshape(-1, 2),columns=['header', 'seq_nt'])

        #print(file_name, '/ chain:', chain, '/dataset:', dataset, '/ antigen:', antigen, )

        df['mice']=mice
        df['chain']=chain
        df['dataset']=dataset
        df['antigen']=antigen
        df['type']='query'
        df.loc[0, 'type']='ref'
        df['ref_nt']=df.loc[0, 'seq_nt']
        df['ref_aa']=''.join(translate(df.loc[0, 'seq_nt'])[0])

        for i in df.index:
            translation=translate(df.loc[i, 'seq_nt'])
            df.loc[i, 'seq_aa']=''.join(translation[0])
            df.loc[i, 'frameshift']=translation[1]
            df.loc[i, 'stopcodon']=translation[2]
            df.loc[i, 'len_nt']=len(df.loc[i, 'seq_nt'])
            df.loc[i, 'len_aa']=len(translation[0])

        df['header'] = df['header'].apply(lambda x: x.lstrip('>'))
        dfs=pd.concat([dfs,df])
    return(dfs)

def ins_dels_miss(df):
    df.reset_index(inplace=True, drop=True)

    df[['nt_ins', 'nt_dels', 'nt_miss']]=0
    for i in df.index:
        query=df.loc[i, 'seq_nt']
        ref=df.loc[i, 'ref_nt']
        for p in range(0, len(ref)):

            if query[p] != ref[p]:
                if ref[p] == '-': df.loc[i, 'nt_ins']+=1
                elif query[p] == '-': df.loc[i, 'nt_dels']+=1
                else: df.loc[i, 'nt_miss']+=1

    df[['aa_ins', 'aa_dels', 'aa_miss']]=0
    for i in df.index:
        query=df.loc[i, 'seq_aa']
        ref=df.loc[i, 'ref_aa']
        for p in range(0, len(ref)):

            if query[p] != ref[p]:
                if ref[p] == del_sign: df.loc[i, 'aa_ins']+=1
                elif query[p] == del_sign: df.loc[i, 'aa_dels']+=1
                else: df.loc[i, 'aa_miss']+=1
    return(df)

def expand_aa(df_initial):
    df=df_initial.copy()
    df_aa_empty=pd.DataFrame(columns=['A{}'.format(i) for i in range(0,int(df['len_aa'].max()))], index=df.index)
    df=pd.concat([df, df_aa_empty], axis=1)

    for i in df.index:
        query=df.loc[i, 'seq_aa']
        ref=df.loc[i, 'ref_aa']
        for p in range(0, len(ref)):
            df.loc[i, 'A{}'.format(p)]=query[p]
    return(df)

def expand_nt(df_initial):
    df=df_initial.copy()
    df_aa_empty=pd.DataFrame(columns=['NT{}'.format(i) for i in range(0,int(df['len_nt'].max()))], index=df.index)
    df=pd.concat([df, df_aa_empty], axis=1)

    for i in df.index:
        query=df.loc[i, 'seq_nt']
        ref=df.loc[i, 'ref_nt']
        for p in range(0, len(ref)):
            df.loc[i, 'NT{}'.format(p)]=query[p]
    return(df)

def df_clean_up(df, zero_miss='exclude'): # By default, excludes sequences without nt mismatches. (mutations=insertions+deletions+mismatches)
    if zero_miss=='exclude':
        df=df[df['nt_miss']!=0]
    df=df[df['frameshift']==False]
    df=df[df['stopcodon']==False]
    df=df[df['type']=='query']
    return(df)

def num_nt_miss(df):
    nt_miss_list=[]

    for chain in ['VH', 'VL']:
        df_per_chain=df[df['chain']==chain]

        grouping=df_per_chain.groupby(by=['mice', 'dataset'])
        dfs_chain=pd.DataFrame()

        for grouped, df_chain in grouping:
            ID='_'.join(grouped)
            # print(ID)
            df_chain.rename(columns={'nt_miss':ID}, inplace=True)
            df_chain=df_chain[[ID]].sort_values(by=ID)
            df_chain.reset_index(inplace=True, drop=True)
            dfs_chain=pd.concat([dfs_chain, df_chain], axis=1)
        nt_miss_list.append(dfs_chain)
    return(nt_miss_list)

def num_aa_miss(df):
    aa_miss_list=[]

    for chain in ['VH', 'VL']:
        df_per_chain=df[df['chain']==chain]

        grouping=df_per_chain.groupby(by=['mice', 'dataset'])
        dfs_chain=pd.DataFrame()

        for grouped, df_chain in grouping:
            ID='_'.join(grouped)
            # print(ID)
            df_chain.rename(columns={'aa_miss':ID}, inplace=True)
            df_chain=df_chain[[ID]].sort_values(by=ID)
            df_chain.reset_index(inplace=True, drop=True)
            dfs_chain=pd.concat([dfs_chain, df_chain], axis=1)
        aa_miss_list.append(dfs_chain)
    return(aa_miss_list)

def num_nt_miss_donuts(df):

    df_categories = pd.DataFrame()
    num_Sequences_list=[]
    bins=[0] + list(range(1, 102,5))
    labels=['0', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30', '31-35', '36-40', '41-45', '46-50', '51-55', \
                       '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90', '91-95', '96-100']
    for column in df.columns:
        df_now = pd.DataFrame()
        df_now=pd.cut(df[column], bins=bins, labels=labels, right=False) # includes left and excludes right boundary
        # display(df_now.dropna())
        num_Sequences=df[column].count()
        print(column,':', num_Sequences)
        num_Sequences_list.append(num_Sequences)
        # df_categories[column]=list(df_now.groupby(df_now).count()/num_Sequences)
        df_categories[column]=list(df_now.groupby(df_now).count())
    df_categories=df_categories.set_axis(labels)
    df_categories.loc['num_seqs'] = num_Sequences_list
    return(df_categories)

def seq_logo_prep(df, seq_type):
    # create Logo object
    df.reset_index(inplace=True, drop=True)
    df_size=len(df)
    # print(df_size)
    # display(df)

    columns_dic={'nt':nts_list, 'aa':aas_list}
    Ref='ref_{}'.format(seq_type)

    ref_seq=df.loc[0, Ref]
    df_counts=pd.DataFrame(0, index=[i for i in range(0, len(ref_seq))], columns=columns_dic[seq_type])

    for seq, p in itertools.product(*[range(0, df_size), range(0, len(ref_seq))]):
        ref=ref_seq[p]
        query=df.loc[seq, 'seq_{}'.format(seq_type)][p]
        if query==ref: continue
        df_counts.loc[p, query]+=1

    df_logo=df_counts/df_size
    return(df_logo)

def seq_logo_plot(df, output_folder, ID):

    plt.rc('font', size=20)
    fig, mxx = plt.subplots(figsize=(20,2))

    nn_logo = logomaker.Logo(df, color_scheme='chemistry', font_name='Times new roman', ax=mxx)

    # style using Logo methods
    nn_logo.style_spines(visible=True) # Up and right borders
    # nn_logo.style_spines(spines=['left'], visible=True, bounds=[0, .75])
    nn_logo.style_spines(spines=['left'], visible=True, bounds=[0, 1.0])

    # style using Axes methods
    # nn_logo.ax.set_xlim([20, 115])
    nn_logo.ax.set_xticks([])
    # nn_logo.ax.set_ylim([-.6, .75])
    nn_logo.ax.set_ylim([0, 1.])
    plt.yticks([0.2,0.4,0.6,0.8], weight='bold')

    if ID in ['B18-383_Passenger_VH_aa', 'B18-383_Passenger_VH_nt', \
              'HA-WT_1_1_mix_VL_aa', 'HA-WT_1_1_mix_VL_nt', \
              'HA-WT_1_1_mix_VH_aa', 'HA-WT_1_1_mix_VH_nt']:

        nn_logo.style_spines(spines=['left'], visible=True, bounds=[0, 0.4])
        nn_logo.ax.set_ylim([0, 0.4])
        plt.yticks([0.1,0.2, 0.3], weight='bold')


    # nn_logo.ax.set_yticks([0, .75])
    # nn_logo.ax.set_yticklabels(['0', '0.75'])
    # nn_logo.ax.figure(figsize=(100,30))
    # nn_logo.ax.set_yticks([])
    # plt.yticks([0.2,0.4], weight='bold')
    print(ID)
    plt.savefig('{}/{}.jpg'.format(output_folder, ID, dpi=1200))
    # plt.show()
    plt.close()
    time.sleep(0.1)
    plt.pause(0.0001)

_modes=modes_dic.keys()
_modes

# Section1: preparation
output_folder_prep=set_output_folder('1_prep')

files_List=glob.glob('{}/input/*/*-*.fasta'.format(data_folder))

dfs=import_fasta(files_List)
dfs.reset_index(drop=True, inplace=True)
dfs

df_ha_wt=pd.concat([dfs[dfs['dataset']=='1_1_OVA'], pd.concat([dfs[dfs['dataset']=='1_1_APC'],dfs[dfs['dataset']=='1_1_CGG']])])
df_ha_wt['dataset']='1_1_mix'
df_ha_wt['antigen']='mix'

df_ha_wt

dfs=pd.concat([dfs, df_ha_wt])
dfs.reset_index(inplace=True, drop=True)
dfs #18046 including reference sequences

dfs=ins_dels_miss(dfs)
dfs

dfs.to_csv('{}/dfs_all.tsv'.format(output_folder_prep), sep = '\t', index=False)

dfs

dfs_expanded_nts=expand_nt(dfs)
dfs_expanded_nts.to_csv('{}/dfs_expanded_nts.tsv'.format(output_folder_prep), sep = '\t', index=False)

dfs_expanded_nts

dfs_expanded_nts_cleaned=df_clean_up(dfs_expanded_nts)
dfs_expanded_nts_cleaned.reset_index(inplace=True, drop=True)
dfs_expanded_nts_cleaned.to_csv('{}/dfs_expanded_nts_cleaned.tsv'.format(output_folder_prep), sep = '\t', index=False)
dfs_expanded_nts_cleaned

df_nts_mismatch=dfs_expanded_nts_cleaned.copy()
df_frq_nts_mismatch=pd.DataFrame()

grouping=df_nts_mismatch.groupby(by=['mice', 'dataset', 'chain'])

for grouped, df in grouping:
    suffix='_'.join(grouped)
    df.reset_index(drop=True, inplace=True)
    num_seqs=len(df)
    ref_nt=df.loc[0, 'ref_nt']
    for p in range(0, len(ref_nt)):

        mismatches_list = ['A', 'T', 'C', 'G']
        mismatches_list.remove(ref_nt[p])
        instances = len(df[df['NT{}'.format(p)].apply(lambda x: True if x in mismatches_list else False)])
        df_frq_nts_mismatch.loc[p, suffix] = instances / num_seqs

df_frq_nts_mismatch.to_csv('{}/frq_nts_mismatches.tsv'.format(output_folder_prep), sep = '\t', index=False)

dfs_expanded_aas=expand_aa(dfs)
dfs_expanded_aas.to_csv('{}/dfs_expanded_aas.tsv'.format(output_folder_prep), sep = '\t', index=False)

dfs_expanded_aas

dfs_expanded_aas_cleaned=df_clean_up(dfs_expanded_aas)
dfs_expanded_aas_cleaned.reset_index(inplace=True, drop=True)
dfs_expanded_aas_cleaned.to_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep = '\t', index=False)
dfs_expanded_aas_cleaned

# Section2: inspecting the number of mismatches per sequence in each dataset excluding sequences without mismatches
# used in Figs 1H, 4A, 5G, 6F
output_folder_num_miss=set_output_folder('2_num_miss')

dfs=pd.read_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)
prefix='w_aa_zmis'

dfs

num_nt_miss_list=num_nt_miss(dfs)

set(num_nt_miss_list[0])

len(num_nt_miss_list)

num_nt_miss_list[1]

sorting_by={
             'HA_0_1_OVA':1, 'HA_0_1_APC':2, 'HA_0_1_CGG':3, \
             'HA_0_1_OVA-CTLA4':4, 'HA_0_1_OVA-Isotype':5,\

             'HA_1_1_OVA': 6, 'HA_1_1_APC':7, 'HA_1_1_CGG':8, \
             'HA_1_100_OVA':9, 'HA_1_100_APC':10, 'HA_1_100_CGG':11,\
             'HA_1_1000_OVA':12, 'HA_1_1000_APC':13, 'HA_1_1000_CGG':14,\
             'HA_1_1000_CGG-CTLA4':15, 'HA_1_1000_CGG-Isotype':16,\

             'HA_1_1_mix':17,\

             'B18-383_OVA':18, 'B18-383_APC':19, 'B18-383_CGG':20,\
             'B18-383_OVA-CTLA4':21, 'B18-383_OVA-Isotype':22, }

sorted_list=sorted(set(num_nt_miss_list[1]), key=lambda x: sorting_by.get(x, float('inf')))

sorted_list

num_nt_miss_vl=num_nt_miss_list[1][sorted_list]
num_nt_miss_vl.to_csv('{}/nt_miss_vl.tsv'.format(output_folder_num_miss), sep = '\t', index=False)
num_nt_miss_vl#.count()

sorted_list=sorted(set(num_nt_miss_list[0]), key=lambda x: sorting_by.get(x, float('inf')))
sorted_list

num_nt_miss_vh=num_nt_miss_list[0][sorted_list]
num_nt_miss_vh.to_csv('{}/nt_miss_vh.tsv'.format(output_folder_num_miss), sep = '\t', index=False)
num_nt_miss_vh#.count()

num_aa_miss_list=num_aa_miss(dfs)

set(num_aa_miss_list[0])

len(num_aa_miss_list)

num_aa_miss_list[1]

sorting_by={
             'HA_0_1_OVA':1, 'HA_0_1_APC':2, 'HA_0_1_CGG':3, \
             'HA_0_1_OVA-CTLA4':4, 'HA_0_1_OVA-Isotype':5,\

             'HA_1_1_OVA': 6, 'HA_1_1_APC':7, 'HA_1_1_CGG':8, \
             'HA_1_100_OVA':9, 'HA_1_100_APC':10, 'HA_1_100_CGG':11,\
             'HA_1_1000_OVA':12, 'HA_1_1000_APC':13, 'HA_1_1000_CGG':14,\
             'HA_1_1000_CGG-CTLA4':15, 'HA_1_1000_CGG-Isotype':16,\

             'HA_1_1_mix':17,\

             'B18-383_OVA':18, 'B18-383_APC':19, 'B18-383_CGG':20,\
             'B18-383_OVA-CTLA4':21, 'B18-383_OVA-Isotype':22, }

sorted_list=sorted(set(num_aa_miss_list[1]), key=lambda x: sorting_by.get(x, float('inf')))

sorted_list

num_aa_miss_vl=num_aa_miss_list[1][sorted_list]
num_aa_miss_vl.to_csv('{}/aa_miss_vl.tsv'.format(output_folder_num_miss), sep = '\t', index=False)
num_aa_miss_vl#.count()

sorted_list=sorted(set(num_aa_miss_list[0]), key=lambda x: sorting_by.get(x, float('inf')))
num_aa_miss_vh=num_aa_miss_list[0][sorted_list]
num_aa_miss_vh.to_csv('{}/aa_miss_vh.tsv'.format(output_folder_num_miss), sep = '\t', index=False)
num_aa_miss_vh#.count()

# Section3: donuts
# used in Fig 1G
output_folder_donuts=set_output_folder('3_donuts')

dfs=pd.read_csv('{}/dfs_all.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)
dfs

dfs=df_clean_up(dfs, zero_miss='include')
dfs.reset_index(inplace=True, drop=True)
dfs

num_nt_miss_list=num_nt_miss(dfs)

num_nt_miss_list[1]

list(num_nt_miss_list[0]['HA_1_1_OVA'].dropna())

print('VH:')
donuts_vh=num_nt_miss_donuts(num_nt_miss_list[0])

donuts_vh=num_nt_miss_donuts(num_nt_miss_list[0])
donuts_vh.to_csv('{}/donuts_vh.tsv'.format(output_folder_donuts), sep = '\t', index=False)
donuts_vh

donuts_vh.sum(axis=1)

print('VL:')
donuts_vl=num_nt_miss_donuts(num_nt_miss_list[1])
donuts_vl.to_csv('{}/donuts_vl.tsv'.format(output_folder_donuts), sep = '\t', index=False)
donuts_vl

donuts_vl.sum(axis=1)

# Section4: rs ratios prep
output_folder_rs_prep=set_output_folder('4_prep_rs')

dfs=pd.read_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)

# Method selection
rs_method='R/S_spike_sil1' #add 1 to silent mutations to avoid division by zero error

def align_regions(seq, mice, chain):
    # print('status 0')

    if chain=='VL':
        # print('status 1')
        FR1=seq[0:24] #8
        CDR1=seq[24:60] #12
        FR2=seq[60:111] #17
        CDR2=seq[111:120] #3
        FR3=seq[120:228] #36
        CDR3=seq[228:255] #9
        FR4=seq[255:285] #10

    elif mice=='HA':
        # print('status 2')
        FR1=seq[0:15] #5
        CDR1=seq[15:39] #8
        FR2=seq[39:90] #17
        CDR2=seq[90:114] #8
        FR3=seq[114:228] #38
        CDR3=seq[228:267] #13
        FR4=seq[267:300] #11

    elif mice=='B18-383':
        # print('status 3')
        FR1=seq[0:33] #11
        CDR1=seq[33:57] #8
        FR2=seq[57:108] #17
        CDR2=seq[108:132] #8
        FR3=seq[132:246] #38
        CDR3=seq[246:285] #13
        FR4=seq[285:318] #11

    return ([FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4])

def rs_ratio(df_initial):
    df=df_initial.copy()
    df.reset_index(inplace=True, drop=True)

    df[['replacement_mut', 'silent_mut']]=''

    df[['R-S_check']]='' #separated calculation

    for i in df.index: # each row
        NT_replacement_list, NT_silent_list, NT_r_s = [], [], []

        query_seq=df.loc[i, 'seq_nt']
        ref_seq=df.loc[i, 'ref_nt']

        seq_regions=align_regions(query_seq, df.loc[i, 'mice'], df.loc[i, 'chain'])
        ref_regions=align_regions(ref_seq, df.loc[i, 'mice'], df.loc[i, 'chain'])

        for seq_region, ref_region in zip (seq_regions, ref_regions): #Iterate for each region in one sequence

            NT_replacement, NT_silent=0,0 #How many AA are replaced or not replaced

            ref_codons=[ref_region[i:i+3] for i in range(0, len(ref_region), 3)] #Split to codons
            seq_codons=[seq_region[i:i+3] for i in range(0, len(seq_region), 3)] #Split to codons
            # print(seq)
            # print(ref)

            for codon_seq, codon_ref in zip (seq_codons, ref_codons): #Iterate for each codon in this region

                if codon_seq=='---': #Codon deletion is not count as a mutation
                    continue  #skip this iteration

                elif codon_ref=='---': #Codon insertion is not count as a mutation
                    continue  #skip this iteration

                elif aas_dic[codon_seq]!=aas_dic[codon_ref]: #AA Replacement (Missense)
                    #print('R')
                    # AA_replacement+=1
                    NT_replacement+=sum(c1!=c2 for c1,c2 in zip(codon_seq,codon_ref))

                elif codon_seq!=codon_ref and aas_dic[codon_seq]==aas_dic[codon_ref]: #AA Silent mutation
                    #print('S')
                    # AA_silent+=1
                    NT_silent+=sum(c1!=c2 for c1,c2 in zip(codon_seq,codon_ref))

                elif codon_seq==codon_ref: #No change
                    pass
                    #print('N')

                else:
                    raise Exception("Something is wrong!")

            NT_replacement_list.append(NT_replacement)
            NT_silent_list.append(NT_silent)
            NT_r_s.append(NT_replacement-NT_silent) #separated calculation

        df.at[i, 'replacement_mut']=NT_replacement_list
        df.at[i, 'silent_mut']=NT_silent_list
        df.at[i, 'R-S_check']=NT_r_s #separated calculation

    return(df)

dfs=rs_ratio(dfs)
dfs

rs_method

for i in dfs.index: # each row
    R_list=dfs.loc[i, 'replacement_mut']
    S_list=dfs.loc[i, 'silent_mut']
    for j, region in enumerate(['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']):

        if rs_method == 'R/S_spike_sil1':
            dfs.loc[i, region]=R_list[j]/(S_list[j]+1) # adds one silent mutation to all cases



dfs

dfs[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].max()

dfs[dfs['FR1']==22.0]

datasets = [
    (dfs['dataset'].str.startswith('1_1_')),
    (dfs['dataset'].str.startswith('1_100_')),
    (dfs['dataset'].str.startswith('1_1000_')),
    (dfs['dataset'].str.startswith('0_1_'))
    ]

rank1 = ['1', '2', '3', '4']

dfs['datasets_rank'] = np.select(datasets, rank1)

antigens = [

    (dfs['antigen'] == 'Passenger'),
    (dfs['antigen'] == 'OVA'),
    (dfs['antigen'] == 'APC'),
    (dfs['antigen'] == 'CGG'),
    (dfs['antigen'] == 'mix'),
    (dfs['antigen'] == 'OVA-Isotype'),
    (dfs['antigen'] == 'OVA-CTLA4'),
    (dfs['antigen'] == 'CGG-Isotype'),
    (dfs['antigen'] == 'CGG-CTLA4'),
]

rank2 = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

dfs['antigens_rank'] = np.select(antigens, rank2)

dfs.to_csv('{}/dfs_rs_ratios.tsv'.format(output_folder_rs_prep), sep = '\t', index=False)
dfs

# Section5: rs ratios
# used in Figs 4B, S5A, S5B
output_folder_rs=set_output_folder('5_rs')


# Figure selection
# rs_figure='Isotype-CTLA4'
# rs_figure='other-antigens-VH'
# rs_figure='other-antigens-VL'
# rs_figure='other-antigens'

for rs_figure in ['Isotype-CTLA4', 'other-antigens']:
    dfs=pd.read_csv('{}/dfs_rs_ratios.tsv'.format(output_folder_rs_prep), sep='\t', header=0, low_memory=False)

    IsCT_datasets=['OVA-CTLA4', 'OVA-Isotype', '0_1_OVA-CTLA4', '0_1_OVA-Isotype', '1_1000_CGG-CTLA4', '1_1000_CGG-Isotype']

    if rs_figure=='Isotype-CTLA4':
        dfs=dfs[dfs["dataset"].isin(IsCT_datasets)] # Subset the relevant datasets

    elif rs_figure=='other-antigens-VH':
        dfs=dfs[~dfs["dataset"].isin(IsCT_datasets+['1_1_mix'])] # Subset the relevant datasets
        dfs=dfs[dfs['chain']=='VH']

    elif rs_figure=='other-antigens-VL':
        dfs=dfs[~dfs["dataset"].isin(IsCT_datasets+['1_1_mix'])] # Subset the relevant datasets
        dfs=dfs[dfs['chain']=='VL']

    elif rs_figure=='other-antigens':
        dfs=dfs[~dfs["dataset"].isin(IsCT_datasets+['1_1_mix'])] # Subset the relevant datasets


    dfs[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].max()

    RS_max=dfs[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].max().max()
    print('RS_max =', RS_max)

    dfs[dfs['FR3']==RS_max]

    dfs

    # dfs[dfs['FR3']==15.0]

    dfs.isna().sum()[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']]

    dfs.isna().sum()[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].sum()

    set(dfs['dataset'])



    # dfs.replace({np.nan : RS_max}, inplace=True) # only for the rs_method='R/S_im' method

    dfs

    dfs_mean=dfs.groupby(by=['mice', 'chain', 'antigen', 'dataset', 'antigens_rank', 'datasets_rank'])[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].mean()
    dfs_mean=dfs_mean.sort_index(level=['chain', 'mice', 'antigens_rank', 'datasets_rank'], ascending=[True,False,True,True])
    dfs_mean.reset_index(inplace=True)
    dfs_mean.to_csv('{}/dfs_rs_mean_no_repeat_{}.tsv'.format(output_folder_rs, rs_figure), sep = '\t', index=False)
    dfs_mean

    # raise Exception("Stop!")

    print('For heatmap range, man mean:', dfs_mean[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']].max().max())

    # Distribution of R/S ratios per dataset single chain
    grouping=dfs.groupby(by=['mice', 'dataset', 'chain'])[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']]

    for grouped, df in grouping:
        suffix='_'.join(grouped)
        # print(suffix)
        # display(df)

        df.rename(columns={'FR1':'FR1_{}'.format(suffix), 'FR1':'FR1_{}'.format(suffix), 'FR2':'FR2_{}'.format(suffix),\
                           'FR3':'FR3_{}'.format(suffix), 'FR4':'FR4_{}'.format(suffix), 'CDR1':'CDR1_{}'.format(suffix),\
                           'CDR2':'CDR2_{}'.format(suffix), 'CDR3':'CDR3_{}'.format(suffix),}, inplace=True)

        # df.to_csv('{}/rs_{}.tsv'.format(output_folder, '_'.join(grouped)), sep = '\t', index=False)

    # Distribution of R/S ratios per dataset both chains
    grouping=dfs.groupby(by=['mice', 'dataset'])[['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4']]

    for grouped, df in grouping:
        suffix='{}_VHVL'.format('_'.join(grouped))
        # print(suffix)

        df.rename(columns={'FR1':'FR1_{}'.format(suffix), 'FR1':'FR1_{}'.format(suffix), 'FR2':'FR2_{}'.format(suffix),\
                           'FR3':'FR3_{}'.format(suffix), 'FR4':'FR4_{}'.format(suffix), 'CDR1':'CDR1_{}'.format(suffix),\
                           'CDR2':'CDR2_{}'.format(suffix), 'CDR3':'CDR3_{}'.format(suffix),}, inplace=True)
        # display(df)
        if rs_figure=='other-antigens': continue
        df.to_csv('{}/rs_distribution_{}_{}.tsv'.format(output_folder_rs, rs_figure, suffix), sep = '\t', index=False)

    df
    def repeat(df_initial):
        df=df_initial.copy()
        chain=list(set(df['chain']))[0]
        mice=list(set(df['mice']))[0]

        if chain=='VH':
            if mice=='HA': repeats=[5,8,17,8,38,13,11]
            elif mice=='B18-383': repeats=[11,8,17,8,38,13,11]
        elif chain=='VL':
            repeats=[8,12,17,3,36,9,10]

        regions=[6,7,8,9,10,11,12]

        df_now=df.iloc[:,[0,1,3]]

        for region, repeat in zip(regions, repeats):
            for r in range(0, repeat):
                df_now=pd.concat([df_now, df.iloc[:,region]], axis=1)
        return(df_now)

    dfs_mean

    grouping=dfs_mean.groupby(by=['mice', 'chain'])

    for grouped, df in grouping:
        print(grouped)
        df=repeat(df)
        # display(df)

        df.to_csv('{}/rs_mean_repeated_{}_{}.tsv'.format(output_folder_rs, rs_figure, '_'.join(grouped)), sep = '\t', index=False)

# Section6: Preparing data for scatter plots
output_folder_prep_plots=set_output_folder('6_prep_plots')

modes_dic

for mode in modes_dic.keys():
    dfs=pd.read_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)
    chain='VH'

    mice=modes_dic[mode]['mice']
    datasets=modes_dic[mode]['datasets']

    print('Mode :', mode, mice, chain, datasets)
    # if mode==1 and chain=='VL':
    #     datasets.remove('Passenger')
    #     print('passenger removed from VL dataset list', datasets)
    mice

    dfs=dfs[dfs['mice']==mice]
    dfs=dfs[dfs['chain']==chain]
    dfs=dfs[dfs["dataset"].isin(datasets)] # Subset the relevant datasets
    set(dfs['dataset'])

    # dfs=my.funcs.df_clean_up(dfs)
    # dfs=dfs[dfs['aa_miss']!=0]
    dfs.reset_index(inplace=True, drop=True)
    dfs

    len_aa=int(dfs['len_aa'].values[0])
    len_aa

    ref_aa=dfs['ref_aa'].values[0]
    ref_aa

    df_aa=pd.DataFrame(columns=[i for i in range(0,len_aa)], index=aas_list)
    df_aa.fillna(np.nan, inplace=True)
    df_aa

    def update_df_aa(case, position, diff_list):

        value = { 1:1, 2:2, 3:3, 12: 4, 13: 5, 23: 6, 123: 7, 4:8, 14:9, 24:10, 34:11, 124:12, 234:13, 134:14, 1234:15, 'ref':np.nan}

        if diff_list:
            for aminoacid in diff_list:
                df_aa.loc[aminoacid,position]=value[case]

    aa_dic = dict()  #{'OVA_p0': {'A', 'S'}, 'APC_p0': {'S'}, 'CGG_p0': {'P', 'S'}}
    # dif_dic = dict() #{'OVA_aa0': ['A'], 'APC_aa0': [],  'CGG_aa0': ['P'], 'intersection_aa0': ['S']}

    for i in range(0, len_aa): #106 for B18-383, 101 for HA-uMT
        for dataset in datasets:
            # dataset='OVA'
            df_now=dfs[dfs['dataset']==dataset]
            # print(dataset, 'df length: ', len(df_now))
            aa_dic['{}_p{}'.format(dataset, i)]=set(df_now.loc[:,'A{}'.format(i)])


        A=aa_dic['{}_p{}'.format(datasets[0],i)]
        B=aa_dic['{}_p{}'.format(datasets[1],i)]
        if len(datasets)==2:
            C=set()
            D=set()
        elif len(datasets)==3:
            C=aa_dic['{}_p{}'.format(datasets[2],i)]
            D=set()
        elif len(datasets)==4:
            C=aa_dic['{}_p{}'.format(datasets[2],i)]
            D=aa_dic['{}_p{}'.format(datasets[3],i)]

        # print(A)
        # print(B)
        # print(C)
        # print(D)

        update_df_aa(1, i, list(A.difference(B, C, D)))
        update_df_aa(2, i, list(B.difference(A, C, D)))
        update_df_aa(3, i, list(C.difference(A, B, D)))
        update_df_aa(4, i, list(D.difference(A, B, C)))


        _12=A & B
        _13=A & C
        _23=B & C

        _14=A & D
        _24=B & D
        _34=C & D

        _123=A & B & C
        _124=A & B & D
        _134=A & C & D
        _234=B & C & D

        _1234=A & B & C & D

        update_df_aa(12, i, _12-_123-_124-_1234)
        update_df_aa(13, i, _13-_123-_134-_1234)
        update_df_aa(23, i, _23-_123-_234-_1234)

        update_df_aa(14, i, _14-_124-_134-_1234)
        update_df_aa(24, i, _24-_124-_234-_1234)
        update_df_aa(34, i, _34-_234-_134-_1234)

        update_df_aa(123, i, _123-_1234)
        update_df_aa(124, i, _124-_1234)
        update_df_aa(134, i, _134-_1234)
        update_df_aa(234, i, _234-_1234)

        update_df_aa(1234, i, _1234)
        update_df_aa('ref', i, list(ref_aa[i]))
        # break

    df_aa

    df_aa.to_csv('{}/df_aa_mode{}_{}.tsv'.format(output_folder_prep_plots, mode, chain), sep = '\t', index=True)

# Section7: frquency and privacy index calculations and scatter plots
# used in Figs 4D, 4E, 4G, 5I, 6H, S4D, S5C
output_folder_scatter_plots=set_output_folder('7_scatter_plots')

for mode in modes_dic.keys():

    df_nts=pd.read_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)
    mice=modes_dic[mode]['mice']
    datasets=modes_dic[mode]['datasets']
    datasets
    chain='VH'

    print('Mode :', mode, mice, chain, datasets)

    df_nts

    df_aa=pd.read_csv('{}/df_aa_mode{}_{}.tsv'.format(output_folder_prep_plots, mode, chain), sep='\t', header=0, low_memory=False, index_col=0)
    # df_aa

    df_nts=df_nts[df_nts['mice']==mice]
    df_nts=df_nts[df_nts['chain']==chain]
    df_nts=df_nts[df_nts["dataset"].isin(datasets)] # Subset the relevant datasets
    df_nts

    #Mode 1
    # 5136 AA MM + AA nMM
    # 4810 AA MM

    print(set(df_nts['mice']), set(df_nts['chain']), set(df_nts['dataset']))

    for dataset in datasets:
        print('{}:{}'.format(dataset, len(df_nts[df_nts['dataset']==dataset])))

    # 2:3
    # 3:7
    # 4:15

    all_combinations=list(np.arange(1.0,16.0)) #15

    if len(datasets)==2:
        cases_list_wanted_all={'shareds':[4.0], 'uniques':[1.0, 2.0], 'oxaxcx':[], 'remainders':[]}

    if len(datasets)==3:
        cases_list_wanted_all={'shareds':[7.0], 'uniques':[1.0, 2.0, 3.0], 'oxaxcx':[5.0, 6.0], 'remainders':[4.0]}

    if len(datasets)==4:
        cases_list_wanted_all={'shareds':[15.0], 'uniques':[1.0, 2.0, 3.0, 8.0], 'oxaxcx':[9.0, 10.0, 11.0], 'remainders':[4.0, 5.0, 6.0, 7.0, 12.0, 13.0, 14.0]}

    dfs_cases_dic={}
    dfs_cases_num_dic={}
    for case, list_wanted in cases_list_wanted_all.items():

        dic_num_blocks={}

        list_exclusion=all_combinations.copy()
        df_aa_now=df_aa.copy()

        for i in list_wanted:
            list_exclusion.remove(i)

        for j in list_exclusion:
            df_aa_now.replace({ j : np.nan}, inplace=True)

        df_aa_now=df_aa_now.iloc[:-2,:]

        for n in list_wanted:
            dic_num_blocks[n]=(df_aa_now == n).sum().sum()
        # print(case, dic_num_blocks)
        dfs_cases_dic[case]=df_aa_now
        dfs_cases_num_dic[case]=dic_num_blocks
        # display(df_aa_now)

    datasets

    dfs_cases_num_dic

    pd.DataFrame.from_dict(dfs_cases_num_dic, orient='index').T.to_csv('{}/mode{}_numbers_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)

    dfs_cases_dic.keys()

    dfs_cases_dic['shareds']

    df_pos_shareds = dfs_cases_dic['shareds'].T
    df_pos_uniques = dfs_cases_dic['uniques'].T
    df_pos_oxaxcx = dfs_cases_dic['oxaxcx'].T
    df_pos_remainders = dfs_cases_dic['remainders'].T

    # df_pos_uniques = pd.read_csv('../data/avneesh/affinity_birth/output/F005/mode{}-1_2_3_8.tsv'.format(mode), sep='\t', header=0, low_memory=False, index_col=0)
    # df_pos_all = pd.read_csv('../data/avneesh/affinity_birth/output/F005/mode{}-1_2_3_4_5_6_7_8_9_10_11_12_13_14_15.tsv'.format(mode), sep='\t', header=0, low_memory=False, index_col=0)

    # df_nts = pd.read_csv('../data/avneesh/affinity_birth/output/F002/all_nt_seq_based.tsv', sep='\t', header=0, low_memory=False) #, index_col='codon'
    # df_nts=df_nts[df_nts['frameshift']==False]
    # df_nts=df_nts[df_nts['stopcodon']==False]

    set(df_nts['dataset'])

    # if mode==5:
    #     df_now=pd.concat([df_nts[df_nts['dataset']=='1_1_OVA'], pd.concat([df_nts[df_nts['dataset']=='1_1_APC'],df_nts[df_nts['dataset']=='1_1_CGG']])])
    #     df_now['dataset']='HA-WT'
    #     df_nts=pd.concat([df_nts,df_now])



    mice

    set(df_nts['mice'])

    aa_ref_seq=df_nts['ref_aa'].values[-1]

    len_ref=len(df_nts['ref_aa'].values[-1])

    df_nts

    multiplier_size=22 # To increase the size of shapes on the plot
    # df_nts.to_csv('{}/mode{}_relevant_dfs.tsv'.format(output_folder, mode), sep = '\t', index=True)

    # To translate AAs names to numbers to sort them by chemical properties
    aa_D = dict()
    aa_L=aas_chemistry_list
    for y, aa in enumerate(aa_L):
        aa_D[aa]=20-y
    print(aa_D)

    list_frequency=['{}_frequency'.format(ds) for ds in datasets]
    list_plot_size=['{}_plot_size'.format(ds) for ds in datasets]
    list_selec_uneven=['{}_selec_uneven'.format(ds) for ds in datasets]

    print(list_frequency)
    print(list_plot_size)
    print(list_selec_uneven)

    df_pos_uniques

    df_aa=df_aa.T

    df_nts['unique_num']=0
    df_nts['shared_num']=0

    # Uniques dataframe
    myDict = dict()

    ds_unique_names={1.0: 0, 2.0: 1, 3.0:2, 8.0:3} # For converting indices from 1238 to 0123

    for p, aa in itertools.product(*[df_aa.index,df_aa.columns[:-2]]):
        v = df_aa.loc[p,aa]
        # if v not in [1.0,2.0,3.0,8.0]: continue
        if v not in cases_list_wanted_all['uniques']: continue

        # myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa}
        myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa': aa, 'aa_n' : aa_D[aa]} # Translate AAs names to numbers to sort them by chemical properties

        n = ds_unique_names[v] # Convert indexing 1238 to 0123

        for i in df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])]['header'].index:
            # print(i)
            df_nts.loc[i,'unique_num']+=1

        instances = len(df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])])
        df_size = len(df_nts[df_nts["dataset"]==datasets[n]])
        df_frequency = instances/df_size*100

        myDict['{}-{}'.format(p,aa)]['{}_instances'.format(datasets[n])] = instances
        myDict['{}-{}'.format(p,aa)]['{}_size'.format(datasets[n])] = df_size

        myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])] = df_frequency
        # myDict['{}-{}'.format(p,aa)]['{}_plot_size'.format(datasets[n])] = df_frequency*multiplier_size
        myDict['{}-{}'.format(p,aa)]['{}_selec_uneven'.format(datasets[n])] = 0.1

    df_all_uniques=pd.DataFrame.from_dict(myDict).transpose()
    # display(df_all_uniques)
    df_all_uniques['status']='uniques'

    df_all_uniques['p']=df_all_uniques['p'].astype('int')+1 # Change the zero-indexing to one-indexing
    df_all_uniques.to_csv('{}/mode{}_uniques_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)

    df_all_uniques

    myDict = dict()

    for p, aa in itertools.product(*[df_aa.index, df_aa.columns[:-2]]):
        total_frequency = 0
        v = df_aa.loc[p,aa]
        # if v not in [9.0,10.0,11.0]: continue
        if v not in cases_list_wanted_all['oxaxcx']: continue

        # myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa}
        myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa, 'aa_n' : aa_D[aa]}
        for n in range(0, len(datasets)):
            # print(n)
            instances = len(df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])])
            df_size = len(df_nts[df_nts["dataset"]==datasets[n]])
            df_frequency = instances/df_size*100
            total_frequency+=df_frequency

            myDict['{}-{}'.format(p,aa)]['{}_instances'.format(datasets[n])] = instances
            myDict['{}-{}'.format(p,aa)]['{}_size'.format(datasets[n])] = df_size

            myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])] = df_frequency
            # myDict['{}-{}'.format(p,aa)]['{}_plot_size'.format(datasets[n])] = df_frequency*multiplier_size

        for n in range(0, len(datasets)):
            myDict['{}-{}'.format(p,aa)]['{}_selec_uneven'.format(datasets[n])] = myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])]/total_frequency

    df_oxaxcx=pd.DataFrame()

    if len(myDict)!=0:
        df_oxaxcx=pd.DataFrame.from_dict(myDict).transpose()

        df_oxaxcx['p']=df_oxaxcx['p'].astype('int')+1 # Change the zero-indexing to one-indexing

        df_oxaxcx['status']='oxaxcx'
        df_oxaxcx.to_csv('{}/mode{}_oxaxcx_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)

    df_oxaxcx

    # Remainder Shareds dataframe
    myDict = dict()

    for p, aa in itertools.product(*[df_aa.index, df_aa.columns[:-2]]):
        total_frequency = 0
        v = df_aa.loc[p,aa]
        if v not in cases_list_wanted_all['remainders']: continue

        # myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa}
        myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa, 'aa_n' : aa_D[aa]}
        for n in range(0, len(datasets)):
            # print(n)
            instances = len(df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])])
            df_size = len(df_nts[df_nts["dataset"]==datasets[n]])
            df_frequency = instances/df_size*100
            total_frequency+=df_frequency

            myDict['{}-{}'.format(p,aa)]['{}_instances'.format(datasets[n])] = instances
            myDict['{}-{}'.format(p,aa)]['{}_size'.format(datasets[n])] = df_size

            myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])] = df_frequency
            # myDict['{}-{}'.format(p,aa)]['{}_plot_size'.format(datasets[n])] = df_frequency*multiplier_size

        for n in range(0, len(datasets)):
            myDict['{}-{}'.format(p,aa)]['{}_selec_uneven'.format(datasets[n])] = myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])]/total_frequency

    df_remainders=pd.DataFrame()

    if len(myDict)!=0:

        df_remainders=pd.DataFrame.from_dict(myDict).transpose()

        df_remainders['p']=df_remainders['p'].astype('int')+1 # Change the zero-indexing to one-indexing

        df_remainders['status']='remainders'
        df_remainders.to_csv('{}/mode{}_remainders_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)

    df_remainders

    # Shareds dataframe
    myDict = dict()

    for p, aa in itertools.product(*[df_aa.index, df_aa.columns[:-2]]):
        total_frequency = 0
        v = df_aa.loc[p,aa]
        if v not in cases_list_wanted_all['shareds']: continue

        # myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa}
        myDict['{}-{}'.format(p,aa)] = {'p': p, 'aa' : aa, 'aa_n' : aa_D[aa]}
        for n in range(0, len(datasets)):
            # print(n)

            for i in df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])]['header'].index:
                df_nts.loc[i,'shared_num']+=1
                # print(datasets[n],i)
            instances = len(df_nts[(df_nts['A{}'.format(p)]==aa) & (df_nts['dataset']==datasets[n])])

            df_size = len(df_nts[df_nts["dataset"]==datasets[n]])
            df_frequency = instances/df_size*100
            total_frequency+=df_frequency

            myDict['{}-{}'.format(p,aa)]['{}_instances'.format(datasets[n])] = instances
            myDict['{}-{}'.format(p,aa)]['{}_size'.format(datasets[n])] = df_size

            myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])] = df_frequency
            # myDict['{}-{}'.format(p,aa)]['{}_plot_size'.format(datasets[n])] = df_frequency*multiplier_size

        for n in range(0, len(datasets)):
            myDict['{}-{}'.format(p,aa)]['{}_selec_uneven'.format(datasets[n])] = myDict['{}-{}'.format(p,aa)]['{}_frequency'.format(datasets[n])]/total_frequency
    df_all_shareds=pd.DataFrame.from_dict(myDict).transpose()
    df_all_shareds['status']='shareds'

    df_all_shareds['p']=df_all_shareds['p'].astype('int')+1 # Change the zero-indexing to one-indexing
    df_all_shareds.to_csv('{}/mode{}_shareds_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)
    df_all_shareds

    grouping=df_nts.groupby(by='dataset')
    dfs_unique_num=pd.DataFrame()
    dfs_shared_num=pd.DataFrame()

    for ds, df_now in grouping:
        cols=['{}_unique_num'.format(ds), '{}_shared_num'.format(ds)]
        df_now.rename(columns={'unique_num':cols[0], 'shared_num':cols[1]}, inplace=True)
        df_now=df_now[[cols[0], cols[1]]]

        df_now=df_now.sort_values(by=cols[0], ascending=False)
        df_now.reset_index(inplace=True, drop=True)
        dfs_unique_num=pd.concat([dfs_unique_num, df_now[cols[0]]], axis=1)

        df_now=df_now.sort_values(by=cols[1], ascending=False)
        df_now.reset_index(inplace=True, drop=True)
        dfs_shared_num=pd.concat([dfs_shared_num, df_now[cols[1]]], axis=1)

    list_dfs_unique_num=['{}_unique_num'.format(ds) for ds in datasets]
    dfs_unique_num[list_dfs_unique_num].to_csv('{}/mode{}_unique_num_{}_aa.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=False)

    list_dfs_shared_num=['{}_shared_num'.format(ds) for ds in datasets]
    dfs_shared_num[list_dfs_shared_num].to_csv('{}/mode{}_shared_num_{}_aa.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=False)

    # raise Exception("Stop!")

    n=10
    for c in list_frequency:
        print('top {} of:{}'.format(n, c))
        df_now=df_all_shareds.dropna().sort_values(by=c, ascending=False).iloc[0:n,]
        df_now=df_now[['p', 'aa']+list_frequency]
        # display(df_now)

        df_now.to_csv('{}/mode{}_top_{}_shared_{}_{}.tsv'.format(output_folder_scatter_plots, mode, n, c, chain), sep = '\t', index=False)


    df_all_shareds

    # Calculating min/max ranges for uniform size/intensity plotting
    min_frequency_uq=df_all_uniques[list_frequency].min(axis=0).min()
    max_frequency_uq=df_all_uniques[list_frequency].max(axis=0).max()

    min_selec_uneven_uq=df_all_uniques[list_selec_uneven].min(axis=0).min()
    max_selec_uneven_uq=df_all_uniques[list_selec_uneven].max(axis=0).max()

    # print('uniques :', min_frequency_uq, max_frequency_uq, min_selec_uneven_uq, max_selec_uneven_uq)

    min_frequency_shd=df_all_shareds[list_frequency].min(axis=0).min()
    max_frequency_shd=df_all_shareds[list_frequency].max(axis=0).max()

    min_selec_uneven_shd=df_all_shareds[list_selec_uneven].min(axis=0).min()
    max_selec_uneven_shd=df_all_shareds[list_selec_uneven].max(axis=0).max()

    # print('shareds :', min_frequency_shd, max_frequency_shd, min_selec_uneven_shd, max_selec_uneven_shd)

    min_frequency=min(min_frequency_uq, min_frequency_shd)
    max_frequency=max(max_frequency_uq, max_frequency_shd)

    min_selec_uneven=min(min_selec_uneven_uq, min_selec_uneven_shd)
    max_selec_uneven=max(max_selec_uneven_uq, max_selec_uneven_shd)

    print('# Mode :', mode, mice, chain, datasets,\
          '\n# uniques :', min_frequency_uq, max_frequency_uq, min_selec_uneven_uq, max_selec_uneven_uq,\
          '\n# shareds :', min_frequency_shd, max_frequency_shd, min_selec_uneven_shd, max_selec_uneven_shd,\
          '\n# uq/shd :', min_frequency, max_frequency, min_selec_uneven, max_selec_uneven)

    multi_index = pd.MultiIndex.from_tuples([('min', 'frequency'),('max', 'frequency'),
                                             ('min', 'selec_uneven'),('max', 'selec_uneven')])

    df_min_max = pd.DataFrame(index=multi_index, columns=['uniques','shareds', 'oxaxcx', 'remainders'])

    df_min_max.loc[('min', 'frequency'), 'uniques'] = min_frequency_uq
    df_min_max.loc[('max', 'frequency'), 'uniques'] = max_frequency_uq
    df_min_max.loc[('min', 'selec_uneven'), 'uniques'] = min_selec_uneven_uq
    df_min_max.loc[('max', 'selec_uneven'), 'uniques'] = max_selec_uneven_uq

    df_min_max.loc[('min', 'frequency'), 'shareds'] = min_frequency_shd
    df_min_max.loc[('max', 'frequency'), 'shareds'] = max_frequency_shd
    df_min_max.loc[('min', 'selec_uneven'), 'shareds'] = min_selec_uneven_shd
    df_min_max.loc[('max', 'selec_uneven'), 'shareds'] = max_selec_uneven_shd

    df_min_max.loc[('min', 'frequency'), 'oxaxcx'] = min_frequency_shd
    df_min_max.loc[('max', 'frequency'), 'oxaxcx'] = max_frequency_shd
    df_min_max.loc[('min', 'selec_uneven'), 'oxaxcx'] = min_selec_uneven_shd
    df_min_max.loc[('max', 'selec_uneven'), 'oxaxcx'] = max_selec_uneven_shd

    df_min_max.loc[('min', 'frequency'), 'remainders'] = min_frequency_shd
    df_min_max.loc[('max', 'frequency'), 'remainders'] = max_frequency_shd
    df_min_max.loc[('min', 'selec_uneven'), 'remainders'] = min_selec_uneven_shd
    df_min_max.loc[('max', 'selec_uneven'), 'remainders'] = max_selec_uneven_shd

    df_min_max.to_csv('{}/mode{}_df_min_max_{}.tsv'.format(output_folder_scatter_plots, mode, chain), sep = '\t', index=True)
    df_min_max

    min_frequency=0
    max_frequency=100
    min_selec_uneven=0
    max_selec_uneven=1.25

    intensity_range=np.linspace(min_selec_uneven, max_selec_uneven, 6, endpoint=True)
    intensity_range=[1.5,1,.75,.5,.25,.01]
    size_range=np.linspace(min_frequency, max_frequency, 6, endpoint=True)
    size_range=[150,100,75,50,25,1]
    aa_range=[3,8,10,12,14,16]

    intensity_range

    size_range

    dfs=[df_all_shareds, df_oxaxcx, df_remainders, df_all_uniques]
    statuses=['shareds', 'oxaxcx', 'remainders', 'uniques']

    for p, df_now, status in zip([len_ref+3, len_ref+5, len_ref+7, len_ref+9], dfs, statuses):

        for P, A, S, I in zip([p]*6, aa_range, size_range, [1]*6):
            if status=='uniques': I=0.1
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'p'] = P
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'aa_n'] = A
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'status'] = status

            for c in list_frequency:
                df_now.loc['m{}S_{}I-{}'.format(S,I,P), c] = S
            for c in list_selec_uneven:
                df_now.loc['m{}S_{}I-{}'.format(S,I,P), c] = I
        # display(df_now)
        for P, A, S, I in zip([p+9]*6, aa_range, [max_frequency*2/3]*6, intensity_range):
            if status=='uniques':
                I=0.1
                A=8
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'p'] = P
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'aa_n'] = A
            df_now.loc['m{}S_{}I-{}'.format(S,I,P), 'status'] = status

            for c in list_frequency:
                df_now.loc['m{}S_{}I-{}'.format(S,I,P), c] = S
            for c in list_selec_uneven:
                df_now.loc['m{}S_{}I-{}'.format(S,I,P), c] = I

    # Calculate shape sizes for ploting
    for n, ds_pl in itertools.product(*[range(0, len(datasets)), dfs]):
        ds_pl['{}_plot_size'.format(datasets[n])]=ds_pl['{}_frequency'.format(datasets[n])]*multiplier_size

    # df_all_shareds[df_all_shareds['p']==len_ref+3]

    # df_oxaxcx[df_oxaxcx['p']==len_ref+5]

    # df_remainders[df_remainders['p']==len_ref+7]

    # df_all_uniques[df_all_uniques['p']==len_ref+9]

    plot_dic={'shareds': df_all_shareds, 'uniques': df_all_uniques, 'oxaxcx': df_oxaxcx, 'remainders': df_remainders}

    # 'Orange' is not a valid value for cmap; supported values are 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r',
    # 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r',
    # 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r',
    # 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2',
    # 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot',
    # 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r',
    # 'copper', 'copper_r', 'crest', 'crest_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat',
    # 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r',
    # 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r',
    # 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r',
    # 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r',
    # 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r'

    Legend=False
    for ds in datasets:
        print(ds)
        plt.rc('font', size=30)

        kwargs  =   {'edgecolor' : "gray", # for edge color
                     'linewidth' : 0.3, # line width of spot
                     'linestyle' : '-', # line style of spot
                    }

        cmap = {'shareds':'Reds', 'uniques':'spring', 'oxaxcx':'Greens', 'remainders':'Blues'}

        if ds==datasets[0] or ds==datasets[-1]: fig, ax = plt.subplots(figsize=(45,8.5))
        else: fig, ax = plt.subplots(figsize=(45,8))



        for status, df_plot in plot_dic.items(): #df shareds or uniques
            df_plot = df_plot[['p', 'aa_n', 'status', '{}_plot_size'.format(ds), '{}_selec_uneven'.format(ds)]].dropna() #, '{}_frequency'.format(ds)
            # df.to_csv('{}/mode{}_{}_df.tsv'.format(output_folder, mode, ds), sep = '\t', index=True)

            # display(df_plot)
            # print(len(df_plot))
            # print(status)
            # s=np.arange(1,111)

            # xx=df_plot['{}_plot_size'.format(ds)]
            # xx=[x * 30 for x in xx]
            # print(xx, len(xx))
            ax.scatter(
                df_plot['p'], df_plot['aa_n'],
                s=list(df_plot['{}_plot_size'.format(ds)]), c=df_plot['{}_selec_uneven'.format(ds)], cmap=cmap[status], **kwargs,# color=color,, marker=marker,
                label=status,
            )

        plt.xticks([])

        if ds==datasets[0]:

            plt.tick_params(which='both', top=True, labeltop=True, bottom=False, labelbottom=False)
            # plt.tick_params(which='minor', top=True, labeltop=True, bottom=False, labelbottom=False)
            plt.xticks(ticks=range(1,len_ref+1, 1), labels=list(aa_ref_seq), minor=False)
            plt.xticks(ticks=range(2,len_ref+1, 2), minor=True)
            plt.xticks(fontsize=40)
            plt.xticks(fontsize=40)

        if ds==datasets[-1]:

            plt.tick_params(which='both', top=False, labeltop=False, bottom=True, labelbottom=True)
            # plt.tick_params(which='minor', top=True, labeltop=True, bottom=False, labelbottom=False)
            plt.xticks(ticks=range(0,len_ref+0, 5), labels=range(0,len_ref+0, 5), minor=False)
            plt.xticks(ticks=range(1,len_ref+1, 1), minor=True)
            plt.xticks(fontsize=40)

        # plt.xticks([])
        # if ds=='Passenger':
        #     plt.xticks(range(0,len_ref+1,2))
        #     # plt.minorticks_on(labelbottom=False)
        # if ds=='OVA':
        #     plt.xticks(range(1,len_ref+1), aa_ref_seq)
        #     plt.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

        plt.yticks(range(20,0,-1), aa_L)
        plt.xlim([0.5, len_ref+1])
        if Legend==True: plt.xlim([0.5, len_ref+20])

        plt.tight_layout()
        # plt.title("title")
        # plt.xlabel('Positions (1-{})'.format(len_ref), )
        # plt.xlabel('')
        plt.ylabel('Amino acids', weight='bold')
        plt.ylabel('')

        # plt.grid(color='lightgray', linestyle='-', linewidth=0.5, axis='both')
        plt.savefig('{}/mode{}-{}-{}.jpg'.format(output_folder_scatter_plots, mode, ds, chain, dpi=150))
        plt.close(fig)
        time.sleep(0.1)
        plt.pause(0.0001)


# raise Exception("Stop!")

# Section8: Sequence logos
# used in Figs 4C, 5H, 6G
output_folder_seq_logos=set_output_folder('8_seq_logos')

dfs=pd.read_csv('{}/dfs_expanded_aas_cleaned.tsv'.format(output_folder_prep), sep='\t', header=0, low_memory=False)
dfs

grouping=dfs.groupby(by=['mice', 'dataset', 'chain'])

for grouped, df in grouping:
    ID='_'.join(grouped)+'_aa'
    df_aa=seq_logo_prep(df, 'aa')
    seq_logo_plot(df_aa, output_folder_seq_logos, ID)

    # ID='_'.join(grouped)+'_nt'
    # df_nt=seq_logo_prep(df, 'nt')
    # seq_logo_plot(df_nt, output_folder_seq_logos, ID)

# raise Exception("Stop!")
print('Finished successfully!')
