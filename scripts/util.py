import os 
import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
from mpl_toolkits.axes_grid1 import make_axes_locatable
from operator import itemgetter

def cleanAxes(axis, top=False, right=False, bottom=True, left=True):
    # adapted from http://nbviewer.ipython.org/github/cs109/content/blob/master/lec_03_statistical_graphs.ipynb    
    axis.spines['top'].set_visible(top)
    axis.spines['right'].set_visible(right)
    axis.spines['left'].set_visible(left)
    axis.spines['bottom'].set_visible(bottom)

    #turn off all ticks
    axis.yaxis.set_ticks_position('none')
    axis.xaxis.set_ticks_position('none')

    #now re-enable visibles
    if top:
        axis.xaxis.tick_top()
    if bottom:
        axis.xaxis.tick_bottom()
    if left:
        axis.yaxis.tick_left()
    if right:
        axis.yaxis.tick_right()

        
def make_score_df(screen, score, rep='ave_Rep1_Rep2'):
    '''
    Make dataframe for given score from CRISPR screening Excel file 
    '''

    score_cols = [i for i,c in enumerate(screen.columns) if score in c]
    # extract screen informations from Excel file 
    screen_info = list(zip(
        score_cols,
        screen.iloc[0, score_cols].tolist(), 
        screen.iloc[1, score_cols].tolist(), 
        screen.columns[score_cols].tolist() 
    ))
    
    # select columns for given score
    cols = [c[0] for c in screen_info if c[1] == rep]
    # remove top 3 rows 
    df = screen.iloc[3:,cols]

    # set dataFrame column names from row 1
    df.columns = screen.iloc[1,cols].tolist()
    # add gene_name column
    gene_names = [str(g) for g in screen.iloc[3:,0]]
    df.insert (0, "gene_name", [g for g,t in zip(gene_names,df.transcripts)])
    # select rows with gene names  
    df = df.iloc[[i for i,g in enumerate(df.gene_name) if 'pseudo_' not in g],:].reset_index(drop=True)

    # only keep one transcript name
    df.transcripts = [str(g).split(',')[0] for g in df.transcripts]    
    # rename phenotype score column to the score name
    df.columns = [score + ' score' if x=='average phenotype of strongest 3' else x for x in df.columns]
    
    # print (df.shape)
    return df

def read_genetable_collapsed(PATH, SCORE):
#     if '.xlsx' in PATH:
#         raw = pd.read_excel(PATH)
#     else:
    raw = pd.read_csv(PATH,sep='\t', low_memory=False) 
    df = make_score_df(raw, SCORE).set_index('gene_name').drop_duplicates()
    out = df.loc[:,[f'{SCORE} score','Mann-Whitney p-value']]
    
    return out


def find_top(df,value, value_thr, stat, stat_thr, drop_dup=False,silent=False):
    # Select rows (genes) which has value >= value_thr & stat < stat_thr 
#     if n_line==None:
    up = df.iloc[
        [i for i,l in enumerate(
            np.array([
                np.array(df.loc[:,value] >= value_thr),
                np.array(df.loc[:,stat] < stat_thr)]).all(axis=0)) if l == 1]
            ,:] 
    dn = df.iloc[
        [i for i,l in enumerate(
            np.array([
                np.array(df.loc[:,value] <= -1*(value_thr)),
                np.array(df.loc[:,stat] < stat_thr)]).all(axis=0)) if l == 1]
            ,:] 
        
    if drop_dup==True:
        up = up.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
        dn = dn.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
    
    if not silent:
        print ('up: ', up.shape[0])
        print ('down:', dn.shape[0])

    return up, dn


# def raw_dict(cells):
#     data = dict(((c,{}) for c in cells))
#     return data

def load_data(comparisons=False, screens=False, wd='/rumi/shams/abe/AML/Decitabine-treatment/'):
    '''Read data into Pandas dataframes'''
    cwd = os.getcwd()
    os.chdir(wd)
    
    data = None
    # Teated vs. non-treated complete result files
    if comparisons==True:
        data = dict()
        comparisons = [
            # HL-60 meRIP-seq - logFC
            'meRIP-seq/hl60_delta_mtyl_table.txt' ,
            # HL-60 Ribo-seq - lnTE
            'Ribo-seq/hl60_delta_te_table.txt' ,
            # x6 cell lines RNA experssion - log2FC
            'RNA-seq/exp/delta_exp_table.txt',
            # x6 cell lines RNA stability  - log2FC
            'RNA-seq/stbl/delta_stbl_table.txt',
        ]
        
        # extract experiment name 
        names = [
            c.split('/')[-1].replace('.txt','').replace('_table','').replace('hl60_','') 
            for c in comparisons
        ]
        
        # read data and write into dictionary 
        tables = [pd.read_csv(x, sep = '\t') for x in comparisons]
        for i,x in enumerate(comparisons): 
            data [names[i]] = tables[i]
        data['delta_te'] = data['delta_te'][
            ['gene_id','gene_name','Estimate_treatmentDRUG','fdr_Pr...z.._treatmentDRUG']
        ].set_index('gene_id')
        data['delta_mtyl'] = data['delta_mtyl'][['ensembl','name','logFC','p_value']]
        data['delta_mtyl'].columns = ['gene_id','gene_name','logFC','pval']

    # include CRISPR screening scores 
    if screens==True:
        if data == None:
            data = dict()

        screens = [
        'CRISPRi-screen/hl60_exp1/genetable_collapsed.txt',
        'CRISPRi-screen/molm13_exp/genetable_collapsed.txt',
        'CRISPRi-screen/molm13_exp/genetable_collapsed.txt',
        'CRISPRi-screen/skm1_exp/genetable_collapsed.txt',
        ]
        
        labels = [f.split('/')[1] for f in screens]
        
        for i,label in enumerate(labels):
            for score in ['rho','gamma','tau']:
                data[f'{label}_{score}'] = read_genetable_collapsed (screens[i],score)
    os.chdir(cwd)
    
    return data


def set_Top_Exp(fc_thr, pv_thr, data=None,comp='hl60_72h_only'):
    print ('Subset Top Exp data frame:')
    fc = f'{comp}_log2FC'
    pv  = f'{comp}_pvalue'
    
    out = {}
    out['up'], out['down'] = find_top(
        data['delta_exp'].loc[:,['gene_name',fc,pv]], 
        fc, fc_thr,pv, pv_thr)
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    print (f'({comp})')
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}')

    return out 


def set_Top_Stbl(fc_thr, pv_thr, data=None,comp='hl60_72h'):
    print ('Subset Top Stbl data frame:')
    fc = f'{comp}_log2FC'
    pv  = f'{comp}_pvalue'
    
    out = {}
    out['up'], out['down'] = find_top(
        data['delta_stbl'].loc[:,['gene_name',fc,pv]], 
        fc, fc_thr,pv, pv_thr)
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    print (f'({comp})')
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}')

    return out 


def two_sided_mtyl(delta_mtyl, fcthr=1,pvthr=0.01):
    '''Read meRIP-seq data into two data frames'''
    ### hyper_methylation gene list
    # subset by threshold 
    hyper = delta_mtyl.iloc[
        np.where([(l and p) for l,p in zip(delta_mtyl.logFC >= fcthr,delta_mtyl.p_value < pvthr)])
    ]
    ### hypo_methylation gene list
    # subset by threshold 
    hypo = delta_mtyl.iloc[
        np.where([(l and p) for l,p in zip(delta_mtyl.logFC <= -(fcthr),delta_mtyl.p_value < pvthr)])
    ]
    
    return hyper, hypo


# RNA methylation    
def set_Top_Mtyl(fc_thr,pv_thr,data=None):
    print ('Subset Top Mtyl data frame:')
    if data==None:
        data = load_data(comparisons=True)
    out = {}
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    out['up'], out['down'] = find_top(
        data['delta_mtyl'], 
        'logFC', fc_thr, 'pval', pv_thr, drop_dup=True
    )
    
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr})' )

    return out 


# Translational efficiency
def set_Top_TE(te_thr,fdr_thr,data=None):
    print ('Subset Top TE data frame:')
    if data==None:
        data = load_data(comparisons=True)
    out = {}
    out['threshold'] = [['te_thr',te_thr],['fdr_thr',fdr_thr]]

    out['up'], out['down'] = find_top(
        data['delta_te'], 'Estimate_treatmentDRUG', te_thr, 'fdr_Pr...z.._treatmentDRUG', fdr_thr
    )
    
    print (f'(te_thr={te_thr}, fdr_thr={fdr_thr})' )

    return out 


def merge_screen_data(cells, exps, scores, data=None):
    if data==None:
        data = load_data(screens=True)

    # find uniq gene names 
    genes = []
    for key in data.keys():
        c, e, s = key.split('_')
        if c in cells and e in exps and s in scores:
            genes.append(
                data[key].index.tolist() 
            )
    
    genes = list(set(genes[0]).intersection(*genes[1:]))
    
    # merge data frames 
    dfs = []
    for key in data.keys():
        c, e, s = key.split('_')
        if c in cells and e in exps and s in scores:
            dfs.append(data[key].loc[genes,:].rename(
                columns={
                    f'{s} score': f'{key} score'#.replace(' ','_')
                    , 
                    'Mann-Whitney p-value': f'{key} pvalue'#.replace(' ','_')
                }).astype(float))

    score_df = pd.concat(dfs, axis=1)
    return score_df


# Rho score
def set_Top_Rho(sc_thr,pv_thr,cells='hl60', exps='exp1', data=None):
    print ('Subset Top Rho data frame:')
    if data==None:
        data = load_data(screens=True)
    
    dfs = []
    
    for key in data.keys():
        c, e, score = key.split('_')
        if c in cells and e in exps and score == 'rho':
            print (f'({c} / {e})')
            dfs.append(
                merge_screen_data(c, e, score, data)
            )

    top = []
    for df in dfs: 
        score,pval = df.columns
        up,dn = find_top (df,score,sc_thr,pval,pv_thr,silent=True)
        print (f'(sc_thr={sc_thr}, pv_thr={pv_thr})')
        print ('up: ', up.shape[0])
        print ('down:', dn.shape[0])
        top.append([up,dn])
    
    up_genes = [up.index.tolist() for up,_ in top]
    up_genes = list(set(up_genes[0]).intersection(*up_genes[1:]))
    dn_genes = [dn.index.tolist() for _,dn in top]
    dn_genes = list(set(dn_genes[0]).intersection(*dn_genes[1:]))
    
    up = pd.concat([df.loc[up_genes,:] for df in dfs],axis=1)
    dn = pd.concat([df.loc[dn_genes,:] for df in dfs],axis=1)
    
    out = {}
    out['threshold'] = [['sc_thr',sc_thr],['pv_thr',pv_thr]]
    out['up'], out['down']  = up, dn

    return out 


def plot_corr(df,title,vmin=None,vmax=None,sub=111):
    fig = plt.figure()
    ax = fig.add_subplot(sub)
    
    alpha = df.columns.values
    cax = ax.matshow(df.corr()) #, interpolation='nearest')
    fig.colorbar(cax)
    cax.set_clim(vmin,vmax)
    
    xaxis = np.arange(len(alpha))
    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    ax.set_xticklabels(alpha, rotation=45)
    ax.set_yticklabels(alpha, rotation=45)
    ax.set_title(title, fontsize=15)
    
    
    plt.show()