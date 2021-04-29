import os 
import pandas as pd 
import numpy as np 
from operator import itemgetter
import matplotlib.pyplot as plt


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
    if '.xlsx' in PATH:
        raw = pd.read_excel(PATH)
    else:
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


def raw_dict(cells):
    data = dict(((c,{}) for c in cells))
    return data

def load_data(comparisons=False, screens=False, wd='/rumi/shams/abe/Gilbertlab/Decitabine-treatment/'):
    '''Read data into Pandas dataframes'''
    cwd = os.getcwd()
    os.chdir(wd)
    
    data = None
    # Teated vs. non-treated complete result files
    if comparisons==True:
        comparisons = [
            # HL-60 meRIP-seq - logFC
            'meRIP-seq/hl60_delta_mtyl_table.txt' ,
            # HL-60 Ribo-seq - lnTE
            'Ribo-seq/hl60_delta_te_table.txt' ,
            # HL-60 RNA-seq 
            # RNA experssion - log2FC
            'RNA-seq/hl60-exp/hl60_delta_exp_table.txt',
            # RNA stability  - logFC
            'RNA-seq/hl60-stbl/hl60_delta_stbl_table.txt',
            ## 5 other AML cell lines RNA-seq
            # RNA experssion - log2FC
            'RNA-seq/other-exp/kg1_delta_exp_table.txt', 'RNA-seq/other-exp/molm14_delta_exp_table.txt',
            'RNA-seq/other-exp/ociaml2_delta_exp_table.txt', 'RNA-seq/other-exp/ociaml3_delta_exp_table.txt',
            'RNA-seq/other-exp/thp1_delta_exp_table.txt',
            # RNA stability - logFC
            'RNA-seq/other-stbl/kg1_delta_stbl_table.txt', 'RNA-seq/other-stbl/molm14_delta_stbl_table.txt',
            'RNA-seq/other-stbl/ociaml2_delta_stbl_table.txt', 'RNA-seq/other-stbl/ociaml3_delta_stbl_table.txt',
            'RNA-seq/other-stbl/thp1_delta_stbl_table.txt'
        ]
        # extract cell line name experiment name 
        names = [c.split('/')[-1].replace('.txt','').replace('_table','') for c in comparisons]
        cells = [names[i].split('_')[0] for i,x in enumerate(comparisons)]
        experiments = [names[i].replace(cells[i]+'_', '') for i,x in enumerate(comparisons)]

        data = raw_dict(cells)

        tables = [pd.read_csv(x, sep = '\t') for x in comparisons]
        for i,x in enumerate(comparisons): 
            cel = cells[i]
            exp = experiments[i]
            # read data and write into dictionary 
            data [cel][exp] = tables[i]

        data['hl60']['delta_te'] = data['hl60']['delta_te'][
            ['gene_id','gene_name','Estimate_treatmentDRUG','fdr_Pr...z.._treatmentDRUG']
        ].set_index('gene_id')
        data['hl60']['delta_mtyl'] = data['hl60']['delta_mtyl'][['ensembl','name','logFC','p_value']]
        data['hl60']['delta_mtyl'].columns = ['gene_id','gene_name','logFC','pval']

    # include CRISPR screening scores 
    if screens==True:
        screens = [
        'CRISPRi-screen/hl60_exp1/hl60_DAC_processing_output_genetable_collapsed.xlsx',
        'CRISPRi-screen/hl60_exp2/hl60_DAC_processing_output_genetable_collapsed.txt',
        'CRISPRi-screen/hl60_exp2/hl60_GSK_processing_output_genetable_collapsed.txt',
        'CRISPRi-screen/molm13_exp/molm13_DAC_processing_output_genetable_collapsed.txt',
        'CRISPRi-screen/molm13_exp/molm13_GSK_processing_output_genetable_collapsed.txt'
        ]

        labels = [itemgetter(0,1,3)(f.replace('CRISPRi-screen/','').replace('/','_').split('_')) for f in screens] 
        cells = set([l[0] for l in labels])
        
        if data == None:
            data = raw_dict(cells)
        
        for i,x  in enumerate(labels):
            cell, exp, drug = x
            for score in ['rho','gamma']:
                data[cell]['_'.join([exp, drug, score])] = read_genetable_collapsed (screens[i],score)
    os.chdir(cwd)
    
    return data


def merge_stbl_data(data=None):
    # Extract and merge experssion data for 6 AML cell lines:
    if data==None:
        data = load_data(comparisons=True)
    S_gene_names = data['hl60']['delta_stbl'].set_index('ensembl_id')[['gene_name']]
    S1 = data['hl60']['delta_stbl'].set_index('ensembl_id')[['logFC_120h','P.Value_120h']].rename(columns={
        'logFC_120h':'logFC',
        'P.Value_120h':'pval'}).add_prefix('hl60.')
    S2, S3, S4, S5, S6 = [
        data[cell_line]['delta_stbl'].set_index('ensembl_id').loc[
            S1.index,
            ['logFC','P.Value']].rename(columns={'P.Value':'pval'}).add_prefix(cell_line+'.') 
        for cell_line in data if cell_line != 'hl60']

    stbl_df = pd.concat((S_gene_names, S1,S2,S3,S4,S5,S6),axis=1)
    # stbl_df.to_csv('delta_stability.txt',sep='\t')
    return stbl_df
    
    
def merge_exp_data(data=None):    
    # Extract and merge experssion data for 6 AML cell lines:
    if data==None:
        data = load_data(comparisons=True)
    E_gene_names = data['hl60']['delta_exp'].set_index('gene_id')[['gene_name']]
    E1 = data['hl60']['delta_exp'].set_index('gene_id')[
        ['log2FC_120h','pval_120h']].rename(columns={"log2FC_120h": "log2FC", "pval_120h": "pval"}).add_prefix('hl60.')
    E2, E3, E4, E5, E6 = [
        data[cell_line]['delta_exp'].set_index('gene_id').loc[
            E1.index,
            ['log2FoldChange','pvalue']].rename(columns={'log2FoldChange':'log2FC','pvalue':'pval'}
        ).add_prefix(cell_line+'.') for cell_line in data if cell_line != 'hl60']

    exp_df = pd.concat((E_gene_names, E1,E2,E3,E4,E5,E6),axis=1)
    # exp_df.to_csv('delta_expression.txt',sep='\t')
    return exp_df


def set_Top_Stbl(fc_thr, pv_thr, cell_lines='hl60',data=None):
    print ('Subset Top Stbl data frame:')
    stbl_df = merge_stbl_data(data).set_index('gene_name')
    cell_lines = cell_lines.split(',')
    
    out = []
    for cl in cell_lines:
        dic = {}
        dic['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
        df = stbl_df.filter(like=cl)

        dic['up'], dic['down'] = find_top(df, f'{cl}.logFC', fc_thr,f'{cl}.pval', pv_thr)

        print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}) in {cl} cell line')
        out.append(dic)
    if len(out) == 1: 
        return out[0]
    else: 
        return out


def set_Top_Exp(fc_thr, pv_thr, cell_lines='hl60',data=None):
    
    print ('Subset Top Exp data frame:')
    exp_df = merge_exp_data(data).set_index('gene_name')
    cell_lines = cell_lines.split(',')
    
    out = []
    for cl in cell_lines:
        dic = {}
        dic['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
        df = exp_df.filter(like=cl)

        dic['up'], dic['down'] = find_top(
            df, 
            f'{cl}.log2FC', fc_thr,
            f'{cl}.pval', pv_thr
        )

        print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}) in {cl} cell line')
        out.append(dic)
    if len(out) == 1: 
        return out[0]
    else: 
        return out


# RNA methylation    
def set_Top_Mtyl(fc_thr,pv_thr,data=None):
    print ('Subset Top Mtyl data frame:')
    if data==None:
        data = load_data(comparisons=True)
    out = {}
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    out['up'], out['down'] = find_top(
        data['hl60']['delta_mtyl'], 
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
        data['hl60']['delta_te'], 'Estimate_treatmentDRUG', te_thr, 'fdr_Pr...z.._treatmentDRUG', fdr_thr
    )
    
    print (f'(te_thr={te_thr}, fdr_thr={fdr_thr})' )

    return out 


def merge_screen_data(cell, score, data=None):
    if data==None:
        data = load_data(screens=True)
    # find uniq gene names 
    genes = [data[cell][i].index.tolist() for i in data[cell].keys() if score in i]
    
    genes = set(genes[0]).intersection(*genes[1:])
    # merge data frames 
    dfs = [data[cell][i].loc[genes,:].rename(columns={
        f'{score} score': f'{cell} {i} score'.replace(' ','_'), 'Mann-Whitney p-value': f'{cell} {i} p-value'.replace(' ','_')
    }).astype(float)
           for i in data[cell].keys() if score in i]
    score_df = pd.concat(dfs, axis=1)
    return score_df

# Rho score
def set_Top_Rho(sc_thr,pv_thr,cell_line='hl60', data=None):
    print ('Subset Top Rho data frame:')
    if data==None:
        data = load_data(screens=True)
    
    dfs = []
    if cell_line=='molm13': 
        dfs.append(merge_screen_data('molm13','rho',data=data).filter(like='DAC').astype(float))
    if cell_line=='hl60':
        dfs.append(merge_screen_data('hl60','rho',data=data).filter(like='DAC').filter(like='exp1').astype(float))
#     dfs.append(merge_screen_data('hl60','rho',data=data).filter(like='DAC').filter(like='exp2').astype(float))

    top = []
    for df in dfs: 
        score,pval = df.columns
        up,dn = find_top (df,score,sc_thr,pval,pv_thr,silent=True)
        top.append([up,dn])
    
    up_genes = [up.index.tolist() for up,_ in top]
    up_genes = set(up_genes[0]).intersection(*up_genes[1:])
    dn_genes = [dn.index.tolist() for _,dn in top]
    dn_genes = set(dn_genes[0]).intersection(*dn_genes[1:])
    
    up = pd.concat([df.loc[up_genes,:] for df in dfs],axis=1)
    dn = pd.concat([df.loc[dn_genes,:] for df in dfs],axis=1)
    
    out = {}
    out['threshold'] = [['sc_thr',sc_thr],['pv_thr',pv_thr]]
    out['up'], out['down']  = up, dn
    print ('up: ', up.shape[0])
    print ('down:', dn.shape[0])

    print (f'(sc_thr={sc_thr}, pv_thr={pv_thr}) in {cell_line} cell line')
    return out 


def plot_corr(df):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    alpha = df.columns.values
    cax = ax.matshow(df.corr()) #, interpolation='nearest')
    fig.colorbar(cax)
    
    xaxis = np.arange(len(alpha))
    ax.set_xticks(xaxis)
    ax.set_yticks(xaxis)
    
    ax.set_xticklabels(alpha, rotation=45)
    ax.set_yticklabels(alpha, rotation=45)
    
    plt.show()