import pandas as pd 
import numpy as np 

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


def find_top(
    df,
    value, value_thr, stat, stat_thr,
    n_line=None, drop_dup=False):
    
    # Select rows (genes) which has value >= value_thr & stat < stat_thr 
    if n_line==None:
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
        
    # found in more than n_line cell lines 
    else:
        up = df.iloc[
            [i for i,l in enumerate(
                np.array([
                    np.array(df.loc[:,value] >= value_thr).sum(axis=1) >= n_line,
                    np.array(df.loc[:,stat] < stat_thr).sum(axis=1) >= n_line]
                ).all(axis=0)) if l == 1]
                ,:]
        dn = df.iloc[
            [i for i,l in enumerate(
                np.array([
                    np.array(df.loc[:,value] <= -1*(value_thr)).sum(axis=1) >= n_line,
                    np.array(df.loc[:,stat] < stat_thr).sum(axis=1) >= n_line,]
                ).all(axis=0)) if l == 1]
                ,:] 
        
    if drop_dup==True:
        up = up.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
        dn = dn.sort_values(stat).drop_duplicates(subset='gene_id', keep="last")
    
    print ('up: ', up.shape[0])
    print ('down:', dn.shape[0])

    return up, dn


def load_data():
    # Teated vs. non-treated complete result files
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
    tables = [pd.read_csv(x, sep = '\t') for x in comparisons]
    cells = [names[i].split('_')[0] for i,x in enumerate(comparisons)]
    experiments = [names[i].replace(cells[i]+'_', '') for i,x in enumerate(comparisons)]

    # read data into Pandas data.Frames 
    data = dict(((c,{}) for c in cells))
    for i,x in enumerate(comparisons): 
        cel = cells[i]
        exp = experiments[i]
        # read data and write into dictionary 
        data [cel][exp] = tables[i]

    data['hl60']['delta_te'] = data['hl60']['delta_te'][
        ['gene_id','gene_name','Estimate_treatmentDRUG','fdr_Pr...z.._treatmentDRUG']
    ]
    data['hl60']['delta_mtyl'] = data['hl60']['delta_mtyl'][['ensembl','name','logFC','p_value']]
    data['hl60']['delta_mtyl'].columns = ['gene_id','gene_name','logFC','pval']
    # Let's include CRISPR screening scores 
    screen = pd.read_excel('screen/CRISPRi_HL60_DAC_genetable_collapsed.xlsx')

    data['hl60']['rho'] = make_score_df(screen, 'rho')
    data['hl60']['gamma'] = make_score_df(screen, 'gamma')
    
    return data


def merge_stbl_data():
    # Extract and merge experssion data for 6 AML cell lines:
    data = load_data()
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
    
    
def merge_exp_data():    
    # Extract and merge experssion data for 6 AML cell lines:
    data = load_data()
    E_gene_names = data['hl60']['delta_exp'].set_index('gene_id')[['gene_name']]
    E1 = data['hl60']['delta_exp'].set_index('gene_id')[['log2FC_120h','pval_120h']].add_prefix('hl60.')
    E2, E3, E4, E5, E6 = [
        data[cell_line]['delta_exp'].set_index('gene_id').loc[
            E1.index,
            ['log2FoldChange','pvalue']].rename(columns={'log2FoldChange':'log2FC','pvalue':'pval'}
        ).add_prefix(cell_line+'.') for cell_line in data if cell_line != 'hl60']

    exp_df = pd.concat((E_gene_names, E1,E2,E3,E4,E5,E6),axis=1)
    # exp_df.to_csv('delta_expression.txt',sep='\t')
    return exp_df


def set_Top_Stbl(fc_thr, pv_thr, n_line):
    '''
    Define top genes in delta(RNA-stability) data space: 
    > genes with log2FC >= `fc_thr` and P-Value < `pv_thr` acrross `n_line` number of cell lines. 
    '''
    print ('Subset Top Stbl data frame:')
    stbl_df = merge_stbl_data()
    out = {}
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    out['up'], out['down'] = find_top(
        stbl_df, 
        [c for c in stbl_df.columns if 'logFC' in c], fc_thr,
        [c for c in stbl_df.columns if 'pval' in c], pv_thr,
        n_line=n_line
    )
    
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}) in more than {n_line} cell lines')

    return out 


def set_Top_Exp(fc_thr, pv_thr, n_line):
    '''
    Define top genes in delta(RNA-expression) data space: 
    genes with log2FC >= `fc_thr` and P-Value < `pv_thr` acrross `n_line` number of cell lines. 
    '''
    print ('Subset Top Exp data frame:')
    exp_df = merge_exp_data()
    out = {}
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    out['up'], out['down'] = find_top(
        exp_df, 
        [c for c in exp_df.columns if 'log2FC' in c], fc_thr,
        [c for c in exp_df.columns if 'pval' in c], pv_thr,
        n_line=n_line
    )
    
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr}) in more than {n_line} cell lines')

    return out 


# RNA methylation    
def set_Top_Mtyl(fc_thr,pv_thr):
    print ('Subset Top Mtyl data frame:')
    data = load_data()
    out = {}
    out['threshold'] = [['fc_thr',fc_thr],['pv_thr',pv_thr]]
    out['up'], out['down'] = find_top(
        data['hl60']['delta_mtyl'], 
        'logFC', fc_thr, 'pval', pv_thr, drop_dup=True
    )
    
    print (f'(fc_thr={fc_thr}, pv_thr={pv_thr})' )

    return out 


# Translational efficiency
def set_Top_TE(te_thr,fdr_thr):
    print ('Subset Top TE data frame:')
    data = load_data()
    out = {}
    out['threshold'] = [['te_thr',te_thr],['fdr_thr',fdr_thr]]

    out['up'], out['down'] = find_top(
        data['hl60']['delta_te'], 'Estimate_treatmentDRUG', te_thr, 'fdr_Pr...z.._treatmentDRUG', fdr_thr
    )
    
    print (f'(te_thr={te_thr}, fdr_thr={fdr_thr})' )

    return out 


# Rho score
def set_Top_Rho(sc_thr,pv_thr):
    print ('Subset Top Rho data frame:')
    data = load_data()
    out = {}
    out['threshold'] = [['sc_thr',sc_thr],['pv_thr',pv_thr]]
    out['up'], out['down']  = find_top(
        data['hl60']['rho'], 
        'rho score', sc_thr,'Mann-Whitney p-value', pv_thr)
    
    print (f'(sc_thr={sc_thr}, pv_thr={pv_thr})')
    
    return out 


def make_final_table(genes):
    # load data
    exp_df = merge_exp_data()
    stbl_df = merge_stbl_data()
    data = load_data()
    
    # get intersects 
    def get_intersect_df(data,k=None,intersect_genes=genes):
        # change intersect_genes based on the biological question
        if k is not None: 
            df = data[k]
        else: 
            df = data

        out = df.iloc[[i for i, g in enumerate (df.gene_name) if g in list(intersect_genes)],]
        
        return out

    E = get_intersect_df(exp_df)
    S = get_intersect_df(stbl_df)
    M = get_intersect_df(data['hl60'],'delta_mtyl')
    T = get_intersect_df(data['hl60'],'delta_te')
    R = get_intersect_df(data['hl60'],'rho')
    G = get_intersect_df(data['hl60'],'gamma')

    out = pd.DataFrame(index=genes)

    out = pd.concat([
        out,
        # Expression
        E.reset_index().set_index('gene_name').iloc[:,range(1,13)].add_prefix('Exp.'),
        # Stability 
        S.reset_index().set_index('gene_name').iloc[:,range(1,13)].add_prefix('Stbl.'),
        # Translational Efficiency 
        T.reset_index().set_index('gene_name').iloc[:,[2,3]].add_prefix('TE.'),
        # CRISPR Screen Rho score
        R.reset_index().set_index('gene_name').iloc[:,[2,3]].add_prefix('Rho.'),
        # CRISPR Screen Gamma score
        G.reset_index().set_index('gene_name').iloc[:,[2,3]].add_prefix('Gamma.')

    ],axis=1)
    
    return out 