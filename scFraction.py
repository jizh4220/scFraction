import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

def fraction_accCT_table(adata, goi, threshold=0.5):
    # Define a threshold for gene expression
    
    accessions = adata.obs['accession'].unique()
    celltypes = adata.obs['celltype'].unique()

    # Initialize a DataFrame to store the results
    df = pd.DataFrame(index=pd.MultiIndex.from_product([accessions, celltypes], names=['accession', 'celltype']), columns=goi)

    # Calculate the number of cells for each combination of accession and cell type
    num_cells = adata.obs.groupby(['accession', 'celltype']).size()

    # Traverse through all combinations of accession and celltype
    for accession in accessions:
        for celltype in celltypes:
            # Filter the data for the accession and celltype of interest
            adata_filtered = adata[(adata.obs['accession'] == accession) & (adata.obs['celltype'] == celltype)]

            # Filter the data for the genes of interest
            adata_filtered = adata_filtered[:, adata_filtered.var_names.isin(goi)]

            # Calculate the fraction of cells where the gene is expressed
            # Convert the sparse matrix to a dense one before applying np.where
            expressed = np.where(adata_filtered.X.toarray() >= threshold, 1, 0)
            fraction_expressed = np.sum(expressed, axis=0) / adata_filtered.X.shape[0] if adata_filtered.X.shape[0] > 0 else 0

            # Store the fraction expressed in the DataFrame
            df.loc[(accession, celltype)] = fraction_expressed

    # Assume that df is your DataFrame and adata is your AnnData object
    metadata = adata.obs
    metadata = metadata.groupby(['accession', 'celltype']).first()
    # Merge the DataFrame with the metadata
    combined_df = df.merge(metadata, left_index=True, right_index=True)
    # Add the number of cells to the DataFrame
    combined_df = combined_df.merge(num_cells.rename('num_cells'), left_index=True, right_index=True)
    return combined_df

def gene_ct_fraction_delta(combined_df, goi, groupby='celltype', x_padding=1.5, y_padding=0.8):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    diseases = combined_df[combined_df['disease'].notna()]['disease'].unique()
    ct_len = len(combined_df.index.unique('celltype'))
    df_exp_cell = combined_df[goi].multiply(combined_df['num_cells'], axis=0)

    if len(goi) >= 100:
        fig_width = len(goi) // 5
    else:
        fig_width = 14

    df_disease1 = combined_df[combined_df['disease'] == diseases[0]]
    d1_exp_cell = df_disease1[goi].multiply(df_disease1['num_cells'], axis=0)
    df_disease2 = combined_df[combined_df['disease'] == diseases[1]]
    d2_exp_cell = df_disease2[goi].multiply(df_disease2['num_cells'], axis=0)
    # Calculate the mean expression for each cell type
    df_mean1 = df_disease1[goi].apply(pd.to_numeric, errors='coerce').fillna(0).groupby(level='celltype', observed=True).mean()
    df_mean2 = df_disease2[goi].apply(pd.to_numeric, errors='coerce').fillna(0).groupby(level='celltype', observed=True).mean()
    df_diff = df_mean1 / (df_mean2 + 0.0001)

    exp_mean1 = d1_exp_cell[goi].apply(pd.to_numeric, errors='coerce').fillna(0).groupby(level='celltype', observed=True).mean()
    exp_mean2 = d2_exp_cell[goi].apply(pd.to_numeric, errors='coerce').fillna(0).groupby(level='celltype', observed=True).mean()
    exp_diff = exp_mean1 / (exp_mean2 + 0.0001) * 500

    # Create a figure
    fig, ax = plt.subplots(figsize=(fig_width, 1.2*ct_len))

    # Assuming 'result' is your DataFrame containing 'mean_expression' and 'fraction_expressing' for each gene
    for gene in goi:
        # Calculate the difference in mean expression
        mean_fraction = df_diff[gene]
        # Calculate the fraction of cells expressing the gene in each group
        mean_exp_cellNum = exp_diff[gene]

        scatter = ax.scatter([gene]*len(mean_fraction), np.arange(len(mean_fraction)), s=mean_exp_cellNum, c=mean_fraction, cmap='Reds', vmin=0, vmax=0.8)
        
    ax.set_yticks(np.arange(len(mean_fraction)))  # Set y-ticks to be the cell types
    ax.set_yticklabels(mean_fraction.index, fontsize=12)  # Set y-tick labels to be the cell types
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9) # Set x-tick labels to be goi
    ax.invert_yaxis()  # Invert y-axis

    # Set the x and y axis limits
    ax.set_xlim(-1.5, len(goi) + x_padding)
    ax.set_ylim(-1.5 + y_padding, len(mean_fraction))

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.05)
    # ax.set_title(f'Mean Gene Expression Fraction for {disease_of_interest}, Expression Threshold: {threshold}')
    plt.colorbar(scatter, label='Mean Fraction', ax=ax, cax=cax, location = 'right')
    plt.show()
    return plt

def gene_ct_fraction_plot(combined_df, goi, groupby='celltype', x_padding=1.5, y_padding=0.8):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ct_len = len(combined_df.index.unique('celltype'))
    df_exp_cell = combined_df[goi].multiply(combined_df['num_cells'], axis=0)
    # df_exp_cell = df_exp_cell / (df_exp_cell.max() + 1e-10)
    # Create a new DataFrame to store the results
    # index = pd.MultiIndex.from_product([goi, ['mean_expression', 'fraction_expressing']])
    # # Create a new DataFrame with the multi-index
    # result = pd.DataFrame(columns=index)

    # # Calculate the mean expression and fraction of cells expressing each gene
    # for gene in goi:
    #     # Calculate the mean expression of the gene in each group
    #     mean_fraction = combined_df.groupby(groupby)[gene].mean()
    #     # Calculate the fraction of cells expressing the gene in each group
    #     mean_exp_cellNum = df_exp_cell.groupby(groupby)[gene].mean()
    #     # Store the results in the DataFrame
    #     result[(gene, 'mean_fraction')] = mean_fraction
    #     result[(gene, 'mean_exp_cellNum')] = mean_exp_cellNum

    if len(goi) >= 100:
        fig_width = len(goi) // 5
    else:
        fig_width = 14

    # Create a figure
    fig, ax = plt.subplots(figsize=(fig_width, 1.2*ct_len))

    # Assuming 'result' is your DataFrame containing 'mean_expression' and 'fraction_expressing' for each gene
    for gene in goi:
        mean_fraction = combined_df.groupby(groupby, observed=True)[gene].mean().fillna(0)
        # Calculate the fraction of cells expressing the gene in each group
        mean_exp_cellNum = df_exp_cell.groupby(groupby, observed=True)[gene].mean().fillna(0)
       
        # print(mean_fraction)
        # mean_exp_cellNum = result[(gene, 'mean_exp_cellNum')] / (result[(gene, 'mean_exp_cellNum')].max() + 1e-10)
        
        # print(mean_exp_cellNum)
        # # Replace any NaN values with 0
        # mean_fraction = mean_fraction
        # mean_exp_cellNum = mean_exp_cellNum.fillna(0)

        # The size and color of the dots are determined by the mean expression and fraction of cells expressing the gene
        # scatter = ax.scatter([gene]*len(mean_fraction), np.arange(len(mean_fraction)), s=mean_fraction*100, c=mean_exp_cellNum, cmap='Reds', vmin=0, vmax=1)
        scatter = ax.scatter([gene]*len(mean_fraction), np.arange(len(mean_fraction)), s=mean_exp_cellNum, c=mean_fraction, cmap='Reds', vmin=0, vmax=0.5)
        
    ax.set_yticks(np.arange(len(mean_fraction)))  # Set y-ticks to be the cell types
    ax.set_yticklabels(mean_fraction.index, fontsize=12)  # Set y-tick labels to be the cell types
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9) # Set x-tick labels to be goi
    ax.invert_yaxis()  # Invert y-axis

    # Set the x and y axis limits
    ax.set_xlim(-1.5, len(goi) + x_padding)
    ax.set_ylim(-1.5 + y_padding, len(mean_fraction))

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="1%", pad=0.05)
    # ax.set_title(f'Mean Gene Expression Fraction for {disease_of_interest}, Expression Threshold: {threshold}')
    plt.colorbar(scatter, label='Mean Fraction', ax=ax, cax=cax, location = 'right')
    plt.show()
    return plt

def refine_metanames(meta):
    rename_dict = {
        "orig.ident": "accession",
        "gse_alias": "experiment",
        "DISEASE": "disease",
        "TISSUE": "tissue",
        "celltype": "celltype"
        # "macro_celltype": "celltype"
    }
    
    # Only rename columns that exist in the dataframe and the new name does not exist
    rename_dict = {k: v for k, v in rename_dict.items() if k in meta.columns and v not in meta.columns}
    meta.rename(columns=rename_dict, inplace=True)
    return meta


def main(adata_path, goi_path):
    # Load the data
    adata = sc.read_h5ad(adata_path)

    # Load genes of interest from the text file
    with open(goi_path, 'r') as f:
        goi = [line.strip() for line in f]
    # make sure you have to accession, experiment, celltype, and disease in your adata
    clean_meta = refine_metanames(adata.obs)

    # Check if 'accession' and 'experiment' exist in clean_meta
    try:
        assert 'accession' in clean_meta.columns, "'accession' not found in clean_meta columns"
        assert 'experiment' in clean_meta.columns, "'experiment' not found in clean_meta columns"
        assert 'celltype' in clean_meta.columns, "'celltype' not found in clean_meta columns"
        assert 'disease' in clean_meta.columns, "'disease' not found in clean_meta columns"
    except AssertionError as e:
        print(e)
        return

    adata.obs = clean_meta
    # where we get the fraction table for downstream visualization
    frac_table = fraction_accCT_table(adata, goi=goi, threshold=0.5)

    ps = gene_ct_fraction_plot(frac_table, goi)
    # Display all plots
    for fig in ps:
        fig.show()


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])