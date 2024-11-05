import pandas as pd
import matplotlib.pyplot as plt


def filter_sc(filename):
    # read the csv file, retain the gene names as indices
    sc_df = pd.read_csv(filename, index_col=0)

    # sum the initial gene expression counts
    initial_count = sc_df.sum(axis=0)

    # create a "mask" data frame that converts all values > 0 to True, otherwise False
    mask_by_gex = (sc_df > 0)
    # Sum genes by columns, only keep genes that have >= 5 True values
    filter_by_gex = mask_by_gex.sum(axis=1) >= 5
    # drop any genes expressed in less than 5 cells
    filtered_by_gex = sc_df[filter_by_gex]

    # create a "mask" data frame that converts all values > 0 to True, otherwise False
    mask_by_cell = (filtered_by_gex > 0)
    # sum cells by row, only keep cells with >= 200 True values
    filter_by_cell = mask_by_cell.sum(axis=0) >= 200
    filtered_sc_df = filtered_by_gex.loc[:, filter_by_cell]

    # sum the filtered gene expression counts
    filtered_count = filtered_sc_df.sum(axis=0)

    return plot_data(initial_count, filtered_count)


def plot_data(initial_count, filtered_count):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
    axes[0].hist(initial_count, color='Green', edgecolor='black')
    axes[0].set_title('Detected Genes Per Cell Before Filtering')

    axes[1].hist(filtered_count, color='Pink', edgecolor='black')
    axes[1].set_title('Detected Genes Per Cell After Filtering')
    for ax in axes:
        ax.set_xlabel('Cells')
        ax.set_ylabel('Gene Expression Counts')

    plt.tight_layout()
    plt.savefig('quality_control_and_filtering.png')


# example implementation
path_to_sc_rna_expression_data = 'scRNA_expression.csv'
filter_sc(path_to_sc_rna_expression_data)