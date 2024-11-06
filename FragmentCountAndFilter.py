import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def count_frag_filter_cells(filename, min_peak_count):
    # read the csv file, retain the genomic regions as indices
    sc_df = pd.read_csv(filename, index_col=0)

    # create a "mask" data frame that converts all values > 0 to True, otherwise False
    # Sum genes by columns, only keep regions that have >= 10 True values
    # drop any genes expressed in less than 10 cells
    filter_by_reg = sc_df.sum(axis=1) >= 10
    filtered_by_reg = sc_df[filter_by_reg]

    # filter and remove all columns which sum to less than 500
    for column in filtered_by_reg.columns.tolist():
        col_sum = filtered_by_reg[column].sum()
        if col_sum < min_peak_count:
            filtered_by_reg.drop(columns=column)

    filtered_sc_df = filtered_by_reg

    # sum the filtered gene expression counts
    peak_counts = filtered_sc_df.sum(axis=0)
    peak_counts = peak_counts.sort_values(axis=0, ascending=False)
    top_5_peak_counts = peak_counts.head(5)
    plot_frag_counts(peak_counts)
    return top_5_peak_counts


def plot_frag_counts(peak_count_df):
    sns.distplot(a=peak_count_df, color='purple')
    #plt.show()
    plt.savefig('fragment_count_and_filtering.png')


# example implementation
path_to_sc_ATAC_fragment_counts = 'scATAC_fragment_counts.csv'
min_peak_count = 500
count_frag_filter_cells(path_to_sc_ATAC_fragment_counts, min_peak_count)


