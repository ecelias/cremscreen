import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import numpy as np
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans


filename = 'merged_atac_rna_counts.csv'
sums_filename = 'merged_atac_rna_count_sums.csv'
merged_data = pd.read_csv(filename, index_col=0)
merged_data_sums = pd.read_csv(filename, index_col=0)


def plot_pca(pca_df, x_axis='PC1', y_axis='PC2', dot_alpha=0.3, dot_size=10, use_name_to_save='plot',
             plot_points=None, use_marginals=0, remove_legend=0):

    plt.figure(figsize=(20, 20))
    g = sns.JointGrid(data=pca_df, x=x_axis, y=y_axis, height=8, ratio=5)

    # plot scatterplot on top
    g.plot_joint(sns.scatterplot, s=dot_size, edgecolor='black', alpha=dot_alpha)

    if use_marginals:
        g.plot_marginals(sns.kdeplot)

    if remove_legend:
        # remove the lefend from the joint plot
        g.ax_joint.legend_.remove()

    g.ax_joint.set_xlim(-40, 40)
    g.ax_joint.set_ylim(-40, 40)

    # adjust font size for ax_joint
    g.ax_joint.set_xlabel(x_axis, fontsize=18)
    g.ax_joint.set_ylabel(y_axis, fontsize=18)
    g.ax_joint.tick_params(axis='both', which='major', labelsize=14)

    if plot_points is not None:
        additional_x = plot_points['PC-1'].values
        additional_y = plot_points['PC-2'].values
        # plot the additional points on top of the existing JointGrid plot
        g.ax_joint.plot(additional_x, additional_y, 'ko', markersize=3)
        #'ko' means black color and circle markers

    plt.tight_layout()
    plt.savefig(use_name_to_save, transparent=True, dpi=300)
    #plt.show()


def plot_umap(umap_df, cluster_labels='clusters', dot_alpha=0.3, dot_size=10, plot_points=None,
              use_name_to_save='plot', remove_legend=0, bw_adjust=1):

    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=umap_df[:, 0], y=umap_df[:, 1], hue=cluster_labels, palette='viridis', s=60)
    plt.title('UMAP projection colored by K-means clusters')
    plt.legend(title='Cluster', loc='best', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(use_name_to_save, transparent=True, dpi=300)
    #plt.show()


def scale_data(df):
    # Extract features and drop the cell identifier column
    cell_ids = df.index  # Keep cell identifiers for labeling

    # remove any features with missing values or where all values equal 0
    df_cleaner = df.drop(columns=[col for col in df if (df[col] == 0).all()])
    df_cleaned = df_cleaner.dropna(axis=1)

    # Normalize and scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df_cleaned)

    # get the number of features in scaled dataset
    columns = df_cleaned.columns.values.tolist()
    index = df_cleaned.index.values.tolist()
    num_cols = len(columns)
    num_samples = len(index)

    return scaled_data, num_cols, num_samples


def pca_umap(df):
    df_scaled, num_cols, num_samples = scale_data(df)
    n_components = min(num_samples, num_cols)

    pca = PCA(n_components=n_components)
    principal_comps = pca.fit_transform(df_scaled)

    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)

    pc_index = 1
    for i in cumulative_explained_variance:
        print(f"PC {pc_index} = {i}")
        pc_index += 1

    umap_reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
    umap_data = umap_reducer.fit_transform(principal_comps)

    pc1 = principal_comps[:, 0]
    pc2 = principal_comps[:, 1]

    kmeans = KMeans(n_clusters=2, random_state=42)
    clusters = kmeans.fit_predict(umap_data)
    M = pd.DataFrame({'PC1': pc1, 'PC2': pc2, 'clusters': clusters})

    df['Cluster'] = clusters

    merged_data_with_clusters = df.to_csv('merged_data_with_clusters.csv')

    plot_pca(M, use_marginals=1, dot_size=30, dot_alpha=1, use_name_to_save="multiome_pca.png")
    plot_umap(umap_data, cluster_labels=clusters, dot_size=30, dot_alpha=1, use_name_to_save="multiome_umap.png")


print(pca_umap(merged_data))
