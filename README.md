<h1><strong>Bioinformatics Research Assistant Screening Test</strong></h1>
<h3>Applicant: Elizabeth Elias</h3>

<p>To the hiring team at the Center for Regenerative Medicine, <br>
Thank your for this opportunity to continue interviewing for this position. I'm super excited for more updates regarding the role. 
<strong>Please treat this README as the report requested by the screening test instructions.</strong>
<br> 
<h4>Task 1</h4>
Plotted with a histogram to capture the distribution of the range of gene detection captured
across the cells in the datasets. By binning cells based on the number of detected genes, 
the side-by-side histogram is able to effectively capture changes in distribution before and after
quality control and filtering was performed.
<h4>Task 2</h4>
Similar to task 1, a histogram is able to capture the distribution in the range of gene detection
captured by this dataset and the additional kernel layer density plot aids in identifying patterns
in the distribution of fragment counts by allowing the audience to observe any peaks, skews or 
multimodel behavior in the distribution. 
<h4>Task 3</h4>
Datasets were merged using dyplr's inner_join method as an inner join will only keep CellID's that are present
in both datasets, which mitigates the need for an additional filtering step if one were to use the built in merge
method. The top 5 highest regions were calculated with row sums and a scatter plot is able to effectively capture
the correlation between two variables, in this case total RNA and ATAC counts. 
<h4>Task 4</h4>
Data was scaled using scikit-learns standard scaler functionality. K means clustering was used as
after scaling, the data is well suited for comparisons which rely on Euclidean distances. 
<h4>Task 5</h4>
DESeq2 was chosen to perform the differential expression analysis as it automatically normalizes data
and provides built-in statistical testing. It is specifically designed for RNAseq data and provides the 
fold changes and p-values needed to identify the most differentially expressed genes. Presenting this data
with a volcano plot, which shows statistical significance versus fold change, highlights the key
findings in this data in an easily interpretable manner. 

<p>A note for Task 2: <ul>
<li>The <code>count_frag_filter_cells</code> function to filter the ATAC fragment counts accepts a file name
and a value for the minimum peak count for each cell in the dataset. With the example data provided, there
were no such cells that had more than 500 peaks for a given cell type. </li></ul> 

<p><strong>Python libraries used:</strong> matplotlib, scikit-learn, seaborn, pandas, uml-learn, numpy
<p><strong>R libraries used:</strong> dyplr, readr, tidyr, ggplot2, data.table, DESeq2, EnhancedVolcano
