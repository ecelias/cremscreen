<h1><strong>Bioinformatics Research Assistant Screening Test</strong></h1>
<h3>Applicant: Elizabeth Elias</h3>

<p>To the hiring team at the Center for Regenerative Medicine, <br>
Thank your for this opportunity to continue interviewing for this position. I'm super excited for more updates regarding the role. 
<br> <br>
<p>All scripts can be run simultaneously from <code>RunScripts.py</code>
<p>Some notes for the scripts in this repo: <ul>
<li>All plots will be saved as .png files. A commented-out line can be uncommented to show all plots as the 
scripts run in each of the functions. If you make this change, <code>RunScripts.py</code> should not be used as it
uses the <code>call</code> function, which is blocking.</li>
<li>The <code>count_frag_filter_cells</code> function to filter the ATAC fragment counts accepts a file name
and a value for the minimum peak count for each cell in the dataset. With the example data provided, there
were no such cells that had more than 500 peaks for a given cell type. </li></ul> 

<p><strong>Python libraries used:</strong> matplotlib, scikit-learn, seaborn, pandas
<p><strong>R libraries used:</strong> dyplr, readr, tidyr, ggplot2, data.table
