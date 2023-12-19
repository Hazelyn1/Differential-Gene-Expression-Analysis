# Differential Gene Expression Analysis
R code that performs various analyses of a GEO microarray dataset to investigate differential gene expression. In this case, the following dataset was used:
https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5306

Starts by loading the GEO data into R using the GEOquery package and saving that data to a variable.

Analysis begins with identifying outliers using visual methods:
- correlation plot
- heirarchical clustering dendrogram
- CV vs. mean plot
- average correlation plot
And identified outliers must be removed before continuing

Next, genes with low expression levels are filtered out using the edgeR package. The expression level of each gene is assessed by its "counts per million" or CPM value, which is given by:

𝐶𝑃𝑀 = (# 𝑟𝑒𝑎𝑑𝑠 𝑚𝑎𝑝𝑝𝑒𝑑 𝑡𝑜 𝑔𝑒𝑛𝑒 × 10^6)/(𝑡𝑜𝑡𝑎𝑙 # 𝑚𝑎𝑝𝑝𝑒𝑑 𝑟𝑒𝑎𝑑𝑠)


