# Differential Gene Expression Analysis
R code that performs various analyses of a GEO microarray dataset to investigate differential gene expression. In this case, the following dataset was used:
https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5306

Starts by loading the GEO data into R using the GEOquery package and saving that data to a variable.

Analysis begins with identifying outliers using visual methods:
- correlation plot
- heirarchical clustering dendrogram
- CV vs. mean plot
- average correlation plot

Identified outliers must be removed before continuing. (Note: in the dataset used above, no outliers were identified so code for outlier removal is not provided)

Next, genes with low expression levels are filtered out using the edgeR package. The expression level of each gene is assessed by its "counts per million" or CPM value, which is given by:

𝐶𝑃𝑀 = (# 𝑟𝑒𝑎𝑑𝑠 𝑚𝑎𝑝𝑝𝑒𝑑 𝑡𝑜 𝑔𝑒𝑛𝑒 × 10^6)/(𝑡𝑜𝑡𝑎𝑙 # 𝑚𝑎𝑝𝑝𝑒𝑑 𝑟𝑒𝑎𝑑𝑠)

Genes are only kept if they have a CPM value greater than 0.5, indicating they have a sufficiently high expression. These highly expressed genes are then subsetted to their own matrix. This matrix is now the dataset used for the analysis. 

Next, feature selection is done using the two-sample student’s t-test. In the case of the dataset linked above, the two samples are metastatic breast cancer samples and primary non-metastatic breast cancer samples. The resulting p-values are stored in a separate variable and a histogram is plotted to visually represent unadjusted p-values.

The p-values from the two-sample student's t-test are then adjusted using the Bonferroni method. This correction is done in two different ways. The first way uses the "adjust" function on the unadjusted p-values using "bonferroni" for the "method" argument, and the second way uses an alpha value of 0.05 and divides 0.05 by the number of rows in the datast following low-expression filtering.

In this case, the results of the second method were chosen for futher manipulation. From those results, genes were selected only if they had a p-values less than 0.05, meaning that the genes that fall below this threshold are statistically significant and the null hypothesis is rejected for those genes. The p-values for the selected genes and the names of the selected genes are assigned to their own variables, and a histogram of the p-value distribution of the selected genes is plotted.

These selected genes and their expresssion values are subsetted from the low-expression filtered dataset from above and assigned to a new dataframe. This new data frame is now the dataset used for the remainder of the analyses.

With the p-value selected genes now in their own data frame, dimensionality reduction is performed using classical multidimensional scaling (MDS). The "cmdscale" function is used and the results are plotted in a scatter plot.

Next, Linear Discriminant Analysis (LDA) is performed for classification of the two samples (metastatic and non-metastatic in this case). LDA is run on all of the p-value selected genes from above using the "lda" function and the "predict" function is used to predict the labels of the samples. These predictions can be visualized via a confusion matrix, which shows in a true/false format how accurately the sample labels were assigned. For example, if I know that I have 18 metastatic samples and 20 non-metastatic samples but the confusion matrix shows that the predict function assigned 16 of the samples as metastatic and 22 as non-metastatic, then I know two metastatic samples were mis-classified. 

The results of the prediction are then plotted, with each sample being assigned a different color (i.e. metastatic samples are blue points, non-metastatic samples are black points).

Lastly, the top 5 highly expressed and bottom 5 lowly-expressed genes are identified from the filtered p-value gene set (same as above). This begins with scaling the LDA results and sorting them in increasing order, where the top five genes are the most lowly expressed and the bottom five are the most highly expressed in the dataset.  





