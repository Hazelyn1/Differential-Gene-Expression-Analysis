#Hazelyn Cates
#8/8/23
#AS.410.671.81.SU23

#Dataset used:
#https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5306 “Human epidermal growth factor receptor 2-positive breast cancer brain metastases”

#Please see end of powerpoint presentnation for all references
#............................................................................................................................................................

#Using GEOquery to load in the data:
#Code:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install("GEOquery") 
library(GEOquery)

data_file <- getGEO(GEO = 'GDS5306') #replace with dataset ID of choice

exprs_data <- data_file@dataTable@table

dim(exprs_data)
[1] 61359    40
#Columns 1 and 2 contain the ID reference and the gene identifier, respectively. Columns 3-20 contain the metastatic samples, columns 21-40 contain the #nonmetastatic samples  

#Retrieving gene names:
exprs_data_names <- exprs_data$IDENTIFIER 
............................................................................................................................................................

#Outlier analysis: 
#Starting with correlation plot:
#First, get correlation matrix sing Pearson’s correlation:
exprs_data_cor <- cor( exprs_data[, 3:40], use = "pairwise.complete.obs", method = "pearson")

#Now plot correlation plot:
corrplot(exprs_data_cor, method = 'color', mar=c(0,0,1,0), tl.cex = 0.5)
title("Expression levels of brain metastatic breast cancer samples\n and primary breast cancer samples", line = 1.7)
...........

#Second, a hierarchical clustering dendrogram:
exprs_data_t <- t(exprs_data[, 3:40]) #transpose data frame
exprs_data_dist <- dist(exprs_data_t, method = "euclidean") #calculate distance b/w data points
exprs_data_clust <- hclust(exprs_data_dist, method = "single") #form clusters within the data

#Plot the dendrogram:
plot(exprs_data_clust, main = "Dendrogram of Metastatic and non-metastatic HER2+ breast cancer samples", xlab = "Distance")
...........

#Third, a CV vs. mean plot:
exprs_data_mean <- apply(log2(exprs_data[, 3:40]), 2, mean) #calculate means across the columns of the expression  data
exprs_data_stddev <- sqrt(apply(log2(exprs_data[, 3:40]), 2, var)) #calculate standard deviation of the samples by column
exprs_data_CV <- exprs_data_stddev/exprs_data_mean #calculate coefficient of variation

#Cast the mean and CV variables to data frames:
exprs_data_meandf <- as.data.frame(exprs_data_mean)
exprs_data_CVdf <- as.data.frame(exprs_data_CV)

#Combine into a single data frame:
exprs_data_df <- data.frame(exprs_data_meandf, exprs_data_CVdf) 

#Rename columns to keep track:
colnames(exprs_data_df) <- c("Mean", "CV")

#Plotting the data frame containing the mean and CV values, using “ggrepel” within ggplot to avoid label overlap: 
exprs_data_lab <- colnames(exprs_data[3:40]) #point labels

library(ggplot2)
library(ggrepel)

ggplot(exprs_data_df, aes(x = Mean, y = CV, label = exprs_data_lab)) +geom_point(color = "lightblue") + geom_text_repel() + labs(title = "HER2+ metastatic and nonmetastatic breast cancer samples\nCV vs. Mean") + theme(plot.title = element_text(hjust = 0.5))

...........
#Lastly, an average correlation plot:
#First, get mean values from correlation results found above:
exprs_data_avg <- apply(exprs_data_cor, 1, mean) #going by row

#Plotting:
plot(c(1, length(exprs_data_avg)), range(exprs_data_avg), type = "n", xlab = "", ylab = "Average r", main = "Average correlation of metastatic and nonmetastatic HER2+ breast cancer samples", axes= F)

#Adding the points:
points(exprs_data_avg, col = "purple", pch = 16, cex = 1)

#Adding the axes:
axis(1, at = c(1:length(exprs_data_avg)), labels = exprs_data_lab, las = 2, cex.lab = 0.4, cex.axis = 0.6)
axis(2)

#Adding vertical separator lines to help visualization:
abline(v = seq(0.5, 62.5, 1,), col = "black")
............................................................................................................................................................

Filtering out genes w/ low expression: 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR") 

library(edgeR)

#Calculating counts-per-million  for all the genes in each sample:
exprs_data_CPM <- cpm(exprs_data[, 3:40])

#The idea is to only keep genes that have a CPM above a certain threshold. In this case, the threshold is going to be CPM = 0.5.
CPM_threshold <- exprs_data_CPM > 0.5 

#To keep track of the genes, assign the row names the gene names:
rownames(CPM_threshold) <- exprs_data_names

#The samples that have a “TRUE” means they have a CPM > 0.5 and therefore have sufficiently high expression. All the others (“FALSE”) are considered to be low expression and can therefore be filtered out.
#This is done as follows:
#First by only including the genes that have at least nine “TRUE” across their samples (columns), meaning that for that gene, it was at least sufficiently expressed in nine (~1/4) of the 38 samples:
#By getting the row sums of “CPM_threshold”, we know how many of the rows have at least x amount of “TRUE” values. 
table(rowSums(CPM_threshold))

#Now given that, we can choose to only keep those rows who have at least (>=9) “TRUE” values:
CPM_keep <- rowSums(CPM_threshold) >= 9
summary(CPM_keep)

#Only want to keep the 58129 “TRUE” genes, so put ONLY those samples into their own matrix:
exprs_data_keep <- exprs_data[CPM_keep, ]

#So now, for further downstream analysis, I’m going to be using the “exprs_data_keep” dataset, since it has the low-expression genes filtered out. And just like usual, start at column 3 for the expression data

#Keep track of which genes these correspond to:
keep_genes <- list() 
for(i in 1:length(CPM_keep)){
     if(CPM_keep[i] == "TRUE"){
         keep_genes <- append(keep_genes, names(CPM_keep[i]))
     }
}
............................................................................................................................................................

#Next, perform a statistical test for feature selection: 
#Using a two-sample student’s t-test for all the genes remaining following the above filtering (now using the “exprs_data_keep” variable) and extracting the p-values:
#Log base 2-transform the data:
exprs_data_keep_log2 <- log2(exprs_data_keep[, 3:40])

s1_samples <- names(exprs_data_keep[3:21]) #metastatic samples
s2_samples <- names(exprs_data_keep[22:40]) #primary breast cancer samples

#Perform two-sample Student’s t-test:
t_test_all_genes <- function(x,s1,s2) {
    x1 <- x[s1]
    x2 <- x[s2]
    x1 <- as.numeric(x1)
    x2 <- as.numeric(x2)
    t_output <- t.test(x1,x2, alternative="two.sided",var.equal=T)
    output <- as.numeric(t_output$p.value)
    return(output)
}
#Calling the function:
exprs_data_keep_pv <- apply(exprs_data_keep_log2, 1, t_test_all_genes, s1_samples, s2_samples)

#Plotting a histogram of the unadjusted p-values:
hist(exprs_data_keep_pv, col = "lightblue", xlab = "p-values", main = "P-value distribution between metastatic and non-metastatic primary HER2+ breast cancer samples")
............................................................................................................................................................

#Next, adjust the p-values for multiplicity using the Bonferroni method: 
exprs_data_keep_pv_adj <- p.adjust(exprs_data_keep_pv, method = "bonferroni")

#Bonferroni correction done a different way, using an alpha value = 0.05: 
exprs_data_keep_BFpv <- 0.05 / nrow(exprs_data_keep)

#Finidng if there are any p-values that are less than the Bonferroni correction value:
BF_pvals <- sum(exprs_data_keep_pv < exprs_data_keep_BFpv)

#Choosing genes that only have a p-value less than 0.05, meaning that the genes that fall below this threshold are statistically significant and the null hypothesis is rejected for these genes: 
pv_0.05 <- sum(exprs_data_keep_pv < 0.05)

#Use cbind to keep track of which gene names correspond to which p-value:
exprs_data_keep_pv_df <- cbind(exprs_data_keep[1:2], exprs_data_keep_pv)

#Create lists that will hold the p-values < 0.05 and the respective gene:
p_vals <- list()
pv_genes <- list()
exprs_data_keep_names <- rownames(exprs_data_keep)

for(i in 1:length(exprs_data_keep_pv)){
    if(exprs_data_keep_pv[i] < 0.05){
p_vals <- append(p_vals, exprs_data_keep_pv[i])
		pv_genes <- append(pv_genes, exprs_data_keep_names[i])
    }
}
pv_genes <- as.numeric(pv_genes)

#Histogram of the p-value distribution of those 6,190 genes:
hist(as.numeric(unlist(p_vals)), col = "lightblue", xlab = "p-values", main = "P-value distribution of genes with a p-value < 0.05")
............................................................................................................................................................

#Subset the 6,190 genes from my “exprs_data_keep” variable, the latter which has the 58000+ genes. So, I need to grab out all 38 columns from each row where the genes match up.
exprs_data_keep_names <- rownames(exprs_data_keep)

#use the row numbers from the “pv_genes” variable instead of the gene name to extract the correct rows in “exprs_data_keep”
genes_data <- data.frame() #to hold the subset of genes and corresponding expression data for the 6,190 genes determine by the p-value analysis

for(i in 1:length(keep_genes)){ #58000+ iterations
   for(j in 1:length(pv_genes)){ #6190 iterations
      if(pv_genes[j] == exprs_data_keep_names[i]){
            gene_row <- exprs_data_keep_names[i]
            #print(gene_row)
            #print(i)
            #print(pv_genes[j])
            #print(exprs_data_keep[gene_row, ])
            #print(exprs_data_keep[(exprs_data_keep_names == i), ])
            genes_data <- rbind(genes_data, exprs_data_keep[gene_row, ])
            #print(genes_data)
       }  
  }   
}
............................................................................................................................................................

#Next, performing dimensionality reduction using classical multidimensional scaling (MDS): 
#”genes_data” is the new data frame that contains the 6,190 genes
genes_data_dist <- dist(t(genes_data[, 3:40])) 

#Metric scaling:
genes_data_loc <- cmdscale(genes_data_dist)

#Plotting results:
plot(genes_data_loc, type = "n", main = "Classical MDS")
 
#Adding points:
points(genes_data_loc[1:19, 1], genes_data_loc[1:19, 2], pch = 16, col = "blue", cex = 1.5) #metastatic results
points(genes_data_loc[20:38, 1], genes_data_loc[20:38, 2], pch = 16, col = "green", cex = 1.5) #non-metastatic results

#Adding legend:
legend("topleft", legend = c("metastatic samples", "primary breast cancer samples"), col = c("blue", "green"), pch = 16, cex = 0.5)
............................................................................................................................................................

#Next, classification via LDA:
genes_data_clas1 <- rep("metastatic", 18)
genes_data_clas2 <- rep("non-metastatic", 20)
genes_data_clas <- append(genes_data_clas1, genes_data_clas2)

genes_data_t <- t(genes_data[, 3:40]) #transpose matrix
genes_data_t <- data.frame(genes_data_clas, genes_data_t) #cast to data frame

#Run LDA on ALL genes:
genes_data_lda <- lda(genes_data_clas~., genes_data_t)

#Prediction:
genes_data_predict <- predict(genes_data_lda, genes_data_t)

#Accessing resulting confusion matrix:
table(genes_data_predict$class,genes_data_clas)

#Plotting all genes following LDA: 
plot(genes_data_predict$x, pch = 16, col = 1, ylab = "Discriminant function", xlab = "Score", main = "Linear discriminant function for metastatic and non-metastatic\nHER2+ breast tumor samples", axes = F, bg=as.numeric(factor(class_names)))
points(genes_data_predict$x[1:18, ], col = "blue", pch = 16)
axis(1, at = c(1:38), names(genes_data[3:40]), las=2, cex.axis=0.4)
axis(2)
legend("topleft", col = c("blue", "black"), pch = 16, cex = 0.6, legend = c("Metastatic", "Non-metastatic"))
............................................................................................................................................................

#Getting the 10 discriminant genes:
#First, obtain scaling data from lda results:
genes_data_lda_dis <- genes_data_lda$scaling

#Sort the data in increasing order:
genes_data_lda_sorted <- sort(genes_data_lda_dis, decreasing = F)

#Getting the top five negative discriminant values:
genes_data_min <- genes_data_lda_sorted[1:5]

#Getting the top five positive discriminant values:
genes_data_max <- genes_data_lda_sorted[6185:6190]

#Match up the scaling values to the correct gene, which at this point is just a number:
rownames(genes_data_lda_dis) <- pv_genes
genes_data_lda_names <- rownames(genes_data_lda_dis) #assign the corresponding gene number from "pv_genes"

lda_genes_min <- list() #empty list to hold the matching gene row number for the top 5 negative discriminant values
lda_genes_max <- list() #empty list to hold the matching gene row number fo the top 5 positive discriminant values

for(i in 1:length(genes_data_min)){ #this iterates 5 times
	for(j in 1:nrow(genes_data_lda_dis)){ #this iterates 6,190 times
		if(genes_data_min[i] == genes_data_lda_dis[j]){
			lda_genes_min <- apppend(lda_genes_min, genes_data_lda_names[j])
		}
		if(genes_data_max[i] == genes_data_lda_dis[j]){
			lda_genes_max <- append(lda_genes_max, genes_data_lda_names[j])
		}
	}
}

#Now extract the actual gene name from the data frame that contains the expression data of the filtered 6,190 genes:
dis_genes <- list() #empty list to hold the gene names of all 10 genes

for(i in 1:length(genes_data_min)){ #iterates 5 times
	for(j in 1:nrow(genes_data)){ #iterates 6,190 times
		if(lda_genes_min[i] == genes_data_names[j]){
			temp <- genes_data_names[j]
			dis_genes <- append(dis_genes, genes_data[temp, 2])
		}
		if(lda_genes_max[i] == genes_data_names[j]){
			temp <- genes_data_names[j]
			dis_genes <- append(dis_genes, genes_data[temp, 2])
		}
	}
}





