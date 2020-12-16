
#############################################################################
# Kernel Maximum Mean Discrepancy 
# "INFERRING DIFFERENTIALLY EXPRESSED PATHWAYS BY USING KERNEL MAXIMUM MEAN DISCREPANCY-BASED TEST
# Esteban Vegas, Ferran Reverter and Josep Maria Oller
# Statistical Department, University of Barcelona
#############################################################################


#############################################################################
# Nutrimouse data set
#############################################################################


if(!require(mixOmics, quietly =T))  install.packages("mixOmics")
require(mixOmics, quietly =T)

if(!require(methods, quietly =T))  install.packages("mixOmics")
require(mixOmics, quietly =T)


data(nutrimouse)
#help(nutrimouse)
#str(nutrimouse)

mydata.gene <- nutrimouse$gene    # gene expressions data set
mydata.FA <- nutrimouse$lipid  # fatty acids      data set

gene.names <- names(mydata.gene)
FA.names <- names(mydata.FA)

ind.diet <- nutrimouse$diet
ind.genotype <- nutrimouse$genotype

# There are 2 factors: genotype with 2-levels factor and diet 5-levels factor.

# Change from genotype factor to genotype list 
ind.gen <- list()
for (i in levels(ind.genotype)) {
  ind.gen[[i]] <- which(ind.genotype==i)
}

# Change from diet factor to diet list 
ind.diet1 <- list()
for (i in levels(ind.diet)) {
  ind.diet1[[i]] <- which(ind.diet==i)
}

#############################################################################
# Subset of genes and fatty acids involved in fatty acids catabolism pathway
#############################################################################


# 1) Selection of genes involved in fatty acids(FA) catabolism pathway

gene.names.sel.1 <- c("ACBP","AOX","BIEN","CPT2","CYP4A10",
                     "HPNCL","L.FABP","PECI","PMDCI", "THIOL",
                     "mHMGCoAS", "CACP","Tpalpha", "Tpbeta", "CYP4A14","ACOTH")

pos.genes <- match(x=c(gene.names.sel.1), table=gene.names)


# 2) Selection of fatty acids involved in fatty acids catabolism pathway
pos.FA <- match(x=FA.names[19:21], table=FA.names) 



# 3) Subset of genes or fatty acids involved in fatty acids catabolism pathway
mydata.gene.sel <- mydata.gene[,pos.genes]
mydata.FA.sel <-mydata.FA[,pos.FA]



#############################################################################
# Assign the x variable to the sample values of the first condition.
# Assign the y variable to the sample values of the second condition.
# Only Gene expresions data set
#############################################################################

x <- as.matrix(mydata.gene.sel[ind.gen$wt,])
y <- as.matrix(mydata.gene.sel[ind.gen$ppar,])



#############################################################################
# Heatmaps
# 
#############################################################################

# my heatmap function

myheatmap <-function(z,my.xlab=NULL, my.ylab=NULL, my.main){
  rc <- cm.colors(nrow(z))
  cc <- cm.colors(ncol(z))
  hv <- heatmap(z, col = cm.colors(256), scale = "column",
                Rowv=NA, Colv=NA,
                RowSideColors = rc,
                ColSideColors = cc, 
                margins = c(5,10),
                xlab = my.xlab, ylab =  my.ylab,
                main = my.main)
}


myheatmap(rbind(x,y), my.main="Gene expression")


#############################################################################
# Hotelling test
# 
#############################################################################

if(!require("Hotelling", quietly =T))  install.packages("Hotelling")
require("Hotelling", quietly =T)

(hotelling.test(x, y))


source("http://bioconductor.org/biocLite.R")
biocLite("GSAR")
library("GSAR")


result <- KStest(object=cbind(t(x),t(y)), group=c(rep(1,20),rep(2,20))) 
result$p.value

result <- WWtest(object=cbind(t(x),t(y)), group=c(rep(1,20),rep(2,20))) 
result$p.value

#############################################################################
# Kernel Maximum Mean Discrepancy: asymptotic distribution
# 
#############################################################################

if(!require("kernlab", quietly =T))  install.packages("kernlab")
require("kernlab", quietly =T)

source("Functions_kmmd_sent.R")

time1 <- proc.time()
#results.kmmd.1 <- by (Tr_C_Y.filt_1[sel_row,],ffactor, kmmd.my1)
results.nl1 <- PermTest.knl1(x,y,R=150)
proc.time() - time1


# Empirical distribution of kernel MDD of fatty acids catabolism pathway under H0

hist(c(results.nl1[["pstats"]],results.nl1[["stat"]]), col="azure3",freq=FALSE, 
     main= "Empirical distribution under H0", xlab="Values")

x0=results.nl1[["stat"]]

arrows(x0=x0, y0=1.5, x1=x0, y1= 0, col="coral3", lty=1, lwd=2)
text(x=results.nl1[["stat"]], y=1.5, pos=3, labels="stat")




#############################################################################
# Assign the x variable to the sample values of the first condition.
# Assign the y variable to the sample values of the second condition.
# Gene expresions data set and Fatty acid data set
#############################################################################


# Gene expresions data set (1)
x1 <- as.matrix(mydata.gene.sel[ind.gen$wt,])
y1 <- as.matrix(mydata.gene.sel[ind.gen$ppar,])

# Fatty acids data set (2)
x2 <- as.matrix(mydata.FA.sel[ind.gen$wt,])
y2 <- as.matrix(mydata.FA.sel[ind.gen$ppar,])

# data set joined  by condition
x <- cbind(x1,x2)
y <- cbind(y1,y2)

#############################################################################
# New results
# 
#############################################################################


# heatmap
myheatmap(rbind(x,y), my.main="Gene expression and Fatty acids")


# Hotelling test

(hotelling.test(x, y))



result <- KStest(object=cbind(t(x),t(y)), group=c(rep(1,20),rep(2,20))) 
result$p.value

result <- WWtest(object=cbind(t(x),t(y)), group=c(rep(1,20),rep(2,20))) 
result$p.value




# Kernel Maximum Mean Discrepancy: asymptotic distribution

time1 <- proc.time()
results.nl1 <- PermTest.knl1(x,y,R=2499)
proc.time() - time1
results.nl1$p.v  #p-value


# Empirical distribution of kernel MDD of fatty acids catabolism pathway under H0

hist(c(results.nl1[["pstats"]],results.nl1[["stat"]]), col="azure3",freq=FALSE, 
     main= "Empirical distribution under H0", xlab="Values")

x0=results.nl1[["stat"]]

arrows(x0=x0, y0=1.5, x1=x0, y1= 0, col="coral3", lty=1, lwd=2)
text(x=results.nl1[["stat"]], y=1.5, pos=3, labels="stat")


#############################################################################
# Fatty acids catabolism pathway: sun vs fish diet
#
#############################################################################


#############################################################################
# Assign the x variable to the sample values of the first condition.
# Assign the y variable to the sample values of the second condition.
# Only Gene expresions data set
#############################################################################

x <- as.matrix(mydata.gene.sel[ind.diet1$sun,]) 
y <- as.matrix(mydata.gene.sel[ind.diet1$fish,])



#############################################################################
# New results
# 
#############################################################################


# heatmap
myheatmap(rbind(x,y), my.main="Gene expression")


# Hotelling test

#(hotelling.test(x, y))  # error number of samples < number of variables



result <- KStest(object=cbind(t(x),t(y)), group=c(rep(1,8),rep(2,8))) 
result$p.value

result <- WWtest(object=cbind(t(x),t(y)), group=c(rep(1,8),rep(2,8))) 
result$p.value

# Kernel Maximum Mean Discrepancy: asymptotic distribution

time1 <- proc.time()
results.nl1 <- PermTest.knl1(x,y,R=2499)
proc.time() - time1

results.nl1$p.v  #p-value


#############################################################################
# Assign the x variable to the sample values of the first condition.
# Assign the y variable to the sample values of the second condition.
# Gene expresions data set and Fatty acid data set
#############################################################################

# Gene expresions data set (1)
x1 <- as.matrix(mydata.gene.sel[ind.diet1$sun,]) 
y1 <- as.matrix(mydata.gene.sel[ind.diet1$fish,])

# Fatty acids data set (2)
x2 <- as.matrix(mydata.FA.sel[ind.diet1$sun,]) 
y2 <- as.matrix(mydata.FA.sel[ind.diet1$fish,])

# data set joined  by condition
x <- cbind(x1,x2)
y <- cbind(y1,y2)



# heatmap
myheatmap(rbind(x,y), my.main="Gene expression and Fatty acids")


# Hotelling test

#(hotelling.test(x, y))  # error number of samples < number of variables


result <- KStest(object=cbind(t(x),t(y)), group=c(rep(1,8),rep(2,8))) 
result$p.value

result <- WWtest(object=cbind(t(x),t(y)), group=c(rep(1,8),rep(2,8))) 
result$p.value


# Kernel Maximum Mean Discrepancy: asymptotic distribution

time1 <- proc.time()
results.nl1 <- PermTest.knl1(x,y,R=2499)
proc.time() - time1

results.nl1$p.v  #p-value


