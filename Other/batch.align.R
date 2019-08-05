# align.batches

align.batches <- function(x,            # data set as a single data.frame/table
                          batch.col,    # columan that defines batches
                          method){      # method of batch alignment
  ## Test parameters
      x <- cell.dat
      batch.col <- GroupName
      method <-

}



dat <- Subset(dat, lg$n2gate)
datr <- gaussNorm(dat, "CD8")$flowset
if(require(flowViz)){
  d1 <- densityplot(~CD8, dat, main="original", filter=curv1Filter("CD8"))
  d2 <- densityplot(~CD8, datr, main="normalized", filter=curv1Filter("CD8"))
  plot(d1, split=c(1,1,2,1))
  plot(d2, split=c(2,1,2,1), newpage=FALSE)
}



library(flowStats)
library(flowCore)
library(flowViz)
library(Biobase)



data(ITN)
dat <- transform(ITN, "CD4"=asinh(CD4), "CD3"=asinh(CD3), "CD8"=asinh(CD8))
lg <- lymphGate(dat, channels=c("CD3", "SSC"),
                preselection="CD4",scale=1.5)
dat <- Subset(dat, lg$n2gate)
datr <- gaussNorm(dat, "CD8")$flowset

if(require(flowViz)){
  d1 <- densityplot(~CD8, dat, main="original", filter=curv1Filter("CD8"))
  d2 <- densityplot(~CD8, datr, main="normalized", filter=curv1Filter("CD8"))
  plot(d1, split=c(1,1,2,1))
  plot(d2, split=c(2,1,2,1), newpage=FALSE)
}






data(GvHD)

## subsetting to flowSet
set <- GvHD[1:4]
GvHD[1:4,1:2]
sel <- sampleNames(GvHD)[1:2]
GvHD[sel, "FSC-H"]
GvHD[sampleNames(GvHD) == sel[1], colnames(GvHD[1]) == "SSC-H"]





## Check data and data column names
head(cell.dat)
dimnames(cell.dat)[[2]]

## Create FCS file metadata - column names with descriptions
metadata <- data.frame(name=dimnames(cell.dat)[[2]], desc=paste('column',dimnames(cell.dat)[[2]],'from dataset'))

## Create FCS file metadata - ranges, min, and max settings -- by default, they are commented out (adjust ranges manually in FlowJo)
#metadata$range <- apply(apply(data,2,range),2,diff) # throws an error because of word entry -- hopefully is ok
#metadata$minRange <- apply(data,2,min)
#metadata$maxRange <- apply(data,2,max)

## Create flowframe with data
data.ff <- new("flowFrame",
               exprs=as.matrix(cell.dat), # in order to create a flow frame, data needs to be read as matrix
               parameters=AnnotatedDataFrame(metadata))

head(exprs(data.ff))

data.fs <- as(list("A"=data.ff),"flowSet")
data.fs


x <- gaussNorm(data.fs, "SA.Bio.Ly6G")$flowset
x







## subsetting to flowSet
set <- data.fs[[1]]
GvHD[1:4,1:2]
sel <- sampleNames(GvHD)[1:2]
GvHD[sel, "FSC-H"]
GvHD[sampleNames(GvHD) == sel[1], colnames(GvHD[1]) == "SSC-H"]




data(ITN)
dat <- transform(ITN, "CD4"=asinh(CD4), "CD3"=asinh(CD3), "CD8"=asinh(CD8))
lg <- lymphGate(dat, channels=c("CD3", "SSC"),
                preselection="CD4",scale=1.5)

dat <- Subset(dat, lg$n2gate)
datr <- gaussNorm(dat, "CD8")$flowset

if(require(flowViz)){
  d1 <- densityplot(~CD8, dat, main="original", filter=curv1Filter("CD8"))
  d2 <- densityplot(~CD8, datr, main="normalized", filter=curv1Filter("CD8"))
  plot(d1, split=c(1,1,2,1))
  plot(d2, split=c(2,1,2,1), newpage=FALSE)
}





###########

test <- as(list("A"=data.ff),"flowSet")

new('flowSet',
    frames = test, # environment with flowFrames
    phenoData = AnnotatedDataFrame(metadata), # object of class AnnotatedDataFrame
    colnames = names(data.ff) # object of class character
)

list.of.dfs <- list()
list.of.dfs[[1]] <- data.ff


flowSet(list.of.dfs, phenoData = AnnotatedDataFrame(metadata))

GvHD



flowSet(GvHD[[1]], GvHD[[2]])
pd <- phenoData(GvHD)[1:2,]
flowSet(s5a01=GvHD[[1]], s5a02=GvHD[[2]],phenoData=pd)



x <- flowSet(data.ff)


datr <- gaussNorm(dat, "FSC.A")$flowset

