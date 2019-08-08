







### Add sample names/numbers

#CytoTools::annotate.files()

## Review the column names
as.matrix(names(data.list))

## Create empty lists
sample.names = list()
sample.assign = list()

## Create sample names and numbers
sample.names <- as.list(c("01_Mock_01",
                          "02_Mock_02",
                          "03_Mock_03",
                          "04_Mock_04",
                          "05_Mock_05",
                          "06_Mock_06",
                          "07_WNV_01",
                          "08_WNV_02",
                          "09_WNV_03",
                          "10_WNV_04",
                          "11_WNV_05",
                          "12_WNV_06"))

sample.names

## Sample assignments
as.matrix(unique(names(data.list)))

sample.assign[[1]] <- 1
sample.assign[[2]] <- 2
sample.assign[[3]] <- 3
sample.assign[[4]] <- 4
sample.assign[[5]] <- 5
sample.assign[[6]] <- 6

sample.assign[[7]] <- 7
sample.assign[[8]] <- 8
sample.assign[[9]] <- 9
sample.assign[[10]] <- 10
sample.assign[[11]] <- 11
sample.assign[[12]] <- 12

sample.names
sample.assign

head(data.list)


##
Spectre::annotate.files(x = data.list, # data
                        col.name = "Sample", # what do you want to call the new column
                        d = sample.names, # what data are you adding (list)
                        assignments = sample.assign, # which files does each keyword get added to
                        add.num = TRUE) # do you also want a 'number' column

head(data.list)



### Add group names

## Review the column names
as.matrix(names(data.list))

## Create empty lists
group.names = list()
group.assign = list()

## Create group names
group.names[[1]] <- "Mock"
group.names[[2]] <- "WNV"

## Specify which files (by number) that belong to each group
group.assign[[1]] <- c(1:6)
group.assign[[2]] <- c(7:12)

group.names
group.assign

## Perform group name embedding and review
Spectre::annotate.files(x = data.list, # data
                        col.name = "Group", # what do you want to call the new column
                        d = group.names, # what data are you adding (list)
                        assignments = group.assign, # which files does each keyword get added to
                        add.num = TRUE) # do you also want a 'number' column

head(data.list)

head(data.list)
head(data.list[[1]])
head(data.list[[7]])
head(data.list[[12]])

### Add batch details

## Review the column names
as.matrix(names(data.list))

## Create empty lists
batch.nums = list()
batch.assign = list()

## Create batch numbers
batch.nums <- as.list(c("01",
                        "02"))

batch.nums

## Sample assignments
as.matrix(unique(names(data.list)))

batch.assign[[1]] <- c(1,3,5,7,9,11)
batch.assign[[2]] <- c(2,4,6,8,10,12)

batch.nums
batch.assign


##
Spectre::annotate.files(x = data.list, # data
                        col.name = "Batch", # what do you want to call the new column
                        d = batch.nums, # what data are you adding (list)
                        assignments = batch.assign, # which files does each keyword get added to
                        add.num = FALSE) # do you also want a 'number' column

head(data.list[[1]])
head(data.list[[12]])

