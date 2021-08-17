### Load libraries

library('Spectre')

Spectre::package.check(type = 'spatial')
Spectre::package.load(type = 'spatial')

### Set PrimaryDirectory

dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### File list

files <- list.files(getwd(), '.csv')
as.matrix(files)

### Table

dt <- data.table("FileName" = files)

### Write

setwd(PrimaryDirectory)
dir.create("FileTable")
setwd("FileTable")

fwrite(dt, 'FileTable.csv')