# Thomas Ashhurst
# 2017-12-06

# Objective: 
    # tSNE dataset/plot -- training data set
    # Use a classifier to assign NEW cells to either
    # a) a specific set of tSNE x/y coordinates
    # b) to a specific (pre-determined) cluster -- could be phenograph cluster, or manually gated 'cluster' with the clutser no added to full dataset

# HERE:
    # https://www.datacamp.com/community/tutorials/machine-learning-in-r#seven
    # Using the iris dataset to learn
    # Classifier is supposed to assign the data points to one of the flower species, based on similarity to existing data points
    # Uses KNN


###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 

    ### 1.1 - Install packages (if not already installed)
      if(!require('ggvis')) {install.packages('ggvis')} # to visualise
      if(!require('dplyr')) {install.packages('dplyr')} # helps ggvis pipelines
      if(!require('class')) {install.packages('class')} # package to run KNN classifier
      if(!require('gmodels')) {install.packages('gmodels')} # evalute the results of the classifier
      if(!require('ggplot2')) {install.packages('ggplot2')} # evalute the results of the classifier
    
    ### 1.2 Load packages 
        library('ggvis')
        library('dplyr')
        library('class')
        library('gmodels')
        library('ggplot2')

    ### 1.3 - Set working directory and assign as 'PrimaryDirectory'
        #setwd("/Users/thomasashhurst/Google Drive/2017 Ashhurst work/1. ASHHURST/01. R&D - Technology/[Script projects]/Classifier/V2.7 modulating k with larger dataset/")
        # Set your working directory here (e.g. "/Users/Tom/Desktop/") -- press tab when selected after the '/' to see options
        
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        #getwd()
        #PrimaryDirectory <- getwd()
        #PrimaryDirectory
        
        getwd()                                                         # Check your working directory has changed correctly
        PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'
        PrimaryDirectory

###################################################### 2. PREPARE DATA ###################################################### 

    
    ### 2.1 Load data ###
    
        # Set up data
        #data_start <- iris
        #data_start
        
        list.files(path=PrimaryDirectory, pattern = ".csv")
        data_start <- read.csv("Output_Sample_Data_with_tSNE.csv")    
        data_new_large <- read.csv("Output_Sample_Data_before_downsample.csv")
        
        head(data_start)
        head(data_new_large)
        
        names(data_start)
        names(data_new_large)
        
        as.matrix(names(data_start))
        as.matrix(names(data_new_large))
        
        ## add columns of 0 for tSNE1, tSNE2, and phenograph cluster to the data_new_large dataset
        
        tSNE1 <- c(0)
        tSNE2 <- c(0)
        phenograph_cluster <- c(0)
        
        data_new_large <- cbind(data_new_large, tSNE1, tSNE2, phenograph_cluster)
        head(data_new_large)
           
        as.matrix(names(data_start))
        as.matrix(names(data_new_large))
        
        head(data_start)
        head(data_new_large)
        
        rm(tSNE1)
        rm(tSNE2)
        rm(phenograph_cluster)
        

        
    ### 2.2 Review data ###    
        
        head(data_start)
        names(data_start)
        str(data_start)
        
        # Division of `Species`
        table(data_start$tSNE1) 
        
        # Percentual division of `Species`
        round(prop.table(table(data_start$tSNE1)) * 100, digits = 1)
        
        # Plot data
        #data_start %>% ggvis(~Sepal.Length, ~Sepal.Width, fill = ~tSNE1) %>% layer_points()
        #data_start %>% ggvis(~Petal.Length, ~Petal.Width, fill = ~tSNE1) %>% layer_points()
    
              #data_start %>% ggvis(~tSNE1, ~tSNE2, fill = ~phenograph_cluster) %>% layer_points()
    
    ### 2.3. Assign 'data'
        
        data <- data_start # pass from data_start to data
        
        
    ### 2.3 OPTIONAL - Normalise necessary columns - SKIP if NOT NEEDED ###  

        # Normalise values of cellular parameters -- the ones used for the classifier
    
        
        ## BUILD SOMETHING HERE -- assign which columns to use
        
        
        
        # Build your own `normalize()` function
        normalise <- function(x) {
          num <- x - min(x)
          denom <- max(x) - min(x)
          return (num/denom)
          }
        
        # Choose the columns to 'normalize' from the the data, then normalise
        as.matrix(names(data_start))
        
              #cols_to_normalise <- c(2:62)
        
              #data_norm <- as.data.frame(lapply(data_start[cols_to_normalise], normalise)) # only normalises columns that should be normalised, once specified
        data_norm <- as.data.frame(data)
        
        data_new_large <- as.data.frame(data_new_large)
        
        # Summarize `iris_norm`
        names(data_norm)
        summary(data_norm)
        
        # data_normalised <- as.data.frame(lapply(data, normalize))
        
        # Replace equivalent columns in data_start with columns from data_norm
        # Add 'normalised' columns to the starting dataset -- replace the original columns
        # That way we can keep tSNE values, cluster numbers etc; run the classifier, and add the results
        
        head(data)
        head(data_norm)
        
        as.matrix(names(data))
        as.matrix(names(data_norm))
        
        #data[cols_to_normalise] <- data_norm[c(1:61)] # data will be original columns, data_norm will be like 1:4 -- DOUBLE CHECK!!!! If you mess this up, then it messes up everything else
        #head(data)
   
        #as.matrix(names(data_start))
        #as.matrix(names(data))
        
        #check <- data.frame(names(data_start), names(data))   # attempts to use 'as.numeric' here for data_predict messes up the result
        
        # Specify column names for `merge`
        #names(check) <- c("Start", "After_norm")
        
        # Inspect `merge` 
        #check
        
        #head(data)
        
        # data now normalised, with all columns intact
        
        
###################################################### 3. SETUP TRAINING AND TEST DATASETS ######################################################         

        
    ### 3.1. Create function for random assignment of data

        # set seed for random number generator -- divide up the dataset EVENLY
        set.seed(1234)
        
        # sample data
        ind <- sample(2, nrow(data), replace=TRUE, prob=c(0.67, 0.33)) # Creates a new column, assigns a 1 to 67% of the datapoints, and 2 to 33%
    
        ## Here we split DATA into DATA.TRAINING and DATA.TEST ##
        
    ### 3.2 Creating the TRAINING SET
            
            # Review column numbers
            head(data) # check for any NAs
            as.matrix(names(data)) 
            
            # Compose pre-training set
            data.pre.training <- data[c(1:66)] #[ind==1, c(1:66)]  # include ALL columns, except any that have 'N/A' -- SampleID sometimes does this
            dim(data.pre.training)
            head(data.pre.training)
            
            # Review column numbers
            as.matrix(names(data.pre.training))
            
            # Choose columns for classifier -- # here select only columns used for the classifier
            cols_for_classifier <- c(8,10,11,12,13,14,15,16,18,19,21,23,24,25,33,34,35,36,37,38,39,48,52,53,55,59,61,62)
            
            # Take ONLY columns used for classifier
            data.training <- data.pre.training[cols_for_classifier] 
            
            # Inspect training set
            head(data.training)
            dim(data.training)
            data.training  #row numbers will be 'random'
        
      ### 3.3. Create the TEST SET
            
            # Review column numbers
            as.matrix(names(data_new_large))
            
            # Compose test set
            data.pre.test <- data_new_large[c(1:66)] #[ind==2, c(1:66)] # include ALL columns, except any that have 'N/A' -- SampleID sometimes does this
            head(data.pre.test)
            
            # Review column numbers
            as.matrix(names(data.pre.test))
            
            # Columns for classifier already chosen above
            
            # Take ONLY columns used for classifier
            data.test <- data.pre.test[cols_for_classifier] # here select only columns used for the classifier
            
            # Inspect training set
            head(data.test)
            dim(data.test)
            data.test #row numbers will be 'random'

            
            dim(data.training)
            dim(data.test)

            #https://discuss.analyticsvidhya.com/t/how-to-count-the-missing-value-in-r/2949/3
            
            #If You need NA count of all --- 
            table(is.na(data.training)) 
            
            #If you need NA count Column wise -- 
            sapply(data.training, function(x) sum(is.na(x)))
            
            #If you need NA count Row wise --- 
            rowSums(is.na(data.training))
            
            
            
            
            #If You need NA count of all --- 
            table(is.na(data.test)) 
            
            #If you need NA count Column wise -- 
            sapply(data.test, function(x) sum(is.na(x)))
            
            #If you need NA count Row wise --- 
            rowSums(is.na(data.test))
            
            
###################################################### 4. SPECIFY TARGET VARIABLE AND RUN CLASSIFIER ######################################################  

    ### SPECIFY K
            
            kvalue <- 3 # default is 3     # try 1, 2, 4, 8, 16, 32
            
            as.matrix(names(data))
            as.matrix(names(data_new_large))
            
            
            col_for_tSNE1             <- 64
            col_for_tSNE2             <- 65
            col_for_PhenographCluster <- 66
            
            
            
    ### STEP 4.1. First run (tSNE1) 
           
        ## Specify target variable
            as.matrix(names(data))
            
            # Compose `data` training labels
            data.trainLabels <- data[col_for_tSNE1] # OK because set.seed was applied to ind
            
            # could also use data.trainLabels data.trainLabels <- data[7], but this includes the original 'row number' -- not sure if that would fit yet
            
            # Inspect result
            print(data.trainLabels)
            names(data.trainLabels)
            
            # Compose `data` test labels
            data.testLabels <- data_new_large[col_for_tSNE1] # OK because set.seed was applied to ind
            
            # Inspect result
            print(data.testLabels)
            names(data.trainLabels)
            
        ## Run KNN
        
            # Build the model
            data_predict <- class:::knn(train = data.training, # had to add class::: because knn is clearly used by multiple packages
                                        test = data.test, 
                                        cl = data.trainLabels[,1], 
                                        k=kvalue)     # Note: the k parameter is often an odd number to avoid ties in the voting scores.
            
            # Inspect `iris_pred`
            data_predict
            
            # See the 'real' and 'predicted' labels for dataset
            as.matrix(data.testLabels)
            as.matrix(data_predict)
        
        ## Manage result
            
            # Put `data.testLabels` in a data frame
            dataTestLabels <- data.frame(data.testLabels)
            
            # Merge `data_predict` and `data.testLabels` 
            merge <- data.frame(data_predict, data.testLabels) # attempts to use 'as.numeric' here for data_predict messes up the result
            
            # Specify column names for `merge`
            names(merge) <- c("Predicted", "Observed")
            
            # Inspect `merge` 
            merge
            
            # Save
            write.csv(x = merge, file = "merge_tSNE1.csv")
            
            # Add 'Predicted' to data.pre.test
            head(data.pre.test)
            
            data.pre.test$tSNE1_new <- merge$Predicted
            data.pre.test$tSNE1_original <- merge$Observed
            head(data.pre.test)
            # data.pre.test now contains the 'new' tSNE1 value
            
    ### STEP 4.2. Second run (tSNE2)
        
            ## Specify target variable
                as.matrix(names(data))
                
                # Compose `data` training labels
                data.trainLabels <- data[col_for_tSNE2] # OK because set.seed was applied to ind
                
                # could also use data.trainLabels data.trainLabels <- data[7], but this includes the original 'row number' -- not sure if that would fit yet
                
                # Inspect result
                print(data.trainLabels)
                
                # Compose `data` test labels
                data.testLabels <- data_new_large[col_for_tSNE2] # OK because set.seed was applied to ind
                
                # Inspect result
                print(data.testLabels)
                
            ## Run KNN
                
                # Build the model
                data_predict <- class:::knn(train = data.training, # had to add class::: because knn is clearly used by multiple packages
                                            test = data.test, 
                                            cl = data.trainLabels[,1], 
                                            k=kvalue)     # Note: the k parameter is often an odd number to avoid ties in the voting scores.
                
                # Inspect `iris_pred`
                data_predict
                
                # See the 'real' and 'predicted' labels for dataset
                as.matrix(data.testLabels)
                as.matrix(data_predict)
                
            ## Manage result
                
                # Put `data.testLabels` in a data frame
                dataTestLabels <- data.frame(data.testLabels)
                
                # Merge `data_predict` and `data.testLabels` 
                merge <- data.frame(data_predict, data.testLabels)   # attempts to use 'as.numeric' here for data_predict messes up the result
                
                # Specify column names for `merge`
                names(merge) <- c("Predicted", "Observed")
                
                # Inspect `merge` 
                merge
                
                # Save
                write.csv(x = merge, file = "merge_tSNE2.csv")
                
                # Add 'Predicted' to data.pre.test
                head(data.pre.test)
                
                data.pre.test$tSNE2_new <- merge$Predicted
                data.pre.test$tSNE2_original <- merge$Observed
                head(data.pre.test)
                # data.pre.test now contains the 'new' tSNE2 value

        
    ### STEP 4.3. Third run (Phenograph?)        
                
                ## Specify target variable
                    as.matrix(names(data))
                    
                    # Compose `data` training labels
                    data.trainLabels <- data[col_for_PhenographCluster] # OK because set.seed was applied to ind
                    
                    # could also use data.trainLabels data.trainLabels <- data[7], but this includes the original 'row number' -- not sure if that would fit yet
                    
                    # Inspect result
                    print(data.trainLabels)
                    
                    # Compose `data` test labels
                    data.testLabels <- data_new_large[col_for_PhenographCluster] # OK because set.seed was applied to ind
                    
                    # Inspect result
                    print(data.testLabels)
                
                ## Run KNN
                    
                    # Build the model
                    data_predict <- class:::knn(train = data.training, # had to add class::: because knn is clearly used by multiple packages
                                                test = data.test, 
                                                cl = data.trainLabels[,1], 
                                                k=kvalue)     # Note: the k parameter is often an odd number to avoid ties in the voting scores.
                    
                    # Inspect `iris_pred`
                    data_predict
                    
                    # See the 'real' and 'predicted' labels for dataset
                    as.matrix(data.testLabels)
                    as.matrix(data_predict)
                    
                ## Manage result
                    
                    # Put `data.testLabels` in a data frame
                    dataTestLabels <- data.frame(data.testLabels)
                    
                    # Merge `data_predict` and `data.testLabels` 
                    merge <- data.frame(data_predict, data.testLabels) # attempts to use 'as.numeric' here for data_predict messes up the result
                    
                    # Specify column names for `merge`
                    names(merge) <- c("Predicted", "Observed")
                    
                    # Inspect `merge` 
                    merge
                    
                    # Save
                    write.csv(x = merge, file = "merge_phenograph.csv")
                    
                    # Add 'Predicted' to data.pre.test
                    head(data.pre.test)
                    
                    data.pre.test$phenograph_new <- merge$Predicted
                    data.pre.test$phenograph_original <- merge$Observed
                    head(data.pre.test)
                    # data.pre.test now contains the 'new' phenograph value

    ### STEP 4.4. Output data
                
                data_output <- as.data.frame(data.pre.test)
                write.csv(x = data_output, file = "data_output.csv")
                
                names(data_output)

                
###################################################### 5. EVALUATE CLASSIFIER ######################################################          
    
    
    ### 5.1. Evalute model
    
        # Plot 'merge'
        data_output %>% ggvis(~tSNE1_original, ~tSNE1_new) %>% layer_points() 

              # plots here are weird because the outcome value is a 'class', not numeric -- once exported to .csv, not an issue

        # gmodels -- works with classes, not continuous numeric really
        CrossTable(x = data.testLabels, 
                   y = data_predict, 
                   prop.chisq=FALSE
                    )


      ### 5.2. Re-import data (from .csv) to visualise
        
        data_viz <- read.csv("data_output.csv")
        head(data_viz)
        
        data_viz %>% ggvis(~tSNE1_original, ~tSNE1_new) %>% layer_points()
        
        data_viz %>% ggvis(~tSNE1_original, ~tSNE2_original, fill = ~phenograph_cluster) %>% layer_points()
        data_viz %>% ggvis(~tSNE1_new, ~tSNE2_new, fill = ~phenograph_cluster) %>% layer_points()
        
        data_viz %>% ggvis(~tSNE1_original, ~tSNE2_original, fill = ~tSNE1_original) %>% layer_points()
        data_viz %>% ggvis(~tSNE1_new, ~tSNE2_new, fill = ~tSNE1_original) %>% layer_points()
        
        
        
        data_viz %>% ggvis(~phenograph_original, ~phenograph_new, fill = ~phenograph_cluster) %>% layer_points()
        
        
        
        ## V2.4 - removing the noisy or rubbish channels definately improved accuracy
        
        
        ## V2.5 -- large dataset seems to have worked well -- lso remove any markers t hat weren't used for tSNE
        

        
                
        

###################################################### OTHER ######################################################         
        
        # ggplot
        ggplot(data = data_output, aes(x = tSNE1, y = tSNE2)) #????
        

        
        
        # Original tSNE1
        data_viz %>% ggvis(~tSNE1, ~tSNE2, fill = ~Event.No) %>% layer_points()
        
        # New tSNE1
        data_viz %>% ggvis(~tSNE1_new, ~tSNE2, fill = ~Event.No) %>% layer_points()
        
        # Old vs New tSNE1
        data_viz %>% ggvis(~tSNE1, ~tSNE1_new, fill = ~Event.No) %>% layer_points()
        
        
        
        
        #If You need NA count of all --- 
        table(is.na(data_viz)) 
        
        #If you need NA count Column wise -- 
        sapply(data.training, function(x) sum(is.na(x)))
        
        #If you need NA count Row wise --- 
        rowSums(is.na(data.training))
        
        
        
    
    
    
    # ggplot
    ggplot(data = data_output,
           aes(tSNE1, tSNE2)
      )
    
    warnings()
    
    ? ggplot
    
    
    if(!require('extrafont')) {install.packages('extrafont')} # evalute the results of the classifier
    library(extrafont)
    
    # need only do this once!
    font_import(pattern="[A/a]rial", prompt=FALSE)
    
    install.packages("ggplot2")
    library(ggplot2)
    
    
    #https://stackoverflow.com/questions/10581440/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon
    #https://stackoverflow.com/questions/34220883/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon

    
    #https://stackoverflow.com/questions/34220883/error-in-grid-calll-textbounds-as-graphicsannotxlabel-xx-xy-polygon
    # --> might have something to do with mis-matched 'row numbers)
    
