##########################################################################################################
#### DRAFT Cytometry Analysis Pipeline for large and compleX data (CAPX) v3.0 - using the Spectre R package
#### Part 2 - Plotting, data exploration, cluster annotation
##########################################################################################################

    # Thomas Myles Ashhurst, Felix Marsh-Wakefield
    # 2019-08-02
    # Workflow: https://sydneycytometry.org.au/capx
    # Spectre R package: https://sydneycytometry.org.au/spectre

##########################################################################################################
#### 1. Install packages, load packages, and set working directory
#########################################################################################################

    ### 1.1. Load 'Spectre' package (using devtools)
        if(!require('devtools')) {install.packages('devtools')}
        library('devtools')

        #if(!require('sydneycytometry/spectre')) {install_github("sydneycytometry/spectre")}
        library("Spectre")

        # We recommend not to update packages that are dependencies of Spectre

    ### 1.2. Install packages
        if(!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if(!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if(!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if(!require("scales")){install.packages("scales")} # for re-scaling if necessary
        if(!require("RColorBrewer")){install.packages("RColorBrewer")} # for re-scaling if necessary
        if(!require("gridExtra")){install.packages("gridExtra")} # for re-scaling if necessary

        ## Plotting
        if (!require("ggplot2")){install.packages("ggplot2")} # for plotting tSNE graphs
        if (!require("colorRamps")){install.packages("colorRamps")} # for colour scheme management
        if (!require("ggthemes")){install.packages("ggthemes")} # for plot themes
        if (!require("scales")){install.packages("scales")} # for re-scaling if necessary
        if (!require("RColorBrewer")){install.packages("RColorBrewer")} # for re-scaling if necessary
        if (!require("gridExtra")){install.packages("gridExtra")} # for re-scaling if necessary

    ### 1.3. Load packages from library
        library('ggplot2')
        library('scales')
        library('colorRamps')
        library('ggthemes')
        library('RColorBrewer')
        library("gridExtra")

        library('ggplot2')
        library('scales')
        library('colorRamps')
        library('ggthemes')
        library('RColorBrewer')
        library("gridExtra")

    ## 1.4. Set working directory

        ## Set working directory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

        # ## Can set manually using these lines, if desired
        # PrimaryDirectory <- "/Users/Tom/Google Drive (t.ashhurst@centenary.org.au)/_Sydney Cytometry/2019_Synced/GitHub/Public github/Spectre/Workflow scripts/Output_CAPX/Output-DimRedData/"
        # setwd(PrimaryDirectory)

        ## Create output directory
        dir.create("Output_CAPX_exploration", showWarnings = FALSE)
        setwd("Output_CAPX_exploration")
        OutputDirectory <- getwd()
        setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Read in data for dot plots
#########################################################################################################

    ### Read in cell.dat

        ## cell.dat
        setwd("./Output_CAPX/Output-ClusteredData/")
        cell.dat_Directory <- getwd()
        list.files(cell.dat_Directory, ".csv")

        cell.dat <- fread("Demo_all_data.csv")

        setwd(PrimaryDirectory)

        ## cell.dat.sub
        setwd("./Output_CAPX/Output-DimRedData/")
        cell.dat.sub_Directory <- getwd()
        list.files(path = cell.dat.sub_Directory, ".csv")

        cell.dat.sub <- fread("DimRed_Demo_all_data.csv")

        setwd(PrimaryDirectory)


##########################################################################################################
#### Basic plotting and data examination
##########################################################################################################

        ### Plot single UMAP MFI plot
        names(cell.dat.sub)

        p1 <- Spectre::colour.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "BV605.Ly6C",
                                   title = paste0("All samples", " - ", "BV605.Ly6C"),
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat.sub
        )
        p1

        ### Plot single UMAP factor plot -- clusters
        names(cell.dat.sub)

        p2 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "FlowSOM_metacluster",
                                   title = paste0("All samples", " - ", "Clusters"),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        p2 # TRUE is fine, FALSE also returns 'NULL' when running this line IF using option for centroid.labels as FALSE

        p2 <- Spectre::labelled.factor.plot(d = cell.dat.sub,
                                            x.axis = "UMAP_42_X",
                                            y.axis = "UMAP_42_Y",
                                            col.axis = "FlowSOM_metacluster",
                                            title = paste0("All samples", " - ", "Clusters"),
                                            dot.size = 1,
                                            align.xy.by = cell.dat.sub,
                                            align.col.by = cell.dat)
        p2


        ### Plot single UMAP factor plot -- samples
        names(cell.dat.sub)

        p3 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "Sample",
                                   title = "Samples",
                                   dot.size = 0.5) # assumes numeric

        p3


        ### Plot single UMAP factor plot -- groups
        names(cell.dat.sub)

        p4 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "Group",
                                   title = "Groups",
                                   dot.size = 0.5) # assumes numeric

        p4

        #ggsave(filename = "p1.jpeg", plot = p1, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p2.jpeg", plot = p2, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p3.jpeg", plot = p3, path = OutputDirectory, width = 9, height = 7)
        #ggsave(filename = "p4.jpeg", plot = p4, path = OutputDirectory, width = 9, height = 7)

        ### Create tiles and save

        gp <- grid.arrange(grobs = list(p1, p2, p3, p4), ncol=2, nrow=2) #top = "Main Title"
        gp

        ggsave(filename = "Grid.jpeg", plot = gp, path = OutputDirectory, width = 18, height = 14)

        ### --> FUNCTION WITH XMAX/MIN etc
        ### --> seperate function for COL MAX/MIN etc





##########################################################################################################
#### 2. Read in data
#########################################################################################################

    ## Read cell.dat
    list.files(PrimaryDirectory, ".csv")
    cell.dat.sub <- fread("DimRed_Demo_all_data.csv")
    head(plot.dat)

    ## Also read the full dataset (not subsampled) for reference
    setwd(PrimaryDirectory)
    setwd("../Output-ClusteredData/")
    list.files(getwd(), ".csv")

    cell.dat <- fread("Demo_all_data.csv")
    head(all.dat)

    setwd(PrimaryDirectory)

    ## Define the sample and group column names
    as.matrix(names(cell.dat))

    samp.col <- "Sample"
    grp.col <- "Group"
    clust.col <- "FlowSOM_metacluster"

    cell.dat[[samp.col]]
    cell.dat[[grp.col]]
    cell.dat[[clust.col]]

    ## COMPARE cell.dat and cell.dat.sub -- anything missing??

            all.samples <- unique(cell.dat[[samp.col]])
            all.groups <- unique(cell.dat[[grp.col]])
            all.clusters <- unique(cell.dat[[clust.col]])

            ## Others if necessary
            batch.col <- "Batch"
            cell.dat[[batch.col]]

            all.batches <- unique(cell.dat[[batch.col]])


##########################################################################################################
#### Standard plots
##########################################################################################################

    ### Plotting aligned colour plots -- CD45, CD117, Ly6C

        as.matrix(names(cell.dat.sub))

        p1 <- Spectre::colour.plot(d = cell.dat.sub, # plot.dat
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "AF700.CD45",
                                   title = paste0("All samples", " - ", "CD45"),
                                   colours = "spectral",
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        p2 <- Spectre::colour.plot(d = cell.dat.sub, # plot.dat
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "BV605.Ly6C",
                                   title = paste0("All samples", " - ", "Ly6C"),
                                   colours = "spectral",
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        p3 <- Spectre::colour.plot(d = cell.dat.sub, # plot.dat
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = "BB515.CD117",
                                   title = paste0("All samples", " - ", "CD117"),
                                   colours = "spectral",
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        f1 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = clust.col,
                                   title = paste0("All samples", " - ", clust.col),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        f2 <- Spectre::factor.plot(d = cell.dat.sub,
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = samp.col,
                                   title = paste0("All samples", " - ", samp.col),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)

        f3 <- Spectre::factor.plot(d = cell.dat.sub, # plot.dat
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = grp.col,
                                   title = paste0("All samples", " - ", grp.col),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)


        gp <- grid.arrange(grobs = list(p1, p2, p3, f1, f2, f3), ncol=3, nrow=2) #top = "Main Title"

        ggsave(filename = "Grid.jpeg", plot = gp, path = OutputDirectory, width = 27, height = 14)

        ## GRID FUNCTION

        as.matrix(names(cell.dat.sub))
        cell.dat.sub[[samp.col]]

        Spectre::factor.plot(d = subset(cell.dat.sub, cell.dat.sub[[samp.col]] == "01_Mock_01"), # plot.dat
                                   x.axis = "UMAP_42_X",
                                   y.axis = "UMAP_42_Y",
                                   col.axis = grp.col,
                                   title = paste0("All samples", " - ", grp.col),
                                   dot.size = 1,
                                   align.xy.by = cell.dat.sub,
                                   align.col.by = cell.dat)





##########################################################################################################
#### Heatmaps
#########################################################################################################

        ### Read in sumtable data

        ## cell freq.
        setwd("./Output_CAPX/Output-SumTables/Output_CellNums/")
        sum.tables.freq_Directory <- getwd()
        list.files(sum.tables.freq_Directory, ".csv")
        setwd(PrimaryDirectory)

        ## MFI per marker
        setwd("./Output_CAPX/Output-SumTables/Output_MFI_per_marker/")
        sum.tables.permarker_Directory <- getwd()
        list.files(sum.tables.permarker_Directory, ".csv")
        setwd(PrimaryDirectory)

        ## MFI per sample
        setwd("./Output_CAPX/Output-SumTables/Output_MFI_per_sample/")
        sum.tables.persample_Directory <- getwd()
        list.files(sum.tables.persample_Directory, ".csv")
        setwd(PrimaryDirectory)










        ## Clusters vs samples -- number of  cells
        setwd(PrimaryDirectory)
        setwd("../Output-SumTables/Output_CellNums/")
        list.files(pattern = ".csv")
        A <- read.csv("SumTable_CellsPerTissue.csv")
        head(A)

        r.height <- 5.5

        dev.off()
        pdf("heatmap.pdf", width = (r.height * 3.268765), height = r.height)
        make.heatmap(x = A,
                     sample.col = "Sample",
                     annot.cols = c(1:3),
                     plot.title = "Cells per tissue",
                     row.sep = 6,
                     y.margin = 4,
                     x.margin = 8)
        dev.off()



        ## Clusters vs markers -- MFI

        ##
        setwd(PrimaryDirectory)
        setwd("../Output-SumTables/Output_MFI_per_sample/")

        list.files(pattern = ".csv")
        B <- read.csv("SumTable_MFI_cluster_x_marker-all_data.csv")
        head(B)
        as.matrix(names(B))

        ###
        r.height <- 8

        dev.off()
        pdf("heatmap2.pdf", width = (r.height*0.8888), height = r.height)
        make.heatmap(x = B,
                     sample.col = "FlowSOM_metacluster",
                     annot.cols = c(1:2,4,8,11,15,16,21,32),
                     plot.title = "All data - MFI",
                     y.margin = 8,
                     x.margin = 4)
        dev.off()



        ###
        Spectre::make.pheatmap(x = B,
                                sample.col = "FlowSOM_metacluster",
                                annot.cols = c(1:2,4,8,11,15,16,21,32),
                                plot.title = "All data - MFI",
                                cell.size = 15
                                )



##########################################################################################################
#### 4. Annotate clusters
#########################################################################################################









##########################################################################################################
#### 7. 'Positive expression' thresholds
#########################################################################################################






##########################################################################################################
#### Make new sumtables
#########################################################################################################










    #library(gtable)

    #g1 <- ggplotGrob(p1)
    #g2 <- ggplotGrob(p2)
    #g <- rbind(g1, g2, size = "first")
    #g$widths <- unit.pmax(g1$widths, g2$widths)

    #g

    #grid.newpage()
    #grid.draw(g)



    gp <- grid.arrange(grobs = list(p1, p2, p3, p4), ncol=2, nrow=2) #top = "Main Title"
    gp

    ggsave(filename = "Grid.jpeg", plot = gp, path = OutputDirectory, width = 18, height = 14)


    ### --> FUNCTION WITH XMAX/MIN etc

    ### --> seperate function for COL MAX/MIN etc



## Group tiles

# another function
# find the X and Y max, and the COL max


str(plot.dat)
names (plot.dat)


plot(plot.dat$FileNo, plot.dat$`Yb176Di :: 176Yb_Ly6C`)


marker
Xmax <-
Xmin <-

Ymax <-
Ymin <-


  scale_colour_gradientn(colours = jet.colors(50),
                         limits = c(quantile(combined.df[[i]], probs = c(0.01)), #0.03-01 seems to work well
                                    quantile(combined.df[[i]], probs = c(0.995))), #0.97-995 seems to work well
                         oob=squish) +

d1 <- subset(plot.dat, plot.dat[, grp.col] == "Air")
gp1 <- Spectre::colour.plot(d = d1,
                            x.axis = "UMAP_42_X",
                            y.axis = "UMAP_42_Y",
                            col.axis = "FlowSOM_metacluster",
                            title = "Air -- Clusters",
                            colours = "spectral",
                            dot.size = 0.5)

d2 <- subset(plot.dat, plot.dat[, grp.col] == "Smoke")
gp2 <- Spectre::colour.plot(d = d2,
                            x.axis = "UMAP_42_X",
                            y.axis = "UMAP_42_Y",
                            col.axis = "FlowSOM_metacluster",
                            title = "Smoke -- Clusters",
                            colours = "spectral",
                            dot.size = 0.5)

d3 <- subset(plot.dat, plot.dat[, grp.col] == "Smoke+FT")
gp3 <- Spectre::colour.plot(d = d3,
                            x.axis = "UMAP_42_X",
                            y.axis = "UMAP_42_Y",
                            col.axis = "FlowSOM_metacluster",
                            title = "Smoke + FT -- Clusters",
                            colours = "spectral",
                            dot.size = 0.5)

gp.grid <- grid.arrange(grobs = list(gp1, gp2, gp3), ncol=3, nrow=1) #top = "Main Title"















##########################################################################################################
#### Automatically generate coloured tSNE/UMAP plots and save to disk
##########################################################################################################

## Loop for colour.plots

plots <- names(cell.dat.sub)
# plots <- c("BV605.Ly6C", "PE.CD115")

### RATHER -- JUST SUBSET TO NUMERIC

setwd(PrimaryDirectory)
setwd(OutputDirectory)
dir.create("Output-ColourPlots")
setwd("Output-ColourPlots")

## Plot all data
for(a in plots){
  p <- Spectre::colour.plot(d = cell.dat.sub, x.axis = "UMAP_42_X", y.axis = "UMAP_42_Y",
                            col.axis = a, title = a, colours = "spectral", dot.size = 1)
  ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
}

## Plot by sample
#for(i in unique(cell.dat.sub[["FileName"]])){
#  # subset
#}
## Plot by group

setwd(PrimaryDirectory)


## Loop for factor.plots

head(cell.dat.sub)
factors <- c("FlowSOM_metacluster")

setwd(OutputDirectory)

for(a in factors){
  p <- Spectre::factor.plot(d = cell.dat.sub,
                            x.axis = "UMAP_42_X",
                            y.axis = "UMAP_42_Y",
                            col.axis = "FlowSOM_metacluster",
                            title = "Cluster",
                            dot.size = 1,
                            add.labels = TRUE) # assumes numeric

  ggsave(p, filename = paste0("All_samples_", a, ".png"), width = 9, height = 7)
}

setwd(PrimaryDirectory)








##########################################################################################





test <- function(d){
  ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
    geom_point(size = point.size) +
    scale_colour_gradientn(colours = colour.scheme(50),
                           limits = c(quantile(d[[col.axis]], probs = c(min.threshold)), #0.03-01 seems to work well
                                      quantile(d[[col.axis]], probs = c(max.threshold))), #0.97-995 seems to work well
                           oob=squish) +
    labs(colour = col.axis)+
    xlab(x.axis)+
    ylab(y.axis)+
    ggtitle(plot.title) +
    theme(
      panel.background = element_rect(fill = "white", colour = "black", size = 0.5) # change 'colour' to black for informative axis
      #axis.line=element_blank(),
      #axis.text.x=element_blank(),
      #axis.text.y=element_blank(),
      #axis.ticks=element_blank(),
      #axis.title.x=element_blank(),
      #axis.title.y=element_blank(),
      #panel.grid.major = element_blank(),
      #panel.background=element_blank(),
      #panel.border=element_blank(),
      #panel.grid.minor=element_blank(),
      #plot.background=element_blank(),
      #legend.position = "right",
      #legend.text=element_text(size=15), # large = 30 # small = 8
      #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
      #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
      #legend.title=element_blank()
      #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
    )
}


test(cell.dat.sub)




#######

jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#spectral.list <- list(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))
spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
spectral.list <- rev(spectral.list)

spectral <- colorRampPalette(c(spectral.list))

YlOrRd.list <- rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(100))
YltoRd.cols <- colorRampPalette(c(YlOrRd.list))

BuPu.list <- rev(colorRampPalette(brewer.pal(9,"BuPu"))(100))
BuPu.cols <- colorRampPalette(c(BuPu.list))


# Custom, from spectral, Red to Blue -->
custom <- colorRampPalette(rev(c("#7F0000", "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B",
                                 "#FFFFBF",
                                 "#E6F598",
                                 "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#00007F",
                                 "#170a45")))

colour.palette.list <- (colorRampPalette(brewer.pal(9, "YlGnBu"))(31)) # 256
colour.palette <- colorRampPalette(c(colour.palette.list))

fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black","#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2")))
fold.palette <- colorRampPalette(rev(c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","black")))
#fold.palette <- colorRampPalette(c(fold.palette.list))


colour.scheme <- jet
#colour.scheme <- spectral
#colour.scheme <- YltoRd.cols
#colour.scheme <- BuPu.cols
#colour.scheme <- custom
#colour.scheme <- colour.palette
#colour.scheme <- fold.palette


head(cell.dat.sub)
filenames <- unique(cell.dat.sub$FileName)


## Plot, coloured by samples
ggplot(data = cell.dat.sub, aes(x = UMAP_42_X, y = UMAP_42_Y, colour = as.factor(cell.dat.sub$GroupName))) +
  geom_point(size = 1)+ # 2 for large # 0.5 for small
  #scale_colour_gradientn(colours = colour.scheme(50)) +
  #scale_colour_manual(name = "FileName", values = c(colour.scheme(length(filenames)))) +
  ggtitle("Samples") +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5)#, # change 'colour' to black for informative axis
    #axis.line=element_blank(),
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    #axis.title.x=element_blank(),
    #axis.title.y=element_blank(),
    #panel.grid.major = element_blank(),
    #panel.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.minor=element_blank(),
    #plot.background=element_blank(),
    #legend.position = "right",
    #legend.text=element_text(size=15), # large = 30 # small = 8
    #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
    #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
    #legend.title=element_blank(),
    #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
  )




head(cell.dat.sub)

## Plot, colour by MFI

d <- cell.dat.sub

x.axis <- "UMAP_42_X"
y.axis <- "UMAP_42_Y"
col.axis <- "BV711.SCA.1"

min.threshold <- 0.01 # 0.01
max.threshold <- 0.995 # 0.995 default

plot.title <- "MFI"
point.size <- 0.5 # 2 for large # 0.5 for small


ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
  geom_point(size = point.size) +
  scale_colour_gradientn(colours = colour.scheme(50),
                         limits = c(quantile(d[[col.axis]], probs = c(min.threshold)), #0.03-01 seems to work well
                                    quantile(d[[col.axis]], probs = c(max.threshold))), #0.97-995 seems to work well
                         oob=squish) +
  labs(colour = col.axis)+
  xlab(x.axis)+
  ylab(y.axis)+
  ggtitle(plot.title) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", size = 0.5) # change 'colour' to black for informative axis
    #axis.line=element_blank(),
    #axis.text.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks=element_blank(),
    #axis.title.x=element_blank(),
    #axis.title.y=element_blank(),
    #panel.grid.major = element_blank(),
    #panel.background=element_blank(),
    #panel.border=element_blank(),
    #panel.grid.minor=element_blank(),
    #plot.background=element_blank(),
    #legend.position = "right",
    #legend.text=element_text(size=15), # large = 30 # small = 8
    #legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
    #legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
    #legend.title=element_blank()
    #plot.title = element_text(color="Black", face="bold", size=15, hjust=0) # size 70 for large, # 18 for small
  )




