# plotting

    ### Load packages
    library('Spectre')

    library('plyr')
    library('data.table')
    library('rstudioapi')
    library('flowViz')
    library('flowCore')
    library('Biobase')
    library('FlowSOM')
    library('Rtsne')
    library('umap')
    library('ggplot2')
    library('scales')
    library('colorRamps')
    library('ggthemes')
    library('RColorBrewer')
    library("gridExtra")

    library("gridExtra")

    ## In order for this to work, a) rstudioapi must be installed and b) the location of this .r script must be in your desired working directory
    dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
    getwd()
    PrimaryDirectory <- getwd()
    PrimaryDirectory

    ## Use this to manually set the working directory
    #setwd("/Users/Tom/Desktop/Experiment")                          # Set your working directory here (e.g. "/Users/Tom/Desktop/") -- press tab when selected after the '/' to see options
    #getwd()                                                         # Check your working directory has changed correctly
    #PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'
    #PrimaryDirectory

    list.files(path=PrimaryDirectory, pattern = ".csv")     # see a list of CSV files


    plot.dat <- as.data.frame(fread(file = "TA270.2.csv"))


    ## Define the sample and group column names
    samp.col <- "SampleName"
    grp.col <- "GroupName"

    ## Check to ensure the correct name has been specified above
    plot.dat[[samp.col]]
    plot.dat[[grp.col]]

    ## Make a sample table
    make.sample.table(x = plot.dat,
                      sample.col.name = samp.col,
                      include.groups = TRUE,
                      group.col.name = grp.col)

    # Check results
    all.sample.names
    all.group.names

    sample.table



    ###





    # https://ggplot2.tidyverse.org/reference/lims.html
    #ggplot(small, aes(mpg, wt, colour = factor(cyl))) +
    #  geom_point() +
    #  lims(colour = c("4", "6", "8"))
    lims(colour = c("4", "6", "8"))





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




