## factor.plot

factor.plot <- function(d,
                        x.axis, # "UMAP1"
                        y.axis, # "UMAP2"
                        col.axis, # "BV605.Ly6C"
                        title,
                        dot.size){

    print(ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = as.factor(d[[col.axis]]))) +
            geom_point(size = dot.size)+ # 2 for large # 0.5 for small
            #scale_colour_gradientn(colours = colour.scheme(50)) +
            #scale_colour_manual(name = "FileName", values = c(colour.scheme(length(filenames)))) +
            ggtitle(title) +
            #labs(colour = col.axis)+
            xlab(x.axis)+
            ylab(y.axis)+
            theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                  #axis.line=element_blank(),
                  #axis.text.x=element_blank(),
                  #axis.text.y=element_blank(),
                  #axis.ticks=element_blank(),
                  axis.title.x=element_text(color="Black", face="bold", size=18),
                  axis.title.y=element_text(color="Black", face="bold", size=18),
                  #panel.grid.major = element_blank(),
                  #panel.background=element_blank(),
                  #panel.border=element_blank(),
                  #panel.grid.minor=element_blank(),
                  #plot.background=element_blank(),
                  #legend.position = "right",
                  legend.text=element_text(size=12), # large = 30 # small = 8
                  legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                  legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                  legend.title=element_blank(),
                  plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
            )
          )
    }


factor.plot.labelled <- function(d,
                                  x.axis, # "UMAP1"
                                  y.axis, # "UMAP2"
                                  col.axis, # "BV605.Ly6C"
                                  title,
                                  dot.size){

  # PREP Data frame with centroids, one entry per centroid
  centroidsDf <- data.frame(
    centroidX = tapply(d[[x.axis]], d[[col.axis]], median), # median
    centroidY = tapply(d[[y.axis]], d[[col.axis]], median),
    centroidCol = tapply(d[[col.axis]], d[[col.axis]], median))

  print(ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = as.factor(d[[col.axis]]))) +
          geom_point(size = dot.size)+ # 2 for large # 0.5 for small
          #scale_colour_gradientn(colours = colour.scheme(50)) +
          #scale_colour_manual(name = "FileName", values = c(colour.scheme(length(filenames)))) +
          ggtitle(title) +
          #labs(colour = col.axis)+
          xlab(x.axis)+
          ylab(y.axis)+
          theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                #axis.line=element_blank(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                axis.title.x=element_text(color="Black", face="bold", size=18),
                axis.title.y=element_text(color="Black", face="bold", size=18),
                #panel.grid.major = element_blank(),
                #panel.background=element_blank(),
                #panel.border=element_blank(),
                #panel.grid.minor=element_blank(),
                #plot.background=element_blank(),
                #legend.position = "right",
                legend.text=element_text(size=12), # large = 30 # small = 8
                legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                legend.title=element_blank(),
                plot.title = element_text(color="Black", face="bold", size=22, hjust=0) # size 70 for large, # 18 for small
          ) +

          geom_point(data = centroidsDf,                     # separate data.frame
                     aes(x = centroidX, y = centroidY),
                     col = "black",                          # notice "col" and "shape" are
                     #shape = 1,
                     size = 2) +

          geom_label(data = centroidsDf,
                     hjust = 0,
                     nudge_x = 0.7,
                     aes(x = centroidX, y = centroidY, label=centroidCol, alpha = 0.5),
                     col = "black",
                     fontface = "bold") +

          guides(alpha = "none")
  )
}


# In GG -- layer of centroids

