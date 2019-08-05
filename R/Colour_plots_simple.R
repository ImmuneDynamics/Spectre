# Spectre::colour.plot

colour.plot <- function(d,
                        x.axis, # "UMAP1"
                        y.axis, # "UMAP2"
                        col.axis, # "BV605.Ly6C"
                        title,
                        colours,
                        dot.size
                        ){

  ## Colour setup
      if(colours == "jet"){
        jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        colour.scheme <- jet
      }

      if(colours == "spectral"){
        spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
        spectral.list <- rev(spectral.list)
        spectral <- colorRampPalette(c(spectral.list))
        colour.scheme <- spectral
      }

  ## Generate and show coloured plot
  print(ggplot(data = d, aes(x = d[[x.axis]], y = d[[y.axis]], colour = d[[col.axis]])) +
    geom_point(size = dot.size) +
    scale_colour_gradientn(colours = colour.scheme(50),
                           limits = c(quantile(d[[col.axis]], probs = c(0.01)), #0.03-01 seems to work well
                                      quantile(d[[col.axis]], probs = c(0.995))), #0.97-995 seems to work well
                           oob=squish) +
    ggtitle(title) +
    labs(colour = col.axis)+
    xlab(x.axis)+
    ylab(y.axis)+
      theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5)#, # change 'colour' to black for informative axis
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
  )

  ## ggsave?

}



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
          labs(colour = col.axis)+
          xlab(x.axis)+
          ylab(y.axis)+
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
  )
  }

