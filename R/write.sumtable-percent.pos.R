#' write.sumtable.percent.pos
#'
#' @export

write.sumtable.percent.pos <- function(x,
                                       sample.name,
                                       cluster.name,
                                       Markers,
                                       Cutoff
                                       )
{

  ### Test data
    # library(Spectre)
    # library(data.table)
    #
    # x <- as.data.table(demo.clustered)
    #
    # Spectre::embed.columns(x = x,
    #                        type = "data.frame",
    #                        base.name = "FlowSOM_metacluster",
    #                        col.name = "PopName",
    #                        match.to = c(1:40),
    #                        new.cols = c(rep("Pop1", 10),
    #                                     rep("Pop2", 10),
    #                                     rep("Pop3", 10),
    #                                     rep("Pop4", 10))
    #                        )
    #
    # x <- cbind(x, embed.res)
    # x
    #
    # sample.name <- "Sample"
    # cluster.name <- "FlowSOM_metacluster"
    # Markers <- c("BV711.SCA.1", "APC.BrdU", "PECy7.CD80")
    # Cutoff <- c(500, 700, 600)


  ###
      x <- as.data.table(x)

      clusters <- unique(x[[cluster.name]])
      clusters <- sort(clusters)
      clusters

      samples <- unique(x[[sample.name]])
      samples <- sort(samples)
      samples

      MarkerCutoffs <- data.frame(Markers, Cutoff)
      MarkerCutoffs

  ###
  for(i in c(1:nrow(MarkerCutoffs))){ ## EACH CSV IS A MARKER
    # i <- 1
    active.marker <- MarkerCutoffs[i,1] #,1 is name, ,2 is cutoff
    active.marker <- as.character(active.marker)
    active.marker

    active.marker.cutoff <- MarkerCutoffs[i,2]
    active.marker.cutoff

    df <- data.frame(matrix(ncol = (1 + length(clusters)), nrow = 0))
    colnames(df) <- c("Sample", as.vector(clusters))
    df

    for(a in samples){ ## EACH ROW IS A SAMPLE
      #a <- "01_Mock_01"
      active.sample <- a

      x.samp <- x[x[[sample.name]] == a]
      x.samp

      # sample.res <- data.table(Samples = samples)
      # sample.res
      #
      # temp <- data.table(matrix(nrow = 12, ncol = length(clusters)))
      # temp
      # names(temp) <- as.character(clusters)
      # temp
      #
      # sample.res <- cbind(sample.res, temp)
      # sample.res

      all.cluster.res <- list()

      for(o in clusters){  ## EACH COLUMN IS THE RESULTS OF A CLUSTER
        #o <- 1
        active.cluster <- o

        active.sample.cluster <- x.samp[x.samp[[cluster.name]] == o,] # Sample == a, Cluster == o
        total <- nrow(active.sample.cluster) # Number of total events of this cluster from this sample

        pos <- active.sample.cluster[active.sample.cluster[[active.marker]] >active.marker.cutoff ,]
        num.pos <- nrow(pos)

        freq <- num.pos / total
        freq <- freq*100
        freq

        res <- data.frame(Cluster = o, Freq = freq)
        res

        all.cluster.res <- rbind(all.cluster.res, res)
      }

      active.sample
      all.cluster.res

      df[nrow(df) + 1,] = c(active.sample, as.vector(all.cluster.res[,2]))
      df

    }

    df
    write.csv(x = df, file = paste0(active.marker, ".csv"))
  }
}


