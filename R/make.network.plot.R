make.network.plot <- function(dat, 
                              timepoint.col,
                              timpoints,
                              cluster.col,
                              marker.cols,
                              node.size = 13,
                              standard.colours = 'spectral',
                              save.to.disk = TRUE) {
  require(Spectre)
  require(gtools)
  require(tidygraph)
  require(ggraph)
  require(RColorBrewer)
  require(viridis)
  
  # Function to find the position of round brackets in cluster id
  get.idx.roundclsbracket <- function(cl) {
    tail(gregexpr("\\)", cl)[[1]], n=1)
  }
  
  # To store transitions
  edge.df <- data.frame(from=character(),
                        to=character())
  
  #### Find all the transitions ####
  for (tp.idx in c(2:length(timepoints))) {
    prev.tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx-1], ]
    prev.tp.clust <- mixedsort(unique(prev.tp.dat[[cluster.col]]))
    
    complex.clusters <- prev.tp.clust[lapply(prev.tp.clust, get.idx.roundclsbracket) > -1]
    simple.clusters <- setdiff(prev.tp.clust, complex.clusters)
    
    # Order the vector based on the length of each element. 
    # Need to use rev as whatever is in it order it in order it in increasing manner. STUPID!
    simple.clusters <- rev(simple.clusters[order(sapply(simple.clusters,length))])
    
    curr.tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx], ]
    curr.tp.clust <- mixedsort(unique(curr.tp.dat[[cluster.col]]))
    
    for (cl in curr.tp.clust) {
      clean.cl <- cl
      ## Match the complex clusters first
      for (cls in complex.clusters) {
        cls.found <- grepl(cls, clean.cl, fixed = TRUE)
        if (cls.found) {
          df <- data.frame(paste0(tp.idx-1,'_', cls), paste0(tp.idx, '_', cl))
          names(df) <- c("from","to")
          edge.df <- rbind(edge.df, df)
          clean.cl <- gsub(cls, "", clean.cl)
        }
      }
      ## Then simple clusters
      for (cls in simple.clusters) {
        cls.found <- grepl(cls, clean.cl, fixed = TRUE)
        if (cls.found) {
          df <- data.frame(paste0(tp.idx-1,'_',cls), paste0(tp.idx,'_',cl))
          names(df) <- c("from","to")
          edge.df <- rbind(edge.df, df)
          clean.cl <- gsub(cls, "", clean.cl)
        }
      }
    }
  }
  
  #### Then get extra details on the nodes ####
  all.clust <- lapply(c(1:length(timepoints)), function(tp.idx) {
    tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx], ]
    tp.clust <- mixedsort(unique(tp.dat[[cluster.col]]))
    tp.clust.name <- sapply(tp.clust, function(c) {
      return(paste0(tp.idx, '_', c))
    })
    return(tp.clust.name)
  })
  all.clust <- unlist(all.clust)
  
  cluster.ids <- sapply(all.clust, function(c) {
    strsplit(c, '_')[[1]][2]
  })
  timepoints.col <- sapply(all.clust, function(c) {
    strsplit(c, '_')[[1]][1]
  })
  node.dat <- data.frame(nodeId = all.clust,
                         clusterId = cluster.ids,
                         timepoints = timepoints.col)
  node.dat$id <- seq.int(nrow(node.dat))
  
  #### Convert edges so it's to and from the numeric id of the node ####
  names(edge.df) <- c('source', 'destination')
  edges <- edge.df %>%
    left_join(node.dat, by = c("source" = "nodeId")) %>%
    rename(from = id)
  edges <- edges %>%
    left_join(node.dat, by = c("destination" = "nodeId")) %>%
    rename(to = id)
  # filter out other node information
  keep.cols <- c('from', 'to')
  edges <- data.table(edges)
  edges <- edges[,..keep.cols]
  
  node.ids <- as.vector(node.dat$nodeId)
  
  #### Compute each cluster's mean expression ####
  for (idx in marker.cols) {
    marker <- names(dat)[idx]
    marker.mean <- sapply(as.vector(node.dat$nodeId), function(id) {
      id_split <- strsplit(id, '_')[[1]]
      timepoint <- timepoints[as.numeric(id_split[1])]
      cluster_id <- id_split[2]
      sub.dat <- dat[dat[[timepoint.col]] == timepoint & dat[[cluster.col]] == cluster_id,]
      mean(sub.dat[[marker]])
    })
    node.dat[[marker]] <- marker.mean
  }
  
  #### Add an extra option which will colour just the node using the very first cluster ####
  cluster.tp1 <- node.dat[node.dat$timepoints == 1,]$clusterId
  cluster.colours <- sapply(c(1:nrow(node.dat)), function(i) {
    row <- node.dat[i,]
    cluster.remain <- row$clusterId %in% cluster.tp1

    if (cluster.remain) {
      return(as.character(row$clusterId))
    } else {
      return("NA")
    }
  })
  node.dat$origin <- cluster.colours
  
  #### Start plotting ####
  img.height <- 15
  img.width <- 15
  
  routes_tidy <- tbl_graph(nodes = as.data.frame(node.dat),
                           edges = as.data.frame(edges),
                           directed = TRUE)
  
  #### Draw plot coloured by timepoints ####
  ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
    geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
    geom_node_point(aes(colour = timepoints), size=node.size) +
    scale_colour_manual(values = rev(
      colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(timepoints))),
      breaks = as.character(timepoints) +
    theme(text = element_text(size=20),
          panel.background = element_rect(fill = 'white'),
          aspect.ratio = 1) + 
    labs(color='DPO bin')
  
  ggsave(paste0("network_colBy_timepoints.pdf"),
         width = img.width,
         height = img.height,
         dpi = 1000,
         limitsize = FALSE)
  
  #### Colour plot based on the origin cluster ####
  ## +1 for the colour brewer because we have NA
  ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
    geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
    geom_node_point(aes(colour = origin), size=node.size) +
    scale_colour_manual(values = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(length(cluster.tp1) + 1)),
                        breaks = as.character(unique(node.dat$origin))) +
    theme(text = element_text(size=20),
          panel.background = element_rect(fill = 'white'),
          aspect.ratio = 1) +
    labs(color='Cluster origin')

  ggsave(paste0("network_colBy_origin.pdf"),
         width = img.width,
         height = img.height,
         dpi = 1000,
         limitsize = FALSE)
  
  #### Colour plot by markers ####
  
  if(standard.colours == "spectral"){
    spectral.list <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(50))
    colour.palette <- scale_colour_gradientn(colours = spectral.list)
  } else if(standard.colours == "inferno"){
    colour.palette <- scale_colour_viridis_c(option='inferno')
  }
  
  
  
  
  for (idx in marker.cols) {
    marker <- names(dat)[idx]
    
    ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
      geom_edge_link(arrow = arrow(length = unit(1, 'mm')), end_cap = circle(2, 'mm')) +
      geom_node_point(aes_string(colour = marker), size=node.size) +
      colour.palette +
      theme(text = element_text(size=20),
            panel.background = element_rect(fill = 'white'))
    
    ggsave(paste0("network_colBy_", marker, '.pdf'),
           width = img.width,
           height = img.height,
           dpi = 1000,
           limitsize = FALSE)
  }
  
}

