#' Draw network diagram
#'
#' Method to draw a network diagram.
#' This rely on a clustering performed by ChronoClust.
#' 
#' @param dat NO DEFAULT. Data.frame. Data to be clustered.
#' @param timepoint.col NO DEFAULT. Column name which represents the time point of each cell (data point) in dat.
#' @param timepoints NO DEFAULT. The time points (in order).
#' @param cluster.col NO DEFAULT. Column denoting the cluster id.
#' @param marker.cols NO DEFAULT. Vector of column names denoting the markers to plot.
#' @param node.size DEFAULT = 13. Size of the node of the diagram.
#' @param arrow.length DEFAULT = 3. Length of the arrow connecting nodes.
#' @param arrow.head.gap DEFAULT = 4. The gap between head of the arrow and node.
#' @param standard.colours DEFAULT = "Spectral". Colour scheme for the markers. Spectral or Inferno or Viridis.
#' 
#'@usage
#'make.network.plot(dat, timepoint.col, timepoints, cluster.col, marker.cols,
#'node.size = 13, arrow.length = 3, arrow.head.gap = 4, standard.colours = 'Spectral')
#'
#'
#' @author Givanna Putri, \email{givanna.haryonoputri@@sydney.edu.au}
#' @export

make.network.plot <- function(dat, 
                              timepoint.col,
                              timepoints,
                              cluster.col,
                              marker.cols,
                              node.size = 'auto',
                              arrow.length = 3,
                              arrow.head.gap = 4,
                              standard.colours = 'Spectral') {
  require(Spectre)
  require(gtools)
  require(tidygraph)
  require(ggraph)
  require(RColorBrewer)
  require(viridis)
  require(stringr)
  
  ## for testing
  # timepoint.col = 'Group'
  # timepoints = c('Mock', 'WNV-01', 'WNV-02', 'WNV-03', 'WNV-04', 'WNV-05')
  # cluster.col = "ChronoClust_cluster_lineage"
  # #cluster.col = "cluster_id"
  # marker.cols = cluster.cols.nos
  # node.size = 4
  # arrow.length = 1
  # arrow.head.gap = 2
  # standard.colours = 'Spectral'
  
  ## Have to make sure column header have no hyphen
  if (TRUE %in% stringr::str_detect(colnames(dat), '-')) {
    message("Some column headers have hyphen (-) in it. Please rename them first!.")
    message("No plots are created.")
    return()
  }
  
  message("Calculating edges")
  # To store transitions
  edge.df <- data.frame(from=character(),
                        to=character())
  
  #### Find all the transitions ####
  for (tp.idx in c(2:length(timepoints))) {
    prev.tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx-1], ]
    prev.tp.clust <- mixedsort(unique(prev.tp.dat[[cluster.col]]))
    
    complex.clusters <- prev.tp.clust[lapply(prev.tp.clust, get.idx.roundclsbracket) > -1]
    # simple.clusters  <- setdiff(prev.tp.clust, complex.clusters)
    simple.n.pipe.clusters <- setdiff(prev.tp.clust, complex.clusters)
    pipe.clusters <- simple.n.pipe.clusters[sapply(simple.n.pipe.clusters, get.idx.pipe) > -1]
    simple.clusters <- setdiff(simple.n.pipe.clusters, pipe.clusters)
    
    # Order the vector based on the length of each element. 
    simple.clusters <- simple.clusters[order(nchar(simple.clusters), simple.clusters, decreasing = TRUE)]
    
    curr.tp.dat <- dat[dat[[timepoint.col]] == timepoints[tp.idx], ]
    curr.tp.clust <- mixedsort(unique(curr.tp.dat[[cluster.col]]))
    
    # this filter out clusters that only exist in current time point.
    # these are new clusters and have no predecessors.
    curr.tp.clust.existing <- lapply(curr.tp.clust, function(cl) {
      if (grepl(",", cl, fixed = TRUE) || grepl("|", cl, fixed = TRUE) || cl %in% simple.clusters) {
        return(cl)
      } 
    })
    curr.tp.clust.existing <- unlist(curr.tp.clust.existing)
    
    for (cl in curr.tp.clust.existing) {
      clean.cl <- cl
      ## Match the complex clusters first
      for (cls in complex.clusters) {
        cls.found <- grepl(cls, clean.cl, fixed = TRUE)
        if (cls.found) {
          df <- data.frame(paste0(tp.idx-1,'_', cls), paste0(tp.idx, '_', cl))
          names(df) <- c("from","to")
          edge.df <- rbind(edge.df, df)
          clean.cl <- gsub(cls, "", clean.cl, fixed = TRUE)
        }
      }
      
      ## Match the split simple clusters first
      for (cls in pipe.clusters) {
        cls.found <- grepl(cls, clean.cl, fixed = TRUE)
        if (cls.found) {
          df <- data.frame(paste0(tp.idx-1,'_', cls), paste0(tp.idx, '_', cl))
          names(df) <- c("from","to")
          edge.df <- rbind(edge.df, df)
          clean.cl <- gsub(cls, "", clean.cl, fixed = TRUE)
        }
      }
      
      ## Then simple clusters
      for (cls in simple.clusters) {
        cls.found <- grepl(cls, clean.cl, fixed = TRUE)
        if (cls.found) {
          df <- data.frame(paste0(tp.idx-1,'_',cls), paste0(tp.idx,'_',cl))
          names(df) <- c("from","to")
          edge.df <- rbind(edge.df, df)
          clean.cl <- gsub(cls, "", clean.cl, fixed = TRUE)
        }
      }
    }
  }
  
  message("Computing node details")
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
  timepoints.idx <- sapply(all.clust, function(c) {
    strsplit(c, '_')[[1]][1]
  })
  timepoints.actual <- sapply(timepoints.idx, function(c) {
    timepoints[as.numeric(c)]
  })
  node.dat <- data.frame(nodeId = all.clust,
                         clusterId = cluster.ids,
                         timepointIdx = timepoints.idx,
                         timepoint = timepoints.actual)
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
  
  message("Calculating marker's average per node")
  
  #### Compute each cluster's mean expression ####
  for (idx in marker.cols) {
    marker <- names(dat)[idx]
    marker.mean <- sapply(as.vector(node.dat$nodeId), function(id) {
      id_split <- strsplit(id, '_')[[1]]
      timepoint <- timepoints[as.numeric(id_split[1])]
      cl_id <- id_split[2]
      sub.dat <- dat[dat[[timepoint.col]] == timepoint & dat[[cluster.col]] == cl_id,]
      mean(sub.dat[[marker]])
    })
    node.dat[[marker]] <- marker.mean
  }
  
  #### Add an extra option which will colour just the node using the very first cluster ####
  # We only want clusters which is alphabetical, the very first cluster
  cluster.origins <- get.cluster.origins(node.dat)
  node.dat$origin <- cluster.origins
  
  #colnames(node.dat) <- gsub(" ","_",colnames(node.dat))
  
  if (node.size == 'auto') {
    node.sizes <- get.cluster.proportions(node.dat = node.dat,
                                          dat = dat,
                                          timepoint.col = timepoint.col,
                                          timepoints = timepoints,
                                          cluster.col = cluster.col)
    
    node.dat$ProportionOfCells <- node.sizes
  }
  
  
  #### Start plotting ####
  img.height <- 15
  img.width <- 15
  
  routes_tidy <- tbl_graph(nodes = as.data.frame(node.dat),
                           edges = as.data.frame(edges),
                           directed = TRUE)
  
  geompoint.prop <- get.geompoint(node.size, 'timepoint')
  
  #### Draw plot coloured by timepoints ####
  
  colour.palette <- get.colour.pallete(standard.colours,
                                       n.col = length(timepoints),
                                       factor = as.character(mixedsort(unique(node.dat$timepoint))))
  
  
  
  ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
    geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')), end_cap = circle(arrow.head.gap, 'mm')) +
    geompoint.prop +
    colour.palette +
    theme(text = element_text(size=20),
          panel.background = element_rect(fill = 'white'),
          aspect.ratio = 1) + 
    labs(color='Timepoint')
  
  ggsave(paste0("network_colBy_timepoints.pdf"),
         width = img.width,
         height = img.height,
         dpi = 1000,
         limitsize = FALSE)
  
  #### Colour plot based on the origin cluster ####
  ## +1 for the colour brewer because we have NA
  colour.palette <- get.colour.pallete(standard.colours,
                                       n.col = length(unique(cluster.origins)) + 1,
                                       factor = as.character(mixedsort(unique(node.dat$origin))))
  geompoint.prop <- get.geompoint(node.size, 'origin')
  
  ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
    geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')), end_cap = circle(arrow.head.gap, 'mm')) +
    geompoint.prop +
    colour.palette +
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
  colour.palette <- get.colour.pallete(standard.colours)
  
  for (idx in marker.cols) {
    marker <- names(dat)[idx]
    
    geompoint.prop <- get.geompoint(node.size, marker)
    
    ggraph(routes_tidy, layout = 'kk', maxiter = 10000) +
      geom_edge_link(arrow = arrow(length = unit(arrow.length, 'mm')), 
                     end_cap = circle(arrow.head.gap, 'mm')) +
      geompoint.prop +
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

# Function to find the position of round brackets in cluster id
get.idx.roundclsbracket <- function(cl) {
  tail(gregexpr("\\)", cl)[[1]], n=1)
}

get.idx.pipe <- function(cl) {
  tail(gregexpr("\\|", cl)[[1]], n=1)
}

get.colour.pallete <- function(standard.colours, n.col=50, factor=NULL) {
  if(tolower(standard.colours) == "spectral"){
    spectral.list <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(n.col))
    if (is.null(factor)) {
      colour.palette <- scale_colour_gradientn(colours = spectral.list)
    } else {
      colour.palette <- scale_colour_manual(values = spectral.list,
                                            breaks = factor)
    }
  } else if(tolower(standard.colours) == "inferno"){
    if (is.null(factor)) {
      colour.palette <- scale_color_viridis(option='B')
    } else {
      colour.palette <- scale_color_viridis(discrete = TRUE, option = "B")
    }
  } else if(tolower(standard.colours) == "viridis"){
    if (is.null(factor)) {
      colour.palette <- scale_color_viridis(option='D')
    } else {
      colour.palette <- scale_color_viridis(discrete = TRUE, option = "D")
    }
    
  }
  return(colour.palette)
}

get.cluster.origins <- function(node.dat) {
  cluster.origins <- sapply(as.vector(unique(node.dat$clusterId)), function(cl.id) {
    alphabet.only <- grepl('^[A-Z]+$', cl.id)
    if (alphabet.only) {
      return(cl.id)
    }
    return("NotOrigin")
  })
  # remove the cluster which is not alphabet only (as above function assign it to null)
  cluster.origins <- cluster.origins[cluster.origins != 'NotOrigin']
  
  # cluster.origins <- node.dat[node.dat$timepoints == 1,]$clusterId
  cluster.colours <- sapply(c(1:nrow(node.dat)), function(i) {
    row <- node.dat[i,]
    cluster.remain <- row$clusterId %in% cluster.origins
    
    if (cluster.remain) {
      return(as.character(row$clusterId))
    } else {
      return("NA")
    }
  })
}

get.cluster.proportions <- function(node.dat, dat, timepoint.col, timepoints, cluster.col) {
  
  # count cell per time point
  cell.count.per.tp <- dat[, .N, by=dat[[timepoint.col]]]
  
  proportions <- sapply(c(1:nrow(node.dat)), function(i) {
    node.row <- node.dat[i,]
    sub.dat <- dat[dat[[timepoint.col]] == node.row$timepoint &
                     dat[[cluster.col]] == node.row$clusterId, ]
    cell.cnt <- cell.count.per.tp[cell.count.per.tp$dat == node.row$timepoint, N]
    return(nrow(sub.dat)/cell.cnt)
  })
  
  
  return(proportions)
}

get.geompoint <- function(node.size, col.by) {
  if (node.size == 'auto') {
    geompoint <- geom_node_point(aes_string(colour = col.by, size='ProportionOfCells')) 
  }
  else {
    geompoint <- geom_node_point(aes_string(colour = col.by), size=node.size)
  }
  return(geompoint)
}

