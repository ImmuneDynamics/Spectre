#' Make an enhanced volcano plot
#'
#' @export

make.volcano.plot <- function(dat.p,
                              dat.fc,
                              vars,
                              title, # ctrl
                              subtitle = NULL,

                              xlim = c(-6, 8),
                              ylim = c(0,4),

                              path = getwd(),
                              width = 7,
                              height = 8,
                              
                              pCutoff = 0.05,
                              FCcutoff = 0.26,
                              
                              hline = c(0.01, 0.05),
                              hlineCol = c('grey75', 'grey25'),
                              
                              pointSize = 3.0,
                              labSize = 3.0,
                              
                              col=c('black', 'black', 'blue', 'red3'),
                              colAlpha = 1,
                              
                              ...){

  ### Packages
      # if (!requireNamespace('BiocManager', quietly = TRUE))
      #   install.packages('BiocManager')
      #
      # BiocManager::install('EnhancedVolcano')

      require(EnhancedVolcano)

  ### Demo

      # dat.p <- test.res[3,3:length(names(test.res))]
      # dat.fc <- test.res[5,3:length(names(test.res))]
      #
      # vars <- names(dat.p)
      #
      # name.1 = "COVID acute (ICU)"
      # name.2 = "COVID conv (hosp)"

  ### Setup

  x <- rbind(dat.p, dat.fc)
  x <- t(x)
  x <- as.data.frame(x)
  rownames(x) <- vars
  colnames(x) <- c("p-value", "log2FoldChange")

  x[,1] <- as.numeric(x[,1])
  x[,2] <- as.numeric(x[,2])

  ### Plot
  v <- EnhancedVolcano(x,
                       lab = rownames(x),
                       x = 'log2FoldChange',
                       y = 'p-value',
                       pCutoff = pCutoff,
                       FCcutoff = FCcutoff,
                       xlim = xlim,
                       #xlim = c(-2,2),
                       ylim = ylim,
                       title = title,
                       subtitle = subtitle,

                       hline = hline,
                       hlineCol = hlineCol,

                       pointSize = pointSize,
                       labSize = labSize,

                       col=col,
                       colAlpha = colAlpha
  )
  ### Save
  ggsave(filename = paste0("Volcano plot - ", title, ".png"),
         plot = v,
         path = path,
         width = width,
         height = height,
         limitsize = FALSE)

  ### Print
  print(v)
}
