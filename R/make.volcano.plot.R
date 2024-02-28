#' Make an enhanced volcano plot
#' 
#' This function creates a volcano plot using 'EnhancedVolcano'.
#' 
#' @param dat.p NO DEFAULT. 
#' @param dat.fc NO DEFAULT. 
#' @param vars NO DEFAULT. 
#' @param title NO DEFAULT. 
#' @param subtitle DEFAULT = NULL. 
#' @param xlim DEFAULT = c(-6, 8). x-axis limits.
#' @param ylim DEFAULT = c(0, 4). y-axis limits.
#' @param path DEFAULT = getwd(). The location to save plots.
#' @param width DEFAULT = 7.
#' @param height DEFAULT = 8.
#' @param pCutoff DEFAULT = 0.05. Statistical significance cut-off.
#' @param FCcutoff DEFAULT = 0.26. Absolute log2 fold-change cut-off.
#' @param hline DEFAULT = c(0.01, 0.05). Horizontal line(s).
#' @param hlineCol DEFAULT = c("grey75", "grey25"). Horizontal line colour(s).
#' @param pointSize DEFAULT = 3.0.
#' @param labSize DEFAULT = 3.0. Label size.
#' @param col DEFAULT = c("black", "black", "blue", "red3").
#' @param colAlpha DEFAULT = 1.
#' 
#' @usage make.volcano.plot(dat.p, dat.fc, vars, title, subtitle, xlim, ylim, path, width, height, pCutoff, FCcutoff, hline, hlineCol, pointSize, labSize, col, colAlpha, ...)
#' 
#' @examples
#' # Prepare data
#' sum.dat <- Spectre::demo.sum
#' colnames(sum.dat)
#' to.plot <- names(sum.dat)[c(4:15)]
#' 
#' # Calculate statistics
#' stat.dat <- Spectre::create.stats(dat = sum.dat, 
#'                                   use.cols = to.plot, 
#'                                   sample.col = 'Sample', 
#'                                   group.col = 'Group', 
#'                                   comparisons = list(c('Mock', 'WNV')
#'                                   )
#' 
#' # Make volcano plot
#' Spectre::make.volcano.plot(dat.p = stat.dat[2,..to.plot], 
#'                            dat.fc = stat.dat[1,..to.plot], 
#'                            vars = to.plot, title = 'Volcano'
#'                            )
#' 
#' @author
#' Thomas M Ashhurst, \email{thomas.ashhurst@@sydney.edu.au}
#' Felix Marsh-Wakefield, \email{felix.marsh-wakefield@@sydney.edu.au}
#' 
#' @export

make.volcano.plot <- function(dat.p,
                              dat.fc,
                              vars,
                              title, # ctrl
                              subtitle = NULL,
                              xlim = c(-6, 8),
                              ylim = c(0, 4),
                              path = getwd(),
                              width = 7,
                              height = 8,
                              pCutoff = 0.05,
                              FCcutoff = 0.26,
                              hline = c(0.01, 0.05),
                              hlineCol = c("grey75", "grey25"),
                              pointSize = 3.0,
                              labSize = 3.0,
                              col = c("black", "black", "blue", "red3"),
                              colAlpha = 1,
                              ...) {

  ### Demo

  # dat.p <- dat[3,3:length(names(dat))]
  # dat.fc <- dat[5,3:length(names(dat))]
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

  x[, 1] <- as.numeric(x[, 1])
  x[, 2] <- as.numeric(x[, 2])

  ### Plot
  v <- EnhancedVolcano::EnhancedVolcano(x,
    lab = rownames(x),
    x = "log2FoldChange",
    y = "p-value",
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    xlim = xlim,
    # xlim = c(-2,2),
    ylim = ylim,
    title = title,
    subtitle = subtitle,
    hline = hline,
    hlineCol = hlineCol,
    pointSize = pointSize,
    labSize = labSize,
    col = col,
    colAlpha = colAlpha
  )
  ### Save
  ggplot2::ggsave(
    filename = paste0("Volcano plot - ", title, ".png"),
    plot = v,
    path = path,
    width = width,
    height = height,
    limitsize = FALSE
  )

  ### Print
  print(v)
}
