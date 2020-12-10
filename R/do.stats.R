#' Statistical analysis
#'
#' @import data.table
#'
#' @export

do.stats <- function(dat,
                     use.cols,
                     sample.col,
                     grp.col,

                     comparisons, # make a tool to create a factorial comparison design -- for now just specify manually

                     variance.test = "kruskal.test", ## add ANOVA
                     pairwise.test = "wilcox.test", ## Add t-test
                     corrections = "fdr"){

  ### Test data
      # dat <- merged.dat
      # names(merged.dat)
      # use.cols <- names(dat)[c(2:147)]
      # sample.col <- "Sample"
      # grp.col <- "DiseaseSeverity"
      # variance.test = "kruskal.test"
      # pairwise.test = "wilcox.text"
      # corrections = 'fdr'
      # comparisons = list(c("COVID_ward", "COVID_ICU"))

  ### Setup

      grps <- unique(dat[[grp.col]])
      ngrps <- c(length(unique(dat[[grp.col]])))

  ### Stat assessment

      res.list <- list()

      for(i in c(1:length(use.cols))){
        # i <- 1
        nm <- use.cols[[i]]
        strt <- dat[[grp.col]]
        temp <- dat[[nm]]
        temp <- data.table("Group" = strt, "Value" = temp)

        k.res <- kruskal.test(x = c(temp[[2]]), g = c(temp[[1]]))
        k.res <- as.data.frame(k.res$p.value)

        colnames(k.res) <- nm
        rownames(k.res) <- "Kruskal"

        ## Comparisons
        w.res.list <- list()

        for(a in comparisons){
          # a <- comparisons[[1]]

          grp1nme <- a[[1]]
          grp2nme <- a[[2]]

          grp1 <- temp[temp[["Group"]] == grp1nme,]
          grp2 <- temp[temp[["Group"]] == grp2nme,]

          w.res <- wilcox.test(grp1[[2]], grp2[[2]])
          w.res.list[[paste0(a[[1]], " to ", a[[2]])]] <- w.res$p.value
        }

        res.w <- cbind(w.res.list)
        colnames(res.w) <- nm

        res <- rbind(k.res, res.w)
        res

        res.list[[nm]] <- res

      }

      new.rows <- do.call(cbind, res.list)
      n <- rownames(new.rows)

      new.rows <- as.data.table(new.rows)

      final.res <- data.table("Comparison" = n,
                              "Type" = rep("p-value", nrow(new.rows)),
                              new.rows)

      names(final.res)[c(1:2)] <- c("Comparison", "Type")
      final.res

      if(!is.null(corrections)){

        nms.save <- names(final.res)
        init.cols <- final.res[,c(1:2)]

        for(a in c(1:(length(names(init.cols))))){
          init.cols[[a]] <- gsub("p-value", paste0("p-value_", corrections), init.cols[[a]])
        }

        adjust.res.list <- list()

        for(i in c(1:nrow(final.res))){
          temp <- final.res[i,c(3:length(names(final.res))), with = FALSE]
          temp <- unlist(temp)
          temp <- unname(temp)

          p.adj <- p.adjust(temp,
                            method = corrections,
                            n = length(temp))

          #plot(temp, p.fdr)
          adjust.res.list[[i]] <- p.adj
        }

        adjust.res <- do.call(rbind, adjust.res.list)
        adjust.res <- as.data.table(adjust.res)

        adjust.res <- cbind(init.cols, adjust.res)
        names(adjust.res) <- nms.save

        final.res <- adjust.res
      }


      #final.res <- list(dat, new.rows)
      #final.res <- rbindlist(final.res, fill = TRUE)

      # tu <- dat[[grp.col]]
      # tu <- as.vector(tu)
      #
      # final.res[[col.append]] <- c(tu, n)

  ### Magnitude

      mag.list <- list()

      for(i in c(1:length(use.cols))){
        # i <- 1
        nm <- use.cols[[i]]
        strt <- dat[[grp.col]]
        temp <- dat[[nm]]
        temp <- data.table("Group" = strt, "Value" = temp)

        mg.comp.list <- list()

        for(a in comparisons){

          grp1nme <- a[[1]]
          grp2nme <- a[[2]]

          grp1 <- temp[temp[["Group"]] == grp1nme,]
          grp2 <- temp[temp[["Group"]] == grp2nme,]

          ctrl.mean <- colMeans(grp1[,2,with = FALSE], na.rm = TRUE)
          comp.mean <- colMeans(grp2[,2,with = FALSE], na.rm = TRUE)

          fold.res <- (comp.mean/ctrl.mean)
          fold.lg2 <- log(fold.res, 2)

          mg.comp.list[[paste0(a[[1]], " to ", a[[2]])]] <- fold.lg2
        }

        mg.comp.list

        mg.res <- cbind(mg.comp.list)
        colnames(mg.res) <- nm

        mag.list[[nm]] <- mg.res
      }

      mag.rows <- do.call(cbind, mag.list)
      n.mag <- rownames(mag.rows)

      n.mag <- as.data.table(n.mag)
      mags <- data.table("Comparison" = n.mag,
                         "Type" = rep("FClog2", nrow(n.mag)),
                         mag.rows)

      names(mags)[c(1:2)] <- c("Comparison", "Type")

  ### Multiple comparisons corrections

      # Coming soon

  ### Wrap up

      rtn <- rbindlist(list(mags, final.res))
      rtn

      # ordered.labs <- list()
      #
      # for(o in c(1:length(comparisons))){
      #   #o <- 1
      #   x <- comparisons[[o]][1]
      #   y <- comparisons[[o]][2]
      #   lab <- paste0(x, " to ", y)
      #   ordered.labs[[o]] <- lab
      # }
      #
      # cbind(ordered.labs)


      # cm <- rbindlist(comparisons)
      # cm <- do.call(cbind, comparisons)

      rtn.k <- rtn[rtn[["Comparison"]] == "Kruskal",]
      rtn.other <- rtn[rtn[["Comparison"]] != "Kruskal",]

      rtn.fc <- rtn.other[rtn.other[["Type"]] == "FClog2",]

      if(is.null(corrections)){
        rtn.p <- rtn.other[rtn.other[["Type"]] == "p-value",]
      }

      if(!is.null(corrections)){
        rtn.p <- rtn.other[rtn.other[["Type"]] == paste0("p-value_", corrections),]
      }

      rtn.fc$Comparison.numbers <- c(1:nrow(rtn.fc))
      rtn.p$Comparison.numbers <- c(1:nrow(rtn.p))

      re.merge <- rbind(rtn.fc, rtn.p)
      setorderv(re.merge, "Comparison")
      setorderv(re.merge, "Comparison.numbers")

      re.merge

      re.merge$Comparison.numbers <- NULL

      all.dat <- rbind(re.merge, rtn.k)

  return(all.dat)

  ###
  #
  # n <- 5
  # r <- 2
  #s
  # n.comp <- factorial(n)/(factorial(r)*factorial(n-r))
  #
  # n.feat <- ncol()
  #
  # p.adjust(p = , method = , n = (n.comp + n.feat))

  ###
}

