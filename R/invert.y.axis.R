### invert.y.axis

invert.y.axis <- function(x,
                          y.axis.name
){x[[y.axis.name]] <- -1*abs(x[[y.axis.name]])}


### arcsinh.transform

#arcsinh.transform <- function(x,
#                              s){
#
#      col.names.dl <- names(data_subset)
#      col.names.SCALE <- col.names.dl[col.nos.scale]
#      data_subset[, col.names.SCALE] <- asinh(data_subset[, col.names.SCALE] / asinh.scale)
#      head(data_subset)
#      summary(data_subset)
#
#      }
