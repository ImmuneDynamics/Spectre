#' do.convert.mask.to.outline - convert mask rasters into polygons to create cell outlines
#'
#' This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.
#'
#' @export

do.convert.mask.to.outline <- function(dat){
  #dat < - spatial.list[["TAXXX 20190705 Liver 243788_s1_p1_r1_a1_ac"]][["TAXXX 20190705 Liver 243788_s1_p1_r1_a1_ac_ilastik_s2_Probabilities_mask.tiff"]]

  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  #message("Converting raster to polygons. Please be patient, this may take some time.")
  pp <- rasterToPolygons(dat, dissolve=TRUE)

  #outline <- fortify(pp)
  outline <- pp
  return(outline)
}
