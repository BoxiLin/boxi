#' PP plot GWAS
#' 
#' Manhattan plot
#'
#' @param dt GWAS summary statisics
#' @param p p values column name
#' 
#' @examples # png(paste("plot/", model,maf, ".png", sep=""),
#' png(paste("plot/test",maf, ".png", sep=""),
#'    height = 480, width = 1000, type="cairo")
#' boxi_manhattan(dt1, p = "Score.pval",
#'               check = 1,
#'               maf = maf, thr = 0.01,chr = "chr", bp = "pos", snp = "variant.id")
#' dev.off()
boxi_manhattan <- function(dt, p = "P", thr = 0.001,maf = 0, check = 1,freq = "MAF",
                           title = "", chr = "CHR", bp = "BP", snp = "SNP"){
  p_val <- dplyr::pull(dt, p)
  p_val<- p_val[!is.na(p_val)]
  frq <- dplyr::pull(dt, freq)
  dt_plot <- dt %>% dplyr::filter(p_val < thr) %>% dplyr::filter(freq>maf)
  man<- function(){qqman::manhattan(dt_plot, chr = chr, bp = bp, p = p, snp = snp,ylim = c(1,10),
                                    col = c("blue4", "orange3"),suggestiveline = FALSE,
                                    main = title)}
  if (check) {
    pp <- p_val[frq>maf]
    histo<- function() {hist(pp, freq = 0, main = "",xlab = "p values")}
    p_pp = boxi::qqunif.plot(pp)
    bottom_row <- cowplot::plot_grid(histo, p_pp, labels = c('Histogram of p-values', 'PP-plot'), label_size = 12,ncol = 1)
    plot_grid(man, bottom_row,
              label_size = 12, ncol = 2,rel_widths = c(3, 1))
  } else {man()}
}