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
    bottom_row <- plot_grid(histo, p_pp, labels = c('Histogram of p-values', 'PP-plot'), label_size = 12,ncol = 1)
    plot_grid(man, bottom_row,
              label_size = 12, ncol = 2,rel_widths = c(3, 1))
  } else {man()}
}




#### Manhattan #######
models <- c("pa", "pac", "pa_int", "pac_int", "pa_iint", "pac_iint")

model <- "pac_residual_all_sib"

dt1 <-  readRDS(paste("gwas_summary/", model,".rds", sep="")) %>% dplyr::filter(P_HWE>1e-50)
dt1 <-  data.table::fread("../../summary_stat/GWAS_discovery_replication/pa_residual_eur_sib.txt") %>% dplyr::filter(P_HWE>1e-50)
dt1$freq<-dt1$MAC/(2*dt1$n.obs)

dt1 <-data.table::fread("../fench_replication/qc/out/french_eur_ethnics_gwas.txt") %>%
  dplyr::filter(P.HWE>1e-50)


maf = 0.01
# png(paste("plot/", model,maf, ".png", sep=""),
png(paste("plot/french_chronic_pa_eur_maf",maf, ".png", sep=""),
    height = 480, width = 1000, type="cairo")
lambda = median(qchisq(dt1$P[!is.na(dt1$P)], df=1, lower.tail=FALSE)) / 0.456
boxi_manhattan(dplyr::filter(dt1, !is.na(P)), p = "P",
               check = 1,
               maf = maf, thr = 0.01,chr = "CHROM", bp = "POS", snp = "ID",
               title = paste("GC Lambda:", round(lambda,3)))
dev.off()

