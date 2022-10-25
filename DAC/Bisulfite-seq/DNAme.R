library(methylKit)


plot_hyperhypo <- function (res, meth.diff.cutoff, qvalue.cutoff, title='', 
                            x_lim=FALSE,y_lim=NA,alpha = 5/10, size = 3,hypo_test = FALSE, mu = 0) {
    if (x_lim==FALSE){
        x_min = res$meth.diff %>% min
        x_max = res$meth.diff %>% max
    } else if (length(x_lim)==2) {
        x_min = x_lim[1]
        x_max = x_lim[2]
    }
    colors_values = c("grey80", "#FF4500") # https://www.rapidtables.com/web/color/RGB_Color.html

    res$sig <- as.factor(res$qvalue < qvalue.cutoff & abs(res$meth.diff) > meth.diff.cutoff)
    relevel(res$sig, ref=TRUE)
    
    vol <- res %>% ggplot(aes(x = meth.diff, y = -log10(qvalue))) + 
        geom_point(
            data = res %>% filter(
                abs(meth.diff) >= meth.diff.cutoff, 
                qvalue < qvalue.cutoff), 
            aes(x = meth.diff, y = -log10(qvalue)), 
            size = size, alpha = alpha, shape = 21, 
            color = colors_values[2],
            fill = colors_values[2]) + 

        geom_point(
            data = res %>% filter(
                (abs(meth.diff) < meth.diff.cutoff) | 
                (qvalue >= qvalue.cutoff)), 
            aes(x = meth.diff, y = -log10(qvalue)), 
            size = 1, alpha = alpha, 
            shape = 21, color = colors_values[1], 
            fill = colors_values[1]) + 
        
        geom_rug(
            data = res %>% filter(
                abs(meth.diff) >= meth.diff.cutoff, 
                qvalue < qvalue.cutoff), 
            alpha = alpha,
            sides = "b", colour= colors_values[2]) + 

        xlim(c(x_min, x_max)) + 
        ylim(c(0, y_lim)) + 
        # ylim(c(0, -log10(res %>% select(qvalue) %>% min ))) + 
        geom_hline(yintercept = -log10(qvalue.cutoff), linetype = "dashed", alpha = 5/10) + 
        ggtitle(title) 
        if (hypo_test){
            w <- wilcox.test(res[res$qvalue < qvalue.cutoff,]$meth.diff, mu=mu, alternative = "less")
            t <- t.test(res[res$qvalue < qvalue.cutoff,]$meth.diff, mu=mu, alternative = "less")

            vol +
            # \nwilcox.test (-log10 p.value): %.2f
            geom_text(aes(
                x_max*0.5,-log10(qvalue.cutoff)/1000,
                label = sprintf(
                    "\n[mu=%i,alter=less]\nt.test (-log10 p.value):%.2f",
                    mu, -log10(t$p.value)),
                vjust = x_min*0.05
            )) -> vol
        }
        
    return(vol + theme(legend.position = "none") + theme_Publication())
}


prepCounts <- function(me){
    # Filter the summarized counts by coverage
    me_filt <- methylKit::filterByCoverage(me,
                          lo.count=10,
                          lo.perc=NULL,
                          hi.count=NULL,
                          hi.perc=99.9)

    # Perform simple normalization
    me_filt_norm <- methylKit::normalizeCoverage(me_filt, method = "median")

    # Merge the samples again
    return(methylKit::unite(me_filt_norm, destrand=FALSE))
}


getDiffList <- function(myDiff,difference.thr=15,qvalue.thr=0.05){
    # get hyper methylated bases and order by qvalue
    myDiff15p.hyper <- methylKit::getMethylDiff(
        myDiff,
        difference=difference.thr,
        qvalue=qvalue.thr,
        type="hyper")
    myDiff15p.hyper <- myDiff15p.hyper[order(myDiff15p.hyper$qvalue),]

    # get hypo methylated bases and order by qvalue
    myDiff15p.hypo <- methylKit::getMethylDiff(
        myDiff,
        difference=difference.thr,
        qvalue=qvalue.thr,
        type="hypo")
    myDiff15p.hypo <- myDiff15p.hypo[order(myDiff15p.hypo$qvalue),]

    # get all differentially methylated bases and order by qvalue
    myDiff15p <- methylKit::getMethylDiff(
        myDiff,
        difference=difference.thr,
        qvalue=qvalue.thr
    )
    myDiff15p <- myDiff15p[order(myDiff15p$qvalue),]
    row.names(myDiff15p) <- 1:nrow(myDiff15p)

    DiffList = list()
    DiffList$hyper = myDiff15p.hyper
    DiffList$hypo = myDiff15p.hypo
    DiffList$all = myDiff15p
    
    return(DiffList)
}


getPromoterDiff <- function(myDiff,myDiff.anot, tss.window=2500){
    # make sure you globally defined `gene.obj` and `tx2name`
    myDiff_dist_tss <- cbind(
        myDiff,
        # getAssociationWithTSS and filter into given window
        genomation::getAssociationWithTSS(myDiff.anot)) %>% 
            dplyr::filter(abs(dist.to.feature) < tss.window) %>% 
            dplyr::select(!target.row)
    myDiff_dist_tss$feature.name = tx2name[myDiff_dist_tss$feature.name,]$gene_name
    
    return(myDiff_dist_tss)
    
}


annotateBedFromDb <- function(gRanges = NULL, db = NULL, window = 5000){
    # https://rdrr.io/bioc/CompGO/man/annotateBedFromDb.html
    genes = transcriptsByOverlaps(ranges = bed, x = db, maxgap=window, columns=c('tx_id', 'gene_id'))
    return(genes)
}