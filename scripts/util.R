## https://medium.com/analytics-vidhya/ggplot2-themes-for-publication-ready-plots-including-dark-themes-9cd65cc5a7e3
## https://rpubs.com/Koundy/71792

theme_Publication <- function(base_size=18,legend.position = "bottom"){
    #, base_family="helvetica") {
    suppressMessages(suppressWarnings(library (ggthemes)))
      (theme_foundation(base_size=base_size)#, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(family="mono"),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               # panel.grid.major = element_line(colour="#f0f0f0"),
               # panel.grid.minor = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = legend.position,
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
}


plot_scatter <- function(df,x_lab,y_lab){
    m = df[,c(x_lab,y_lab)]
    colnames(m) <- c('x','y')
    
    dens <- kde2d(df[,x_lab], df[,y_lab])
    gr <- data.frame(with(dens, expand.grid(x,y)), as.vector(dens$z))

    names(gr) <- c("xgr","ygr",'zgr')

    colfunc <- colorRampPalette(c(rep("slateblue1", each=10), rep("slateblue3", each=20), rep("slateblue4", each=1)))

    s<-cor.test(df[,x_lab], df[,y_lab], method="pearson")
    # https://stackoverflow.com/questions/25367943/extract-verbatim-p-value-from-cor-test
    pv = ifelse(s$p.value==0,-log(2.2e-16,base=10),-log(s$p.value,base=10))
    ann <- sprintf("Rho = %.3g\n-log10(P) = %.4g", s$estimate, pv)

    mod <- loess(zgr~xgr*ygr, data=gr)

    m$pointdens <- predict(mod, newdata=data.frame(xgr=x, ygr=y))

    m %>% ggplot(aes(x=x,y=y, color=pointdens)) + 
        geom_point(size=1, alpha=0.4,show.legend = FALSE) +
        geom_smooth(method=lm,level=0.9999999, size=1, se=FALSE) +
        scale_colour_gradientn(colours=colfunc(30),guide = "none") +
        xlab(x_lab) + ylab(y_lab) +
        annotate("text", x=0,y=0, label=ann,size=8)+
        theme_bw(30) + 
        theme(panel.background = element_rect(colour = "black")) + 
        theme(panel.background = element_rect(colour = "black"), panel.grid.minor = element_blank()) +
        labs(colour="") -> p
    return (p)
}


scale_fill_Publication <- function(...){
    suppressMessages(suppressWarnings(library (scales)))
    discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}


scale_colour_Publication <- function(...){
    suppressMessages(suppressWarnings(library (scales)))
    discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}


plot_Save <- function (p, name_it,w=NA, h=NA){
    ggsave(paste(name_it,'png',sep='.'), plot = p, device = 'png', width = w, height=h, dpi = 300)
    ggsave(paste(name_it,'pdf',sep='.'), plot = p, device = 'pdf', width = w, height=h, dpi = 300)
}
