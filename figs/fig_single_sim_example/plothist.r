#!/usr/bin/Rscript --vanilla

library("lattice")
library("grid")
library("colorRamps")
source("/home/bram/R/src/bramlib.r")

line.par <- 5002

heights <- c(0.2, 1, 0.05,1,0.05,1,0.05,1,0.3)
widths <- c(0.2,1,0.2)

filename <- "sim_conflict_25_2_2014_124253_108543000"

pred.file <- "../fig_resol/summary_n4.csv"

dat <- read.table(file=filename,header=T,sep=";",nrow=5000)
pred.dat <- read.table(file=pred.file,header=T,sep=";")

pred.dat.sub <- pred.dat[pred.dat$c1==0.5 & pred.dat$c2==0.65,]

lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)

tick.cex <- 0.5
label.cex <- 0.8
lwd <- 0.3

panel.xh <- function(x, y, dat, yvar, colors,varvar=NA,pred.vars,pred.colors, ...)
{
    for (pred_i in 1:length(pred.vars))
    {
        panel.abline(a=pred.vars[[pred_i]],b=0,col=pred.colors[[pred_i]])
    }


    for (i in 1:length(yvar))
    {
        if (!is.na(varvar))
        {
            panel.polygon(x=c(dat$generation,rev(dat$generation)),
                            y=c(dat[,yvar[[i]]]+sqrt(dat[,varvar[[i]]]),rev(dat[,yvar[[i]]]+sqrt(dat[,varvar[[i]]]))),
                            fill=colors[[i]],
                            border="transparent",
                            alpha=0.5)
        }

        panel.xyplot(x=dat$generation,
                        y=dat[,yvar[[i]]],
                        type="l",
                        col=colors[[i]],
                        lwd=1.0)
    }
}
    

init.plot(filename="sim_example",width=350,height=500)
    pushViewport(viewport(layout=lo))

    xp <- xyplot(d0 ~ generation,
                        data=dat,
                        ylim=c(0,1),
                        dat=dat,
                        yvar=c("d0","d1"),
                        colors=c("blue","red"),
                        pred.vars=pred.dat.sub[,c("d1","d2")],
                        pred.colors=c("black","green"),
                        panel=panel.xh)

    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.xh",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent",lwd=lwd))

    upViewport()

    xp <- xyplot(d0 ~ generation,
                        data=dat,
                        ylim=c(0,1),
                        dat=dat,
                        yvar=c("s0","s1"),
                        colors=c("blue","red"),
                        pred.vars=pred.dat.sub[,c("s1","s2")],
                        pred.colors=c("black","green"),
                        panel=panel.xh)

    pushViewport(viewport(layout.pos.row=4,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.xh",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent",lwd=lwd))

    upViewport()


    xp <- xyplot(d0 ~ generation,
                        data=dat,
                        ylim=c(0,1),
                        dat=dat,
                        yvar=c("dsignal","dnosignal"),
                        colors=c("blue","red"),
                        pred.vars=pred.dat.sub[,c("dS","dNS")],
                        pred.colors=c("black","green"),
                        panel=panel.xh)

    pushViewport(viewport(layout.pos.row=6,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.xh",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent",lwd=lwd))

    upViewport()

    pushViewport(viewport(layout.pos.row=3,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=F,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text="")
    upViewport()
    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        labels=T,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("signal envt 1", italic(s)[1],sep="")))
    upViewport()

upViewport()
dev.off()
