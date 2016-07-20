#!/usr/bin/Rscript --vanilla

library("lattice")
library("grid")
library("colorRamps")
source("/home/bram/R/src/bramlib.r")

line.par <- 10002

heights <- c(0.2, 1, 0.05,1,0.05,1,0.05,1,0.05,1,0.05,1,0.3)
widths <- c(0.1,1,0.2)

lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)

tick.cex <- 0.5
label.cex <- 0.8
main.tck <- -0.4

file.outcomes <- "/home/bram/Projects/conflict_dispersal/figs/fig_resol/summary_n4.csv"

outcomes <- read.table(file.outcomes,header=T,sep=";")
outcomes <- outcomes[outcomes$p1 == 0.5,]


flist <- commandArgs(trailingOnly=T)[[1]]

panel.lvp <- function(x,y,z,a_ab,...)
{
    panel.levelplot(x=x,y=y,z=z,...)
    panel.abline(a=a_ab,b=0,col="red",lty=2)
}

for (file in flist)
{
    print(file)
    data <- read.table(file,sep=";",header=T)
    print(gsub("(.*)dist$","\\1",file))
    par.data <- read.table(gsub("(.*)dist$","\\1",file),sep=";",header=F,skip=line.par)
    names(par.data) <- c("name","val")

    current.outcome <- outcomes[outcomes$c1 == par.data[par.data$name == "c0","val"] &
                                    outcomes$c2 == par.data[par.data$name == "c1","val"],]

    if (nrow(current.outcome) > 1)
    {
        print("falen.")
    }

    pdf(file=paste("branchplot_",basename(file),".pdf",sep=""),width=10,height=10)
    pushViewport(viewport(layout=lo))
    xp <- levelplot(log10(s0+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=1,s=seq(0,1,0.01),v=1),
                        a_ab=current.outcome$s1,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=3,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=F,
                        tck=main.tck,
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
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("signal envt 1", italic(s)[1],sep="")))
    upViewport()
    xp <- levelplot(log10(s1+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=240/360,s=seq(0,1,0.01),v=1),
                        a_ab=current.outcome$s2,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=4,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=5,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        tck=main.tck,
                        labels=F,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text="")
    upViewport()
    pushViewport(viewport(layout.pos.row=4,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        tck=main.tck,
                        labels=T,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("signal envt 2", italic(s)[2],sep="")))
    upViewport()




    xp <- levelplot(log10(dS+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=199/360,s=seq(0,1,0.01),v=1),
                        a_ab=current.outcome$dS,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=6,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=7,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=F,
                        cex=tick.cex,
                        tck=main.tck,
                        labelcex=label.cex,
                        text="")
    upViewport()
    pushViewport(viewport(layout.pos.row=6,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        labels=T,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("", italic(d)["S"],sep="")))
    upViewport()

    xp <- levelplot(log10(dNS+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=130/360,s=seq(0,1,0.01),v=1),
                        a_ab=current.outcome$dNS,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=8,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=8,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        labels=T,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("", italic(d)["NS"],sep="")))
    upViewport()

    xp <- levelplot(log10(d0+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=37/360,s=seq(0,1,0.01),v=1),
,
                        a_ab=current.outcome$d1,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=10,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=11,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=F,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text="")
    upViewport()
    pushViewport(viewport(layout.pos.row=10,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        labels=T,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("", italic(d)[1],sep="")))
    upViewport()



    xp <- levelplot(log10(d1+1) ~ generation * bin,
                        data=data,
                        col.regions=hsv(h=306/360,s=seq(0,1,0.01),v=1),
                        a_ab=current.outcome$d2,
                        panel=panel.lvp)

    pushViewport(viewport(layout.pos.row=12,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            do.call("panel.lvp",trellis.panelArgs(xp,1))
            grid.rect(gp=gpar(fill="transparent"))

    upViewport()
    pushViewport(viewport(layout.pos.row=13,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=T,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text="generation, t")
    upViewport()
    pushViewport(viewport(layout.pos.row=12,
                            layout.pos.col=1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$y.limits,
                        side="right",
                        labels=T,
                        tck=main.tck,
                        cex=tick.cex,
                        labelcex=label.cex,
                        text=expression(paste("", italic(d)[2],sep="")))
    upViewport()
    upViewport()
    dev.off()
}
