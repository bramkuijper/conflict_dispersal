library("MASS")
library("lattice")
source("~/R/src/bramlib.r")

if (!exists("dat"))
{
    dat <- read.table("summary_rufus.csv",sep=";",header=T)
    str(dat)
}

heights <- c(0.2, 1, 0.3)
widths <- c(0.3,1,0.3)

line.lwd <- 0.75

lwd <- 0.5
tick.cex <- 0.7
label.cex <-1.0 
legend.cex <- 0.8
plot.tck <- -0.25
ind.label.cex <- 1.0


lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)

init.plot(filename="fig_sim_resol",width=650,height=650,font="helvetica",type="pdf")
pushViewport(viewport(layout=lo))
k <- trellis.par.get("box.3d")

k$lwd <- 0.5

trellis.par.set("box.3d",k)

xp <- wireframe(d0 ~ c0 * c1,
                aspect=c(1,0.5),
                data=dat,
                zlim=c(0,1),
                xlim=c(0,1),
                ylim=c(0,1),
                cex=0.5,
                zoom=1.25,
                zlab="",
                xlab="",
                ylab="",
                lwd=0.1,
                drape=F,
                shade=T,
                #             shade.colors.palette = function(irr, ref, height, w = .5) {
                #                      grey(w * irr + (1 - w) * (1 - (1-ref)^.4)) },
                #                alpha=0.0,
                col="transparent",
                #screen=list(z=60,x=-20,y=-50),
                screen=list(z=-20,x=-55,y=-20),
                # par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                scales=list(x=list(cex=tick.cex,at=pretty(range(dat$c1),n=6)),
                        y=list(cex=tick.cex,at=pretty(range(dat$c0),n=6)),
                        z=list(cex=tick.cex,at=pretty(range(dat$d1),n=6))
                        ),
                #                col.regions=grey(seq(0,1,1/100)),
                alpha.regions=rep(1.0,times=100),
                default.scales=list(arrows=F,lwd=0.5))

    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.wireframe",trellis.panelArgs(xp,1))

        #        grid.text(x=-0.2,y=1.05,label=label,gp=gpar(cex=ind.label.cex*1.1))
        grid.text(x=0.29,y=0.03,label=expression(paste("dispersal costs patch 1, ",italic(c)[1])),rot=-20, gp=gpar(cex=label.cex))
        grid.text(x=0.96,y=0.2,label=expression(paste("dispersal costs patch 2, ",italic(c)[2])),rot=52,gp=gpar(cex=label.cex))
        grid.text(x=-0.15,y=0.35,label=expression(paste("dispersal from patch 1, ",italic(d)[1])),rot=105,gp=gpar(cex=label.cex))
        grid.text(x=.5,y=.97,label=expression(paste("dispersal when offspring rely on maternal signal")),rot=0,gp=gpar(cex=ind.label.cex))
    upViewport()


    upViewport()

exit.plot() 
