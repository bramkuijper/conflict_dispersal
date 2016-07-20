library("MASS")
library("lattice")
source("~/R/src/bramlib.r")


if (!exists("dat"))
{
    dat <- read.table("summary_overview_p.csv",sep=";",header=T)

    dat <- dat[dat$p1 %in% c(0.1, 0.25) & dat$n == 5,]

    if (nrow(dat) < 1)
    {
        stop("no data.")
    }
}

# calculate rm of all values

heights <- c(0.2, 1, 0.3,1,0.2)
widths <- c(0.5,1,0.5,1,0.05)

line.lwd <- 0.75

lwd <- 0.5
tick.cex <- 0.5
sub.label.cex <- 0.7
label.cex <- 0.8
legend.cex <- 0.8
plot.tck <- -0.25
ind.label.cex <- 1.0


lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)

init.plot(filename="fig_overview_p",width=650,height=575,font="helvetica",type="pdf")
pushViewport(viewport(layout=lo))
k <- trellis.par.get("box.3d")

k$lwd <- 0.5

trellis.par.set("box.3d",k)



make.ticks <- function(row,col, xp)
{

        pushViewport(viewport(layout.pos.row=row,
                                layout.pos.col=col-1,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            single.axis(range=xp$y.limits,
                            side="right",
                            tck=plot.tck,
                            labels=T,
                            cex=tick.cex,
                            labelcex=label.cex,
                            nsub=5,
                            x.text.off=0.7,
                            text=expression(paste("dispersal costs patch 2, ",italic(c)[2])))
        upViewport()
        pushViewport(viewport(layout.pos.row=row,
                                layout.pos.col=col+1,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            single.axis(range=xp$y.limits,
                            side="left",
                            tck=plot.tck,
                            labels=F,
                            cex=tick.cex,
                            labelcex=label.cex,
                            nsub=5,
                            text="")
        upViewport()
        pushViewport(viewport(layout.pos.row=row-1,
                                layout.pos.col=col,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            single.axis(range=xp$x.limits,
                            side="bottom",
                            tck=plot.tck,
                            labels=F,
                            cex=tick.cex,
                            labelcex=label.cex,
                            nsub=5,
                            text="")
        upViewport()
        pushViewport(viewport(layout.pos.row=row+1,
                                layout.pos.col=col,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            single.axis(range=xp$x.limits,
                            side="top",
                            tck=plot.tck,
                            labels=T,
                            cex=tick.cex,
                            labelcex=label.cex,
                            nsub=5,
                            text=expression(paste("dispersal costs patch 1, ",italic(c)[1])))
        upViewport()
}

dat.sub <- dat[dat$p1 == 0.1,]
xp <- wireframe(d1 ~ c1 * c2,
                aspect=c(1,0.5),
                data=dat.sub,
                zlim=c(0,1),
                xlim=c(0,1),
                ylim=c(0,1),
                #cex=0.5,
                zoom=1.25,
                zlab="",
                xlab="",
                ylab="",
             shade.colors.palette = function(irr, ref, height, w = .5) {
                      grey(w * irr + (1 - w) * (1 - (1-ref)^.4)) },
                #                alpha=0.0,
                lwd=0.1,
                drape=T,
                shade=T,
                #col="transparent",
                #screen=list(z=60,x=-20,y=-50),
                screen=list(z=-20,x=-55,y=-20),
                # par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                scales=list(x=list(cex=tick.cex,at=pretty(range(dat$c1),n=6)),
                        y=list(cex=tick.cex,at=pretty(range(dat$c2),n=6)),
                        z=list(cex=tick.cex,at=pretty(range(dat$d1),n=6))
                        ),
                col.regions=blue2green2red(300),
                alpha.regions=rep(1.0,times=300),
                default.scales=list(arrows=F,lwd=0.5))


    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.wireframe",trellis.panelArgs(xp,1))

        label <- "A"
        grid.text(x=-0.2,y=.97,label=label,gp=gpar(cex=ind.label.cex*1.1))
        grid.text(x=0.29,y=0.03,label=expression(paste("dispersal costs patch 1, ",italic(c)[1])),rot=-20, gp=gpar(cex=label.cex))
        grid.text(x=0.96,y=0.2,label=expression(paste("dispersal costs patch 2, ",italic(c)[2])),rot=52,gp=gpar(cex=label.cex))
        grid.text(x=-0.15,y=0.35,label=expression(paste("dispersal from patch 1, ",italic(d)[1])),rot=105,gp=gpar(cex=label.cex))
        grid.text(x=.5,y=1.05,label=expression(paste("dispersal when offspring rely on maternal signal")),rot=0,gp=gpar(cex=ind.label.cex))
        grid.text(x=-0.4,y=0.97,label=expression(paste(bolditalic(p)[bold("1")],bold(" = 0.1"))),rot=0, gp=gpar(cex=label.cex*1.1))
    upViewport()


    # function to draw smoothed contours
    panel.cp <- function(x,y,z,the.data,texts=list(a=list(label="fully informative",x=0.1,y=0.9),
                    b=list(label="partially informative",x=0.15,y=0.7,rot=20),
                    c=list(label="uninformative",x=0.38,y=0.5),
                    d=list(label="partially informative",x=0.28,y=0.02,rot=40),
                    e=list(label="fully informative",x=0.8,y=0.2,rot=80)),...)
    {
        # grey background color of the square
        panel.polygon(x=c(0,0,1,1),y=c(0,1,1,0),
                        border="transparent",
                        col="grey90")

        c1 <- sort(unique(the.data[,"c1"]))
        c2 <- sort(unique(the.data[,"c2"]))

        matrix.contour <- matrix(data=the.data[,"siginfo"], 
                                    nrow=length(c1), 
                                    ncol=length(c2),
                                    byrow=T)

        # contour for 0.98
        cnt <- contourLines(x=c1,
                            y=c2,
                            matrix.contour,
                            levels=.98)

        contour.part.info <- as.data.frame(cnt[[1]])

        # plot the first one for which y > x
        xy1 <- smooth.spline(x=contour.part.info[,"x"],
                        y=contour.part.info[,"y"],
                        spar=0.6
                        )

        xy2 <- as.data.frame(cbind(x=xy1$x,y=xy1$y))
        xy2 <- xy2[xy2$y <= 1.0,]
        str(xy2)
        print(nrow(xy2))

        panel.polygon(x=c(xy2$x,max(xy2$x),min(xy2$x)),
                       y=c(xy2$y,1.0,1.0),
                       border="transparent",
                       col="grey70")

        panel.xyplot(x=xy2$x,
                        y=xy2$y,
                        type="l",
                        lwd=0.5,
                        col="black")

        # now the one for which x<y 

        contour.part.info <- as.data.frame(cnt[[2]])

        xy1 <- smooth.spline(
                        y=contour.part.info[,"x"],
                        x=contour.part.info[,"y"],
                        cv=T
                        )

        xy2 <- as.data.frame(cbind(x=xy1$y,y=xy1$x))
        xy2 <- xy2[xy2$y <= 1.0,]

        panel.polygon(x=c(xy2$x,max(xy2$x),min(xy2$x)),
                       y=c(xy2$y,0.0,0.0),
                       border="transparent",
                       col="grey70")

        panel.xyplot(x=xy2$x,
                        y=xy2$y,
                        type="l",
                        lwd=0.5,
                        col="black")
 

        ##############################################################

        # contour for .001
        cnt <- contourLines(x=c1,
                            y=c2,
                            matrix.contour,
                            levels=.001)

        contour.no.info <- as.data.frame(cnt[[1]])

        #first part of contour where x > y
        contour.no.info.1 <- contour.no.info[contour.no.info$x <= contour.no.info$y,]
        contour.no.info.2 <- contour.no.info[contour.no.info$x >= contour.no.info$y,]

        xy1 <- smooth.spline(x=contour.no.info.1[,"x"],
                        y=contour.no.info.1[,"y"]
                        )

        xy2 <- smooth.spline(x=contour.no.info.2[,"x"],
                        y=contour.no.info.2[,"y"]
                        )

        panel.polygon(x=c(xy1$x,rev(xy2$x)),
                        y=c(xy1$y,rev(xy2$y)),
                        type="l",
                        lwd=0.5,
                        col="white"
                        )

        # now print all the texts
        for (i in 1:length(texts))
        {
            grid.text(x=texts[[i]]$x,
                    y=texts[[i]]$y,
                    just="left",
                    label=texts[[i]]$label,
                    rot=ifelse(!is.null(texts[[i]]$rot),texts[[i]]$rot,0),
                    gp=gpar(cex=sub.label.cex))
        }
    }

        # contourplot depicting information content in the maternal
        # signal

        xp <- contourplot(siginfo ~ c1 * c2,
                        aspect=c(1.5,0.8),
                        labels=F,
                        data=dat.sub,
                        the.data=dat.sub,
                        panel=panel.cp,
                        at=c(0.001,0.98),
                        zlim=c(0,1),
                        xlim=c(0.0,1.0),
                        ylim=c(0.0,1.0),
                        #cex=0.5,
                        zoom=1.0,
                        zlab="",
                        xlab="",
                        ylab="",
                        alpha=0.0,
                        lwd=0.1,
                        drape=T,
                        pretty=T,
                        #region=T,
                        col="grey40",
                        #                shade=T,
                        #screen=list(z=60,x=-20,y=-50),
                        # par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                        scales=list(x=list(cex=tick.cex,at=pretty(range(dat$c1),n=6)),
                                y=list(cex=tick.cex,at=pretty(range(dat$c2),n=6)),
                                z=list(cex=tick.cex,at=pretty(range(dat$d1),n=6))
                                ),
                        #                        col.regions=c("red","green","blue"),
                        alpha.regions=rep(.9,times=300)[100:300],
                        default.scales=list(arrows=F,lwd=0.5))

        pushViewport(viewport(layout.pos.row=2,
                                layout.pos.col=4,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            do.call("panel.cp",trellis.panelArgs(xp,1))

            grid.rect(gp=gpar(fill="transparent",lwd=.5))

            label <- "B"
            grid.text(x=-0.1,y=.97,label=label,gp=gpar(cex=ind.label.cex))
        grid.text(x=.5,y=1.05,label=expression(paste("maternal signal: information content")),rot=0,gp=gpar(cex=ind.label.cex))
        upViewport()

        make.ticks(2, 4, xp)
        #################################################### p = 0.25 ####################################################


dat.sub <- dat[dat$p1 == 0.25,]

xp <- wireframe(d1 ~ c1 * c2,
                aspect=c(1,0.5),
                data=dat.sub,
                zlim=c(0,1),
                xlim=c(0,1),
                ylim=c(0,1),
                #cex=0.5,
                zoom=1.25,
                zlab="",
                xlab="",
                ylab="",
             shade.colors.palette = function(irr, ref, height, w = .5) {
                      grey(w * irr + (1 - w) * (1 - (1-ref)^.4)) },
                #                alpha=0.0,
                lwd=0.1,
                drape=T,
                shade=T,
                #col="transparent",
                #screen=list(z=60,x=-20,y=-50),
                screen=list(z=-20,x=-55,y=-20),
                # par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                scales=list(x=list(cex=tick.cex,at=pretty(range(dat$c1),n=6)),
                        y=list(cex=tick.cex,at=pretty(range(dat$c2),n=6)),
                        z=list(cex=tick.cex,at=pretty(range(dat$d1),n=6))
                        ),
                col.regions=blue2green2red(300),
                alpha.regions=rep(1.0,times=300),
                default.scales=list(arrows=F,lwd=0.5))


    pushViewport(viewport(layout.pos.row=4,
                            layout.pos.col=2,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.wireframe",trellis.panelArgs(xp,1))

        label <- "C"
        grid.text(x=-0.2,y=.97,label=label,gp=gpar(cex=ind.label.cex*1.1))
        grid.text(x=0.29,y=0.03,label=expression(paste("dispersal costs patch 1, ",italic(c)[1])),rot=-20, gp=gpar(cex=label.cex))
        grid.text(x=0.96,y=0.2,label=expression(paste("dispersal costs patch 2, ",italic(c)[2])),rot=52,gp=gpar(cex=label.cex))
        grid.text(x=-0.15,y=0.35,label=expression(paste("dispersal from patch 1, ",italic(d)[1])),rot=105,gp=gpar(cex=label.cex))
        grid.text(x=-0.38,y=0.97,label=expression(paste(bolditalic(p)[bold("1")],bold(" = 0.25"))),rot=0, gp=gpar(cex=label.cex*1.1))
        #        grid.text(x=.5,y=.97,label=expression(paste("dispersal when offspring rely on maternal signal")),rot=0,gp=gpar(cex=ind.label.cex))
    upViewport()

        texts=list(a=list(label="fully informative",x=0.1,y=0.9),
                    b=list(label="partially informative",x=0.15,y=0.65,rot=20),
                    c=list(label="uninformative",x=0.38,y=0.5),
                    d=list(label="partially informative",x=0.3,y=0.02,rot=40),
                    e=list(label="fully informative",x=0.8,y=0.2,rot=80))
        xp <- contourplot(siginfo ~ c1 * c2,
                        aspect=c(1.5,0.8),
                        labels=F,
                        data=dat.sub,
                        texts=texts,
                        the.data=dat.sub,
                        panel=panel.cp,
                        at=c(0.001,0.98),
                        zlim=c(0,1),
                        xlim=c(0.0,1.001),
                        ylim=c(0.0,1.001),
                        #cex=0.5,
                        zoom=1.0,
                        zlab="",
                        xlab="",
                        ylab="",
                        alpha=0.0,
                        lwd=0.1,
                        drape=T,
                        pretty=T,
                        #region=T,
                        col="grey40",
                        #                shade=T,
                        #screen=list(z=60,x=-20,y=-50),
                        # par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                        scales=list(x=list(cex=tick.cex,at=pretty(range(dat$c1),n=6)),
                                y=list(cex=tick.cex,at=pretty(range(dat$c2),n=6)),
                                z=list(cex=tick.cex,at=pretty(range(dat$d1),n=6))
                                ),
                        #                        col.regions=c("red","green","blue"),
                        alpha.regions=rep(.9,times=300)[100:300],
                        default.scales=list(arrows=F,lwd=0.5))

        pushViewport(viewport(layout.pos.row=4,
                                layout.pos.col=4,
                                xscale=xp$x.limits,
                                yscale=xp$y.limits
                                ))

            do.call("panel.cp",trellis.panelArgs(xp,1))

            grid.rect(gp=gpar(fill="transparent",lwd=.5))

            label <- "D"
            grid.text(x=-0.1,y=.97,label=label,gp=gpar(cex=ind.label.cex))
        upViewport()
        make.ticks(4, 4, xp)

upViewport()
exit.plot() 
