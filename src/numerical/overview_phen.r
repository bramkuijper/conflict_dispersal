library("lattice")
library("grid")
library("colorRamps")
source("/home/bram/R/src/bramlib.r")

rect.lwd <- 0.75
annotcex <- 0.8
annotcex.2 <- 0.85
annotcex.3 <- 0.9
annotcex.4 <- 1
legendcex <- 0.8
tickcex <- 0.7
maincex <- 0.9
labelcex <- 1.1
circle.x.offset <- 0.1
circle.y.offset <- 0.9
circle.rad <- 0.05
circle.small.rad <- 0.03

color <- F

intermed.val <- function(val)
{
    p1 <- as.numeric(val["p1"])
    p2 <- as.numeric(val["p2"])

    smval <- choices[choices$a==ifelse(p1 <= 0.05,1,ifelse(p1>=0.95,2,3)) & choices$b==ifelse(p2 <= 0.05,1,ifelse(p2>=0.95,2,3)),"out"]

    return(smval)
}

if (color)
{
    main.colors <- c("#3bb5ff","yellow","blue","#ff6900","#76d06b","darkred","green","#b800bb","red")
} else
{
    main.colors <- c("transparent","grey40","grey63","transparent","grey80","grey30","grey45","grey90","black")
}

entro.expression <- expression(-(p1*((1/2)*fa1+fa2+(1/2)*fb1+fb0)*log(p1*((1/2)*fa1+fa2+(1/2)*fb1+fb0))+
                                (1-p1)*((1/2)*fa1+fa2+(1/2)*fb1+fb0)*log((1-p1)*((1/2)*fa1+fa2+(1/2)*fb1+fb0))+
                                p2*((1/2)*fa1+fa0+(1/2)*fb1+fb2)*log(p2*((1/2)*fa1+fa0+(1/2)*fb1+fb2))+
                                (1-p2)*((1/2)*fa1+fa0+(1/2)*fb1+fb2)*log((1-p2)*((1/2)*fa1+fa0+(1/2)*fb1+fb2))))


hatch.squares <- function(x,y,width,npoints.x=5,npoints.y=5,height,pch=20,cex=0.01)
{
    xy <- as.data.frame(cbind(x=x,y=y))

    for (iter in 1:nrow(xy))
    {
        points <- point.square(
                                    points.x=seq(xy[iter,"x"]-0.5*width,
                                                    xy[iter,"x"]+0.5*width,
                                                    width/npoints.x),
                                    points.y=seq(xy[iter,"y"]-0.5*height,
                                                    xy[iter,"y"]+0.5*height,
                                                    height/npoints.y))

        grid.points(x=points$x,
                        y=points$y,
                        gp=gpar(cex=cex),
                        pch=pch)
    }
}

calc.entropy <- function(row)
{
    if (NA %in% row)
    {
        return(0)
    } else
    {
        entro <- with(as.list(row), { eval(entro.expression) })

        if (is.na(entro))
        {
            entro <- 0
        }
        return(entro)
    }
}

choices <- expand.grid(a=c(1,2,3),
                        b=c(1,2,3))

choices[,"out"] <- seq(1,nrow(choices))

if (!exists("mumu"))
{
    mumu <- NULL
}

draw.p1p2.contour <- function(data.sub, k)
{
    data.order <- data.sub[order(data.sub$x2,data.sub$y2),]

    x2.un <- sort(unique(data.order$x2))
    y2.un <- sort(unique(data.order$y2))

    mat.dat <- matrix(data=(data.order$p1 - data.order$p2 + 0.01), 
                    nrow=length(x2.un),
                    ncol=length(y2.un),
                    byrow=T)

    contour <- contourLines(x=x2.un,
                            y=y2.un,
                            mat.dat,
                            levels=0)

    if (length(contour) < 1)
    {
        print(paste("no contour for k=",data.order$k[[1]]))
        return()
    }

    print(length(contour))

    for (i in 1:length(contour))
    {
        ct <-contour[[i]]
        panel.xyplot(x=unlist(ct["x"]),
                        y=unlist(ct["y"]),
                        type="l",
                        col="black",
                        lty="23",
                        lwd=2.5)
        panel.xyplot(x=unlist(ct["x"]),
                        y=unlist(ct["y"]),
                        type="l",
                        col="white",
                        lty="23",
                        lwd=1.5)

        points <- as.data.frame(cbind(x=unlist(ct["x"]),
                                    y=unlist(ct["y"])))

                    #        if (k==0.1)
                    #        {
                    #            points <- points[points$x > 0.42 & points$x < 0.58 & points$y < -0.21,]
                    #        } else if (k==0.5)
                    #        {
                    #            points <- points[points$x > 0.36 & points$x < 0.64 & points$y < -0.2,]
                    #        } else if (k==0.95)
                    #        {
                    #            points <- points[points$x > 0.3 & points$x < 0.7 & points$y < -0.5,]
                    #        }
                    #
                    #
                    #        str(points)
                    #
                    #        panel.xyplot(x=points$x,
                    #                        y=points$y,
                    #                        type="l",
                    #                        col="white",
                    #                        lty="24",
                    #                        lwd=2)
    }
}

plot.overview <- function(file="overview_phen", data)
{
    the.data <- read.table(data, sep=";", header=T)
    
    the.data <- NULL
    if (length(data) > 1)
    {
        for (i in 1:length(data))
        {
            data.sub <- read.table(data[[i]], sep=";", header=T)

            if (!("file" %in% names(data.sub)))
            {
                data.sub[,"file"] <- NA
            }

            the.data <- rbind(the.data,data.sub)
        }
    } else
    {
        the.data <- read.table(data, sep=";", header=T)
    }

    #    the.data <- the.data[the.data$time >= 450000,


    # calculate heights and widths
    widths <- c(0.3,rep(c(1,0.1),times=length(unique(the.data$k))))

    widths[[length(widths)]] <- 0.9

    heights <- c(0.4,0.9,0.2,1,0.2,1,0.2)

    init.plot(file=file, type="pdf", width=850, height=758,font="helvetica")

        layout <- grid.layout(widths=widths,
                                heights=heights,
                                ncol=length(widths),
                                nrow=length(heights))

        # the first levelplot
        pushViewport(viewport(layout=layout))

        lev.panel <- function(x,y,z,at,data.sub,...)
        {
            panel.levelplot(x,y,z,at=at,...)
            #            panel.abline(v=0.5,lwd=1.5,lty="48",col="white")
            draw.p1p2.contour(data.sub, k= data.sub$k[[1]])
        }

        if (is.null(mumu))
        {
            print("mumu null")
            the.data[,"intermed"] <- apply(the.data, 1, intermed.val)

            print(unique(the.data[,"intermed"]))


            #            the.data[,"x2"] <- round(the.data[,"s1"] / (the.data[,"s1"] + the.data[,"s2"]),digits=2)
            #            the.data[,"y2"] <- round((log10(the.data[,"s1"]) + log10(the.data[,"s2"]) )/2,digits=2)

            #            the.data <- the.data[the.data$x2 >= 0.2 & the.data$x2 <= 0.8,]
            #            the.data <- the.data[the.data$y2 != -0.38 & the.data$y2 != -0.58,]
        
            mumu <<- the.data
        }
        else
        {
            the.data <- mumu
        }

        at.range=seq(range(the.data[,"intermed"])[[1]]-1,range(the.data[,"intermed"])[[2]]+1,1)

        by(the.data, list(the.data$k), function(x)
        {
            lvp <- levelplot(intermed ~ x2 * y2,
                                    data=x,
                                    panel=lev.panel,
                                    data.sub=x,
                                    at=at.range,
                                    #                                    col.regions=blue2green2red(length(at.range)))
                                    # pure, pure, p1 mixed, pure, pure,
                                    col.regions=main.colors)#blue2green2red(3))

#            zlim.p <- c(0,1)
#            cld1 <- wireframe(p1 ~ x2 * y2,
#                                    data=x,
#                                    zoom=1.1,
#                                    zlim=zlim.p,
#                                    cex=0.5,
#                                    lwd=0.3,
#                                    col="grey30",
#                                    screen=list(z=-35,x=-60,y=-10),
#                                    drape=F,
#                                    xlab="",#list(label=expression(paste(italic(f)[italic(e)[2]])),cex=labelcex,rot=30),
#                                    ylab="",#list(label=expression(paste(bar(italic(s)))),cex=labelcex,rot=-40),
#                                    zlab="",#list(label=expression(paste(italic(p)[1],"      ",sep="")),cex=labelcex),
#                                    scales=list(x=list(cex=tickcex,at=pretty(range(x$x2),n=6)),
#                                            y=list(cex=tickcex,at=pretty(range(x$y2),n=6)),
#                                    z=list(cex=tickcex,at=pretty(range(zlim.p),n=6),labels=paste(as.character(pretty(range(zlim.p),n=6)),"  ",sep=""))
#                                    ),
#                            default.scales=list(arrows=F))
#
#            cld2 <- wireframe(p2 ~ x2 * y2,
#                                    data=x,
#                                    zoom=1.1,
#                                    zlim=zlim.p,
#                                    cex=0.1,
#                                    lwd=0.3,
#                                    col="grey30",
#                                    drape=F,
#                                    screen=list(z=-35,x=-60,y=-10),
#                                    #                                    col.regions=gray((10:85)/100),
#                                    xlab="",
#                                    ylab="",
#                                    zlab="",
#                                    scales=list(x=list(cex=tickcex,at=pretty(range(x$x2),n=6)),
#                                            y=list(cex=tickcex,at=pretty(range(x$y2),n=6)),
#                                    z=list(cex=tickcex,at=pretty(range(zlim.p),n=6),
#                                            labels=paste(as.character(pretty(range(zlim.p),n=6)),"  ",sep=""))
#                                            ),
#                                    default.scales=list(arrows=F))
#
                            #            cld.entro <- wireframe(entropy ~ x2 * y2,
                            #                                    data=x,
                            #                                    zoom=1.1,
                            #                                    cex=0.5,
                            #                                    lwd=0.3,
                            #                                    col="black",
                            #                                    drape=F,
                            #                                    zlim=c(-0,1.5),
                            #                                    screen=list(z=-35,x=-60,y=-10),
                            #                                    #                                    col.regions=gray((10:85)/100),
                            #                                    xlab=list(label=expression(paste(italic(f)[italic(e)[2]])),cex=labelcex,rot=30),
                            #                                    ylab=list(label=expression(paste(bar(italic(s)))),cex=labelcex,rot=-40),
                            #                                    zlab=list(label=expression(paste(italic(H),"  ",sep="")),cex=labelcex,rot=0),
                            #                                    scales=list(x=list(cex=tickcex,at=pretty(range(x$x2),n=6)),
                            #                                            y=list(cex=tickcex,at=pretty(range(x$y2),n=6)),
                            #                                            z=list(cex=tickcex,at=pretty(c(0,1.5),n=6),labels=paste(as.character(pretty(c(0,1.5),n=6)),"  ",sep=""))
                            #                                            ),
                            #                                    default.scales=list(arrows=F))
                            #
            colno <- match(x$k[[1]],sort(unique(the.data$k)))*2

            labelno <- match(x$k[[1]],sort(unique(the.data$k)))

            labels <- c("A","B","C")
            labels.2 <- c("D","E","F")
            labels.3 <- c("G","H","I")
            labels.4 <- c("J","K","L")
            
            pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=colno,
                            xscale=lvp$x.limits,
                            yscale=lvp$y.limits
                            ))

                do.call("lev.panel",trellis.panelArgs(lvp,1))

                grid.rect(gp=gpar(fill="transparent", lwd=rect.lwd, col="black"))

                mixture.coords.y <- c(0.6, 0.5, 0.3)
                z1.coords.y <- c(0.55, 0.55, 0.7)
                z1.coords.x <- c(0.2, 0.15, 0.15)

                z2.coords.x <- c(0.85, 0.85, 0.85)
                z2.coords.y <- c(0.55, 0.55, 0.7)
                
                bactp1.coords.x <- c(0.45, 0.42, 0.42)
                bactp1.coords.y <- c(0.7, 0.7, 0.7)
                
                bactp2.coords.x <- c(0.55, 0.57, 0.57)
                bactp2.coords.y <- c(0.7, 0.7, 0.7)
                
                opp.coords.x <- c(0.5, 0.5, 0.5)
                opp.coords.y <- c(0.9, 0.9, 0.9)

                if (labelno %in% c(2,3))
                {
                    # text for the region where nongenetic inheritance occurs
                    grid.text(x=0.5,just="centre",y=mixture.coords.y[[labelno]]-0.13,
                            label=expression(paste(bold("low-fidelity\ninheritance"))),
                            gp=gpar(cex=annotcex,fontface="bold",col="white",lineheight=annotcex))
                }
                
                if (labelno == 1)
                {
                    # text for the region where nongenetic inheritance occurs
                    grid.text(x=0.8,just="centre",y=0.65,
                            label=expression(paste(bold("low-fidelity\ninheritance"))),
                            gp=gpar(cex=annotcex,fontface="bold",col="black",lineheight=annotcex))

                    grid.lines(x=c(0.8, 0.55),
                                y=c(0.62,0.58),
                                gp=gpar(lwd=1),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))

                }


                if (labelno %in% c(1,2,3))
                {
                    #                    grid.text(x=0.5,y=mixture.coords.y[[labelno]]-0.23,
                    #                            just=c("centre","centre"),
                    #                            label=expression(paste(bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("1")],"+",bolditalic(z)[bold("2")])),
                    #                            gp=gpar(cex=annotcex.2,fontface="bold",col="white"))
                    #                    grid.text(x=0.5,y=mixture.coords.y[[labelno]]-0.28,
                    #                            just=c("centre","centre"),
                    #                            label=expression(paste(bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("1")],"+",bolditalic(z)[bold("2")])),
                    #                            gp=gpar(cex=annotcex.2,fontface="bold",col="white"))

                    grid.text(x=z1.coords.x[[labelno]],just=c("centre","centre"),y=z1.coords.y[[labelno]]-0.07,
                            label=expression(paste(bold(z)[bold("1")],bold(" only"),sep="")),
                            gp=gpar(cex=annotcex.4,fontface="bold",col="black",lineheight=annotcex))

                    #                    grid.text(x=z1.coords.x[[labelno]],y=z1.coords.y[[labelno]]-0.14,just=c("left","centre"),label=expression(paste(bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("1")])),gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                    #                    grid.text(x=z1.coords.x[[labelno]],y=z1.coords.y[[labelno]]-0.20,just=c("left","centre"),label=expression(paste(bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("1")])),gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                    #                   
                    grid.text(x=z2.coords.x[[labelno]],just=c("centre","centre"),y=z2.coords.y[[labelno]]-0.07,
                            label=expression(paste(bold(z)[bold("2")],bold(" only"),sep="")),
                            gp=gpar(cex=annotcex.4,fontface="bold",col="black",lineheight=annotcex))
                    #                    # text for the region where both produce z2
                    #                    grid.text(x=z2.coords.x[[labelno]],y=z2.coords.y[[labelno]]-0.07,just=c("centre","centre"),label=expression(paste(bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("2")])),gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                    #                    grid.text(x=z2.coords.x[[labelno]],y=z2.coords.y[[labelno]]-0.14,just=c("centre","centre"),label=expression(paste(bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("2")])),gp=gpar(cex=annotcex.2,fontface="bold",col="black"))

                    grid.text(x=0.82,just="centre",y=0.85,
                            label=expression(paste(bold("one heritable\nphenotype"))),
                            gp=gpar(cex=annotcex*0.95,fontface="bold",col="black",lineheight=annotcex))
                    
                    grid.lines(x=c(0.8, 0.55),
                                y=c(0.82,0.74),
                                gp=gpar(lwd=1),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))
                    
                    grid.lines(x=c(0.8, 0.45),
                                y=c(0.82,0.74),
                                gp=gpar(lwd=1),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))
                }

                if (labelno %in% c(1,2))
                {

                    grid.text(x=0.5,just=c("centre","centre"),y=c(0.16,0.05)[[labelno]],
                            label=expression(paste(bold("high-fidelity\ninheritance"),sep="")),
                            gp=gpar(cex=annotcex,fontface="bold",col="white",lineheight=annotcex))
                }

                if (labelno == 1000)
                {

                   
                    # text for the region where bacterial persistence occurs 
                    grid.text(x=0.4,
                                y=1.11,
                                just=c("right","centre"),
                                label=expression(
                                        paste(
                                                bold("4: "),bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("1")],"+",bolditalic(z)[bold("2")],sep=""
                                                    )),
                                gp=gpar(cex=annotcex.2,fontface="bold",col="black")
                                )

                    grid.text(x=0.4,
                            y=1.05,
                            just=c("right","centre"),
                            label=expression(paste(bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("1")])),
                            gp=gpar(cex=annotcex.2,fontface="bold",col="black"))

                    grid.lines(x=c(0.42, 0.44, 0.44),
                                y=c(1.08, 1.08, 0.85),
                                gp=gpar(lwd=1.5),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))
                    
                    # text for the region where bacterial persistence occurs 
                    grid.text(x=0.6,
                                y=1.11,
                                just=c("left","centre"),
                                label=expression(
                                        paste(
                                                bold("5: "),bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("2")])),
                                gp=gpar(cex=annotcex.2,fontface="bold",col="black")
                                )

                    grid.text(x=0.6,
                            y=1.05,
                            just=c("left","centre"),
                            label=expression(paste(phantom("5: "),bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("1")],"+",bolditalic(z)[bold("2")])),
                            gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                    
                    grid.lines(x=c(0.58, 0.56, 0.56),
                                y=c(1.08, 1.08, 0.85),
                                gp=gpar(lwd=1.5),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))
                    
                    # text for the region where both play opposite
                    grid.text(x=0.5,
                                y=1.23,
                                just=c("centre","centre"),
                                label=expression(
                                        paste(
                                                bold("6: "),bolditalic(z)[bold("1")]%->%phantom(z),bolditalic(z)[bold("2")])),
                                gp=gpar(cex=annotcex.2,fontface="bold",col="black")
                                )

                    grid.text(x=0.5,
                            y=1.18,
                            just=c("centre","centre"),
                            label=expression(paste(phantom("6: "),bolditalic(z)[bold("2")]%->%phantom(z),bolditalic(z)[bold("1")])),
                            gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                    
                    grid.lines(x=c(0.5, 0.5),
                                y=c(1.13, 0.95),
                                gp=gpar(lwd=1.5),
                                arrow=arrow(type="closed",length=unit(0.01,"npc")))
                }
                grid.text(x=0.5,y=mixture.coords.y[[labelno]],just=c("centre","centre"),label="2",gp=gpar(cex=annotcex.2,fontface="bold",col="white"))
                grid.text(x=z1.coords.x[[labelno]],y=z1.coords.y[[labelno]],just=c("centre","centre"),label="1",gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                grid.text(x=z2.coords.x[[labelno]],y=z2.coords.y[[labelno]],just=c("centre","centre"),label="3",gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                
                grid.text(x=bactp1.coords.x[[labelno]],
                            y=bactp1.coords.y[[labelno]],
                            just=c("centre","centre"),
                            label="4",gp=gpar(cex=annotcex.2,fontface="bold",col="white"))
                
                grid.text(x=bactp2.coords.x[[labelno]],y=bactp2.coords.y[[labelno]],just=c("centre","centre"),label="5",gp=gpar(cex=annotcex.2,fontface="bold",col="white"))
                grid.text(x=opp.coords.x[[labelno]],y=opp.coords.y[[labelno]],just=c("centre","centre"),label="6",gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                
                x.coor <- c(0.32, 0.29, 0.23) 
                y.coor <- c(0.1, 0.1, 0.1) 
                grid.text(x=x.coor[[labelno]],y=y.coor[[labelno]],just=c("centre","centre"),label="8",gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                
                x.coor <- c(0.68, 0.71, 0.77) 
                y.coor <- c(0.1, 0.1, 0.1) 
                grid.text(x=x.coor[[labelno]],y=y.coor[[labelno]],just=c("centre","centre"),label="9",gp=gpar(cex=annotcex.2,fontface="bold",col="black"))
                
                if (labelno < 3)
                {
                    co7 <- c(0.3,0.17)
                
                    grid.text(x=0.5,
                                y=co7[[labelno]],
                                just=c("centre","centre"),
                                label="7",gp=gpar(cex=annotcex.2,fontface="bold",col="white"))
                }

                grid.circle(x=circle.x.offset,
                                y=circle.y.offset,
                                r=circle.rad,
                                gp=gpar(lwd=0.5))
                grid.text(x=circle.x.offset,
                                y=circle.y.offset,
                                label=labels[[labelno]])

            upViewport(1)
            pushViewport(viewport(layout.pos.row=1,
                            layout.pos.col=colno,
                            xscale=lvp$x.limits,
                            yscale=lvp$y.limits
                            ))

                disp.text.y <- 0.2

                if (labelno == 2)
                {
                    disp.text.y <- 0.2
                }            

                approxvals <- c(0.83,0.51,0.37)

                grid.text(x=0.5,y=disp.text.y,
                                label=eval(substitute(expression(paste(italic(k)," = ",k.val,", dispersal, ",italic(d)%~~% approx.val,sep="")),list(k.val=x$k[[1]],approx.val=approxvals[[labelno]]))),
                                gp=gpar(cex=maincex))

            upViewport(1)


            pushViewport(viewport(layout.pos.row=3,
                            layout.pos.col=colno,
                            xscale=lvp$x.limits,
                            yscale=lvp$y.limits
                            ))

                single.axis(lvp$x.limits,side="top",nmain=8,nsub=4,cex=tickcex, labels=T, text=expression(paste("frequency environment 2, ",italic(f)[italic(e)[2]],sep="")),tck=0.5, y.text.off=0.3,labelcex=labelcex)
            upViewport(1)
            
            pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=colno-1,
                            xscale=lvp$x.limits,
                            yscale=lvp$y.limits
                            ))

                labels <- colno==2
                text <- ifelse(labels,expression(paste("mean log transition rate, ",bar(italic(s)),sep="")),"")

                single.axis(lvp$y.limits,side="right",nmain=8,nsub=4,cex=tickcex, labels=labels, text=text,tck=0.5,x.text.off=0.3,labelcex=labelcex)
            upViewport(1)

            pushViewport(viewport(layout.pos.row=4,
                                layout.pos.col=colno,
                                xscale=cld1$x.limits,
                                yscale=cld1$y.limits))

                do.call("panel.wireframe",trellis.panelArgs(cld1,1))

            if (x$k[[1]] == 0.1)
            {
                grid.text(x=-0.06,y=0.25,rot=100,just=c("left","bottom"), label=expression(paste("proportion ",italic(z)[1]," offspring",sep="")), gp=gpar(cex=legendcex))
                grid.text(x=-0.015,y=0.35,rot=100,just=c("left","bottom"), label=expression(paste("by ",italic(z)[1]," parent, ",italic(p)[1],sep="")), gp=gpar(cex=legendcex))
            }

            grid.text(x=0,y=0.15,rot=-31,just=c("left","bottom"), label=expression(paste("frequency environment 2, ",italic(f)[italic(e)[2]])), gp=gpar(cex=legendcex))
            grid.text(x=0.61,y=-0.17,rot=38,just=c("left","bottom"), label=expression(paste("mean log transition rate, ",bar(italic(s)))), gp=gpar(cex=legendcex))

            grid.circle(x=circle.x.offset,
                            y=circle.y.offset+0.05,
                            r=circle.rad,
                            gp=gpar(lwd=0.5))
            grid.text(x=circle.x.offset,
                            y=circle.y.offset+0.05,
                            label=labels.2[[labelno]])
            upViewport(1)

            
            pushViewport(viewport(layout.pos.row=6,
                                layout.pos.col=colno,
                                xscale=cld2$x.limits,
                                yscale=cld2$y.limits))
            
                do.call("panel.wireframe",trellis.panelArgs(cld2,1))
            
            if (x$k[[1]] == 0.1)
            {
                grid.text(x=-0.06,y=0.25,rot=100,just=c("left","bottom"), label=expression(paste("proportion ",italic(z)[1]," offspring",sep="")), gp=gpar(cex=legendcex))
                grid.text(x=-0.015,y=0.35,rot=100,just=c("left","bottom"), label=expression(paste("by ",italic(z)[2]," parent, ",italic(p)[2],sep="")), gp=gpar(cex=legendcex))
            }

            grid.text(x=0,y=0.15,rot=-31,just=c("left","bottom"), label=expression(paste("frequency environment 2, ",italic(f)[italic(e)[2]])), gp=gpar(cex=legendcex))
            grid.text(x=0.61,y=-0.17,rot=38,just=c("left","bottom"), label=expression(paste("mean log transition rate, ",bar(italic(s)))), gp=gpar(cex=legendcex))
            
            grid.circle(x=circle.x.offset,
                            y=circle.y.offset+0.05,
                            r=circle.rad,
                            gp=gpar(lwd=0.5))
            grid.text(x=circle.x.offset,
                            y=circle.y.offset+0.05,
                            label=labels.3[[labelno]])
            upViewport(1)
            
            #            pushViewport(viewport(layout.pos.row=8,
            #                                layout.pos.col=colno,
            #                                xscale=cld.entro$x.limits,
            #                                yscale=cld.entro$y.limits))
            #            
            #                do.call("panel.wireframe",trellis.panelArgs(cld.entro,1))
            #            
            #            grid.circle(x=circle.x.offset,
            #                            y=circle.y.offset+0.05,
            #                            r=circle.rad,
            #                            gp=gpar(lwd=0.5))
            #            grid.text(x=circle.x.offset,
            #                            y=circle.y.offset+0.05,
            #                            label=labels.4[[labelno]])
            #            upViewport(1)
        })
            
        pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=length(widths),
                            ))
            x.pos.text <- 0.09
            grid.text(x=x.pos.text,y=0.95,just=c("left","centre"),label="1",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.95,just=c("left","center"),expression(paste("both parents produce only ",italic(z)[1],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.9,just=c("left","center"),expression(paste(italic(p)[1]," = ",italic(p)[2]," = 1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=0.8,just=c("left","centre"),label="2",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.8,just=c("left","centre"),expression(paste("both parents produce mixture",sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.75,just=c("left","centre"),expression(paste("of ",italic(z)[1]," and ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.7,just=c("left","centre"),expression(paste("0<",italic(p)[1],", ",italic(p)[2],"<1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=0.6,just=c("left","centre"),label="3",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.6,just=c("left","center"),expression(paste("both parents produce only ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.55,just=c("left","centre"),expression(paste(italic(p)[1]," = ",italic(p)[2]," = 0",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=0.45,just=c("left","centre"),label="4",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.45,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce mixture",sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.4,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce only ",italic(z)[1],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.35,just=c("left","centre"),expression(paste("0<",italic(p)[1],"<1, ",italic(p)[2]," = 1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=0.25,just=c("left","centre"),label="5",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.2,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce only ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.25,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce mixture",sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.15,just=c("left","centre"),expression(paste(italic(p)[1]," = 0, 0<",italic(p)[2],"<1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=0.05,just=c("left","centre"),label="6",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=0.05,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce only ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=0.0,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce only ",italic(z)[1],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.05,just=c("left","centre"),expression(paste(italic(p)[1]," = 0, ",italic(p)[2]," = 1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=-0.15,just=c("left","centre"),label="7",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=-0.15,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce only ",italic(z)[1],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.20,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce only ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.25,just=c("left","centre"),expression(paste(italic(p)[1]," = 1, ",italic(p)[2]," = 0",sep="")),gp=gpar(cex=annotcex))

            grid.text(x=x.pos.text,y=-0.35,just=c("left","centre"),label="8",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=-0.35,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce only ",italic(z)[1],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.40,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce mixture",sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.45,just=c("left","centre"),expression(paste(italic(p)[1]," = 1, 0<",italic(p)[2],"<1",sep="")),gp=gpar(cex=annotcex))
            
            grid.text(x=x.pos.text,y=-0.55,just=c("left","centre"),label="9",gp=gpar(cex=annotcex.2,fontface="bold"))

            grid.text(x=0.15,y=-0.55,just=c("left","centre"),expression(paste(italic(z)[1]," parents produce mixture",sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.60,just=c("left","centre"),expression(paste(italic(z)[2]," parents produce only ",italic(z)[2],sep="")),gp=gpar(cex=annotcex))
            grid.text(x=0.15,y=-0.65,just=c("left","centre"),expression(paste("0<",italic(p)[1],"<1, ",italic(p)[2],"=0",sep="")),gp=gpar(cex=annotcex))

            grid.lines(x=c(0.1,0.2),y=c(-0.75,-0.75),gp=gpar(lwd=5,col="black",lty="24"))
            grid.lines(x=c(0.1,0.2),y=c(-0.75,-0.75),gp=gpar(lwd=1,col="white",lty="24"))
            grid.text(x=0.22,y=-0.75,just=c("left","centre"),label="negative parental effects",gp=gpar(cex=annotcex,lineheight=annotcex))
            grid.text(x=0.22,y=-0.81,just=c("left","centre"),label=expression(paste(italic(p)[2],">",italic(p)[1],sep="")),gp=gpar(cex=annotcex,lineheight=annotcex))
        upViewport(1)


        upViewport(1)

    exit.plot()
}

#plot.overview(data="../results/summary_pPhen_m_total_stopruns.csv")
#plot.overview(data=c("xpphen_fulltime_plus_sims"))
#plot.overview(data=c("summary_pphen_new_setup_e_0001.csv"))
#plot.overview(data=c("summary_pphen_full_range.csv"))
#plot.overview(data=c("summary_pphenerror.dat"))
#plot.overview(data=c("tah"))
plot.overview(data=c("mortsel_p1p2_summary.csv"))#c("summary_pphen_seems_okay_k01.csv","summary_pphen_k05_now_really_ok.csv","summary_pphen_now_okay_k095.csv"))
#plot.overview(data="summary_xpphen_full_completed.csv")
