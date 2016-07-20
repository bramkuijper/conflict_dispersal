# figure that prints the course of evolution
# when only m or b, or both evolve

library("lattice")
library("grid")
source("/home/bram/R/src/bramlib.r")

library("rootSolve")

heights <- c(0.2, 1, 0.3,1,0.3)
widths <- c(0.3,1,0.3,1,0.3)

line.lwd <- 0.75

lwd <- 0.5
tick.cex <- 0.5
label.cex <- 0.8
legend.cex <- 0.8
plot.tck <- -0.25
ind.label.cex <- 0.8 

Power <- function(x,y)
{
    return(x^y)
}

dmat1 <- expression(
-(((1 - d1)*(p1 + p2))/
       (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)*
         (-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
            (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
           ((-1 + n)*Power(-1 + d2,2)*p2)/
            (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2))))) + 
    ((1 - (p1 + p2)/
          (n*(-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
               (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
              ((-1 + n)*Power(-1 + d2,2)*p2)/
               (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2)))))*
       (1/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2) + 
         (1 - c1)*(-(p1/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2)) - 
            p2/(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2)))))/2.)

dmat2 <- expression(
-(((1 - d2)*(p1 + p2))/
       (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2)*
         (-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
            (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
           ((-1 + n)*Power(-1 + d2,2)*p2)/
            (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2))))) + 
    ((((-1 + c2)*p1)/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2) + 
         (1 + (-1 + c2)*p2)/(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2)))*
       (1 - (p1 + p2)/
          (n*(-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
               (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
              ((-1 + n)*Power(-1 + d2,2)*p2)/
               (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2))))))/2.)

doff1 <- expression((1/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2) - 
    ((1 - d1)*(p1 + p2))/
     (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)*
       (-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
          (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
         ((-1 + n)*Power(-1 + d2,2)*p2)/
          (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2)))) + 
    (1 - c1)*(-(p1/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2)) - 
       p2/(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2)))))

doff2 <- expression(1/(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2)) - 
    ((1 - d2)*(p1 + p2))/
     (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2)*
       (-1 + ((-1 + n)*Power(-1 + d1,2)*p1)/
          (n*Power(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2,2)) + 
         ((-1 + n)*Power(-1 + d2,2)*p2)/
          (n*Power(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2),2)))) + 
    (1 - c2)*(-(p1/(-1 + d1*(1 + (-1 + c1)*p1) + (-1 + c2)*d2*p2)) - 
       p2/(-1 + (-1 + c1)*d1*p1 + d2*(1 + (-1 + c2)*p2))))

if (!exists("dat"))
{
    dat <- as.data.frame(expand.grid(p=.5,n=4,d1=NA,d2=NA,d1.off=NA,d2.off=NA,c1=c(0.001,seq(0.01,0.99,0.01),0.999),c2=c(0.001,seq(0.01,0.99,0.01),0.999)))

    # function to provide to multiRoot
    mat.sol <- function(x,pars) {
        n <- as.numeric(pars["n"])
        c1 <- as.numeric(pars["c1"])
        c2 <- as.numeric(pars["c2"])
        p1 <- as.numeric(pars["p"])
        p2 <- 1 - p1

        d1 <- x[1]
        d2 <- x[2]

        F1 <- eval(dmat1)
        F2 <- eval(dmat2)

        c(F1,F2)
    }
    off.sol <- function(x,pars) {
        n <- as.numeric(pars["n"])
        c1 <- as.numeric(pars["c1"])
        c2 <- as.numeric(pars["c2"])
        p1 <- as.numeric(pars["p"])
        p2 <- 1 - p1

        d1 <- x[1]
        d2 <- x[2]

        F1 <- eval(doff1)
        F2 <- eval(doff2)

        c(F1,F2)
    }

    # loop through the dataset
    for (i in 1:nrow(dat))
    {
        pars <- list(n=dat[i,"n"],
                        p=dat[i,"p"],
                        c1=dat[i,"c1"],
                        c2=dat[i,"c2"])
                        
        ss <- multiroot(f=mat.sol, start=c(0.1,0.25),pars=pars)

        dd <- ss$root
        dd[[1]] <- ifelse(dd[[1]] < 0.001,0.001,ifelse(dd[[1]] > 0.999,0.999,dd[[1]]))
        dd[[2]] <- ifelse(dd[[2]] < 0.001,0.001,ifelse(dd[[2]] > 0.999,0.999,dd[[2]]))

        ss <- multiroot(f=off.sol, start=c(0.1,0.25),pars=pars)
        dd.off <- ss$root
        dd.off[[1]] <- ifelse(dd.off[[1]] < 0.001,0.001,ifelse(dd.off[[1]] > 0.999,0.999,dd.off[[1]]))
        dd.off[[2]] <- ifelse(dd.off[[2]] < 0.001,0.001,ifelse(dd.off[[2]] > 0.999,0.999,dd.off[[2]]))

        # add to data frame
        dat[i,c("d1","d2","d1.off","d2.off")] <- c(dd,dd.off)
    }
}



lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)



block <- function(row, col, dataset, yvar,
        colors,alpha,label="A",label.main,ylab,xlab,ylim,
        x.labels=T,
        abline=NA,
        lwds,ltys)
{
    xp <- wireframe(as.formula(paste(paste(yvar,collapse="+")," ~ c1 * c2")),
                    yvar=yvar,
                    aspect=c(1,0.5),
                    data=dataset,
                    zlim=c(0,1),
                    xlim=c(0,1),
                    ylim=c(0,1),
                    cex=0.5,
                    zoom=1.2,
                    col="grey80",
                    alpha=0.1,
                    lwd=0.5,
                    drape=T,
                    screen=list(x=-50,y=-30),
                    par.settings=list(box.3d=list(lwd=0.3,col="transparent")),
                    scales=list(x=list(cex=tick.cex,at=pretty(range(dataset$c1),n=6)),
                            y=list(cex=tick.cex,at=pretty(range(dataset$c2),n=6)),
                            z=list(cex=tick.cex,at=pretty(range(dataset$d1),n=6))
                            ),
                    col.regions=colors,
                    alpha.regions=alpha,
                    default.scales=list(arrows=F,lwd=0.5),
                    the.dat=dataset)


    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.wireframe",trellis.panelArgs(xp,1))

        grid.text(x=0.5,y=1.1,label=label.main,gp=gpar(cex=label.cex))

        grid.text(x=0.07,y=0.95,label=label,gp=gpar(cex=ind.label.cex))
    upViewport()
    
    pushViewport(viewport(layout.pos.row=row+1,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        x.at <- c()
        if (x.labels)
        {
            x.labels=c("1","10","100","1000","10000","100000")
            x.at=c(0,1,2,3,4,5)
        }
    upViewport()
    
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col-1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            expr <- ylab
    upViewport()
}

 
init.plot(filename="fig_battleground",width=650,height=320,font="times")
pushViewport(viewport(layout=lo))
k <- trellis.par.get("box.3d")

k$lwd <- 0.5

trellis.par.set("box.3d",k)

str(dat)

    block(row=2,col=2,
            ylim=c(-0.02,1.0),
            dataset=dat,
            yvar=c("d1"),
            colors=c(rep("white",times=9),rep("grey70",times=120)),
            alpha=rep(0.5,times=120),
            label=expression(italic(a)),
            lwds=list(list(1.0,1.0,0.5,0.5)),
            ltys=list(list(1,1,"12","12")),
            ylab=expression(paste("log mean fitness, log(",bar(italic(W))[italic(t)],")",sep="")),
            label.main="",
            x.labels=F,
            xlab="")

    block(row=2,col=4,
            ylim=c(-1,21),
            dataset=dat,
            alpha=rep(0.5,times=120),
            yvar=c("d1.off"),
            colors=c(rep("white",times=9),rep("grey10",times=120)),
            label=expression(italic(b)),
            lwds=list(list(1.0,1.0,0.5,0.5)),
            ltys=list(list(1,1,"12","12")),
            label.main="",
            x.labels=F,
            ylab=expression(paste("mean phenotype, ",bar(italic(z))[italic(t)],sep="")),
            xlab="")

    block(row=4,col=2,
            ylim=c(-1,21),
            dataset=dat,
            alpha=rep(0.5,times=120),
            yvar=c("d1","d1.off"),
            colors=c(rep("white",times=9),rep("grey10",times=120)),
            label=expression(italic(b)),
            lwds=list(list(1.0,1.0,0.5,0.5)),
            ltys=list(list(1,1,"12","12")),
            label.main="",
            x.labels=F,
            ylab=expression(paste("mean phenotype, ",bar(italic(z))[italic(t)],sep="")),
            xlab="")
    
legendje <- function(x,y,width=0.05,height=0.05,label,
                    lty=1,
                    lwd=1.0,fill.color="blue")
{
    pushViewport(viewport(x=x,y=y,width=width,height=height))
#        grid.rect(x=unit(0.2,"npc"),width=unit(0.4,"npc"),height=unit(0.5,"npc"),gp=gpar(fill=fill.color,alpha=0.5,lwd=0.5,col="transparent"))
        grid.lines(x=unit(c(0,.4),"npc"),y=unit(c(0.5,0.5),"npc"),gp=gpar(col=fill.color,lty=lty,lwd=lwd))
       # grid.points(x=unit(.2,"npc"),y=unit(0.5,"npc"),gp=gpar(col=fill.color,pch=21,cex=0.1,fill=fill.color))
        grid.text(x=unit(0.5,"npc"),y=0.5,just="left",label=label,gp=gpar(cex=legend.cex))
    upViewport()
}

exit.plot()


