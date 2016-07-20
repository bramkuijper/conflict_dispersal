
library("lattice")
library("grid")
library("colorRamps")
source("/home/bram/R/src/bramlib.r")

# p=0.4, d=0.7
data <- read.table("iter_conflict_22_3_2016_162539_228799606",sep=";",header=T)

#data <- read.table("iter_conflict_22_3_2016_162539_228800502",sep=";",header=T)


heights <- c(0.2, 1, 0,1,0,1,0.2)
widths <- c(0.5,1,0.1,1,0.05)


line.lwd <- 0.75

lwd <- 0.5
tick.cex <- 0.5
label.cex <- 0.8
sub.label.cex <- 0.67
legend.cex <- 0.8
plot.tck <- -0.25
ind.label.cex <- 1.0


lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)


block <- function(row, col, 
                    formula,
                    zlim=c(0,1),
                    ind.label="A",
                    screen=list(z=90,x=-60,y=-10),
                    data)
{
    xp <- wireframe(as.formula(formula),
                        data=data,
                        lwd=0.1,
                        scales=list(x=list(cex=tick.cex),
                                    y=list(cex=tick.cex),
                                    z=list(cex=tick.cex)),
                        ylim=c(0,1),
                        xlim=c(0,1),
                        zlim=zlim,
                        zoom=1.1,
                        screen=screen,
                        aspect=c(1,0.7),
                        xlab=list(cex=label.cex,label=expression(paste(italic(c)[1]))),
                        ylab=list(cex=label.cex,label=expression(paste(italic(c)[2]))),
                        default.scales=list(arrows=F)
                        )

    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.wireframe",trellis.panelArgs(xp,1))

        grid.text(x=0.1,y=0.9,label=ind.label)

    upViewport(1)
    

}


init.plot(filename="all_wireframes",width=650,height=800,font="helvetica",type="pdf")
pushViewport(viewport(layout=lo))

data[data$qS == data$qNS,"siginfo"] <- 0

block(2,2, formula="s1 ~ c1 * c2",data=data)
block(2,4,ind.label="B", formula="s2 ~ c1 * c2",data=data)
block(4,2,ind.label="C", formula="qS ~ c1 * c2",data=data)
block(4,4,ind.label="D", formula="qNS ~ c1 * c2",data=data)
block(6,2,ind.label="E", formula="siginfo ~ c1 * c2",data=data,zlim=c(0,1.3), screen=list(z=60,x=-60,y=-10))

upViewport(1)

exit.plot()
