library("lattice")
library("MASS")

dat <- read.table("summary_overview_p.csv",header=T,sep=";")

panel.cp <- function(x,y,z,the.data,...)
{
    c1vals <- sort(unique(the.data[,"c1"]))
    c2vals <- sort(unique(the.data[,"c2"]))

    #    the.data <- the.data[the.data$siginfo == 0.001,]


    matrix.contour <- matrix(data=the.data[,"siginfo"], 
                                nrow=length(c1vals), 
                                ncol=length(c2vals),
                                byrow=F)

    cnt <- contourLines(x=c1vals,
                        y=c2vals,
                        matrix.contour,
                        levels=.001)

    contour.mut <- as.data.frame(cnt[[1]])


    #    panel.polygon(x=c(contour.mut[,"x"],add.x.coords),
    #                    y=c(contour.mut[,"y"],add.y.coords),
    #                    border="black",
    #                    alpha=0.5,
    #                    col="black")

    print(length(contour.mut[,"x"]))
    #                npoints=15,
    #                type="l", 
    #                col="black", 
    #                lty=1,
    #                lwd=1.0)


    #first part of contour where x > y
    contour.mut.1 <- contour.mut[contour.mut$x <= contour.mut$y,]
    contour.mut.2 <- contour.mut[contour.mut$x >= contour.mut$y,]

    xy1 <- smooth.spline(x=contour.mut.1[,"x"],
                    y=contour.mut.1[,"y"]
                    )

    xy2 <- smooth.spline(x=contour.mut.2[,"x"],
                    y=contour.mut.2[,"y"]
                    )

    panel.polygon(x=c(xy1$x,rev(xy2$x)),
                    y=c(xy1$y,rev(xy2$y)),
                    type="l"
                    )
}

sub.dat <- dat[dat$p1 == 0.1,]

print(contourplot(siginfo ~ c1 * c2, 
                data=sub.dat,
                the.data=sub.dat,
                panel=panel.cp,
                at=c(0.001,0.98)))
