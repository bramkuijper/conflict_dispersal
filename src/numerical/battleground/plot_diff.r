#!/usr/bin/env Rscript
# plot difference between offspring and maternal strategies

args <- commandArgs(trailingOnly = T)

library("lattice")
library("colorRamps")

the.data <- read.table(args[[1]],sep=";",header=T)

# data is in long format (offspring and mother are on separate rows)
# whereas for substracting stuff it would be great to have it in wide format
# in which offspring values are next to the maternal ones
# see http://www.ats.ucla.edu/stat/r/faq/reshape.htm 

the.data.rs <- reshape(
        the.data,
        direction="wide",
        timevar="type", # the value that varies within a given parameter combi: offspring or mother
        idvar=c("p1","d","n","c1","c2","sigma12","sigma21") # the parameters that should not vary within a subject
        )

the.data.rs.2 <- the.data.rs

the.data.rs["diff"] <- the.data.rs$s1.mother - the.data.rs$s1.offspring
the.data.rs["var"] <- "s1"

the.data.rs.2["diff"] <- the.data.rs$s2.mother - the.data.rs$s2.offspring
the.data.rs.2["var"] <- "s2"

the.data.rs.diff <- rbind(the.data.rs.2,the.data.rs)


pdf(file="diff_s1s2.pdf")
print(wireframe(diff ~ c1 * c2 | var
                ,data=the.data.rs.diff
                ,default.scales=list(arrows=F)
                ,main=paste("d=",the.data.rs[1,"d"],", s12=",the.data.rs[1,"sigma12"],", s21=",the.data.rs[1,"sigma21"],sep="")
                #,col.regions=matlab.like
                )
        )
dev.off()

pdf(file="overlay.pdf")
print(wireframe(s1.mother + s1.offspring ~ c1 * c2
                ,data=the.data.rs
                ,alpha=0.8,
                ,main=paste("d=",the.data.rs[1,"d"],", s12=",the.data.rs[1,"sigma12"],", s21=",the.data.rs[1,"sigma21"],sep="")
                ,default.scales=list(arrows=F)
                #,col.regions=matlab.like
                )
        )
dev.off()
