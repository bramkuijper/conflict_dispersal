#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

mom <- args[[1]]
off <- args[[2]]

library("lattice")



data.offspring <- read.table(off,sep=";",header=T,nrow=1000)
data.mother <- read.table(mom,sep=";",header=T,nrow=1000)


pdf(paste("plot_conflict",basename(mom),".pdf"))
print(xyplot(data.mother$s1 + data.mother$s2 + data.offspring$s1 + data.offspring$s2 + data.offspring$p1 ~ data.mother$generation,
                ,type="l"
                ,ylim=c(0,1),
                ,lwd=2
                ,auto.key=T))
dev.off()

#print(xyplot(data.mother$vars1 + data.mother$vars2 + data.offspring$vars1 + data.offspring$vars2  ~ data.mother$generation,
#                ,type="l"
#                ,lwd=2
#                ,auto.key=T))
