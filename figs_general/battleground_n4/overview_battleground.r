library("lattice")
library("colorRamps")

if(!exists("the.data"))
{
    the.data <- read.table("total_output_n4.csv",sep=";",header=T)
}

# reshape the data so that we get maternal and offspring strategies 
# next to one another and we can substract
the.data.rs <- reshape(
        the.data,
        direction="wide",
        timevar="type", # the value that varies within a given parameter combi: offspring or mother
        idvar=c("p1","d","n","c1","c2","sigma12","sigma21") # the parameters that should not vary within a subject
        )

the.data.rs["diff.s1"] <- the.data.rs$s1.mother - the.data.rs$s1.offspring
the.data.rs["diff.s2"] <- the.data.rs$s2.mother - the.data.rs$s2.offspring

pdf("overview_s1.pdf",width=20,height=20)
print(levelplot(diff.s1 ~ c1 * c2 | d * sigma12 * sigma21
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,col.regions=matlab.like
                ,xlab=expression(paste("Cost of maladaptation in envt 1")),
                ,ylab=expression(paste("Cost of maladaptation in envt 2")),
                ,main=expression(paste("Difference between prop z1 offspring produced by mother and offspring in envt 1")),
                ,col="grey90"
                ,data=the.data.rs
                )
        )
dev.off()

pdf("overview_s2.pdf",width=20,height=20)
print(levelplot(diff.s2 ~ c1 * c2 | d * sigma12 * sigma21
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,col.regions=matlab.like
                ,xlab=expression(paste("Cost of maladaptation in envt 1")),
                ,ylab=expression(paste("Cost of maladaptation in envt 2")),
                ,main=expression(paste("difference between prop z1 offspring produced by mother and offspring in envt 2")),
                ,col="grey90"
                ,data=the.data.rs
                )
        )
dev.off()
