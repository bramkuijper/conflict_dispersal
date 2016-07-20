library("lattice")

if(!exists("the.data"))
{
    the.data <- read.table("total_output.csv",sep=";",header=T)
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

pdf("overview.pdf",width=20,height=20)
print(wireframe(diff.s1 ~ c1 * c2 | d * sigma12 * sigma21,
                strip=function(strip.levels,...) { strip.default(strip.levels=T,...) },
                lwd=0.3,
                col="grey90",
                data=the.data.rs
                )
        )
dev.off()
