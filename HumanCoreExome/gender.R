        args=commandArgs()[7]

        # Extracting key user specified parameter from 'arg's
        datafile=substr(args, 2, nchar(args))

        args2=commandArgs()[8]
        cutoff=as.integer(substr(args2, 2, nchar(args2)))

        print(datafile)
        print(cutoff)


data.raw <- read.table(file = datafile, header=T)

postscript("gender.ps", paper="letter", horizontal=T)
N <- sum((data.raw$SEX==1 & data.raw$N_ySNP < cutoff) | (data.raw$SEX==2 & data.raw$N_ySNP>cutoff) )
text <- paste(sep="", "Gender Checking in JDRF DN Samples (", N, " Samples Mislabeled)")
plot(data.raw$N_ySNP[data.raw$SEX==2], data.raw$xHeterozygosity[data.raw$SEX==2], type="p",
            col = "red", xlim=c(0,max(data.raw$N_ySNP)), ylim=c(0,max(data.raw$xHeterozygosity)),
            main = text,
            xlab="# Y-Chr SNPs", ylab = "X-Chr Heterozygosity")
points(data.raw$N_ySNP[data.raw$SEX==1], data.raw$xHeterozygosity[data.raw$SEX==1], col = "blue")
points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>cutoff], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>cutoff], col = "red")
points(data.raw$N_ySNP[data.raw$SEX==0], data.raw$xHeterozygosity[data.raw$SEX==0], col = "black")
abline(v=cutoff,col="purple", lty=2)
abline(v=cutoff*2/3, col="purple")
abline(v=cutoff*4/3, col="purple")
abline(h=0.1, col="purple")

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
legend=c("Female", "Male", "Unknown"), col=c("red", "blue", "black"),
text.col = c("red", "blue", "black"), pch = 19, cex = 1.2)

dev.off()



