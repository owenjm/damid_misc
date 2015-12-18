#!/usr/bin/env Rscript 

X11()

plotr2 <- function(x,y,title="",pchar=20, varx=deparse(substitute(x)), vary=deparse(substitute(y)), cormethod='pearson') {
    lsum <- summary(lm(y ~ x))
    r2 <- lsum$r.squared
    intercept <- lsum$coefficients[1,1]
    xfact <- lsum$coefficients[2,1]
    
    equation <- cat(paste(vary,"=",xfact,varx,"+",intercept))
    
    corr <- cor(x,y,use="complete",method=cormethod)
    
    plot(x,y,pch=pchar,col=rgb(0,0,0,0.1),xlab=varx,ylab=vary)
    
    rp = vector('expression',3)
    rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
		  list(MYVALUE = format(r2,digits=3)))[2]
    rp[2] = substitute(expression(cor == MYOTHERVALUE), 
		  list(MYOTHERVALUE = format(corr, digits = 3)))[2]
    rp[3] = equation
    
    legend("topleft",legend=rp,bty='n')
    title(title)
}

read.gff <- function (x,name="score") {
  temp.data <- read.table(x,row.names=NULL)[,c(1,4,5,6)]
  names(temp.data) <- c("chr","start","end",name)
  return(temp.data)
}

args <- commandArgs(trailingOnly = TRUE)

file.1 <- args[1]
file.2 <- args[2]

name.1 <- gsub(".gff","",file.1)
name.1 <- gsub(".*/","",file.1)

name.2 <- gsub(".gff","",file.2)
name.2 <- gsub(".*/","",file.2)

cat(paste("Reading",file.1,"...\n"))
n1 <- read.gff(file.1,"n1")
n1 <- as.data.frame(lapply(n1,function(x) if(is.character(x)|is.factor(x)) gsub("chr","",x) else x))

cat(paste("Reading",file.2,"...\n"))
n2 <- read.gff(file.2,"n2")
n2 <- as.data.frame(lapply(n2,function(x) if(is.character(x)|is.factor(x)) gsub("chr","",x) else x))

cat("Merging data ...\n")
exp.f <- merge(n1,n2,by=c("chr","start","end"))

cat("Plotting ...\n")
plotr2(exp.f$n1, exp.f$n2, title="Correlation plot", varx=paste(name.1,"log2(DAM-fusion/DAM)"), vary=paste(name.2,"log2(DAM-fusion/DAM)"), pchar=".")

cat("Writing files ...\n")
dev.copy(png, paste(name.1,"-",name.2,"correlation.plot.png",sep="."), width=800,height=800);dev.off()
dev.copy(pdf,paste(name.1,"-",name.2,"correlation.plot.pdf",sep="."),width=10,height=10);dev.off()

cat("All done.\n\n")