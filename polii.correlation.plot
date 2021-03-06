#!/usr/bin/env Rscript

# Copyright © 2014-15, Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

X11()

greptidy <- function (x) {
  y <- x[grep("^(?!CR|CG|snoR|tRNA|scaR|snRNA|snmRNA|5srRNA|His-Psi|His|5S|mir-)",x,value=F,perl=T)]
  return(y)
}

plotr2.iii <- function(x,y,gene.names,title="",pchar=20, varx=deparse(substitute(x)), vary=deparse(substitute(y)),highlight="",plot.names=F,highlight.names="",grep.tidy=T,fontsize=1) {
    r2 <- summary(lm(y ~ x), na.action="na.exclude")$r.squared
    corr <- cor(x,y,use="complete",method="spearman")
    
    plot(x,y,pch=pchar,col=rgb(0,0,0,0.1),xlab=varx,ylab=vary)
    
   
    datf <- data.frame(name=gene.names,x,y)
    hlgenes <- highlight
	
    points(datf$x[datf$name %in% hlgenes], datf$y[datf$name %in% hlgenes], col=rgb(1,0,0,0.2), pch=20 )
    
    if (grep.tidy) {
      hlgenes.txt <- greptidy(highlight.names)
    } else {
    	hlgenes.txt <- highlight.names
    }
    
    rp = vector('expression',2)
    rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
		list(MYVALUE = format(r2,digits=3)))[2]
    rp[2] = substitute(expression(cor == MYOTHERVALUE), 
		list(MYOTHERVALUE = format(corr, digits = 3)))[2]
    
    legend("topleft",legend=rp,bty='n')
    title(title)
    
    if (plot.names == T) {
      outl <- datf[ datf$name %in% hlgenes.txt,]
      par(cex=fontsize)
      
	  if (length(outl$name) > 0) {
		text(outl[,2],outl[,3],labels=outl$name,adj=c(-0.2/fontsize,0.5),cex=0.7,col=rgb(0,0,0,0.7))
	  }
    }
  
    abline(h=0,v=0)
}

### Read CLI options
input.args <- commandArgs(trailingOnly = TRUE)

in.files <- vector()
read.ops <- function (x) {
  for (op in x) {
	if (any(grepl("^--",op))) {
		op <- gsub("^--","",op)
		y <- unlist(strsplit(op,"="))
		
		if (y[1] == "help") {
		  cat(paste("polii.correlation.plot\nGenerates expression correlation plot from .details.csv files generated by polii.gene.call.\n\nUsage: polii.correlation.plot [options] --in1=[polii.n1.details.csv] --in2=[polii.n2.details.csv] --name1=[short_name_of_file_1] --name2=[short_name_of_file_2]\n\n", sep=""))
		  
		  cat("Options:\n")
		  for (n in names(op.args)) {
			cat(paste("  ",n,"=",op.args[[n]],"\n",sep=""))
		  }
		  cat("\n")
		  quit("no",1)
		}
		
		if (!is.null(op.args[[ y[1] ]])) {
		  op.args[[ y[1] ]] <<- y[2]
		} else {
		  cat("Error: Option",y[1],"not recognised ...\n")
		  quit("no",1)
		}
	} else {
		in.files <<- c(in.files,op)
	}
  }
}

write.ops <- function () {
  out.df <- data.frame()
  for (n in names(op.args)) {
	v <<- as.character(op.args[[n]])
	df.line <- data.frame(
	  option=n,
	  value=v
	)
	out.df <- rbind(out.df, df.line)
  }
  write.table(out.df,"input.args.single.txt",row.names=F)
}

op.args <- list(
  "in1" = "",
  "in2" = "",
  "name1" = "file1",
  "name2" = "file2"
)

read.ops(input.args)

if (length(op.args[["in1"]]) == 0 || length(op.args[["in2"]]) == 0) {
	cat("polii.correlation.plot\nGenerates expression correlation plot from .details.csv files generated by polii.gene.call.\n\nUsage: polii.correlation.plot [options] --in1=[polii.n1.details.csv] --in2=[polii.n2.details.csv]\n\n")
	quit("no",1)
}

write.ops()

n1 <- read.table(op.args[["in1"]],header=T)
names(n1) <- c("name","n1","gatcs","FDR.n1")
n2 <- read.table(op.args[["in2"]],header=T)
names(n2) <- c("name","n2","gatcs","FDR.n2")


curr.path <- paste("comparison-",op.args[["name1"]],"-",op.args[["name2"]],sep="")
dir.create(curr.path)
setwd(curr.path)


exp.f <- merge(n1,n2,by=c("name","gatcs"))
exp.f[exp.f == Inf] <- NA
exp.f[exp.f == -Inf] <- NA
exp.f <- na.omit(exp.f)

norm.factor <- vector();
for (gene in c("Gapdh1","RpL32","Act5C","alphaTub84B")) {
	norm.factor <- c(norm.factor, exp.f$n1[exp.f$name == gene] / exp.f$n2[exp.f$name == gene])
}
n <- mean(norm.factor)

cat(paste("Norm factor was",n,"\n"))

exp.f.norm <- exp.f
exp.f.norm$n1 <- exp.f.norm$n1 / n


outliers <- greptidy(subset(exp.f, (FDR.n1 < 0.01 | FDR.n2 < 0.01) & (abs(n1 - n2) > 1), name)[[1]])
outliers.norm <- greptidy(subset(exp.f.norm, (FDR.n1 < 0.01 | FDR.n2 < 0.01) & (abs(n1 - n2) > 1), name)[[1]])


## un-norm
plotr2.iii(exp.f$n1, exp.f$n2, exp.f$name, title="Average gene expression", varx=paste(op.args[["name1"]],"log2(DAM-fusion/DAM)"), vary=paste(op.args[["name2"]],"log2(DAM-PolII/DAM)"), plot.names=T, highlight=subset(exp.f, (FDR.n1 < 0.01 | FDR.n2 < 0.01), name)[[1]], highlight.names=outliers,fontsize = 0.7, grep.tidy=F)
mtext("FDR < 0.01 highlighted")

dev.copy(png, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.png",sep=""), width=800,height=800); dev.off()
dev.copy(pdf, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.pdf",sep=""), width=10, height=10 ); dev.off()
dev.copy(pdf, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.hires.pdf",sep=""), width=25, height=25 ); dev.off()

## norm
plotr2.iii(exp.f.norm$n1, exp.f.norm$n2, exp.f.norm$name, title="Normalised average gene expression", varx=paste(op.args[["name1"]],"log2(DAM-fusion/DAM)"), vary=paste(op.args[["name2"]],"log2(DAM-PolII/DAM)"), plot.names=T, highlight=subset(exp.f.norm, (FDR.n1 < 0.01 | FDR.n2 < 0.01), name)[[1]], highlight.names=outliers.norm,fontsize = 0.7, grep.tidy=F)
mtext("FDR < 0.01 highlighted")

dev.copy(png, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.norm.png",sep=""), width=800,height=800); dev.off()
dev.copy(pdf, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.norm.pdf",sep=""), width=10, height=10 ); dev.off()
dev.copy(pdf, paste(op.args[["name1"]],"-vs-",op.args[["name2"]],".correlation.plot.norm.hires.pdf",sep=""), width=25, height=25 ); dev.off()


## write tables
data.out <- function (exp.f, norm=F) {
	exp.f$diff <- exp.f$n1 - exp.f$n2
	
	normtext <- ifelse (norm, ".norm", "")
	
	exp.out <- exp.f
	
	exp.out <- subset(exp.out, (FDR.n1 < 0.01 | FDR.n2 < 0.01))
	
	names(exp.out)[3] <- op.args[["name1"]]
	names(exp.out)[4] <- paste("FDR",op.args[["name1"]],sep=".")
	names(exp.out)[5] <- op.args[["name2"]]
	names(exp.out)[6] <- paste("FDR",op.args[["name2"]],sep=".")
	
	write.table(exp.out,paste(op.args[["name1"]],"-",op.args[["name2"]],normtext,".expression.csv",sep=""), row.names=F, col.names=T, quote=T, sep=",")
}

data.out(exp.f)
data.out(exp.f.norm,norm=T)

setwd("..")

