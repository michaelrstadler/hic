# This is a file of old R code that I didn't end up using but has some interesting stuff that may be useful


# Grabs the diagonal of a HiC file and rotates it so the diagonal is horizontal, also just takes upper half
HiC.diagonalize <- function(x, width){
	x.horiz <- diag(x)
	x.width <- length(x[1,])
	x.height <- length(x[,1])
	for (i in 1:width){
		new.row <- c(rep(0,i), x[row(x) == col(x) + i])
		x.horiz <- rbind(new.row, x.horiz)
	}
	return(x.horiz)
}

#Zeros out everything but some given number of diagonal rows. So a width of 3 would leave the diagonal and the +1 and +2 diagonals, make everything else zero
HiC.zeroOffDiag <- function(x, width){
	x.rows <- length(x[,1])
	x.cols <- length(x[1,])
	x1 <- x
	cols.to.zero.above <- width + 2 # leave two diagonals above the diagonal
	cols.to.zero.below <- 1 # leave two diagonals below the diagonal
	for (i in 1:x.rows){
		if (cols.to.zero.above <= x.cols){
			x1[i,cols.to.zero.above:x.cols] <- 0
		}
		cols.to.zero.above <- cols.to.zero.above + 1
		
		if (i > width + 1){
			x1[i,1:cols.to.zero.below] <- 0
			cols.to.zero.below <- cols.to.zero.below + 1
		}		
	}
	return(x1)
}

#Zeros the diagonal and additional (+1, +2, etc. diagonal) diagonals as specified by width; leaves everything else.
HiC.zeroDiag <- function(x, width){
	x.rows <- length(x[,1])
	x.cols <- length(x[1,])
	x1 <- x
	for (i in 1:x.rows){
		start <- max(1, i - width)
		stop <- min(x.cols, i + width)
		x1[i, start:stop] <- 0
	}	
	return(x1)
}

# Makes little blocks of 'doped' values along the diagonal. Makes boxes of size n by n (n defined by size variable) that have all counts increased by multiplying by supplied factor. Boxes alternate, with an altered box followed by an unaltered box.
diagonal.dope <- function(x, size, factor){
	x1 <- x
	for (i in seq(1, nrow(x), 2*size)){
		for (j in 0:(size - 1)){
			for (k in 0:(size - 1)){
				x1[i + j, i + k] <- factor * x[i + j, i + k]
			}
		}
	}
	return(x1)
}

# Takes two files of Hi-C binned matrices, produces a simple difference map for the "cool region" at left end of 3R. I also set it up to do autonaming for hte output file with the date string at the beginning. Not sure why.
cool.region.difference.heatmap <- function(filename1, filename2, name1, name2){
	process <- function(filename){
		x <- read.matrix(filename)
		x <- grab.chr(x, '3R')
		x <- HiC.matrix.vanilla.normalize(x)
		x <- HiC.matrix.scale.log(x)
		x <- x[55:390,55:390]
		return(x)
	}
	plot.subtraction <- function(y1, y2, y1.name, y2.name){
		date.string <- gsub('-','',Sys.Date())
		y3 <- y1 - y2
		#y3 <- histogram.equalize(y3)
		jpeg(paste(date.string, '_', y1.name, '_minus_',y2.name,'.jpeg',sep=''),4000,4000)
		heatmap.natural(y3)
		dev.off()
	}
	x1 <- process(filename1)
	x2 <- process(filename2)
	plot.subtraction(x1, x2, name1, name2)	
	plot.subtraction(x2, x1, name2, name1)	
}

#numbers are for 1300 bp bins
polytene.plot.heat <- function(filename, outname, topcompress=0.995, bottomcompress=0.75, startbin=565, endbin=645, LOG=TRUE){
	x <- read.matrix(filename)
	x <- x[startbin:endbin,startbin:endbin]
	x[is.na(x)] <- 0
	x <- HiC.matrix.compress(x,topcompress, bottomcompress)
	if (LOG){
		x <- HiC.matrix.scale.log(x)
	}
	jpeg(outname, 4000,4000)
	heatmap.natural(x)
	dev.off()
}


ROC.draw <- function(model, dep.variable){
	predict <- predict(model, type='response')
	ROCRpred <- prediction(predict, dep.variable)
	ROCRperf <- performance(ROCRpred, 'tpr','fpr')
	plot(ROCRperf, colorize = FALSE, text.adj = c(-0.2,1.7))
}