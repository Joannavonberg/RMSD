options(stringsasFactors = FALSE)

# unit cell dimensions
a <- 79.3
b <- 79.3
c <- 38.2

rmsd <- c()

RMSD <- function(timeframe){
     tot <- 0
     chainA <- timeframe[,1]
     tot <- sum(apply(timeframe, 2, MeanDifference, chainA = chainA))
     tot <- sqrt(tot/length(timeframe[1,]))
     tot 
}

MeanDifference <- function(col, chainA){
     sum(abs(col - chainA)^2)/length(col)
}



for(n in 1:100){
      tmp <- scan(sprintf("x%.0f.txt",n))
      x <- data.frame(matrix(tmp, ncol = 8))
      colnames(x) <- c(LETTERS[1:8])		
      tmp <- scan(sprintf("y%.0f.txt",n))
      y <- data.frame(matrix(tmp, ncol = 8))
      colnames(y) <- c(LETTERS[1:8])
      tmp <- scan(sprintf("z%.0f.txt",n))
      z <- data.frame(matrix(tmp, ncol = 8))
      colnames(z) <- c(LETTERS[1:8]) 

      x2 <- x
      y2 <- y
      z2 <- z

      #x <- x%%a
      #y <- y%%b
      #z <- z%%c

      x2$A <- x$A%%a
      y2$A <- y$A%%b
      z2$A <- z$A%%c

      x2$B <- (-x$B + a)%%a		# -X+1
      y2$B <- (-y$B + b)%%b	  	# -Y+1
      z2$B <- (z$B - 0.5*c)%%c	  	# Z+1/2

      x2$C <- (y$C - 0.5*b)%%a		# -Y+1/2
      y2$C <- (-x$C + 0.5*a)%%b   		# X+1/2
      z2$C <- (z$C + 0.25*c)%%c		# Z-1/4

      x2$D <- (-y$D + 0.5*a)%%a		# Y+1/2
      y2$D <- (x$D - 0.5*b)%%b 		 # -X+1/2
      z2$D <- (z$D - 0.25*c)%%c		# Z+1/4

      x2$E <- (-x$E + 0.5*a)%%a	  	# -X+1/2
      y2$E <- (y$E - 0.5*b)%%b		# Y+1/2
      z2$E <- (-z$E + 0.75*c)%%c		# -Z+3/4

      x2$F <- (x$F - 0.5*a)%%a		# X+1/2
      y2$F <- (-y$F + 0.5*b)%%b		# -Y+1/2
      z2$F <- (-z$F + 1.25*c)%%c	  	# -Z+5/4

      x2$G <- (y$G)%%b			# Y
      y2$G <- x$G%%a			# X
      z2$G <- (-z$G + c)%%c		# -Z+1

      x2$H <- (-y$H + a)%%a		# -Y+1
      y2$H <- (-x$H + b)%%b		# -X+1
      z2$H <- (-z$H - 0.5*c)%%c		# -Z+1/2

      rmsd <- c(rmsd, mean(RMSD(x2), RMSD(y2), RMSD(z2)))
}

png("rmsd_300K_NVT_cryodim.png")
plot(rmsd)
dev.off()