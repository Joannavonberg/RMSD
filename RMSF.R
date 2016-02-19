install.packages("rPython") 
library(rPython)

options(stringsasFactors = FALSE)

# unit cell dimensions
#a <- 78.64	# NB!!! these are cryo-dimensions!
#b <- 79.64
#c <- 37.06

a <- 79.3	# NB!!! these are RT-dimensions!
b <- 79.3 
c <- 38.2 

rmsd <- c()

# reference values are needed, these are the positions of chain A at timestep 0

RMSD <- function(timeframe, ref, ucp){
     tot <- 0
     tot <- sum(apply(timeframe, 2, MeanDifference, ref = ref, ucp = ucp))
     tot <- sqrt(tot/length(timeframe[1,]))
     tot 
}

MeanDifference <- function(col, ref, ucp){
     #mean(abs(col - ref)^2)
     total <- 0
     for(n in 1:length(col)){
     	 total <- total + min( (ucp - max(ref[n], col[n]) + min(ref[n], col[n])) , abs(ref[n] - col[n]) )^2
	 n <- n + 1
     }
     total <- total/length(col)
}

rmsf <- function(row, ucp, ref){
     total <- 0
     for(e in row){
     	   total <- total + min( ucp - max(ref, e) + min(ref, e) , abs(ref - e) )^2
     }
     total <- sqrt(total/8)
     total
}

dif <- function(row, ucp, ref){
    total <- 0
    n <- 1
    for (e in row){
    	total <- total + (abs(e - ref[n])%%ucp)^2
	n <- n + 1
    }
    total <- sqrt(total/8)
    total
}

x2$RMSF <- apply(x2, 1, rmsf, ucp = a, ref = refx)
y2$RMSF <- apply(y2, 1, rmsf, ucp = b, ref = refy)
z2$RMSF <- apply(z2, 1, rmsf, ucp = c, ref = refz)

png("frame45_z2_rmsf.png")
plot(z2$RMSF, pch = 19, cex = 0.5, xlab = "atoms", ylab = "RMSF")
lines(z2$RMSF)
dev.off()

x2$RMSF <- NULL
y2$RMSF <- NULL
z2$RMSF <- NULL

for(n in 1:101){
      tmp <- scan(sprintf("x%.0f.txt",n))
      x <- data.frame(matrix(tmp, ncol = 8))
      colnames(x) <- c(LETTERS[1:8])		
      tmp <- scan(sprintf("y%.0f.txt",n))
      y <- data.frame(matrix(tmp, ncol = 8))
      colnames(y) <- c(LETTERS[1:8])
      tmp <- scan(sprintf("z%.0f.txt",n))
      z <- data.frame(matrix(tmp, ncol = 8))
      colnames(z) <- c(LETTERS[1:8]) 

      if(n == 1){
      	   refx <- x$A
	   refy <- y$A
	   refz <- z$A
      }

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

      rmsd <- c(rmsd, mean(RMSD(x2, ref = refx, ucp = a), RMSD(y2, ref = refy, ucp = b), RMSD(z2, ref = refz, ucp = c)))
}

png("rmsd_300K_NPT.png")
plot(rmsd, pch = 19, cex = 0.5, xlab = "timesteps", ylab = "RMSD")
lines(rmsd)
dev.off()

# to write to a pdb file

x3 <- c()
for (n in 1:8){
    x3 <- c(x3, x2[,n])
}

y3 <- c()
for (n in 1:8){
    y3 <- c(y3, y2[,n])
}

z3 <- c()
for (n in 1:8){
    z3 <- c(z3, z2[,n])
}

python.load('../../scripts/changePDB.py')
python.call('change', x3, y3, z3)