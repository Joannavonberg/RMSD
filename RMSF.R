#install.packages("rPython")
#install.packages(repos=NULL, pkgs="/var/autofs/pauken/home/berg/R/i486-pc-linux-gnu-library/2.15/rPython")
library(rPython)

options(stringsasFactors = FALSE)

# unit cell dimensions
#a <- 78.64	# NB!!! these are cryo-dimensions!
#b <- 78.64
#c <- 37.06

a <- 79.3	# NB!!! these are RT-dimensions!
b <- 79.3 
c <- 38.2 

rmsd <- c()

# reference values are needed, these are the positions of the crystal structure

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

rmsf <- function(row, ucp){
     total <- 0
     for(e in row[1:8]){
     	   total <- total + min( ucp - max(row["ref"], e) + min(row["ref"], e) , abs(row["ref"] - e) )^2
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

load <- function(a, n){
     tmp <- scan(sprintf("%s%.0f.txt", a, n))
     x <- data.frame(matrix(tmp, ncol = 8))
     colnames(x) <- c(LETTERS[1:8])
     x
}

python.load('/work/berg/scripts/old/changePDB.py')

#	CRYO
#refx <- scan("/work/berg/Git/ref/x_cryo_protein_noH.txt")
#refy <- scan("/work/berg/Git/ref/y_cryo_protein_noH.txt")
#refz <- scan("/work/berg/Git/ref/z_cryo_protein_noH.txt")

#	RT
refx <- scan("/work/berg/Git/ref/x_RT_protein_noH.txt")
refy <- scan("/work/berg/Git/ref/y_RT_protein_noH.txt")
refz <- scan("/work/berg/Git/ref/z_RT_protein_noH.txt")

x_trans <- c()
y_trans <- c()
z_trans <- c()

for(n in 1:400){
      x <- load("x", n)
      y <- load("y", n)
      z <- load("z", n)

      x2 <- x
      y2 <- y
      z2 <- z

      #x <- x%%a
      #y <- y%%b
      #z <- z%%c

      #x2$A <- x$A%%a
      #y2$A <- y$A%%b
      #z2$A <- z$A%%c

      x2$B <- (-x$B)			# -X+1
      y2$B <- (-y$B)	  		# -Y+1
      z2$B <- (z$B - 0.5*c)	  	# Z+1/2

      x2$C <- (y$C - 0.5*b)		# -Y+1/2
      y2$C <- (-x$C + 0.5*a)   		# X+1/2
      z2$C <- (z$C + 0.25*c)		# Z-1/4

      x2$D <- (-y$D + 0.5*a)		# Y+1/2
      y2$D <- (x$D - 0.5*b) 		 # -X+1/2
      z2$D <- (z$D - 0.25*c)		# Z+1/4

      x2$E <- (-x$E + 0.5*a)	  	# -X+1/2
      y2$E <- (y$E - 0.5*b)		# Y+1/2
      z2$E <- (-z$E + 0.75*c)		# -Z+3/4

      x2$F <- (x$F - 0.5*a)		# X+1/2
      y2$F <- (-y$F + 0.5*b)		# -Y+1/2
      z2$F <- (-z$F + 0.25*c)	  	# -Z+5/4

      x2$G <- (y$G)			# Y
      y2$G <- x$G			# X
      z2$G <- (-z$G)		# -Z+1

      x2$H <- (-y$H)		# -Y+1
      y2$H <- (-x$H)		# -X+1
      z2$H <- (-z$H + 0.5*c)		# -Z+1/2

      # to put the first (or 505th) atom in the unit cell, and keep molecules whole
      if(n == 1){
      	   for(t in 1:8){
      	   	 x_premodulo <- x2[505,t]
	   	 y_premodulo <- y2[505,t]
	   	 z_premodulo <- z2[505,t]

      	   	 x_trans <- c(x_trans, x_premodulo%%a - x_premodulo)
	   	 y_trans <- c(y_trans, y_premodulo%%b - y_premodulo)
	   	 z_trans <- c(z_trans, z_premodulo%%c - z_premodulo)
	   }
	   refx <- apply(x2, 1, mean)
	   refy <- apply(y2, 1, mean)
	   refz <- apply(z2, 1, mean)
      }

      # to make sure the same translation is used in every timestep, it is only calculated for the first timestep
      for(t in 1:8){
      	    x2[,t] <- x2[,t] + x_trans[t]
	    y2[,t] <- y2[,t] + y_trans[t]
	    z2[,t] <- z2[,t] + z_trans[t]
      }

      rmsd <- c(rmsd, sqrt(RMSD(x2, ref = refx, ucp = a)^2 + RMSD(y2, ref = refy, ucp = b)^2 + RMSD(z2, ref = refz, ucp = c)^2) )

      # to write to a pdb file

      #x3 <- c()
      #for (t in 1:8){
      #	  x3 <- c(x3, x2[,t])
      #}

      #y3 <- c()
      #for (t in 1:8){
      #	  y3 <- c(y3, y2[,t])
      #}

      #z3 <- c()
      #for (t in 1:8){
      #	  z3 <- c(z3, z2[,t])
      #}

      #python.call('change', x3, y3, z3, n, FALSE)
}



png("rmsd_step1asref.png")

plot(rmsd, type = "l", axes=NULL, xlab = "timestep (ns)", ylab = "RMSD (Angstrom)")#, ylim = c(1, 1.6))
axis(1, at=seq(0,400, 100), lab = "timestep (ns)")

dev.off()
