plotps <-
function(data, plot = 1, ps, a, basepollution = NULL, outcome = NULL){
ps <- data[,ps] 
a  <- data[,a]
bp <- data[,basepollution]
out<- data[,outcome]

if(all.equal(unique(a), c(0,1)) == F) stop("a must be coded as values 0, 1.")


if(plot == 1) {
plot(ps, col=a+1, pch=16, ylab="Propensity Score")
}

if(plot == 2){
par(mfrow=c(1,2))
hist(ps[a==1], main="Non-attainment", xlab="PS")
hist(ps[a==0], main="Attainment", xlab="PS")
}
if(plot == 3) {
plot(ps[a==1], bp[a==1], xlim=c(0,1), ylim=range(bp, na.rm=T), pch=16, 
xlab="Propensity Score", ylab="Baseline pollution", col=2)
points(ps[a==0], bp[a==0], pch=16)
legend("topleft", col=c("black", "red"), c("Attainment", "Non-attainment"), pch=16)
}
if(plot == 4){
  plot(ps[a==1 & is.na(out) == F ], 
                   bp[a==1 & is.na(out) == F ], 
xlim=c(0,1), ylim=range(bp, na.rm=T), pch=1, 
xlab="Propensity Score", ylab="Baseline Pollution", col=2, main="Missing excluded")
points(ps[a==0 & is.na(out) == F ], 
 bp[a==0 & is.na(out) == F ], pch=1)
points(ps[a==0 & is.na(out) == T ], 
 bp[a==0 & is.na(out) == T ], pch=16,  col=1)
points(ps[a==1 & is.na(out) == T ], 
 bp[a==1 & is.na(out) == T ], pch=16,  col=2)
legend("topleft", col=c("black", "red", "black", "red"), c("Attainment", "Non-attainment", 
"Missing Attain", "Missing NA"), pch=c(1, 1, 16, 16))
}
}
