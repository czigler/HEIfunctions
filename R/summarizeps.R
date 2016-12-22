summarizeps <-
function(ps, a, data, nextreme){

a <- data[,a]
ps<- data[,ps]

cat("A = 0 \n")
cat("Min/Max: (", round(min(ps[a==0]), digits=4), ", ", round(max(ps[a==0]), digits=4), ") \n\n", sep="")

cat("A = 1 \n")
cat("Min/Max: (", round(min(ps[a==1]), digits=4), ", ", round(max(ps[a==1]), digits=4), ") \n\n", sep="")

cat("Top ", nextreme, " smallest", "\n")
table <- round(cbind(sort(ps[a==0])[1:nextreme], sort(ps[a==1])[1:nextreme]), digits=4)
colnames(table) <- c("A = 0", "A = 1")
row.names(table) <- 1:nextreme
print(table)
cat("\n")

cat("Top ", nextreme, " largest", "\n")
table <- round(cbind(sort(ps[a==0], decreasing=T)[1:nextreme], sort(ps[a==1], decreasing=T)[1:nextreme]), digits=4)
colnames(table) <- c("A = 0", "A = 1")
row.names(table) <- 1:nextreme
print(table)

}
