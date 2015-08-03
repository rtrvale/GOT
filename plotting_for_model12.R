source("got_readin.R") # gets char.names
source("inference_test_model12.R")

#### Plots of posterior distributions ####
##########################################

plot_bar <- function(fit, char){

# make a single barplot for one character from the fitted model "fit".

pred1 <- fit$pred1
N <- nrow(pred1)
par(cex.lab=3, cex.axis=1.9)
plot(0:20, 0:20, ylim=c(0,1), type="n", xaxt="n", yaxt="n",
xlab=char.names[char], ylab="")

for (i in 0:20) rect(i-0.3, 0, i+0.3, 
sum(pred1[,char]==i)/N, col="black")
axis(1, at=seq(0,20,4), paste(seq(0,20,4)), padj=1 )
axis(1, at=seq(2,18,4), paste(seq(2,18,4)), padj=1 )
axis(2, at=c(0, 0.5, 1), paste(c(0, 0.5, 1)), las=1)

par(cex.lab=1, cex.axis=1)
}

plot_24 <- function(fit){

# make two pages of bar plots from a model fitted to a 24x5 data matrix
# Save the results as an A4 png file for reproducibility

opar <- par()

png(file="posteriors_1.png", width=500, height=707, units="mm", res=72)
par(mfrow=c(4,3), mar=c(6.1,3.1,4.1,2.1), mgp = c(5, 1, 0))
for (i in 1:12) plot_bar(fit, i)
dev.off()

dev.new()

png(file="posteriors_2.png", width=500, height=707, units="mm", res=72)
par(mfrow=c(4,3), mar=c(6.1,3.1,4.1,2.1), mgp = c(5, 1, 0))
for (i in 13:24) plot_bar(fit, i)
dev.off()

par(opar)
}

#### Plots of posterior probabilities of zero POV chapters ####

plot_prob0 <- function(fit, error.bars=T, to.png=F){

  pred1 <- fit$pred1
  pred2 <- fit$pred2
  prob0 <- apply(pred1, 2, function(x) sum(x==0))/nrow(pred1)
  prob0.2 <- apply(pred2, 2, function(x) sum(x==0))/nrow(pred2)

  ord <- order(prob0)
  
  if (to.png) png(file="prob0.png", width=210, height=210, units="mm", res=72)
  plot(1:24, prob0[ord], pch=16, bg="black", xlab="", xaxt="n", ylim=c(0,1), ylab="",
  las= 1)
  axis(1, las=2, at=1:24, char.names[ord], cex.axis=0.9)
  if (error.bars){ 
    for (i in 1:24) segments(i, prob0[ord][i]-0.03, i, 
                    min(prob0[ord][i] + 0.03,1))
  }
  for (i in 1:24) points(i, prob0.2[ord][i], col="blue") 

  legend("bottomright", pch=c(16,1), col=c("black", "blue"),
legend=c("Posterior probability 0 POV chapters in book 6",
"Posterior probability 0 POV chapters in book 7") )

  if (to.png) dev.off()
}

#### Plots of posterior intervals for a model fit to show coverage.
#### Used in the validation for M.1 = M[1:9, 1:2]

plot_validation <- function(fit, true, to.png=F){

# fit -- a model fit
# true -- actual values of M[,d+1] in a vector of length nrow(M)

pred1 <- fit$pred1
ord <- order(apply(pred1, 2, median))
lower50 <- colQ(pred1, 0.25)
upper50 <- colQ(pred1, 0.75)
lower80 <- colQ(pred1, 0.1)
upper80 <- colQ(pred1, 0.9)

top.y <- max(pred1)
N <- ncol(pred1)

if (to.png) png(file="validation.png", width=210, height=210, units="mm", res=72)
plot(1:N, lower50[ord], ylim=c(0, top.y), "l", xaxt="n", xlab="", 
ylab="POV Chapters", las=1)
lines(upper50[ord])
lines(lower80[ord], lty=2)
lines(upper80[ord], lty=2)

points(1:N, true[ord], pch=16, col=grey(0.1), cex=1.3)
axis(1, las=2, at=1:N, char.names[1:N][ord])

legend("topleft", lty=c(1,2), legend=c("50% credible bounds from fitted model",
"80% credible bounds from fitted model"))
if (to.png) dev.off()
}

### Function to take quantiles by column

colQ <- function(X, alpha) apply(X, 2, function(x) quantile(x, alpha))

### Function to write matrix for use in Latex file

put_matrix <- function(matrix, char.names){

for (i in 1:length(char.names)){
  cat(char.names[i])
  for (j in 0:22){
    cat("&", sum(matrix[,i]==j))
  }
  cat("\\\\","\n")
}

}

### Function to plot results of the inference test

plot_coverage_pars <- function(OUT){

  alpha <- seq(0.01, 0.99, 0.01)
  hits <- rep(0, length(alpha))
  for (i in 1:99) hits[i] <- count_pars_hits(OUT, alpha[i])
  plot(alpha, hits, xlab="Credibility level", ylab="Actual coverage", ylim=c(0,1),
  las=1)
  abline(0,1)
}

plot_coverage_pred1 <- function(OUT){
  
  alpha <- seq(0.01, 0.99, 0.01)
  hits <- rep(0, length(alpha))
  for (i in 1:99) hits[i] <- count_pred1_hits(OUT, alpha[i])
  plot(alpha, hits, xlab="Credibility level", ylab="Actual coverage", ylim=c(0,1),
  las=1)
  abline(0,1)  

}

plot_coverage_png <- function(OUT){

# Make the two coverage plots into a single png file.

png(file="inference_test.png", width=210, height=105, units="mm", res=72)

par(mfrow=c(1,2))
plot_coverage_pars(OUT)
plot_coverage_pred1(OUT)
par(mfrow=c(1,1))

dev.off()
}
