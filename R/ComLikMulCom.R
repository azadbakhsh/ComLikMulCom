


ComLikMulCom <- function(data, yno , idno , y.type =c("normal", "probit", "quadratic exponential")[1] ,  type = c("many_one", "all_pair")[1], f=1){
  
  source("CLcov.R")
  
  require(mvtnorm)
  
  N   <-  length((data[,idno]))               #total number of measurements
  n   <-  length(unique(data[,idno]))         #number of individuals
  p   <-  ncol(data)-2                        # number of factors
  m   <-  N/n
  
  if("many_one" %in% type){
    
    mt1     <- matrix(0, nrow=(p-1), ncol=p)
    mt1[,f] <- -1
    
    for (i in 1:(p-1)){
      if(mt1[i,i]==0) {        
        mt1[i,i]<- 1}
      else if (mt1[i,i]==-1){
        k <- i
{break}}
    }

for (j in k:(p-1)){
  mt1[j,(j+1)]<- 1}

C     <-    mt1

  }


if("all_pair" %in% type){
  
  C     <-   matrix(0,nrow=p,ncol=p*(p-1)/2)
  comb   <- 	combn(p,2)
  C[cbind(comb[1,],1:(p*(p-1)/2))] <- 1
  C[cbind(comb[2,],1:(p*(p-1)/2))] <- -1
  C     <-    t(C)
}



covcl <- CLcov(data = data, yno = yno, idno = idno, N = N, n = n, p= p , y.type = y.type, type=type)



bstr      <- C%*%covcl$beta_hat
cov_str   <- C%*%covcl$Cov_beta%*%t(C)
se_str    <- sqrt(diag(cov_str))^(-1)
v_mat     <- diag(se_str)
corr_beta <- v_mat %*% cov_str %*% v_mat

T         <-	 bstr*se_str
cov_T 	  <- 	 v_mat%*%cov_str%*%t(v_mat)  


###	computing the m-variate normal quantile to consider it as threshold and performing hypothesis test	

qu <- qmvnorm(0.95, tail ="both.tails", mean = 0,  corr=cov_T)
quantile <- qu$quantile

rm <- abs(T)<quantile
sM     <-    ifelse(rm ==TRUE,"Accept","Reject") 

###Bonferroni

Bon <- 0.05/nrow(C)
qBon <- qnorm(1-Bon/2, lower.tail = TRUE)  

rb <- abs(T)<qBon
sB    <-  ifelse(rb==TRUE,"Accept","Reject")                      

###Dunn-Sidak

dunn <- 1-(1-0.05)^(1/nrow(C))
qDunn <- qnorm(1-dunn/2)  

rd <- abs(T)<qDunn
sD   <-   ifelse(rd==TRUE,"Accept","Reject") 

### Scheffe

qscheffe  <- sqrt(qchisq(0.95,(p-1)))  

rs <- abs(T)< qscheffe
sS     <-  ifelse( rs ==TRUE,"Accept","Reject") 


res <- list(m = m, n = n, p= p, Cov_beta= covcl$Cov_beta , T_stat= T,normal_quantile =quantile , MNQ = sM  , Bon = sB ,  Sidak = sD , Scheffe = sS)
return(res)

}

