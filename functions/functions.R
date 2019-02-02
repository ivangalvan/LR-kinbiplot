Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor,lwd=1)
} 


gametes <- function(x) {
  individual <- cbind(substr(x,1,1),substr(x,2,2))
  thealleles <- apply(individual,1,sample,1)
  return(thealleles)
}

zygote <- function(x,y) {
  child <- paste(x,y,sep="")
  child[child=="BA"] <- "AB"
  return(child)
}


IBSvector = function(x,y){
  
  n00 <- sum(x==0 & y==0,na.rm=TRUE)
  n11 <- sum(x==1 & y==1,na.rm=TRUE)
  n22 <- sum(x==2 & y==2,na.rm=TRUE)
  n01 <- sum(x==0 & y==1,na.rm=TRUE) + sum(x==1 & y==0,na.rm=TRUE)
  n02 <- sum(x==0 & y==2,na.rm=TRUE) + sum(x==2 & y==0,na.rm=TRUE)
  n12 <- sum(x==1 & y==2,na.rm=TRUE) + sum(x==2 & y==1,na.rm=TRUE)
  v <- c(n00=n00,n11=n11,n22=n22,n01=n01,n02=n02,n12=n12)
  return(v)
  
}

generateUN <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," UN pair \n")
    pair <- sample(1:n,2)
    i1 <- pair[1]
    i2 <- pair[2]
    out.mat[i,] <- IBSvector(HardyWeinberg::recode(Y[1:nvar,i1],"A/B"),
                               HardyWeinberg::recode(Y[1:nvar,i2],"A/B"))
    i <- i+1
  }
  return(list(out.mat=out.mat))
}


generateFS <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," FS pair \n")
    pair <- sample(1:n,2)
    i1 <- pair[1]
    i2 <- pair[2]
    kid1  <- zygote(gametes(Y[1:nvar,i1]),gametes(Y[1:nvar,i2]))
    kid1  <- HardyWeinberg::recode(kid1,"A/B")
    kid2  <- zygote(gametes(Y[1:nvar,i1]),gametes(Y[1:nvar,i2]))
    kid2  <- HardyWeinberg::recode(kid2,"A/B")
    out.mat[i,] <- IBSvector(kid1,kid2)
    i <- i+1
  }
  return(list(out.mat=out.mat))
}

generate2ND <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 2ND pair \n")
    indv <- sample(1:n,3)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
    kid1  <- zygote(gametes(Y[1:nvar,i1]),gametes(Y[1:nvar,i2]))
    kid1  <- HardyWeinberg::recode(kid1,"A/B")
    kid2  <- zygote(gametes(Y[1:nvar,i1]),gametes(Y[1:nvar,i3]))
    kid2  <- HardyWeinberg::recode(kid2,"A/B")
    out.mat[i,] <- IBSvector(kid1,kid2)
      i <- i+1
  }
  return(list(out.mat=out.mat))
}


generate3RD <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 3RD pair \n")
    indv <- sample(1:n,4)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
    i4 <- indv[4]
    gr.dad <- Y[1:nvar,i1]
    gr.mum <- Y[1:nvar,i2]
      
    child1 <- zygote(gametes(gr.dad),gametes(gr.mum))
    child2 <- zygote(gametes(gr.dad),gametes(gr.mum))
      
    outsider1 <- Y[1:nvar,i3]
    outsider2 <- Y[1:nvar,i4]
      
    grchild1 <- zygote(gametes(child1),gametes(outsider1))
    grchild1 <- HardyWeinberg::recode(grchild1,alleles="A/B")
      
    grchild2 <- zygote(gametes(child2),gametes(outsider2))
    grchild2 <- HardyWeinberg::recode(grchild2,alleles="A/B")
      
    out.mat[i,] <- IBSvector(grchild1,grchild2)
      
    i <- i+1
    }
  return(list(out.mat=out.mat))
}

generate4TH <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 4TH pair \n")
    indv <- sample(1:n,5)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
    i4 <- indv[4]
    i5 <- indv[5]
    gr.dad <- Y[1:nvar,i1]
    gr.mum <- Y[1:nvar,i2]
      
    child1 <- zygote(gametes(gr.dad),gametes(gr.mum))
    child2 <- zygote(gametes(gr.dad),gametes(gr.mum))
      
    outsider1 <- Y[1:nvar,i3]
    outsider2 <- Y[1:nvar,i4]
    outsider3 <- Y[1:nvar,i5]
      
    grchild1 <- zygote(gametes(child1),gametes(outsider1))
      
    grchild2 <- zygote(gametes(child2),gametes(outsider2))
    grchild2 <- HardyWeinberg::recode(grchild2,alleles="A/B")
      
    gr.grchild1 <- zygote(gametes(grchild1),gametes(outsider3))
    gr.grchild1 <- HardyWeinberg::recode(gr.grchild1,alleles="A/B")
      
      
    out.mat[i,] <- IBSvector(gr.grchild1,grchild2)
      
    i <- i+1
  }
  return(list(out.mat=out.mat))
}


generate5TH <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 5TH pair \n")
    indv <- sample(1:n,6)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
    i4 <- indv[4]
    i5 <- indv[5]
    i6 <- indv[6]
    gr.dad <- Y[1:nvar,i1]
    gr.mum <- Y[1:nvar,i2]
      
    child1 <- zygote(gametes(gr.dad),gametes(gr.mum))
    child2 <- zygote(gametes(gr.dad),gametes(gr.mum))
      
    outsider1 <- Y[1:nvar,i3]
    outsider2 <- Y[1:nvar,i4]
    outsider3 <- Y[1:nvar,i5]
    outsider4 <- Y[1:nvar,i6]
      
      
    grchild1 <- zygote(gametes(child1),gametes(outsider1))
    grchild2 <- zygote(gametes(child2),gametes(outsider2))
      
    gr.grchild1 <- zygote(gametes(grchild1),gametes(outsider3))
    gr.grchild1 <- HardyWeinberg::recode(gr.grchild1,alleles="A/B")
      
    gr.grchild2 <- zygote(gametes(grchild2),gametes(outsider4))
    gr.grchild2 <- HardyWeinberg::recode(gr.grchild2,alleles="A/B")
      
    out.mat[i,] <- IBSvector(gr.grchild1,gr.grchild2)
      
    i <- i+1
    
  }
  return(list(out.mat=out.mat))
}


generate6TH <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 6TH pair \n")
    indv <- sample(1:n,7)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
    i4 <- indv[4]
    i5 <- indv[5]
    i6 <- indv[6]
    i7 <- indv[7]
    gr.dad <- Y[1:nvar,i1]
    gr.mum <- Y[1:nvar,i2]
      
    child1 <- zygote(gametes(gr.dad),gametes(gr.mum))
    child2 <- zygote(gametes(gr.dad),gametes(gr.mum))
      
    outsider1 <- Y[1:nvar,i3]
    outsider2 <- Y[1:nvar,i4]
    outsider3 <- Y[1:nvar,i5]
    outsider4 <- Y[1:nvar,i6]
    outsider5 <- Y[1:nvar,i7]
      
      
    grchild1 <- zygote(gametes(child1),gametes(outsider1))
    grchild2 <- zygote(gametes(child2),gametes(outsider2))
      
    gr.grchild1 <- zygote(gametes(grchild1),gametes(outsider3))
      
    gr.gr.grchild1 <- zygote(gametes(gr.grchild1),gametes(outsider5))
    gr.gr.grchild1 <- HardyWeinberg::recode(gr.gr.grchild1,alleles="A/B")
      
    gr.grchild2 <- zygote(gametes(grchild2),gametes(outsider4))
    gr.grchild2 <- HardyWeinberg::recode(gr.grchild2,alleles="A/B")
      
    out.mat[i,] <- IBSvector(gr.gr.grchild1,gr.grchild2)
      
    i <- i+1
    
  }
  return(list(out.mat=out.mat))
}


generate3.4S <- function(Y,np=20,nvar=nrow(Y)) {
  n <- ncol(Y)
  out.mat <- matrix(NA,ncol=6,nrow=np)
  i <- 1
  while(i <= np) {
    cat(i," 3/4S pair \n")
    indv <- sample(1:n,3)
    i1 <- indv[1]
    i2 <- indv[2]
    i3 <- indv[3]
	  i4 <- indv[3]
    
	  gr.dad <- Y[1:nvar,i1]
    gr.mum <- Y[1:nvar,i2]
      
    child1 <- zygote(gametes(gr.dad),gametes(gr.mum))
    child2 <- zygote(gametes(gr.dad),gametes(gr.mum))
      
    outsider1 <- Y[1:nvar,i3]
      
    grchild1 <- zygote(gametes(child1),gametes(outsider1))
    grchild1 <- HardyWeinberg::recode(grchild1,alleles="A/B")
      
    grchild2 <- zygote(gametes(child2),gametes(outsider1))
    grchild2 <- HardyWeinberg::recode(grchild2,alleles="A/B")
      
    out.mat[i,] <- IBSvector(grchild1,grchild2)
      
    i <- i+1
    
  }
  return(list(out.mat=out.mat))
}


logratiopca <- function(X,Y,relationships){
  
  # Compute Xcclr
  
  Xl = log(X)
  n=nrow(X)
  D=ncol(X)
  ones = rep(1,D)
  I = diag(D)
  Hr = I-(1/D)*(ones)%*%t(ones)
  Xclr = Xl%*%Hr
  I = diag(n)
  ones = rep(1,n)
  Hc = I-(1/n)*(ones)%*%t(ones)
  Xc_clr = Hc%*%Xclr
  
  
  # calculate coordinates
  
  svd.out = svd(Xc_clr)
  Fp = svd.out$u%*%diag(svd.out$d)
  Gs = svd.out$v
  
  # Compute Ycclr
  
  Yl = log(Y)
  Yclr = as.matrix(Yl)%*%Hr
  nsup <- nrow(Yclr)
  m <- colMeans(Xclr)
  Yc_clr = Yclr-rep(1,nsup)%o%m
  
  
  # project
  
  projection = t(solve(t(Gs)%*%Gs)%*%t(Gs)%*%t(Yc_clr))
  
  # GOF PCA
  
  pca.results <- princomp(Xclr,cor=FALSE)
  
  pca1 = round(pca.results$sdev[1]^2/sum(pca.results$sdev^2)*100,2)
  pca2 = round(pca.results$sdev[2]^2/sum(pca.results$sdev^2)*100,2)
  pca3 = round(pca.results$sdev[3]^2/sum(pca.results$sdev^2)*100,2)
  pca4 = round(pca.results$sdev[4]^2/sum(pca.results$sdev^2)*100,2)
  pca5 = round(pca.results$sdev[5]^2/sum(pca.results$sdev^2)*100,2)
  
  # predict
  
  pca_out = as.data.frame(Fp)
  
  pca_out$group = relationships
  
  lda_out = lda(relationships~V1+V2+V3,pca_out)
  
  lda_pred = predict(lda_out,as.data.frame(projection))
  
  
  return(list(empirical = projection,
              pca1=pca1,
              pca2=pca2,
              pca3=pca3,
              pca4=pca4,
              pca5=pca5,
              Fp=Fp,
              Gs=Gs,
              prediction = lda_pred$class))
  
}

