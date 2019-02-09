LR_kinbiplot <- function(data.in = data,
                                data.out = "study",
                                chr = "1-22",
                                geno = 0,
                                hwe = 0.05,
                                MAF = 0.4,
                                LD.window = 50,
                                LD.step = 5,
                                LD.correlation = 0.2,
                                rel.cutoff = 0.025,
                                mz.ibd2 = 0.7,
                                po.ibd1 = 0.7,
                                nsim.rel = 100,
                                relationships = c("UN","6TH","5TH","4TH","3RD","2ND","3.4S","FS"),
                                rel.colors = c("forestgreen","darkorange","deeppink2","cyan2",
                                               "darkgoldenrod2","darkorchid2","black","dodgerblue2"),
                                peel.and.zoom =TRUE){
  
  dir.create("outputs")
  dir.create("outputs/plots")
  
  list.of.packages <- c("data.table","dplyr","HardyWeinberg","compositions","plyr","calibrate","BMhyb")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  require(data.table)
  require(plyr)
  require(dplyr)
  require(HardyWeinberg)
  require(compositions)
  require(MASS)
  require(calibrate)
  require(BMhyb)
  
  source("functions/functions.R")
  
  if(identical(relationships,
               c("UN","6TH","5TH","4TH","3RD","2ND","3.4S","FS"))==F){
    peel.and.zoom=F}
  
  # Apply QC filters: chr, MAF, hwe, geno, LD
  
  system(paste0("plink --bfile  ",data.in," --maf ",MAF," --hwe midp ",hwe," --chr ",chr," --geno ",geno," --make-bed ",
                "--out  outputs/",data.out,"_qc"))
  
  system(paste0("plink --bfile  outputs/",data.out,"_qc -indep-pairwise ",LD.window," ",LD.step," ",LD.correlation,
                " --out outputs/",data.out))
  
  system(paste0("plink --bfile  outputs/",data.out,"_qc  --extract  outputs/",data.out,".prune.in ",
                " --make-bed --out  outputs/",data.out,"_qc2"))
  
  # Detect potential unrelated individuals with rel.cutoff
  
  system(paste0("plink --bfile outputs/",data.out,"_qc2 --rel-cutoff ",rel.cutoff,
                " --out outputs/",data.out))
  
  system(paste0("plink --bfile outputs/",data.out,"_qc2 --keep  outputs/",data.out,".rel.id  ",
                "--genome full --out outputs/",data.out,"_un"))
  
  genome = fread(paste0("outputs/",data.out,"_un.genome"))
  
  genome$p0 = genome$IBS0/(genome$IBS0+genome$IBS1+genome$IBS2)
  
  genome$p2 = genome$IBS2/(genome$IBS0+genome$IBS1+genome$IBS2)
  
  p0_p2 = genome %>% dplyr::select(p0,p2) %>% round(4) %>% unique()
  
  png(paste0("outputs/plots/",data.out,"_unrelated.png"),res=300,width = 6,height = 6,units = "in")
  plot(p0_p2$p0,p0_p2$p2,xlab=expression(p["0"]),ylab=expression(p["2"]),
       xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="blue4")
  title(paste0(data.out," unrelated pairs: rel.cutoff ",rel.cutoff))
  dev.off()
  
  # Remove potential MZ and PO pairs
  
  system(paste0("plink --bfile outputs/",data.out,"_qc2 ",
                "--genome full --out outputs/",data.out))
  
  genome = fread(paste0("outputs/",data.out,".genome"))
  
  mz_pairs = genome %>% filter(Z2>mz.ibd2)
  po_pairs = genome %>% filter(Z1>po.ibd1)
  
  mz_po_pairs = rbind(mz_pairs,po_pairs)
  
  mz_po_pairs %>% dplyr::select(IID1,IID2) %>% 
    write.table(paste0("outputs/MZ_PO_pairs.txt"),row.names = F,col.names = T,quote = F)
  
  k0_k1 = genome %>% dplyr::select(Z0,Z1) %>% round(4) %>% unique()
  
  png(paste0("outputs/plots/",data.out,"_ibd.png"),res=300,width = 6,height = 6,units = "in")
  plot(k0_k1$Z0,k0_k1$Z1,xlab=expression(hat('k')['0']),ylab="",
       xlim=c(0,1),ylim=c(0,1),pch=16,cex=1,col="blue4")
  title(ylab=expression(hat('k')['1']), line=2.5)
  title(paste0("IBD estimation: ",data.out))
  dev.off()
  
  
  # Simulate relationships
  
  system(paste0("plink --bfile outputs/",data.out,"_qc2 --keep  outputs/",data.out,".rel.id ",
                "--recode A-transpose --out outputs/",data.out,"_un"))
  
  unrelated = as.data.frame(fread(paste0("outputs/",data.out,"_un.traw")))
  
  unrelated = unrelated[,-c(1:6)]
  
  unrelated[unrelated==0] <- "AA"
  unrelated[unrelated==1] <- "AB"
  unrelated[unrelated==2] <- "BB"

    # create a list of functions with all the relationships
  
  funcList = list()
  
  for(i in 1:length(relationships)){
    
    funcList[[i]] <- paste0("generate",relationships[i])
  }
  
  # start simulations
  
  sim_data = NULL
  
  for(i in 1:length(relationships)){
    
    cat(paste0(" \n Starting ",nsim.rel," simulations for ",relationships[i]," relationships \n"))
    
    f <- get(funcList[[i]])
    sim = f(unrelated,nsim.rel)
    sim = as.data.frame(sim$out.mat)
    sim$relationship = relationships[i]
    sim_data = rbind(sim_data,sim)
  }
  
  # logratio PCA for the simulations
  
  sim_data$color = mapvalues(sim_data$relationship , from=relationships, 
                          to = rel.colors)
  
  # Compute Xcclr

  X = as.matrix(sim_data[,1:6])
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
  
  pca.results <- princomp(Xclr,cor=FALSE)
  
  pca1 = round(pca.results$sdev[1]^2/sum(pca.results$sdev^2)*100,2)
  pca2 = round(pca.results$sdev[2]^2/sum(pca.results$sdev^2)*100,2)
  
  x_lab = paste0("PC1 (",pca1,"%)")
  y_lab = paste0("PC2 (",pca2,"%)")

  rel_label = relationships
  
  if(sum(relationships %in% "3.4S")==1){
    rel_label[rel_label=="3.4S"] = "3/4S"
    }
  
  png(paste0("outputs/plots/",data.out,"_logratiopca_simulations.png"),
      res=300,width = 7,height = 7,units = "in")
  plot(Fp[,1],
       Fp[,2],col=sim_data$color,asp=1,xlab=x_lab,ylab=y_lab)
  legend("topright",rel_label,col=rel.colors,pch=16)  
  title(paste0("logratio PCA: simulated related pairs"))
    dev.off()
  
    
  # compute confussion matrix of simulated data
  
  Fp = as.data.frame(Fp)
  
  Fp$relationship = sim_data$relationship          
  
  my_lda = lda(relationship ~ V1+V2+V3,data = Fp)  
    
  my_prediction = predict(my_lda,Fp)
  
  conf_mat = table(Fp$relationship,my_prediction$class)
  
  cat("Confusion matrix for simulated relationships (%): \n ")
  print((conf_mat/nsim.rel)*100)
  
  
  # project empirical data
  
  system(paste0("plink --bfile outputs/",data.out,"_qc2 ",
                "--recode A-transpose --out outputs/",data.out))
  
  Y <- fread(paste0("outputs/",data.out,".traw"))
  Y = Y[,-c(1:6)]
  
  fam_file = fread(paste0("outputs/",data.out,"_qc2.fam"))
  
  cat(paste0(" \n Pairwise calculation of the six allele sharing IBS parts for ",ncol(Y)," individuals \n"))
  
  Y_6parts <- NULL
  pair_ids <- NULL
  
  for(i in 1:165) {
    cat(i," \n")
    for(j in 1:165) {
      
      if(i<j){
        
        aux <- IBSvector(t(Y[,i,with=F]),t(Y[,j,with=F]))
        aux2 <- c(fam_file$V2[i],fam_file$V2[j])
        
        Y_6parts <- rbind(Y_6parts,aux)
        pair_ids <- rbind(pair_ids,aux2)
        
      }
    }
  }
  
  rownames(Y_6parts) = NULL
  rownames(pair_ids) = NULL
  
  pair_ids = as.data.frame(pair_ids)
  pair_ids$pair = do.call(paste0,list(pair_ids$V1,"_",pair_ids$V2))
  mz_po_pairs$pair = do.call(paste0,list(mz_po_pairs$IID1,"_",mz_po_pairs$IID2))
  
  Y_6parts = Y_6parts[-which(pair_ids$pair %in% mz_po_pairs$pair),]
  pair_ids = pair_ids[-which(pair_ids$pair %in% mz_po_pairs$pair),]
  
  if(peel.and.zoom){
  
    relationships = c("UN","6TH","5TH","4TH","3RD","2ND","3.4S","FS")
    rel.colors = c("forestgreen","darkorange","deeppink2","cyan2",
                   "darkgoldenrod2","darkorchid2","black","dodgerblue2")
    
    rel_label = relationships
    
    rel_label[rel_label=="3.4S"] = "3/4S"

  # without PO ######
  
  png(paste0("outputs/plots/",data.out,"_logratio_biplot_pca_peel_and_zoom.png"),
      res=300,width = 10,height = 14,units = "in")
  
  par(mfrow=c(3,2))
  
  logpca = logratiopca(X,Y_6parts,relationships = sim_data$color)
  first_prediction = logpca$prediction

  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC2 (",logpca$pca2,"%)")
  
  dim1 = 1
  dim2 = 2
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:7){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  
  Gs <- logpca$Gs
  sf <- 0.5
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  legend("topright",rel_label,col=rel.colors,pch=16)  
  title("A: without PO",cex.main=2)
  
  
  # without PO, FS and 3/4S #####
  
  Y_6parts_2 = Y_6parts[-which(first_prediction %in% "dodgerblue2"),]
  rel.out = which(sim_data$relationship %in% c("FS","3.4S"))
  
  logpca = logratiopca(X[-rel.out,],Y_6parts_2,relationships = sim_data$color[-rel.out])
  
  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC2 (",logpca$pca2,"%)")
  
  dim1 = 1
  dim2 = 2
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:5){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  Gs <- logpca$Gs
  sf <- 0.4
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  title("B: without PO, FS and 3/4S",cex.main=2)
  
  
  # without PO, FS, 3/4S and 2ND ####
  
  Y_6parts_2 = Y_6parts[-which(first_prediction %in% c("dodgerblue2","darkorchid2")),]
  rel.out = which(sim_data$relationship %in% c("FS","3.4S","2ND"))
  
  logpca = logratiopca(X[-rel.out,],Y_6parts_2,relationships = sim_data$color[-rel.out])
  
  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC2 (",logpca$pca2,"%)")
  
  dim1 = 1
  dim2 = 2
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:4){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  Gs <- logpca$Gs
  sf <- 0.2
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  title("C: without PO, FS, 3/4S and 2ND",cex.main=2)
  
  
  # without PO, FS, 3/4S, 2ND and 3RD ####
  
  Y_6parts_2 = Y_6parts[-which(first_prediction %in% c("dodgerblue2",
                                                       "darkorchid2",
                                                       "darkgoldenrod2")),]
  rel.out = which(sim_data$relationship %in% c("FS","3.4S","2ND","3RD"))
  
  logpca = logratiopca(X[-rel.out,],Y_6parts_2,relationships = sim_data$color[-rel.out])
  
  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC2 (",logpca$pca2,"%)")
  
  dim1 = 1
  dim2 = 2
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:3){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  Gs <- logpca$Gs
  sf <- 0.12
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  title("D: without PO, FS, 3/4S, 2ND and 3RD",cex.main=2)
  
  # without PO, FS, 3/4S, 2ND, 3RD and 4TH #####
  
  Y_6parts_2 = Y_6parts[-which(first_prediction %in% c("dodgerblue2",
                                                       "darkorchid2",
                                                       "darkgoldenrod2",
                                                       "cyan2")),]
  
  rel.out = which(sim_data$relationship %in% c("FS","3.4S","2ND","3RD","4TH"))
  
  logpca = logratiopca(X[-rel.out,],Y_6parts_2,relationships = sim_data$color[-rel.out])
  
  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC2 (",logpca$pca2,"%)")
  
  dim1 = 1
  dim2 = 2
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:2){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  
  Gs <- logpca$Gs
  sf <- 0.12
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  title("E: without PO, FS, 3/4S, 2ND, 3RD and 4TH",cex.main=2)
  
  # without PO, FS, 3/4S, 2ND, 3RD and 4TH #####
  
  Y_6parts_2 = Y_6parts[-which(first_prediction %in% c("dodgerblue2",
                                                       "darkorchid2",
                                                       "darkgoldenrod2",
                                                       "cyan2")),]
  
  rel.out = which(sim_data$relationship %in% c("FS","3.4S","2ND","3RD","4TH"))
  
  logpca = logratiopca(X[-rel.out,],Y_6parts_2,relationships = sim_data$color[-rel.out])
  
  x_lab = paste0("PC1 (",logpca$pca1,"%)")
  y_lab = paste0("PC3 (",logpca$pca3,"%)")
  
  dim1 = 1
  dim2 = 3
  
  logpca$empirical[,1] = -logpca$empirical[,1] 
  logpca$Fp[,1] = -logpca$Fp[,1]
  
  plot(1e7,1e7,pch=1,cex=0.8,asp=1,
       xlab =x_lab, ylab=y_lab,
       xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
              max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
       ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
              max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
  
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  
  points(logpca$empirical[,dim1],logpca$empirical[,dim2],
         col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
  
  PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                  ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                  lcolor = "gray")
  
  for(i in 1:2){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
  }
  
  
  Gs <- logpca$Gs
  sf <- 0.12
  Gs <- sf*Gs
  arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
  textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
  title("F: without PO, FS, 3/4S, 2ND, 3RD and 4TH",cex.main=2)
  
  dev.off()
  
  }
  
  if(!peel.and.zoom){
    
    
    png(paste0("outputs/plots/",data.out,"_logratio_biplot_pca.png"),
        res=300,width = 7,height = 7,units = "in")

    logpca = logratiopca(X,Y_6parts,relationships = sim_data$color)
    
    x_lab = paste0("PC1 (",logpca$pca1,"%)")
    y_lab = paste0("PC2 (",logpca$pca2,"%)")
    
    dim1 = 1
    dim2 = 2
    
    logpca$empirical[,1] = -logpca$empirical[,1] 
    logpca$Fp[,1] = -logpca$Fp[,1]
    
    plot(1e7,1e7,pch=1,cex=0.8,asp=1,
         xlab =x_lab, ylab=y_lab,
         xlim=c(min(logpca$Fp[,dim1],logpca$empirical[,dim1]),
                max(logpca$Fp[,dim1],logpca$empirical[,dim1])),
         ylim=c(min(logpca$Fp[,dim2],logpca$empirical[,dim2]),
                max(logpca$Fp[,dim2],logpca$empirical[,dim2])))
    
    abline(v=0,lty=2)
    abline(h=0,lty=2)
    
    points(logpca$empirical[,dim1],logpca$empirical[,dim2],
           col=as.character(logpca$prediction),pch=1,cex=0.8,asp=1)
    
    
    PlotConvexHull(xcoord = logpca$Fp[1:(nsim.rel),dim1], 
                    ycoord = logpca$Fp[1:(nsim.rel),dim2], 
                    lcolor = "gray")
    
    
    for(i in 1:(length(relationships)-1)){
    
    PlotConvexHull(xcoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim1], 
                    ycoord = logpca$Fp[((nsim.rel)*i+1):((nsim.rel)*(i+1)),dim2], 
                    lcolor = "gray")
    }
    
    Gs <- logpca$Gs
    sf <- 0.5
    Gs <- sf*Gs
    arrows(0,0,-Gs[,1],Gs[,2],length=0.1,col="blue")
    textxy(-Gs[,1],Gs[,2], c(expression(k["00"]),expression(k["11"]),expression(k["22"]),                         expression(k["01"]),expression(k["02"]),expression(k["12"])),cex=0.75)
    legend("topright",rel_label,col=rel.colors,pch=16)  

    dev.off()
    
  }
  
  logpca = logratiopca(X,Y_6parts,relationships = sim_data$relationship)

  colnames(pair_ids)[1:2] = c("ID1","ID2")
  
  pair_ids$prediction = logpca$prediction
  
  genome$pair = do.call(paste0,list(genome$IID1,"_",genome$IID2))

  pair_ids = left_join(pair_ids,genome %>% dplyr::select(pair,Z0,Z1,Z2,PI_HAT))
  
  pair_ids %>% 
    write.table(paste0("outputs/",data.out,"_predicted_pairs.txt"),
                row.names = F,col.names = T,quote = F)
  
  cat(paste0(" \n Predicted relationships are saved in: outputs/",data.out,"_predicted_pairs.txt \n"))
  
  cat(paste0(" \n The predicted relationships are: \n"))
  
  print(table(logpca$prediction))
  
  empirical = logpca$empirical
  rownames(empirical) = NULL
  
  return(list(empirical = empirical,
              pca1=logpca$pca1,
              pca2=logpca$pca2,
              pca3=logpca$pca3,
              pca4=logpca$pca4,
              pca5=logpca$pca5,
              Fp=logpca$Fp,
              Gs=logpca$Gs,
              pairs.ids = pair_ids,
              conf.matrix.sim = (conf_mat/nsim.rel)*100))
  
}



