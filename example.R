setwd("github")

source("LR-kinbiplot.R")


# logratio biplot with "UN","6TH","5TH","4TH","3RD","2ND","3.4S","FS" relationships 
# and peel and zoom approach

logpca = LR_kinbiplot(data.in="data/CEU",
                                data.out = "CEU_population",
                                nsim.rel = 20,
                                peel.and.zoom = T)

table(logpca$pairs.ids$prediction)


# logratio biplot with "UN","6TH","5TH","4TH","3RD","2ND","3.4S","FS" relationships 

logpca = LR_kinbiplot(data.in="data/CEU",
                                data.out = "CEU_population",
                                nsim.rel = 20,
                                peel.and.zoom = F)


# logratio biplot with "UN","4TH","3RD","2ND","3.4S","FS" relationships 


logpca = LR_kinbiplot(data.in="data/CEU",
                                data.out = "CEU_population",
                                nsim.rel = 20,
                                peel.and.zoom = F,
                                relationships = c("UN","4TH","3RD","2ND","3.4S","FS"),
                                rel.colors = c("forestgreen","cyan2",
                                              "darkgoldenrod2","darkorchid2","black","dodgerblue2"))

