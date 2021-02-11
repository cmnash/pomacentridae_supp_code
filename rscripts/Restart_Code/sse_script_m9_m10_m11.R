library(hisse)
library(parallel)
library(ape)
library(diversitree)
library(data.table)
library(picante)
library(phytools)

StartingValueTry <-
  function(phy, 
           data, 
           f,
           trans.rate, 
           root.type = "herr_als", 
           tau, 
           eps,
           hidden,
           name) {
    freqs <- table(apply(data[, 2:3], 1,
                         function(x)
                           switch(
                             paste0(x, collapse = ""),
                             "00" = 1,
                             "01" = 2,
                             "10" = 3,
                             "11" = 4
                           )))
    samp.freq.tree <- Ntip(phy) / sum(freqs / f[2:4])
    init.pars <-
      hisse:::starting.point.generator(phy, 4, samp.freq.tree, yule = FALSE)
    
    turnover <- tau
    eps <- eps
    NewStarting <- function(iteration) {
      turn.start <- exp(rnorm(4, log(init.pars[1] + init.pars[5]), 1))
      eps.start <- runif(1, 0, 1)
      trans.start <- exp(rnorm(12, log(init.pars[9])))
      starting.vals <- c(turn.start, rep(eps.start, 4), trans.start)
      print(starting.vals)
      tmp <-
        MuHiSSE(
          phy,
          data,
          f = f,
          turnover = turnover,
          eps = eps,
          trans.rate = trans.rate,
          hidden.states = hidden,
          starting.vals = starting.vals,
          root.type = root.type
        )
      save(tmp,
           file = paste("output/", name, "_", iteration, ".Rsave", sep = ""))
    }
    mclapply(1:100, NewStarting, mc.cores = 28)
  }





#load in phylogeny
phy <- read.nexus("data/DamselChrono330.phy")

#load in traits and extract
states <- read.csv("data/DamTraits3.csv")
diet <- states[,c(1, 2)]
colnames(diet) <- c("name", "diet")
diet <- diet[is.na(diet$diet) == F,]
diet <- diet[is.na(diet$name) == F,]


#convert intermediate from 0 -> 1 and benthic from 1 -> 0
diet$diet <- gsub(0, "I", diet$diet) 
diet$diet <- gsub(1, "P", diet$diet)
diet$diet <- gsub(2, "B", diet$diet)

diet$diet <- gsub("I", "1_1", diet$diet) 
diet$diet <- gsub("B", "1_0", diet$diet)
diet$diet <- gsub("P", "0_1", diet$diet)

diet <- data.frame(species = diet$name, benthic = tstrsplit(diet$diet, "_")[[1]], pelagic = tstrsplit(diet$diet, "_")[[2]])



#add tip states
states_ben <- as.numeric(diet[,2]) - 1
names(states_ben) <- as.character(diet$species)
states_pel <- as.numeric(diet[,3]) - 1
names(states_pel) <- as.character(diet$species)


#match to phylo
pom_benthic <- match.phylo.data(phy,states_ben)
pom_pelagic <- match.phylo.data(phy, states_pel)

#sampling fraction
f <- c(1, 0.68, 0.72, 0.76)

#states
state <- data.frame(species = pom_benthic$phy$tip.label, benthic = pom_benthic$data, pelagic = pom_pelagic$data, row.names = NULL)



#Model 9: CID4
#turnover, exctinction, f
turnover_m9 <- c(0,1,1,1,0,2,2,2,0,3,3,3,0,4,4,4)
extinction.fraction_m9 <- c(0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m9 <- TransMatMakerMuHiSSE(hidden.traits = 3, make.null = T, include.diagonals = F, cat.trans.vary = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m9 <- ParDrop(trans.rate_m9, c(1,2,3,5))
trans.rate.mod_m9[1,] <- 0
trans.rate.mod_m9[5,] <- 0
trans.rate.mod_m9[9,] <- 0
trans.rate.mod_m9[13,] <- 0
trans.rate.mod_m9[,1] <- 0
trans.rate.mod_m9[,5] <- 0
trans.rate.mod_m9[,9] <- 0
trans.rate.mod_m9[,13] <- 0
diag(trans.rate.mod_m9) <- NA



state <- data.frame(species = pom_benthic$phy$tip.label, benthic = pom_benthic$data, pelagic = pom_pelagic$data, row.names = NULL)

#MuHisse_m9 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m9, eps = extinction.fraction_m9, hidden.states = T, trans.rate = trans.rate.mod_m9)



#Model 10: CID5
#turnover, exctinction, f
turnover_m10 <- c(0,1,1,1,0,2,2,2,0,3,3,3,0,4,4,4,0,5,5,5)
extinction.fraction_m10 <- c(0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m10 <- TransMatMakerMuHiSSE(hidden.traits = 4, make.null = T, include.diagonals = F, cat.trans.vary = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m10 <- ParDrop(trans.rate_m10, c(1,2,3,5))
trans.rate.mod_m10[1,] <- 0
trans.rate.mod_m10[5,] <- 0
trans.rate.mod_m10[9,] <- 0
trans.rate.mod_m10[13,] <- 0
trans.rate.mod_m10[17,] <- 0
trans.rate.mod_m10[,1] <- 0
trans.rate.mod_m10[,5] <- 0
trans.rate.mod_m10[,9] <- 0
trans.rate.mod_m10[,17] <- 0
trans.rate.mod_m10[,13] <- 0

diag(trans.rate.mod_m10) <- NA




#MuHisse_m10 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m10, eps = extinction.fraction_m10, hidden.states = T, trans.rate = trans.rate.mod_m10)



#Model 11: CID6
#turnover, exctinction, f
turnover_m11 <- c(0,1,1,1,0,2,2,2,0,3,3,3,0,4,4,4,0,5,5,5,0,6,6,6)
extinction.fraction_m11 <- rep(c(0,1,1,1), 6)
f_m11 = c(0,1,1,1)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m11 <- TransMatMakerMuHiSSE(hidden.traits = 5, make.null = T, include.diagonals = F, cat.trans.vary = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m11 <- ParDrop(trans.rate_m11, c(1,2,3,5))
trans.rate.mod_m11[1,] <- 0
trans.rate.mod_m11[5,] <- 0
trans.rate.mod_m11[9,] <- 0
trans.rate.mod_m11[13,] <- 0
trans.rate.mod_m11[17,] <- 0
trans.rate.mod_m11[21,] <- 0
trans.rate.mod_m11[,1] <- 0
trans.rate.mod_m11[,5] <- 0
trans.rate.mod_m11[,9] <- 0
trans.rate.mod_m11[,17] <- 0
trans.rate.mod_m11[,13] <- 0
trans.rate.mod_m11[,21] <- 0
diag(trans.rate.mod_m11) <- NA




#MuHisse_m11 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m11, eps = extinction.fraction_m11, hidden.states = T, trans.rate = trans.rate.mod_m11)


StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m9, tau = turnover_m9, eps = extinction.fraction_m9, hidden = T, name = "Pom_MuHisse_m9")
StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m10, tau = turnover_m10, eps = extinction.fraction_m10, hidden = T, name = "Pom_MuHisse_m10")
StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m11, tau = turnover_m11, eps = extinction.fraction_m11, hidden = T, name = "Pom_MuHisse_m11")

