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
phy <- read.tree("data/Damsel345.phy")
phy$tip.label[293] <- "Amblypomacentrus_tricinctus"
phy$tip.label[50] <- "Stegastes_rectifraenum"


#load in traits and extract
states <- read.csv("data/S3Table_DamselTraitsRev.csv")
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



#Model 6: MuHISSE: Pelagic Absorbing
#turnover, exctinction, f
turnover_m6 <- c(0,1,2,3,0,4,5,6)
extinction.fraction_m6 <- c(0,1,2,3,0,4,5,6)

#make transition matrix without a hidden state
trans.rate_m6 <- TransMatMakerMuHiSSE(hidden.traits = 1)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m6 <- ParDrop(trans.rate_m6, c(1,2,3,5))
trans.rate.mod_m6[1,5] <- 0
trans.rate.mod_m6[5,1] <- 0
trans.rate.mod_m6[5,6] <- 0
trans.rate.mod_m6[5,7] <- 0
trans.rate.mod_m6[6,5] <- 0
trans.rate.mod_m6[7,5] <- 0
trans.rate.mod_m6[trans.rate.mod_m6 == 3] <- 0
trans.rate.mod_m6[3,4] <- 3
trans.rate.mod_m6[trans.rate.mod_m6 == 8] <- 4
trans.rate.mod_m6[trans.rate.mod_m6 == 10] <- 5
trans.rate.mod_m6[trans.rate.mod_m6 == 11] <- 0
trans.rate.mod_m6[trans.rate.mod_m6 == 12] <- 6
trans.rate.mod_m6[trans.rate.mod_m6 == 13] <- 7
trans.rate.mod_m6[2,6] <- 8
trans.rate.mod_m6[3,7] <- 8
trans.rate.mod_m6[4,8] <- 8


#MuHisse_m6 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m6, eps = extinction.fraction_m6, hidden.states = T, trans.rate = trans.rate.mod_m6)



#Model 7: CID2
#turnover, exctinction, f
turnover_m7 <- c(0,1,1,1,0,2,2,2)
extinction.fraction_m7 <- c(0,1,1,1,0,2,2,2)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m7 <- TransMatMakerMuHiSSE(hidden.traits = 1,make.null = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m7 <- ParDrop(trans.rate_m7, c(1,2,3,5))
trans.rate.mod_m7[1,5] <- 0
trans.rate.mod_m7[5,1] <- 0
trans.rate.mod_m7[5,6] <- 0
trans.rate.mod_m7[5,7] <- 0
trans.rate.mod_m7[6,5] <- 0
trans.rate.mod_m7[7,5] <- 0
trans.rate.mod_m7[2,6] <- 6
trans.rate.mod_m7[3,7] <- 6
trans.rate.mod_m7[4,8] <- 6


#MuHisse_m7 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m7, eps = extinction.fraction_m7, hidden.states = T, trans.rate = trans.rate.mod_m7)




#Model 8: CID3
#turnover, exctinction, f
turnover_m8 <- c(0,1,1,1,0,2,2,2,0,3,3,3)
extinction.fraction_m8 <- c(0,1,1,1,0,2,2,2,0,3,3,3)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m8 <- TransMatMakerMuHiSSE(hidden.traits = 2, make.null = T, include.diagonals = F, cat.trans.vary = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m8 <- ParDrop(trans.rate_m8, c(1,2,3,5))
trans.rate.mod_m8[1,] <- 0
trans.rate.mod_m8[5,] <- 0
trans.rate.mod_m8[9,] <- 0
trans.rate.mod_m8[,1] <- 0
trans.rate.mod_m8[,5] <- 0
trans.rate.mod_m8[,9] <- 0
diag(trans.rate.mod_m8) <- NA




#MuHisse_m8 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m8, eps = extinction.fraction_m8, hidden.states = T, trans.rate = trans.rate.mod_m8)


StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m6, tau = turnover_m6, eps = extinction.fraction_m6, hidden = T, name = "Pom_MuHisse_m6")
StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m7, tau = turnover_m7, eps = extinction.fraction_m7, hidden = T, name = "Pom_MuHisse_m7")
StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m8, tau = turnover_m8, eps = extinction.fraction_m8, hidden = T, name = "Pom_MuHisse_m8")

