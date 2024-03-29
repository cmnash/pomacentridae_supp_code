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


#Model 4: Shifts in diversification across the tree that are dependent on state (with hidden state)
#turnover, exctinction, f
turnover_m4 <- c(0,1,2,3,0,4,5,6)
extinction.fraction_m4 <- c(0,1,2,3,0,4,5,6)

#make transition matrix without a hidden state. This way the transition rates are not impacted by changes in the diversification rate regime -- that is, $q_{MB,A\rightarrow MP,A} = q_{MB,B\rightarrow MP,B}$. 
trans.rate_m4 <- TransMatMakerMuHiSSE(hidden.traits = 1, make.null = F)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m4 <- ParDrop(trans.rate_m4, c(1,2,3,5))
trans.rate.mod_m4[1,5] <- 0
trans.rate.mod_m4[5,1] <- 0
trans.rate.mod_m4[5,6] <- 0
trans.rate.mod_m4[5,7] <- 0
trans.rate.mod_m4[6,5] <- 0
trans.rate.mod_m4[7,5] <- 0
trans.rate.mod_m4[8,6] <- 5
trans.rate.mod_m4[8,7] <- 6
trans.rate.mod_m4[6,8] <- 7
trans.rate.mod_m4[7,8] <- 8
trans.rate.mod_m4[6,2] <- 9
trans.rate.mod_m4[7,3] <- 9
trans.rate.mod_m4[8,4] <- 9
trans.rate.mod_m4[2,6] <- 10
trans.rate.mod_m4[3,7] <- 10
trans.rate.mod_m4[4,8] <- 10


#MuHisse_m4 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m4, eps = extinction.fraction_m4, hidden.states = T, trans.rate = trans.rate.mod_m4)

#Model 5: MuHISSE: Benthic Absorbing
#turnover, exctinction, f
turnover_m5 <- c(0,1,2,3,0,4,5,6)
extinction.fraction_m5 <- c(0,1,2,3,0,4,5,6)

#make transition matrix without a hidden state
trans.rate_m5 <- TransMatMakerMuHiSSE(hidden.traits = 1)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m5 <- ParDrop(trans.rate_m5, c(1,2,3,5))
trans.rate.mod_m5[1,5] <- 0
trans.rate.mod_m5[5,1] <- 0
trans.rate.mod_m5[5,6] <- 0
trans.rate.mod_m5[5,7] <- 0
trans.rate.mod_m5[6,5] <- 0
trans.rate.mod_m5[7,5] <- 0
trans.rate.mod_m5[trans.rate.mod_m5 == 4] <- 0
trans.rate.mod_m5[trans.rate.mod_m5 == 8] <- 4
trans.rate.mod_m5[trans.rate.mod_m5 == 10] <- 5
trans.rate.mod_m5[trans.rate.mod_m5 == 11] <- 6
trans.rate.mod_m5[trans.rate.mod_m5 == 12] <- 0
trans.rate.mod_m5[trans.rate.mod_m5 == 13] <- 7
trans.rate.mod_m5[2,6] <- 8
trans.rate.mod_m5[3,7] <- 8
trans.rate.mod_m5[4,8] <- 8


#MuHisse_m5 <- MuHiSSE(phy = pom_benthic$phy, data = state, f = f, turnover = turnover_m5, eps = extinction.fraction_m5, hidden.states = T, trans.rate = trans.rate.mod_m5)



StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m4, tau = turnover_m4, eps = extinction.fraction_m4, hidden = T, name = "Pom_MuHisse_m4")
StartingValueTry(phy = pom_benthic$phy, data = state, f = f, trans.rate = trans.rate.mod_m5, tau = turnover_m5, eps = extinction.fraction_m5, hidden = T, name = "Pom_MuHisse_m5")

