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
size <- states[,c(1, 5)]
colnames(size) <- c("name", "size")
size <- size[is.na(size$size) == F,]
size <- size[is.na(size$name) == F,]


#convert intermediate from 0 -> 1 and small from 1 -> 0
size$size <- gsub("M", "1_1", size$size) 
size$size <- gsub("S", "1_0", size$size)
size$size <- gsub("L", "0_1", size$size)

size <- data.frame(species = size$name, small = tstrsplit(size$size, "_")[[1]], large = tstrsplit(size$size, "_")[[2]])



#add tip states
states_small <- as.numeric(size[,2]) - 1
names(states_small) <- as.character(size$species)
states_large <- as.numeric(size[,3]) - 1
names(states_large) <- as.character(size$species)


#match to phylo
pom_small <- match.phylo.data(phy,states_small)
pom_large <- match.phylo.data(phy, states_large)

#sampling fraction
f <- c(0,1,1,1)

#states
state <- data.frame(species = pom_small$phy$tip.label, small = pom_small$data, large = pom_large$data, row.names = NULL)



##NO HIDDEN STATE
#Model 1: No effect of size on diversification (homogeneous diversification independent of states)
#turnover
turnover_m1 <- c(0,1,1,1)
extinction.fraction_m1 <- c(0,1,1,1)

#make transition matrix without a hidden state
trans.rate_m1 <- TransMatMakerMuHiSSE(hidden.traits = 0, include.diagonals = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m1 <- ParDrop(trans.rate_m1, c(1,2,3,4,7,10))


#MuHisse_m1 <- MuHiSSE(phy = pom_small$phy, data = state, f = f, turnover = turnover_m1, eps = extinction.fraction_m1, hidden.states = F, trans.rate = trans.rate.mod_m1, condition.on.survival = T)



#Model 2: Effect of state on diversification (MuSSE)
#turnover
turnover_m2 <- c(0,1,2,3)
extinction.fraction_m2 <- c(0,1,2,3)

#make transition matrix without a hidden state
trans.rate_m2 <- TransMatMakerMuHiSSE(hidden.traits = 0, include.diagonals = T)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m2 <- ParDrop(trans.rate_m1, c(1,2,3,4,7,10))


#MuHisse_m2 <- MuHiSSE(phy = pom_small$phy, data = state, f = f, turnover = turnover_m2, eps = extinction.fraction_m2, hidden.states = F, trans.rate = trans.rate.mod_m2)



#Model 3: MuHISSE: state dependent diversification across the tree that are independent of state
#turnover, exctinction, f
turnover_m3 <- c(0,1,2,3,0,4,5,6)
extinction.fraction_m3 <- c(0,1,2,3,0,4,5,6)

#make transition matrix without a hidden state
trans.rate_m3 <- TransMatMakerMuHiSSE(hidden.traits = 1)

#remove transitions to and from the 1st state in the model (00)
trans.rate.mod_m3 <- ParDrop(trans.rate_m3, c(1,2,3,5))
trans.rate.mod_m3[1,5] <- 0
trans.rate.mod_m3[5,1] <- 0
trans.rate.mod_m3[5,6] <- 0
trans.rate.mod_m3[5,7] <- 0
trans.rate.mod_m3[6,5] <- 0
trans.rate.mod_m3[7,5] <- 0
trans.rate.mod_m3[trans.rate.mod_m3 == 8] <- 5
trans.rate.mod_m3[trans.rate.mod_m3 == 10] <- 6
trans.rate.mod_m3[trans.rate.mod_m3 == 11] <- 7
trans.rate.mod_m3[trans.rate.mod_m3 == 12] <- 8
trans.rate.mod_m3[trans.rate.mod_m3 == 13] <- 9

state <- data.frame(species = pom_small$phy$tip.label, small = pom_small$data, large = pom_large$data, row.names = NULL)

#MuHisse_m3 <- MuHiSSE(phy = pom_small$phy, data = state, f = f, turnover = turnover_m3, eps = extinction.fraction_m3, hidden.states = T, trans.rate = trans.rate.mod_m3)




StartingValueTry(phy = pom_small$phy, data = state, f = f, trans.rate = trans.rate.mod_m1, tau = turnover_m1, eps = extinction.fraction_m1, hidden = F, name = "Pom_MuHisse_m1")
StartingValueTry(phy = pom_small$phy, data = state, f = f, trans.rate = trans.rate.mod_m2, tau = turnover_m2, eps = extinction.fraction_m2, hidden = F, name = "Pom_MuHisse_m2")
StartingValueTry(phy = pom_small$phy, data = state, f = f, trans.rate = trans.rate.mod_m3, tau = turnover_m3, eps = extinction.fraction_m3, hidden = T, name = "Pom_MuHisse_m3")