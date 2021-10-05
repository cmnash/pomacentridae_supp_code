#this script contains code for phylogenetic signal, phylogenetic anova, and lineage through time plots

library(geiger)
library(car)
library(tidyverse)
library(picante)
library(data.table)
library(picante)
library(viridis)
library(phytools)

phy <- read.tree("~/data/Damsel345.phy") %>% ladderize()
#phy <- force.ultrametric(phy, method = "nnls")

phy$tip.label[293] <- "Amblypomacentrus_tricinctus"
phy$tip.label[50] <- "Stegastes_rectifraenum"


states <- read.csv("~/data/S3Table_DamselTraitsRev.csv", blank.lines.skip = T)
trait <- states[,c(1, 2)]
colnames(trait) <- c("name", "diet")
trait <- trait[is.na(trait$diet) == F,]
trait <- trait[is.na(trait$name) == F,]


#convert intermediate from 0 -> 1 and benthic from 1 -> 0
trait$diet <- gsub(0, "Int", trait$diet) 
trait$diet <- gsub(1, "Benthic", trait$diet)
trait$diet <- gsub(2, "Pelagic", trait$diet)


#body size
tl <- states[,c(1, 7)]
colnames(tl) <- c("name", "TL")
tl <- tl[is.na(tl$TL) == F,]
tl <- tl[is.na(tl$name) == F,]

## FARMING ##
farm <- states[,c(1, 4)]
colnames(farm) <- c("name", "farm")
farm <- farm[is.na(farm$farm) == F,]
farm <- farm[is.na(farm$name) == F,]

#convert intermediate from 0 -> 1 and benthic from 1 -> 0
farm$farm <- gsub(0, "No Farming", farm$farm) 
farm$farm <- gsub(1, "Farming", farm$farm)

#add tip states
states_farm <- farm$farm
names(states_farm) <- farm$name




## SIZE 3 STATE ##
size3 <- states[,c(1, 5)]
colnames(size3) <- c("name", "size3")
size3 <- size3[is.na(size3$size3) == F,]
size3 <- size3[is.na(size3$name) == F,]


#add tip states
states_size3 <- as.vector(size3[,2])
names(states_size3) <- size3$name

#match to phylo
pom_size3 <- match.phylo.data(phy, states_size3)

#body size
bd <- states[,c(1, 8)]
colnames(bd) <- c("name", "BD")
bd <- bd[is.na(bd$BD) == F,]
bd <- bd[is.na(bd$name) == F,]

#join them together
state <- inner_join(trait, tl) 
state <- inner_join(state, bd)

#add tip states
states_diet <- state$diet
names(states_diet) <- state$name

#add tip states
states_tl <- state$TL
names(states_tl) <- state$name

#add tip states
states_bd <- state$BD
names(states_bd) <- state$name

#match to phylo
pom_diet <- match.phylo.data(phy, states_diet)
pom_tl <- match.phylo.data(phy, states_tl)
pom_bd <- match.phylo.data(phy, states_bd)
pom_farm <- match.phylo.data(phy, states_farm)


#calculate the phylogenetic ANOVA with ecotypes as groups
p_anv_diet_bd <- phylANOVA(tree = pom_diet$phy, x = pom_diet$data, y = pom_bd$data, nsim = 10000, posthoc = T, p.adj = "holm")

p_anv_diet_tl <- phylANOVA(tree = pom_diet$phy, x = pom_diet$data, y = pom_tl$data, nsim = 10000, posthoc = T, p.adj = "holm")

#calculate the phylogenetic ANOVA with ecotypes as groups
p_anv_farm_bd <- phylANOVA(tree = pom_farm$phy, x = pom_farm$data, y = pom_bd$data, nsim = 10000, posthoc = T, p.adj = "holm")

p_anv_farm_tl <- phylANOVA(tree = pom_farm$phy, x = pom_farm$data, y = pom_tl$data, nsim = 10000, posthoc = T, p.adj = "holm")





#leveneTest(y = pom_r$data, group = rdiet)



##Phylogenetic signal

#add tip states
states_diet <- trait$diet
names(states_diet) <- trait$name

#add tip states
states_tl <- tl$TL
names(states_tl) <- tl$name

#add tip states
states_bd <- bd$BD
names(states_bd) <- bd$name

#match to phylo
pom_diet <- match.phylo.data(phy, states_diet)
pom_tl <- match.phylo.data(phy, states_tl)
pom_bd <- match.phylo.data(phy, states_bd)



#compute K and conduct a hypothesis test of the null hypothesis of no signal
pom_tl_k <- phylosig(tree = pom_tl$phy, x = pom_tl$data, test = T )

#estimate lambda using maximum likelihood
pom_tl_lambda <- phylosig(tree = pom_tl$phy, x = pom_tl$data, method = "lambda", test = T )


#compute K and conduct a hypothesis test of the null hypothesis of no signal
pom_bd_k <- phylosig(tree = pom_bd$phy, x = pom_bd$data, test = T )

#estimate lambda using maximum likelihood
pom_bd_lambda <- phylosig(tree = pom_bd$phy, x = pom_bd$data, method = "lambda", test = T) 


#compute K and conduct a hypothesis test of the null hypothesis of no signal
#pom_diet_k <- phylosig(tree = pom_diet$phy, x = pom_diet$data, test = T )

#estimate lambda using maximum likelihood
#pom_diet_lambda <- phylosig(tree = pom_diet$phy, x = pom_diet$data, method = "lambda", test = T) 

#fit discrete lambda model
#create lambda 0 tree (no phylo signal)
pom_phy_lambda0 <- rescale(pom_diet$phy, model = "lambda", 0)

#ARD model
m2 <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)

#fit model with actual phy
fit_pom_ace_ARD_geiger_lambda <-
  fitDiscrete(phy = pom_diet$phy,
              dat = pom_diet$data,
              model = m2,
              type = "discrete",
              transform = "lambda",
              niter = 1000
  )

#fit model with lambda 0 phy
fit_pom_ace_ARD_geiger_lambda0 <-
  fitDiscrete(phy = pom_phy_lambda0,
              dat = pom_diet$data,
              model = m2,
              type = "discrete",
              transform = "lambda",
              niter = 1000
  )


#fit with delta
fit_pom_ace_ARD_geiger_delta <-
  fitDiscrete(phy = pom_diet$phy,
              dat = pom_diet$data,
              model = m2,
              type = "discrete",
              transform = "delta",
              niter = 500
  )

#if delta = 0
pom_phy_delta0 <- rescale(pom_diet$phy, model = "delta", 0)

#fit with delta0
fit_pom_ace_ARD_geiger_delta0 <-
  fitDiscrete(phy = pom_phy_delta0,
              dat = pom_diet$data,
              model = m2,
              type = "discrete",
              transform = "delta",
              niter = 500
  )






#fit discrete lambda model
#create lambda 0 tree (no phylo signal)
pom_phy_lambda0 <- rescale(pom_size5$phy, model = "lambda", 0)

#ARD model
m2 <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)

#fit model with actual phy
fit_pom_ace_ARD_geiger_lambda <-
  fitDiscrete(phy = pom_size5$phy,
              dat = pom_size5$data,
              model ="ARD",
              type = "discrete",
              transform = "lambda",
              niter = 1000
  )

#fit model with lambda 0 phy
fit_pom_ace_ARD_geiger_lambda0 <-
  fitDiscrete(phy = pom_phy_lambda0,
              dat = pom_size5$data,
              model = "ARD",
              type = "discrete",
              transform = "lambda",
              niter = 1000
  )


#fit with delta
fit_pom_ace_ARD_geiger_delta <-
  fitDiscrete(phy = pom_size5$phy,
              dat = pom_size5$data,
              model = "ARD",
              type = "discrete",
              transform = "delta",
              niter = 500
  )

#if delta = 0
pom_phy_delta0 <- rescale(pom_size5$phy, model = "delta", 0)

#fit with delta0
fit_pom_ace_ARD_geiger_delta0 <-
  fitDiscrete(phy = pom_phy_delta0,
              dat = pom_size5$data,
              model = "ARD",
              type = "discrete",
              transform = "delta",
              niter = 500
  )


##lineage through time plot
lt <- ltt(pom_diet$phy,show.tree=TRUE,lwd=3, gamma = T, n = 6)
gt <- gtt(pom_diet$phy, n = 100)
mc <- mccr(lt, rho = sampling.f, nsim = 10000)
pdf("Pom_LTT.pdf", height = 8, width = 12)
ltt(pom_diet$phy,show.tree=TRUE,lwd=3)
dev.off()
