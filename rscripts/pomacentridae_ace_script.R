setwd("~/Box/Pom_Phylo_Paper/")
## this script was used to estimate and plot ancestral states on the phylogeny and examine their transitions
## used in Fig 3 and 4
library(tidyverse)
library(ape)
library(phytools)
library(picante)
library(viridis)
library(geiger)
library(RColorBrewer)
library(fishualize)
library(diversitree)

phy <- read.tree("~/data/Damsel345.phy") %>% ladderize()
phy <- force.ultrametric(phy)
#phy with just pomacentrinae
#pom_phy <- extract.clade(phy, node = 496)

#set clade labels
#labels <- c("Microspathodontidae", "Chrominae","Glyphisodonitdae", "Pomacentrinae")
#nodes <- c(332, 385, 476, 496)

phy$tip.label[293] <- "Amblypomacentrus_tricinctus"
phy$tip.label[50] <- "Stegastes_rectifraenum"


# load traits
states <- read.csv("~/data/S3Table_DamselTraitsRev.csv", blank.lines.skip = T)


## DIET ##
diet <- states[,c(1, 2)]
colnames(diet) <- c("name", "diet")
diet <- diet[is.na(diet$diet) == F,]
diet <- diet[is.na(diet$name) == F,]

#convert values
diet$diet <- gsub(0, "I", diet$diet) 
diet$diet <- gsub(1, "P", diet$diet)
diet$diet <- gsub(2, "B", diet$diet)

#add tip states
states_diet <- diet[,2]
names(states_diet) <- diet$name

#match to phylo
pom_diet <- match.phylo.data(phy, states_diet)

#set up transition matrices that corrrespond to models
#model 1: all rates equal
m1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3)
rownames(m1) <- c("Benthic", "Int", "Pelagic")
colnames(m1) <- c("Benthic", "Int", "Pelagic")

fit_diet_m1 <- ace(pom_diet$data, pom_diet$phy, model = m1, type = "discrete")
fd_1 <- fitDiscrete(pom_diet$phy, dat = pom_diet$data, model = m1)

#model 2: all rates different (constrained)
m2 <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)
rownames(m2) <- c("Benthic", "Int", "Pelagic")
colnames(m2) <- c("Benthic", "Int", "Pelagic")
fit_diet_m2 <- ace(pom_diet$data, pom_diet$phy, model = m2, type = "discrete")
fd_2 <- fitDiscrete(pom_diet$phy, dat = pom_diet$data, model = m2)

#model 3: symmetrical (1 -> 2 = 2 -> 1)
m3 <- matrix(c(0, 1, 0, 1, 0, 2, 0, 2, 0), 3)
rownames(m1) <- c("Benthic", "Int", "Pelagic")
colnames(m1) <- c("Benthic", "Int", "Pelagic")
fit_diet_m3 <- ace(pom_diet$data, pom_diet$phy, model = m3, type = "discrete")
fd_3 <- fitDiscrete(pom_diet$phy, dat = pom_diet$data, model = m3)

#compare AIC scores 
fd_1[[4]]$aicc #ER
fd_2[[4]]$aicc #ARD
fd_3[[4]]$aicc #SYM

#set colors for plotting
cols <-c("maroon","darkorange","deepskyblue2")
names(cols) <- c("B", "I", "P")

pdf("Pomacentridae_Ecotype_ACE_newcols.pdf", height = 15, width = 15)

#par(mar=c(5,5,5,5))#, oma = c(1,1,1,1))

#plot 
plotTree(pom_diet$phy, type = "fan", fsize = 0.3, type = "i", offset = 6, tips = seq(5, 349, 1), asp = 0.5)
#for(i in 1:4) {
  #arc.cladelabels(text=labels[i],node=nodes[i], lab.offset = 1.33, ln.offset = 1.28, cex = 0.7, lwd = 3.5)}


nodelabels(node=1:pom_diet$phy$Nnode+Ntip(pom_diet$phy), pie=fit_diet_m2$lik.anc,piecol=cols[colnames(fit_diet_m2$lik.anc)],cex=0.3)

tiplabels(pie=to.matrix(pom_diet$data,unique(pom_diet$data)),piecol=cols[unique(pom_diet$data)], cex = 0.2)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(pom_diet$phy)),fsize=0.8)

dev.off()

#plot transitions
ftmk_2_diet <- fitMk(tree = pom_diet$phy, x = pom_diet$data, model = m2)
plot(ftmk_2_diet)




## FARMING ##
farm <- states[,c(1, 4)]
colnames(farm) <- c("name", "farm")
farm <- farm[is.na(farm$farm) == F,]
farm <- farm[is.na(farm$name) == F,]

#convert intermediate from 0 -> 1 and benthic from 1 -> 0
farm$farm <- gsub(0, "No Farming", farm$farm) 
farm$farm <- gsub(1, "Farming", farm$farm)

#add tip states
states_farm <- farm[,2]
names(states_farm) <- farm$name

#match to phylo
pom_farm <- match.phylo.data(phy, states_farm)

#set up transition matrices that corrrespond to models
#model 1: all rates equal
fit_farm_m1 <- ace(pom_farm$data, pom_farm$phy, model = "ER", type = "discrete")
fd_1_farm <- fitDiscrete(pom_farm$phy, dat = pom_farm$data, model = "ER")

#model 2: all rates different
fit_farm_m2 <- ace(pom_farm$data, pom_farm$phy, model = "ARD", type = "discrete")
fd_2_farm <- fitDiscrete(pom_farm$phy, dat = pom_farm$data, model = "ARD")

#model 3: symmetrical (1 -> 2 = 2 -> 1)
fit_farm_m3 <- ace(pom_farm$data, pom_farm$phy, model = "SYM", type = "discrete")
fd_3_farm <- fitDiscrete(pom_farm$phy, dat = pom_farm$data, model = "SYM")

#compare AIC scores 
fd_1_farm[[4]]$aicc #ER
fd_2_farm[[4]]$aicc #ARD
fd_3_farm[[4]]$aicc #SYM

#set colors for plotting
cols <- c("darkgoldenrod1", "forestgreen")
names(cols) <- unique(pom_farm$data)

#plot 
pdf("Pomacentridae_Farming_ACE.pdf", height = 11, width = 11)

plotTree(pom_farm$phy, type = "fan", fsize = 0.3, type = "i", offset = 5, tips = seq(5, 349, 1))

nodelabels(node=1:pom_farm$phy$Nnode+Ntip(pom_farm$phy), pie=fit_farm_m2$lik.anc,piecol=cols[colnames(fit_farm_m2$lik.anc)],cex=0.3)

tiplabels(pie=to.matrix(pom_farm$data,unique(pom_farm$data)),piecol=cols[unique(pom_farm$data)], cex = 0.2)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(pom_farm$phy)),fsize=0.8)

dev.off()
#plot transitions
ftmk_2_farm <- fitMk(tree = pom_farm$phy, x = pom_farm$data, model = "ARD")
plot(ftmk_2_farm)






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

#set up transition matrices that corrrespond to models
#model 1: all rates equal
fit_size3_m1 <- ace(pom_size3$data, pom_size3$phy, model = "ER", type = "discrete")
fd_1_size3 <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = "ER")

#model 2: all rates different
#plot transitions
#(diag, ml, sl, lm, diag, sm, ls, ms, diag)
m2_size_small <- matrix(c(0, 1, 0, 2, 0, 0, 3, 4, 0), 3)
m2_size_large <- matrix(c(0, 1, 2, 0, 0, 3, 0, 4, 0), 3)
m2_size_constrained <- matrix(c(0, 1, 0, 2, 0, 3, 0, 4, 0), 3)

fit_size3_m2 <- ace(pom_size3$data, pom_size3$phy, model = "ARD", type = "discrete")
fd_2_size3 <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = "ARD")

#small absorbing
fit_size3_m2_small <- ace(pom_size3$data, pom_size3$phy, model = m2_size_small, type = "discrete")
fd_2_size3_small <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = m2_size_small)

#large absorbing
fit_size3_m2_large <- ace(pom_size3$data, pom_size3$phy, model = m2_size_large, type = "discrete")
fd_2_size3_large <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = m2_size_large)

#no S-L transitions
fit_size3_constrianed <- ace(pom_size3$data, pom_size3$phy, model = m2_size_constrained , type = "discrete")
fd_2_size3_constrained <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = m2_size_constrained )



#model 3: symmetrical (1 -> 2 = 2 -> 1)
fit_size3_m3 <- ace(pom_size3$data, pom_size3$phy, model = "SYM", type = "discrete")
fd_3_size3 <- fitDiscrete(pom_size3$phy, dat = pom_size3$data, model = "SYM")

#compare AIC scores 
fd_1_size3[[4]]$aic
fd_2_size3[[4]]$aic
fd_3_size3[[4]]$aic
fd_2_size3_small[[4]]$aic
fd_2_size3_large[[4]]$aic
fd_2_size3_constrained[[4]]$aic
#set colors for plotting
cols <- fish(n = 3, option = "Hypsypops_rubicundus", end = 0.9)
names(cols) <- c("S", "M", "L")

setwd("../size/")

#plot 
pdf("Pomacentridae_Size3_ACE_constrained.pdf", height = 11, width = 11)

plotTree(pom_size3$phy, type = "fan", fsize = 0.3, type = "i", offset = 5, tips = seq(5, 349, 1))

nodelabels(node=1:pom_size3$phy$Nnode+Ntip(pom_size3$phy), pie=fit_size3_constrianed$lik.anc,piecol=cols[colnames(fit_size3_constrianed$lik.anc)],cex=0.3)

tiplabels(pie=to.matrix(pom_size3$data,unique(pom_size3$data)),piecol=cols[unique(pom_size3$data)], cex = 0.2)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(pom_size3$phy)),fsize=0.8)

dev.off()

ftmk_3_size3 <- fitMk(tree = pom_size3$phy, x = pom_size3$data, model = "ARD")
plot(ftmk_3_size3)


#plotTree(pom_diet$phy, type = "fan", fsize = 0.3, type = "i", label.offset = 5)

# nodelabels(node=1:tree$Nnode+Ntip(tree), pie=fit_diet_m2$lik.anc,piecol=cols,cex=0.3)

# col <- c("red", "blue", "white", "green")

#tiplabels(pie=to.matrix(pom_diet$data,sort(unique(pom_diet$data))),piecol=cols, cex = 0.3)
#tiplabels(pie=to.matrix(t,sort(unique(t))),piecol=col, cex = 0.3)

#tiplabels(pie=pie, piecol=cols, cex = 0.2)


#match to phylo
pom_size3 <- match.phylo.data(phy, states_size3)
pom_size3$data <- gsub("S", 1, pom_size3$data)
pom_size3$data <- gsub("M", 2, pom_size3$data)
pom_size3$data <- gsub("L", 3, pom_size3$data)

pom_size3$data <- as.numeric(pom_size3$data)
names(pom_size3$data) <- pom_size3$phy$tip.label



#alternative liklihood method
lik <- make.mkn(tree = pom_size3$phy, states = pom_size3$data, k = 3)
fitted.ard <- find.mle(func = lik, x.init = fit_size3_m2$rates)
asr.ard <- asr.marginal(lik, pars = fitted.ard$par)
rownames(asr.ard) <- c("S", "M", "L")




plotTree(pom_size3$phy, type = "fan", fsize = 0.3, type = "i", offset = 5, tips = seq(5, 349, 1))
# for(i in 1:4) {
#   arc.cladelabels(text=labels[i],node=nodes[i], lab.offset = 1.42, ln.offset = 1.38)}

nodelabels(node=1:pom_size3$phy$Nnode+Ntip(pom_size3$phy), pie=t(asr.ard),piecol=cols[colnames(t(asr.ard))],cex=0.3)

tiplabels(pie=to.matrix(pom_size3$data,unique(pom_size3$data)),piecol=cols[unique(pom_size3$data)], cex = 0.2)

add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1], y=-max(nodeHeights(pom_size3$phy)),fsize=0.8)


#plot transition plots in same fig
pdf("Pomacentridae_ACE_Transition_Plots.pdf", height = 8, width = 8)

layout(matrix(c(1:4), 2,2, byrow = T))

plot(ftmk_2_diet)
title("Ecotype Transition Rates")
plot(ftmk_2_farm)
title("Farming Transition Rates")
plot(ftmk_3_size3)
title("Size (3 State) Transition Rates")
#plot(ftmk_2_size5)
#title("Size (5 State) Transition Rates")

dev.off()


## mirror tree with ecotype and farming
pdf("Pomacentridae_Ecotype_Farming_MirrorPhy_newcols.pdf", height = 15, width = 10)


layout(matrix(1:3,1,3),widths=c(0.4,0.07,0.4))
plotTree(pom_diet$phy, ftype = "off")

cols <-c("maroon","darkorange","deepskyblue2")
names(cols) <- c("B", "I", "P")


nodelabels(node=1:pom_diet$phy$Nnode+Ntip(pom_diet$phy), pie=fit_diet_m2$lik.anc,piecol=cols[colnames(fit_diet_m2$lik.anc)],cex=0.4)

tiplabels(pie=to.matrix(pom_diet$data,unique(pom_diet$data)),piecol=cols[unique(pom_diet$data)], cex = 0.2)

#add.simmap.legend(colors=cols,prompt=FALSE, x = 5,fsize=0.8)

plot.new()
plot.window(xlim=c(-0.1,0.1),ylim=c(1, length(pom_farm$phy$tip.label)))
par(cex=0.26)

is_tip <- pom_farm$phy$edge[,2] <= length(pom_farm$phy$tip.label)

ordered_tips <- pom_farm$phy$edge[is_tip, 2]

tip_labels <- gsub("_", " ", pom_farm$phy$tip.label[ordered_tips])

graphics::text(rep(0,length(pom_farm$phy$tip.label)), 1:length(pom_farm$phy$tip.label),tip_labels, font = 3)

#set colors for plotting
cols <- c("gold", "forestgreen")
names(cols) <- unique(pom_farm$data)

plotTree(pom_farm$phy, direction = "leftwards",ftype = "off")


nodelabels(node=1:pom_farm$phy$Nnode+Ntip(pom_farm$phy), pie=fit_farm_m2$lik.anc,piecol=cols[colnames(fit_farm_m2$lik.anc)],cex=0.4)

tiplabels(pie=to.matrix(pom_farm$data,unique(pom_farm$data)),piecol=cols[unique(pom_farm$data)], cex = 0.2)

#add.simmap.legend(colors=cols,prompt=FALSE,x = 40,fsize=1)

dev.off()





