## this script processes the MuHiSSE restart data for ML, plotting, generating support regions, and parameter best fits

library(tidyverse)
library(hisse)
library(RColorBrewer)
library(gghisse)
library(viridis)
library(data.table)
options(scipen = 2)

#load adapted functions from utilhisse
source("~/rscripts/muhisse_divplot.R")


#get names of restart files by model
setwd("~/data/MuHiSSE/100restart_output/")

m1.files <- list.files(pattern = "m1_", full.names = TRUE)
m2.files <- list.files(pattern = "m2_", full.names = TRUE)
m3.files <- list.files(pattern = "m3_", full.names = TRUE)
m4.files <- list.files(pattern = "m4_", full.names = TRUE)
m5.files <- list.files(pattern = "m5_", full.names = TRUE)
m6.files <- list.files(pattern = "m6_", full.names = TRUE)
m7.files <- list.files(pattern = "m7_", full.names = TRUE)
m8.files <- list.files(pattern = "m8_", full.names = TRUE)


num <- list()
for (i in 1:8){
  num[[i]] <- gsub(pattern = ".Rsave",replacement = "", x = tstrsplit(eval(parse(text = paste("m", i, ".files", sep = ""))), split = "_")[[4]])
}

#load in model fits by model into list
m1.fits <- purrr::map(m1.files, ~ get(load(.)))
names(m1.fits) <- num[[1]]

m2.fits <- purrr::map(m2.files, ~ get(load(.)))
names(m2.fits) <- num[[2]]

m3.fits <- purrr::map(m3.files, ~ get(load(.)))
names(m3.fits) <- num[[3]]

m4.fits <- purrr::map(m4.files, ~ get(load(.)))
names(m4.fits) <- num[[4]]

m5.fits <- purrr::map(m5.files, ~ get(load(.)))
names(m5.fits) <- num[[5]]

m6.fits <- purrr::map(m6.files, ~ get(load(.)))
names(m6.fits) <- num[[6]]

m7.fits <- purrr::map(m7.files, ~ get(load(.)))
names(m7.fits) <- num[[7]]

m8.fits <- purrr::map(m8.files, ~ get(load(.)))
names(m8.fits) <- num[[8]]




#function to get the likilihood and aic scores for all the fits
get_results <- function(List_names, List_fits){
  res <- tibble(
    num = gsub(
      pattern = ".Rsave",
      replacement = "",
      tstrsplit(List_names, split = "_")[[4]]
    ),
    lik = purrr::map(List_fits, ~ .["loglik"]) %>% unlist,
    aic = purrr::map(List_fits, ~ .["AICc"]) %>% unlist
    
  )
  res <- res %>% arrange(desc(lik)) %>% mutate(daic = aic - min(aic)) %>% mutate(aicw = exp(-0.5*daic)/sum(exp(-0.5*daic)))
  return(res)
}

#apply to all the models
all_fits <- list()
for (i in 1:8){
  all_fits[[i]] <- get_results(eval(parse(text = paste("m", i, ".files", sep = ""))), eval(parse(text = paste("m", i, ".fits", sep = ""))))
}

#sort based on liklihood and aic
mod.best <- bind_rows(all_fits, .id = "ID") %>% mutate(model = case_when(
  ID == 1 ~ "m1",
  ID == 2 ~ "m2",
  ID == 3 ~ "m3",
  ID == 4 ~ "m4",
  ID == 5 ~ "m5",
  ID == 6 ~ "m6",
  ID == 7 ~ "m7",
  ID == 8 ~ "m8"
))


m <- mod.best %>%
  group_by(model) %>%
  arrange(desc(lik)) %>% 
  top_n(n = 1, wt = lik) %>% 
  ungroup() %>% 
  mutate(daic=aic-min(aic)) %>% 
  mutate(aicw=exp(-0.5*daic)/sum(exp(-0.5*daic))) %>% 
  arrange(aic) %>% 
  data.frame

#extract restart with the best liklihood for each model (based on table above)
muhisse_m1_best <- m1.fits$`38`
muhisse_m2_best <- m2.fits$`18`
muhisse_m3_best <- m3.fits$`27`
muhisse_m4_best <- m4.fits$`11`
muhisse_m5_best <- m5.fits$`18`
muhisse_m6_best <- m6.fits$`19`
muhisse_m7_best <- m7.fits$`54`
muhisse_m8_best <- m8.fits$`52`


#convert into a list
mod_fit <-
  list(
    muhisse_m1_best,
    muhisse_m2_best,
    muhisse_m3_best,
    muhisse_m4_best,
    muhisse_m5_best,
    muhisse_m6_best,
    muhisse_m7_best,
    muhisse_m8_best
  )


#estimate the likeliest states using marginal reconstruction algorithm and model averaging
## models are split by whether hidden states are used
res <- list()
for (i in 1:length(mod_fit)){
  mod <- mod_fit[[i]]
  
  if (i < 3){ #all non-hidden state models
  recon <-
    MarginReconMuHiSSE(
      phy = mod$phy,
      data = mod$data,
      f = mod$f,
      pars = mod$solution,
      hidden.states = 1,
      root.type = mod$root.type,
      root.p = mod$root.p,
      aic = mod$AICc,
      n.cores = 4
    )
  } else if (i > 2 & i < 8){ #models with 2 hidden states
    recon <-
      MarginReconMuHiSSE(
        phy = mod$phy,
        data = mod$data,
        f = mod$f,
        pars = mod$solution,
        hidden.states = 2,
        root.type = mod$root.type,
        root.p = mod$root.p,
        aic = mod$AICc,
        n.cores = 4
      )
  } else { #models with more than 2 hidden states (CID)
    recon <-
      MarginReconMuHiSSE(
        phy = mod$phy,
        data = mod$data,
        f = mod$f,
        pars = mod$solution,
        hidden.states =  i - 5,
        root.type = mod$root.type,
        root.p = mod$root.p,
        aic = mod$AICc,
        n.cores = 4
      )
  } 
  
  res[[i]] <- recon
}

#load("~/Pom_Muhisse_recon_ecotype_Jul21_21.RData")

#prepare MuHiSSE marginal reconstruction for plotting
m_proc_eco <- m_process_recon(res1[[7]])

#match ecotype colors
cols_eco <- c("maroon","darkorange","deepskyblue2","white")
cols_eco <- cols_eco[c(1,2,3,4)]


##ridgelines of estimates
pdf("Pomacentridae_Muhisse_ridgelines_netdiv_UPDATED.pdf", height = 8, width = 10)
m_ridgelines_adjust(
  processed_recon = m_proc_eco,
  parameter = 'net.div',
  states_names = c('00','Pelagic','Benthic','Intermediate'),
  state_order <- c('00','Benthic','Intermediate','Pelagic'),
  plot_as_waiting_time = F,
  line_colors = cols_eco,
  fill_colors = cols_eco,
  plot_min = T
) + theme_classic()
dev.off()

#remove underscores from phylo
m_proc_eco$tree_data@phylo$tip.label <- gsub("_", " ",m_proc_eco$tree_data@phylo$tip.label)


##ancestral state estimation under muhisse
col <- viridis(3)
cols <- viridis(3, alpha = 0.2)
mtc <- m_trait_recon_direction(
  processed_recon = m_proc_eco,
  states_of_first_character = c('-','Benthic'),
  states_of_second_character = c('-','Pelagic'),
  cutoff = as.numeric(c('0.2','0.2')),
  colors = c(col[3], cols[3], col[1], cols[1], col[2], col[2], "white"),
  show_tip_labels = TRUE,
  tree_layout = 'fan',
  tree_direction = 'up',
  time_axis_ticks = 6,
  open_angle = 5,
  tip_size = 1,
  size = 0.8,
  right = T, 
  offset = 1.2)
pdf("muhisse_trait_ecotype_UPDATED.pdf", height = 8, width = 8)
rotate_tree(mtc, 5)
dev.off()


##diversification plot


#plot MuHiSSE model-averaged marginal ancestral state estimatoin for diversification rates
mrc <- m_rate_recon_direction(
  processed_recon = m_proc_eco,
  parameter = 'net.div',
  show_tip_labels = T,
  discrete = T,
  breaks = seq(0,0.3,0.01),
  plot_as_waiting_time = F,
  tree_layout = 'fan',
  time_axis_ticks = 6,
  colors = viridis(12),
  open_angle = 5,
  size = 0.8,
  tip_size = 1, 
  right = T
)

pdf("Pomacentridae_MuHiSSE_Diversification_Rate_Phy_netdiv_UPDATED.pdf", height = 10, width = 10)

rotate_tree(mrc, 5)

dev.off()



#Get model averages rates
mar <- GetModelAveRates(res, type = "both")

#infer support region
sr_eco <- SupportRegionMuHiSSE(muhisse.obj = muhisse_m3_best, n.points = 1000, desired.delta = 2, min.number.points = 10, verbose = T)
class(sr)
sr[[1]][,1:9]
sr[[1]][,c(1, 50:57)]

save(sr, file = "~/Box/Pom_Phylo_Paper/Supplementary Code/data/MuHiSSE/Pom_supportregion_m3_aug2_ecotype.RData")


#get model averaged rates
mod_avgrates <- GetModelAveRates(res1, x = "tips")

#transition and diversification rates from MuHiSSE best fit (m4)
m3_rates <- m_collect_rates(muhisse_m3_best, hidden_traits = T, character_states = c("N", "Pelagic", "Benthic", "Intermediate"))

#confidence intervals
con_int <-sr_eco$ci
con_int <- con_int[,con_int[4,] != 0][,-1]
sol <- muhisse_m3_best$solution
sol <- sol[sol != 0]

#clean data for table for best fit
sol_round <- paste(round(sol, 4), " [", round(con_int[1,],4), "-", round(con_int[4,],4), "]", sep = "" )

solu <- tibble(Parameter = colnames(con_int[,con_int[4,] != 0]), Solution = sol_round)

write_csv(solu, "~/Pomacentridae_MuHiSSE_m3best_SR_ectotype.csv")



