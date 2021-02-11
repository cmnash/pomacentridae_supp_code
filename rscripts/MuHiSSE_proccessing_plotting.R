## this script processes the MuHiSSE restart data for ML, plotting, generating support regions, and parameter best fits

library(tidyverse)
library(hisse)
library(RColorBrewer)
library(utilhisse)
library(viridis)
library(data.table)
options(scipen = 2)

#load adapted functions from utilhisse
source("~/Box/Pom_Phylo_Paper/Supplementary Code/rscripts/muhisse_divplot.R")

#get names of restart files by model
setwd("~/Box/Pom_Phylo_Paper/Supplementary Code/data/MuHiSSE/100restart_output/")
m1.files <- list.files(pattern = "m1_", full.names = TRUE)
m2.files <- list.files(pattern = "m2_", full.names = TRUE)
m3.files <- list.files(pattern = "m3_", full.names = TRUE)
m4.files <- list.files(pattern = "m4_", full.names = TRUE)
m5.files <- list.files(pattern = "m5_", full.names = TRUE)
m6.files <- list.files(pattern = "m6_", full.names = TRUE)
m7.files <- list.files(pattern = "m7_", full.names = TRUE)
m8.files <- list.files(pattern = "m8_", full.names = TRUE)
m9.files <- list.files(pattern = "m9_", full.names = TRUE)
m10.files <- list.files(pattern = "m10_", full.names = TRUE)
m11.files <- list.files(pattern = "m11_", full.names = TRUE)
m12.files <- list.files(pattern = "m12_", full.names = TRUE)
m13.files <- list.files(pattern = "m13_", full.names = TRUE)

num <- list()
for (i in 1:13){
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

m9.fits <- purrr::map(m9.files, ~ get(load(.)))
names(m9.fits) <- num[[9]]

m10.fits <- purrr::map(m10.files, ~ get(load(.)))
names(m10.fits) <- num[[10]]

m11.fits <- purrr::map(m11.files, ~ get(load(.)))
names(m11.fits) <- num[[11]]

m12.fits <- purrr::map(m12.files, ~ get(load(.)))
names(m12.fits) <- num[[12]]

m13.fits <- purrr::map(m13.files, ~ get(load(.)))
names(m13.fits) <- num[[13]]



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
for (i in 1:13){
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
  ID == 8 ~ "m8",
  ID == 9 ~ "m9",
  ID == 10 ~ "m10",
  ID == 11 ~ "m11",
  ID == 12 ~ "m12",
  ID == 13 ~ "m13"
))


mod.best %>%
  group_by(model) %>%
  arrange(desc(lik)) %>% 
  top_n(n = 1, wt = lik) %>% 
  ungroup() %>% 
  mutate(daic=aic-min(aic)) %>% 
  mutate(aicw=exp(-0.5*daic)/sum(exp(-0.5*daic))) %>% 
  arrange(aic) %>% 
  data.frame

#extract restart with the best liklihood for each model (based on table above)
muhisse_m1_best <- m1.fits$`19`
muhisse_m2_best <- m2.fits$`16`
muhisse_m3_best <- m3.fits$`71`
muhisse_m4_best <- m4.fits$`26`
muhisse_m5_best <- m5.fits$`66`
muhisse_m6_best <- m6.fits$`50`
muhisse_m7_best <- m7.fits$`69`
muhisse_m8_best <- m8.fits$`1`
muhisse_m9_best <- m9.fits$`94`
muhisse_m10_best <- m10.fits$`8`
muhisse_m11_best <- m11.fits$`56`
muhisse_m12_best <- m12.fits$`12`
muhisse_m13_best <- m13.fits$`96`

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
    muhisse_m8_best,
    muhisse_m9_best,
    muhisse_m10_best,
    muhisse_m11_best,
    muhisse_m12_best,
    muhisse_m13_best
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


#load("../All_recon_Oct30.RData")
#prepare MuHiSSE marginal reconstruction for plotting
m_proc <- m_process_recon(res)

#match ecotype colors
cols_eco <- c(viridis(3),"white")
cols_eco <- cols_eco[c(3,1,2,4)]


##ridgelines of estimates
pdf("Pomacentridae_Muhisse_ridgelines_wcorrectcols.pdf", height = 8, width = 10)
m_ridgelines_adjust(
  processed_recon = m_proc,
  parameter = 'net.div',
  states_names = c('00','Pelagic','Benthic','Intermediate'),
  state_order <- c('00','Benthic','Intermediate','Pelagic'),
  plot_as_waiting_time = F,
  line_colors = cols_eco,
  fill_colors = cols_eco
) + theme_classic()
dev.off()

##ancestral state estimation under muhisse
# col <- viridis(3)
# pdf("muhisse_trait.pdf", height = 8, width = 8)
# m_trait_recon(
#   processed_recon = m_proc,
#   states_of_first_character = c('-','Benthic'),
#   states_of_second_character = c('-','Pelagic'),
#   cutoff = as.numeric(c('0.5','0.5')),
#   colors = c(col[3], col[1], col[2], "white"),
#   show_tip_labels = TRUE,
#   tree_layout = 'fan',
#   tree_direction = 'up',
#   time_axis_ticks = 10,
#   open_angle = 10)
# dev.off()


##diversification plot

#remove underscores from phylo
m_proc$tree_data@phylo$tip.label <- gsub("_", " ",m_proc$tree_data@phylo$tip.label)

#plot MuHiSSE model-averaged marginal ancestral state estimatoin for diversification rates
mrc <- m_rate_recon_direction(
  processed_recon = m_proc,
  parameter = 'net.div',
  show_tip_labels = T,
  discrete = T,
  breaks = seq(0,0.3,0.01),
  plot_as_waiting_time = F,
  tree_layout = 'fan',
  time_axis_ticks = 6,
  colors = viridis(18),
  open_angle = 5,
  size = 0.8,
  tip_size = 1, 
  right = T
)

pdf("Pomacentridae_MuHiSSE_Diversification_Rate_Phy.pdf", height = 10, width = 10)

rotate_tree(mrc, 5)

dev.off()



#Get model averages rates
mar <- GetModelAveRates(res, type = "both")

#infer support region
sr <- SupportRegionMuHiSSE(muhisse.obj = m4.fits$`26`, n.points = 1000, desired.delta = 2, min.number.points = 10, verbose = T)
class(sr)
sr[[1]][,1:9]
sr[[1]][,c(1, 50:57)]

save(sr, file = "../././Data/MuHISSE_Results/muhisse_use/Pom_supportregion_m4_oct30.RData")


#get model averaged rates
mod_avgrates <- GetModelAveRates(res, x = "tips")

#transition and diversification rates from MuHiSSE best fit (m4)
m4_rates <- m_collect_rates(m4.fits$`26`, hidden_traits = T, character_states = c("N", "Pelagic", "Benthic", "Intermediate"))

#confidence intervals
con_int <-sr$ci
con_int <- con_int[,con_int[4,] != 0][,-1]
sol <- muhisse_m4_best$solution
sol <- sol[sol != 0]

#clean data for table for best fit
sol_round <- paste(round(sol, 4), " [", round(con_int[1,],4), "-", round(con_int[4,],4), "]", sep = "" )

solu <- tibble(Parameter = colnames(con_int[,con_int[4,] != 0]), Solution = sol_round)

write_csv(solu, "Pomacentridae_MuHiSSE_m4best_SR.csv")
