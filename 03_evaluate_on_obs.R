#Run these if you haven't ran any previous scripts
#directory where this script is located
basedir <- "/net/xenon/climphys_backedup/maegli/Precipitation/Scripts/pub/"
setwd(basedir)
source("./functions.R")
source("./00_global_vars.R")


#load obs products
load("./data/DJF/GPCP_5d00_DJF_anom_XAX.RData")
load("./data/DJF/GHCN_pr_5d00_DJF_anom_XAX.RData")
load("./data/DJF/GPCC_pr_5d00_DJF_anom_XAX.RData")
load("./data/DJF/CR20_psl_5d00_DJF_anom_XAX.RData")
load("./data/DJF/PREC_pr_5d00_DJF_anom_XAX.RData")
load("./data/DJF/Smith_XAX_DJF.RData")
load("./data/DJF/ERA5_pr_5d00_DJF_anom_XAX.RData")

#load pretrained ridge models
load("./models/mod_DJF_prGHCN1950_ZM_single.RData")
mod_prGHCN1950_sg <- mod[[1]]

load("./models/mod_DJF_prpslGHCN1950_ZM_single.RData")
mod_prpslGHCN1950_sg <- mod[[1]]

load("./models/mod_DJF_prGPCC1950_ZM_single.RData")
mod_prGPCC1950 <- mod[[1]]

load("./models/mod_DJF_prpslGPCC1950_ZM_single.RData")
mod_prpslGPCC1950<- mod[[1]]


#define evaluation periods
eval_period <- 1979:2008

GHCN_ix <- which(GHCN_pr_5d00_DJF_anom_XAX$M$year %in% eval_period)
CR20_ix <- which(CR20_psl_5d00_DJF_anom_XAX$M$year %in% eval_period)
GPCP_ix <- which(GPCP_5d00_DJF_anom_XAX$M$year %in% eval_period)
GPCC_ix <- which(GPCC_pr_5d00_DJF_anom_XAX$M$year %in% eval_period)
PREC_ix <- which(PREC_pr_5d00_DJF_anom_XAX$M$year %in% eval_period)
CRU_ix <-  which(CRU_tas_5d00_DJF_anom_XAX$M$year %in% eval_period)
Smith_ix <- which(Smith_XAX_DJF$M$year %in% eval_period)
ERA5_ix <- which(ERA5_pr_5d00_DJF_anom_XAX$M$year %in% eval_period)

lambda <- "1se"

temp <- list()
for (lat in 1:36){
  
  GHCN_prpsl1950 <- get_pred(GHCN_pr_5d00_DJF_anom_XAX$X[GHCN_ix, -GHCN_na_1950_DJF],
                             CR20_psl_5d00_DJF_anom_XAX$X[CR20_ix, ],
                             mod_prpslGHCN1950_sg[[lat]],
                             lambda,
                             multi = T) %>% c()
  
  GHCN_pr1950 <- get_pred(GHCN_pr_5d00_DJF_anom_XAX$X[GHCN_ix, -GHCN_na_1950_DJF],
                          CR20_psl_5d00_DJF_anom_XAX$X[CR20_ix, ], 
                          mod_prGHCN1950_sg[[lat]], 
                          lambda, 
                          multi = F) %>% c()
  
  
  GPCC_pr1950 <- get_pred(GPCC_pr_5d00_DJF_anom_XAX$X[GPCC_ix, -GPCC_na_1950_DJF],
                            CR20_psl_5d00_DJF_anom_XAX$X[CR20_ix,], 
                            mod_prGPCC1950[[lat]], 
                            lambda, 
                            multi = F) %>% c()
  
  GPCC_prpsl1950 <- get_pred(GPCC_pr_5d00_DJF_anom_XAX$X[GPCC_ix, -GPCC_na_1950_DJF],
                               CR20_psl_5d00_DJF_anom_XAX$X[CR20_ix,], 
                               mod_prpslGPCC1950[[lat]], 
                               lambda, 
                               multi = T) %>% c()
  
  
  Smith <- Smith_XAX_DJF$Y$ZM[lat,Smith_ix]
  PREC <- PREC_pr_5d00_DJF_anom_XAX$Y$ZM[lat,PREC_ix]
  GPCP <- GPCP_5d00_DJF_anom_XAX$Y$ZM[lat,GPCP_ix]
  ERA <- ERA5_pr_5d00_DJF_anom_XAX$Y$ZM[lat, ERA5_ix]
  GHCN_simple <- GHCN_pr_5d00_DJF_anom_XAX$Y$AZMP[lat, GHCN_ix]
  
  temp[[lat]] <-  tibble(
    GHCN_prpsl1950 = GHCN_prpsl1950,
    GHCN_pr1950 = GHCN_pr1950,
    GPCC_prpsl1950 = GPCC_prpsl1950,
    GPCC_pr1950 = GPCC_pr1950,
    GHCN_simple = GHCN_simple,
    Smith = Smith,
    PREC = PREC,
    ERA5 = ERA,
    GPCP = GPCP,
    lat = latitudes[lat],
    year = as.integer(GHCN_pr_5d00_DJF_anom_XAX$M$year[GHCN_ix])) %>%
    pivot_longer(-c(GPCP,year, lat), names_to = "type", values_to = "value") 
  
  
}

GHCN_pred <- bind_rows(temp)

GHCN_pred  %>% 
  filter(lat %in% c(42.5, -47.5, 2.5), type != "GHCN_pr1950", type != "GHCN_prtas1950") %>% 
  ggplot(aes(x = year, y = value, group = type, fill = type))+
  geom_hline(yintercept = 0) +
  geom_line(aes(color = type)) +
  ylab("dPr [mm/month]") +
  xlab("Year")  +
  facet_wrap(~lat,  scales = "free") +
  theme(legend.title = element_blank())


sum_GPCP_DJF <- GHCN_pred %>% filter(year >= 1979, year <= 2008) %>% 
  group_by(lat, type) %>% 
  summarise(cor = cor(GPCP, value),
            rmse = rmse(GPCP, value),
            lambda = lambda)



sum_GPCP_DJF  %>% 
  ggplot(aes(lat, cor, color = type)) +
  geom_line() +
  ylab("Correlation with GPCP") + 
  xlab("Lattitude [°N]") +
  scale_x_continuous(breaks = c(-45, 0, 45)) +
  #scale_color_manual(values = c("#f0ae14","#FD7C53", "#5396FD", "#424242")) +
  theme(legend.title = element_blank())

ggsave(file = "./figures/cor_with_GPCP_DJF.png", height = 4, width = 8)



#==========JJA======================


load("./data/JJA/GPCP_5d00_JJA_anom_XAX.RData")
load("./data/JJA/ERA5_slp_5d00_JJA_anom_XAX.RData")
load("./data/JJA/GHCN_pr_5d00_JJA_anom_XAX.RData")
load("./data/JJA/ERA5_pr_5d00_JJA_anom_XAX.RData")
load("./data/JJA/GPCC_pr_5d00_JJA_anom_XAX.RData")
load("./data/JJA/CR20_psl_5d00_JJA_anom_XAX.RData")
load("./data/JJA/PREC_pr_5d00_JJA_anom_XAX.RData")
load("./data/JJA/Smith_XAX_JJA.RData")



load("./models/mod_JJA_prGHCN1950_ZM_single.RData")
mod_prGHCN1950_JJA <- mod[[1]]

load("./models/mod_JJA_prGPCC1950_ZM_single.RData")
mod_prGPCC1950_JJA <- mod[[1]]

load("./models/mod_JJA_prpslGHCN1950_ZM_single.RData")
mod_prpslGHCN1950_JJA <- mod[[1]]

load("./models/mod_JJA_prpslGPCC1950_ZM_single.RData")
mod_prpslGPCC1950_JJA <- mod[[1]]


eval_period <- 1979:2008

GHCN_ix <- which(GHCN_pr_5d00_JJA_anom_XAX$M$year %in% eval_period)
CR20_ix <- which(CR20_psl_5d00_JJA_anom_XAX$M$year %in% eval_period)
GPCP_ix <- which(GPCP_5d00_JJA_anom_XAX$M$year %in% eval_period)
GPCC_ix <- which(GPCC_pr_5d00_JJA_anom_XAX$M$year %in% eval_period)
PREC_ix <- which(PREC_pr_5d00_JJA_anom_XAX$M$year %in% eval_period)
Smith_ix <- which(Smith_XAX_JJA$M$year %in% eval_period)
ERA_ix <- which(ERA5_pr_5d00_JJA_anom_XAX$M$year %in% eval_period)


lambda <- "1se"

temp <- list()
for (lat in 1:36){
  
  GHCN_pr1950 <- get_pred(GHCN_pr_5d00_JJA_anom_XAX$X[GHCN_ix, -GHCN_na_1950_JJA],
                          CR20_psl_5d00_JJA_anom_XAX$X[CR20_ix, ], 
                          mod_prGHCN1950_JJA[[lat]], 
                          lambda, 
                          multi = F) %>% c()
  
  GPCC_pr1950 <- get_pred(GPCC_pr_5d00_JJA_anom_XAX$X[GPCC_ix, -GPCC_na_1950_JJA],
                          CR20_psl_5d00_JJA_anom_XAX$X[CR20_ix, ], 
                          mod_prGPCC1950_JJA[[lat]], 
                          lambda, 
                          multi = F) %>% c()
  
  GHCN_prpsl1950 <- get_pred(GHCN_pr_5d00_JJA_anom_XAX$X[GHCN_ix ,-GHCN_na_1950_JJA],
                             CR20_psl_5d00_JJA_anom_XAX$X[CR20_ix,], 
                             mod_prpslGHCN1950_JJA[[lat]], 
                             lambda, 
                             multi = T) %>% c()
  
  GPCC_prpsl1950 <- get_pred(GPCC_pr_5d00_JJA_anom_XAX$X[GPCC_ix ,-GPCC_na_1950_JJA],
                             CR20_psl_5d00_JJA_anom_XAX$X[CR20_ix,], 
                             mod_prpslGPCC1950_JJA[[lat]], 
                             lambda, 
                             multi = T) %>% c()
  

  #=======================================
  
  Smith <- Smith_XAX_JJA$Y$ZM[lat, Smith_ix]
  PREC <- PREC_pr_5d00_JJA_anom_XAX$Y$ZM[lat, PREC_ix]
  simple <- GHCN_pr_5d00_JJA_anom_XAX$Y$AZMP[lat, GHCN_ix]
  ERA5 <- ERA5_pr_5d00_JJA_anom_XAX$Y$ZM[lat, ERA_ix]
  GPCP <-  GPCP_5d00_JJA_anom_XAX$Y$ZM[lat, GPCP_ix]
  
  temp[[lat]] <-  tibble(
    GHCN_prpsl1950 = GHCN_prpsl1950,
    GHCN_pr1950 = GHCN_pr1950,
    GPCC_prpsl1950 = GPCC_prpsl1950,
    GPCC_pr1950 = GPCC_pr1950,
    GHCN_simple = as.numeric(simple),
    Smith= Smith,
    PREC = PREC,
    ERA5 = ERA5,
    GPCP = GPCP,
    lat = latitudes[lat],
    year = as.integer(GHCN_pr_5d00_JJA_anom_XAX$M$year[GHCN_ix])) %>%
    pivot_longer(-c(year, lat, GPCP), names_to = "type", values_to = "value")

}

GHCN_pred_JJA <- bind_rows(temp)


GHCN_pred_JJA %>% 
ggplot(aes(x = year, y = value, group = type, fill = type))+
  geom_hline(yintercept = 0) +
  geom_line(aes(color = type)) +
  ylab("dPr [mm/month]") +
  xlab("Year")  +
  facet_wrap(~lat,  scales = "free") +
  theme(legend.title = element_blank())




sum_GPCP_JJA <- GHCN_pred_JJA %>% 
  filter(year >= 1979, year <= 2008) %>% 
  group_by(lat, type) %>% 
  summarise(cor = cor(value, GPCP),
            rmse = rmse(value, GPCP))



sum_GPCP_recon <- tibble(sum_GPCP_JJA, seas = "JJA") %>% bind_rows(tibble(sum_GPCP_DJF, seas = "DJF")) %>% 
  filter(stringr::str_detect(type, "_")) %>% 
  #filter(!stringr::str_detect(type, "simple")) %>% 
  separate(type, c("obs", "mask"), sep = "_")  %>%
  filter(!stringr::str_detect(mask, "GHCN")) %>% 
  mutate(mask = case_when(mask == "pr1950" ~ "Pr",
                          mask == "prpsl1950" ~ "Pr + SLP",
                          mask == "simple" ~ "Simple"))



sum_GPCP_ref <- tibble(sum_GPCP_JJA, seas = "JJA") %>% bind_rows(tibble(sum_GPCP_DJF, seas = "DJF")) %>% 
  filter(type %in%  c("Smith", "PREC", "ERA5")) %>% 
  pivot_longer(c(cor, rmse), names_to = "loss")

ggplot() +
  geom_line(data = sum_GPCP_recon, aes(lat, cor, color = obs), size= 0.8) +
  geom_line(data = sum_GPCP_ref, aes(lat, value, linetype = type), color = "grey50") +
  facet_grid(mask ~ seas, scales = "free") +
  scale_x_continuous(breaks = c(-45, 0, 45)) +
  ylab("Correlation to GPCP") + 
  xlab("Lattitude [°N]") +
  scale_color_manual(values = c("forestgreen", "#f0bb35"), labels = c("GHCN", "GPCC")) +
  pub_theme
  
#supplement figure
sum_GPCP_recon %>% filter(obs == "GHCN") %>% 
  pivot_longer(c(cor, rmse), names_to = "loss") %>% 
ggplot() +
  geom_line(aes(lat, value, color = mask), size= 0.8) +
  geom_line(data = filter(sum_GPCP_ref, type == "ERA5") , aes(lat, value, linetype = type), color = "grey50") +
  facet_grid(loss ~ seas, scales = "free") +
  scale_x_continuous(breaks = c(-45, 0, 45)) +
  ylab("Reconstruction success compared to GPCP") + 
  xlab("Lattitude [°N]") +
  scale_color_manual(values = c("forestgreen", "#f0bb35", "springgreen3"), labels = c("GHCN", "GPCC", "GHCN Simple")) +
  pub_theme

ggsave(file = "./figures/ridge_vs_simple_GHCN.png", height = 6, width = 9)





plot(mod_prGHCN1950_sg[[16]])
