#beta plots

#Run these if you haven't ran any previous scripts
#directory where this script is located
#basedir <- "/net/xenon/climphys_backedup/maegli/Precipitation/Scripts/pub/"
setwd(basedir)
source("./functions.R")
source("./00_global_vars.R")

#==================load and prepare model coeffiecients
ridge_file_DJF_prGHCN <- "./models/mod_DJF_prGHCN1950_ZM_single.RData"
load(ridge_file_DJF_prGHCN)
mod_DJF_prGHCN <- mod[[1]]

ridge_file_DJF_prpslGHCN <- "./models/mod_DJF_prpslGHCN1950_ZM_single.RData"
load(ridge_file_DJF_prpslGHCN)
mod_prpsl1950_single <- mod[[1]]


sds_pr <- LENS_pr_DJF_5d00_XAX_norm$X[,-GHCN_na_1950_DJF] %>% matrixStats::colSds()
sds_psl <- LENS_psl_DJF_5d00_XAX_norm$X %>% matrixStats::colSds()


spdf1_list <- list()
spdf2pr_list <- list()
spdf2psl_list <- list()

for (lat in 1:36){
  
  breakprpsl <- 2592 - length(GHCN_na_1950_DJF)
  n_betas_prpsl <- mod_prpsl1950_single[[lat]]$glmnet.fit$beta %>% nrow()
  
  spdf1_list[[lat]] <- map_betas_norm(mask = GHCN_na_1950_DJF, mod = mod_DJF_prGHCN[[lat]], sds_pr) %>% tibble(lat = latitudes[lat])
  spdf2pr_list[[lat]] <- map_betas_norm(GHCN_na_1950_DJF, mod_prpsl1950_single[[lat]], sds_pr, break_sec = 1:breakprpsl)  %>% tibble(lat = latitudes[lat])
  spdf2psl_list[[lat]] <- map_betas_norm(mask = 0,mod_prpsl1950_single[[lat]], sds_psl, break_sec = ((breakprpsl+1):n_betas_prpsl)) %>% tibble(lat = latitudes[lat])
  
}


prpsl_pr <- bind_rows(spdf2pr_list) %>% mutate(type = "Pr+SLP model, Pr coefficients")
only_pr <- bind_rows(spdf1_list) %>% mutate(type = "Pr model, Pr coefficients")


spdf_plot <- bind_rows(spdf2psl_list) %>% mutate(type = "Pr+SLP model, SLP coefficients") %>% 
  bind_rows(prpsl_pr, only_pr) 


library(sf)
land <- rnaturalearth::ne_download(
  scale = 110,
  type = "land",
  category = "physical",
  returnclass = "sf")


#Pr only coefficients
spdf_plot %>% filter(lat %in% c(47.5) & type ==  "Pr model, Pr coefficients") %>% 
  ggplot() + 
  geom_tile(aes(fill= Betas, x = x, y = y)) + 
  scale_fill_gradient2(na.value="transparent",
                       mid = "white",
                       low = "#061c8c",
                       high = "#b50404") +
  geom_hline(aes(yintercept = lat), linetype='dotted') + 
  geom_sf(data = land, fill = NA, color = "grey50") +
  facet_grid(lat ~type) +
  labs(fill = "Regression Coefficients") +
  theme(axis.ticks = element_blank(),
        axis.title=element_blank(),
        legend.position = "bottom", 
        legend.key.width = unit(1.8, "cm"),
        legend.key.height = unit(3, "mm"),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect("#f5f2f2")) +
  ylab(element_blank()) + xlab(element_blank()) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))


#Pr + SLP coefficients
spdf_plot %>% filter(lat %in% c(47.5) & stringr::str_detect(type, "SLP")) %>% 
  ggplot() + 
  geom_tile(aes(fill= Betas, x = x, y = y)) + 
  scale_fill_gradient2(na.value="transparent",
                       mid = "white",
                       low = "#061c8c",
                       high = "#b50404") +
  geom_hline(aes(yintercept = lat), linetype='dotted') + 
  geom_sf(data = land, fill = NA, color = "grey50") +
  facet_grid(lat ~ type) +
  labs(fill = "Regression Coefficients") +
  theme(axis.ticks = element_blank(),
        axis.title=element_blank(),
        legend.position = "bottom", 
        legend.key.width = unit(1.8, "cm"),
        legend.key.height = unit(3, "mm"),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect("#f5f2f2")) +
  ylab(element_blank()) + xlab(element_blank()) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))


