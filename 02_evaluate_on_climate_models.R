#Run these if you haven't ran any previous scripts
#directory where this script is located
basedir <- "/net/xenon/climphys_backedup/maegli/Precipitation/Scripts/pub/"
setwd(basedir)
source("./functions.R")
source("./00_global_vars.R")

#load pretrained models
load("./models/model_as_truth/mod_DJF_prpsl1950_ZM.RData")
full_mod_prpsl1950 <- mod

load("./models/model_as_truth/mod_DJF_pr1950_ZM.RData")
full_mod_pr1950 <- mod

load("./models/model_as_truth/mod_DJF_prpsl1900_ZM.RData")
full_mod_prpsl1900 <- mod

load("./models/model_as_truth/mod_DJF_pr1900_ZM.RData")
full_mod_pr1900 <- mod

lambda <- "1se" #flag: either lambda = "min" or "1se" ||| 1se is more heavily regularized than min and was used for all the analysis presented
big_temp <- list()

for (i in 1:length(LENS_models)){
  #loop through climate model which is to be excluded
  mod_prpsl1950 <- full_mod_prpsl1950[[i]] 
  mod_pr1950 <- full_mod_pr1950[[i]]
  mod_prpsl1900 <- full_mod_prpsl1900[[i]]
  mod_pr1900 <- full_mod_pr1900[[i]]
  
  test_ix <- test_ix_list2005[[i]]
  temp <- list()
  
  for (lat in 1:36){
    
    pred_prpsl1950 <- get_pred(LENS_pr_DJF_5d00_XAX_norm$X[test_ix, -GHCN_na_1950_DJF],
                               LENS_psl_DJF_5d00_XAX_norm$X[test_ix, ], 
                               mod_prpsl1950[[lat]], 
                               lambda, 
                               multi = T) %>% c()
    
    pred_pr1950 <- get_pred(LENS_pr_DJF_5d00_XAX_norm$X[test_ix, -GHCN_na_1950_DJF],
                            LENS_psl_DJF_5d00_XAX_norm$X[test_ix, ],
                            mod_pr1950[[lat]], 
                            lambda, 
                            multi = F) %>% c()
    
    #pred_prpsl1900 <- get_pred(LENS_pr_DJF_5d00_XAX_norm$X[test_ix, -GHCN_na_1900_DJF],
    #                           LENS_psl_DJF_5d00_XAX_norm$X[test_ix, ], 
    #                           mod_prpsl1900[[lat]], 
    #                           lambda, 
     #                          multi = T) %>% c()
    
   # pred_pr1900 <- get_pred(LENS_pr_DJF_5d00_XAX_norm$X[test_ix, -GHCN_na_1900_DJF],
    #                        LENS_psl_DJF_5d00_XAX_norm$X[test_ix, ],
    #                        mod_pr1900[[lat]], 
     #                       lambda, 
      #                      multi = F) %>% c()
    
    
    naive_ix <- which(lonlatDT5$lat == latitudes[lat] & lonlatDT5$GHCN_DJF == T)
    len_naive_ix <- length(naive_ix)
    
    if(len_naive_ix == 0){ 
      pred_naive <- NA 
    }else if(len_naive_ix == 1){
      pred_naive <- LENS_pr_DJF_5d00_XAX_norm$X[test_ix, naive_ix]
    }else{
      pred_naive <- rowMeans(LENS_pr_DJF_5d00_XAX_norm$X[test_ix, naive_ix])
    }
    
    
    
    ZM <-  LENS_pr_DJF_5d00_XAX_norm$Y$ZM[test_ix, lat]
    temp[[lat]] <-  data.table(prpsl1950 = pred_prpsl1950,
                           pr1950 = pred_pr1950,
                           #prpsl1900 = pred_prpsl1900,
                          # pr1900 = pred_pr1900,
                           lat = latitudes[lat],
                           Simple = as.numeric(pred_naive),
                           mod =  LENS_pr_DJF_5d00_XAX_norm$M$mod[test_ix],
                           year = as.numeric(LENS_pr_DJF_5d00_XAX_norm$M$year[test_ix]),
                           ens = LENS_pr_DJF_5d00_XAX_norm$M$ens.mem[test_ix],
                           ZM = ZM)  %>% 
      melt(id.vars = c("year", "ens", "mod", "lat", "ZM"), variable.name = "type", value.name = "value")
  }
  big_temp[[i]] <- bind_rows(temp)
}

pr <- rbindlist(big_temp)


#time series plot
lat_sel <- c(57.5)
pr_ens1 <-  pr %>% filter(ens == "r10i1p1" & mod == "CESM1-CAM5" & lat == lat_sel & type %in% c("prpsl1950", "pr1950", "Simple"))
pr_ens1 %>% 
  ggplot(aes(x = year, y = value, group = type))+
  geom_line(aes(color = type), size = 0.9) +
  geom_line(aes(year, ZM, color = "ZM"), linetype = 2) +
  ylab("Precip anomalies [mm/month]") +
  xlab("Year")  +
  facet_grid(lat ~., scales = "free")  +
  scale_color_manual(values = c("#FD7C53", "#5396FD", "gainsboro", "grey30"), 
                     labels = c("Masked Pr", "Masked Pr + SLP","Simple", "Zonal Mean Pr (target)")) +
  pub_theme +
  theme(strip.text.y.right = element_text(angle = 0))

ggsave(file = "./figures/ZM_CESM1_ens1_DJF_legend.png", height = 2.5, width = 7)


#calculate correlations and summarise over over mods
#warnings are okay here, sometimes the ridge model only has an intercept and then its not possible to calculate correlations for that latititude
sum <- pr[,.(cor = cor(ZM, value), rmse = rmse(ZM,value)), by = .(ens,mod,lat,type)]

sumsum <-  sum[,.(cor = mean(cor), rmse = mean(rmse)), by = .(mod,lat,type)]

sumsumsum <- sumsum %>% melt(id.vars = c("mod","lat","type"), variable.name = "loss") %>% 
  .[,.(upper = max(value), lower = min(value), med = median(value)), by =.(lat,type,loss)]

#cells_per_lat <- spdf1_list[[lat]] %>% data.table() %>% .[!is.na(Betas), .N, by = y] %>% .[,N_norm := N/72]
sumsumsum[loss == "cor", loss := "Correlation"]
sumsumsum[loss == "rmse", loss := "RMSE"]

#correlation plot (Figure 2)
sumsumsum[type %in% c("Simple", "pr1950", "prpsl1950")] %>% 
  ggplot() +
  geom_line(aes(lat, med, group = type, color = type)) +
  geom_ribbon(aes(x = lat, ymax = upper, ymin = lower, group = type, fill = type), alpha = 0.3) +
  ylab("Reconstruction Success") +
  xlab("Latitude [°N]") +
  scale_x_continuous(breaks = c(-45, 0, 45)) +
  theme(legend.title = element_blank()) +
  scale_color_manual(labels = c("Masked Pr + SLP", "Masked Pr", "Simple"), values = c("#5396FD","#FD7C53", "grey50")) +
  scale_fill_manual( labels = c("Masked Pr + SLP", "Masked Pr", "Simple"), values = c("#5396FD","#FD7C53", "grey50")) +
  facet_wrap(~loss, scales = "free")+
  pub_theme


ggsave(file = "./figures/modelCV_ridgevssimple.png", height = 4, width = 8)
#===============================================
#================JJA============================
#===============================================

ridge_file_prpsl1950_JJA <- "./models/model_as_truth/mod_JJA_prpsl1950_ZM.RData"
ridge_file_pr1950_JJA <-    "./models/model_as_truth/mod_JJA_pr1950_ZM.RData"

load(ridge_file_prpsl1950_JJA)
full_mod_prpsl1950_JJA <- mod

load(ridge_file_pr1950_JJA)
full_mod_pr1950_JJA <- mod

lambda <- "1se"
big_temp_JJA <- list()


latitudes <- seq(-87.5, 87.5 ,5)

for (i in 1:length(LENS_models)){
  mod_prpsl1950_JJA <- full_mod_prpsl1950_JJA[[i]] 
  mod_pr1950_JJA <- full_mod_pr1950_JJA[[i]]
  
  test_ix <- test_ix_list2005[[i]]
  temp <- list()
  
  for (lat in 1:36){
    
    
    pred_prpsl1950 <- get_pred(LENS_pr_JJA_5d00_XAX_norm$X[test_ix, -GHCN_na_1950_JJA],
                               LENS_psl_JJA_5d00_XAX_norm$X[test_ix, ], 
                               mod_prpsl1950_JJA[[lat]], 
                               lambda, 
                               multi = T) %>% c()
    
    pred_pr1950 <- get_pred(LENS_pr_JJA_5d00_XAX_norm$X[test_ix, -GHCN_na_1950_JJA],
                            LENS_psl_JJA_5d00_XAX_norm$X[test_ix, ],
                            mod_pr1950_JJA[[lat]], 
                            lambda, 
                            multi = F) %>% c()
    
    naive_ix <- which(lonlatDT5$lat == latitudes[lat] & lonlatDT5$GHCN_JJA == T)
    len_naive_ix <- length(naive_ix)
    
    if(len_naive_ix == 0){ 
      pred_naive <- NA 
    }else if(len_naive_ix == 1){
      pred_naive <- LENS_pr_JJA_5d00_XAX_norm$X[test_ix, naive_ix]
    }else{
      pred_naive <- rowMeans(LENS_pr_JJA_5d00_XAX_norm$X[test_ix, naive_ix])
    }
    
    ZM <-  LENS_pr_JJA_5d00_XAX_norm$Y$ZM[test_ix, lat]
    temp[[lat]] <-  data.table(prpsl1950 = pred_prpsl1950,
                           pr1950 = pred_pr1950,
                           Simple = as.numeric(pred_naive),
                           lat = latitudes[lat],
                           mod =  LENS_pr_JJA_5d00_XAX_norm$M$mod[test_ix],
                           year = as.numeric(LENS_pr_JJA_5d00_XAX_norm$M$year[test_ix]),
                           ens = LENS_pr_JJA_5d00_XAX_norm$M$ens.mem[test_ix],
                           ZM = ZM) %>% 
      melt(id.vars = c("year", "ens", "mod", "lat", "ZM"), variable.name = "type", value.name = "value")
  }
  big_temp_JJA[[i]] <- bind_rows(temp)
}
ZM_JJA_mod <- rbindlist(big_temp_JJA)


ens1_sel_JJA <- ZM_JJA_mod %>%  filter(ens == "r10i1p1" & mod == "CESM1-CAM5" & lat == -2.5)


sum_ZMJJA <- ZM_JJA_mod[, .(cor = cor(ZM, value), rmse = rmse(ZM,value)), by = .(ens,mod,lat,type)
                        ][, .(cor = mean(cor), rmse = mean(rmse)), by = .(lat,mod,type)]


ens1_sel_JJA %>% 
  ggplot(aes(x = year, y = value, group = type))+
  geom_line(aes(color = type), size = 0.8) +
  geom_line(aes(year, ZM, color = "ZM"), linetype = 2) +
  ylab("Precip anomalies [mm/month]") +
  xlab("Year")  +
  facet_grid(lat ~., scales = "free")  +
  scale_color_manual(values = c("#FD7C53", "#5396FD", "gainsboro", "grey30"),
                     labels = c("Masked Pr", "Masked Pr + SLP","Simple", "Zonal Mean Pr (target)")) +
  pub_theme +
  theme(strip.text.y.right = element_text(angle = 0),
        legend.position = "none")

ggsave(file = "./figures/ZM_CESM1_ens1_JJA.png", height = 2.5, width = 7)


sumsum_ZMJJA <- melt(sum_ZMJJA, id.vars = c("lat", "mod", "type"), variable.name = "loss")[, .(med = median(value), 
                                                                      upper = max(value), 
                                                                      lower = min(value)),
                                                                  by = .(lat, type, loss)]
sumsum_ZMJJA[loss == "cor", loss := "Correlation"]
sumsum_ZMJJA[loss == "rmse", loss := "RMSE"]


#==========DJF + JJA plot ==================
vline <- data.table(lat = c(57.5,-2.5), seas = c("DJF", "JJA"))


sumsum_ZMJJA %>% data.table(seas = "JJA") %>% 
  rbind(data.table(sumsumsum, seas = "DJF")) %>% 
  .[type %in% c("Simple", "pr1950", "prpsl1950")] %>% 
  .[loss == "Correlation"] %>% 
  ggplot() +
  geom_line(aes(lat, med, group = type, color = type)) +
  geom_ribbon(aes(x = lat, ymax = upper, ymin = lower, group = type, fill = type), alpha = 0.3) +
  ylab("Correlation coefficient between reconstruction and target") +
  xlab("Latitude [°N]") +
  scale_x_continuous(breaks = c(-45, 0, 45))  +
  theme(legend.title = element_blank()) +
  scale_color_manual(labels = c("Masked Pr + SLP", "Masked Pr", "Simple"), values = c("#5396FD","#FD7C53", "gainsboro")) +
  scale_fill_manual( labels = c("Masked Pr + SLP", "Masked Pr", "Simple"), values = c("#5396FD","#FD7C53", "gainsboro")) +
  pub_theme +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0)) +
  facet_grid(seas ~., scales = "free")+
  geom_vline(data = vline, aes(xintercept = lat), alpha = 0.8)

ggsave(file = "./figures/avg_cor_bylat_DJF_JJA.png", height = 5, width = 7)



