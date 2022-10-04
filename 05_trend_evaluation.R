#Run these if you haven't ran any previous scripts

#directory where this script is located
basedir <- "/net/xenon/climphys_backedup/maegli/Precipitation/Scripts/pub/"
setwd(basedir)
source("./functions.R")
source("./00_global_vars.R")

#=============================
load("./data/DJF/CMIP5_pr_PIC_DJF_XAX.RData")
load("./data/DJF/CR20_psl_5d00_DJF_anom_XAX.RData")
load("./data/DJF/CMIP5_psl_PIC_DJF_XAX.RData")

load("./data/JJA/CMIP5_pr_PIC_JJA_XAX.RData")
load("./data/JJA/CMIP5_psl_PIC_JJA_XAX.RData")
load("./data/JJA/CR20_psl_5d00_JJA_anom_XAX.RData")


startyear <- 1950
lambda_sel <- "1se"

load("./models/mod_JJA_prGPCC1950_ZM_single.RData")
mod_prGPCC_JJA <-  mod[[1]]
GPCC_JJA_Recon <- FR_comp("JJA", 
                          GPCC_pr_5d00_JJA_anom_XAX, 
                          GPCC_na_1950_JJA, 
                          mod_prGPCC_JJA, 
                          lambda_sel,
                          startyear)


load("./models/mod_JJA_prGHCN1950_ZM_single.RData")
mod_prGHCN_JJA <-  mod[[1]]
GHCN_JJA_Recon <- FR_comp("JJA", 
                          GHCN_pr_5d00_JJA_anom_XAX, 
                          GHCN_na_1950_JJA, 
                          mod_prGHCN_JJA, 
                          lambda_sel,
                          startyear)


load("./models/mod_DJF_prGHCN1950_ZM_single.RData")
mod_prGHCN_DJF <-  mod[[1]]
GHCN_DJF_Recon <- FR_comp("DJF", 
                          GHCN_pr_5d00_DJF_anom_XAX, 
                          GHCN_na_1950_DJF, 
                          mod_prGHCN_DJF, 
                          lambda_sel,
                          startyear)

load("./models/mod_DJF_prGPCC1950_ZM_single.RData")
mod_prGPCC_DJF <-  mod[[1]]
GPCC_DJF_Recon <- FR_comp("DJF", 
                          GPCC_pr_5d00_DJF_anom_XAX, 
                          GPCC_na_1950_DJF, 
                          mod_prGPCC_DJF, 
                          lambda_sel,
                          startyear)


load("./models/mod_DJF_prpslGHCN1950_ZM_single.RData")
mod_prpslGHCN_DJF <-  mod[[1]]
GHCN_DJF_PSLRecon <- FR_comp("DJF", 
                             GHCN_pr_5d00_DJF_anom_XAX, 
                             GHCN_na_1950_DJF, 
                             mod_prpslGHCN_DJF, 
                             lambda_sel,
                             startyear,
                             multi = T,
                             CR20_psl_5d00_DJF_anom_XAX)

load("./models/mod_DJF_prpslGPCC1950_ZM_single.RData")
mod_prpslGPCC_DJF <-  mod[[1]]
GPCC_DJF_PSLRecon <- FR_comp("DJF", 
                             GPCC_pr_5d00_DJF_anom_XAX, 
                             GPCC_na_1950_DJF, 
                             mod_prpslGPCC_DJF, 
                             lambda_sel,
                             startyear,
                             multi = T,
                             CR20_psl_5d00_DJF_anom_XAX)

load("./models/mod_JJA_prpslGHCN1950_ZM_single.RData")
mod_prpslGHCN_JJA <-  mod[[1]]
GHCN_JJA_PSLRecon <- FR_comp("JJA", 
                             GHCN_pr_5d00_JJA_anom_XAX, 
                             GHCN_na_1950_JJA, 
                             mod_prpslGHCN_JJA, 
                             lambda_sel,
                             startyear,
                             multi = T,
                             CR20_psl_5d00_JJA_anom_XAX)


load("./models/mod_JJA_prpslGPCC1950_ZM_single.RData")
mod_prpslGPCC_JJA <-  mod[[1]]
GPCC_JJA_PSLRecon <- FR_comp("JJA", 
                             GPCC_pr_5d00_JJA_anom_XAX, 
                             GPCC_na_1950_JJA, 
                             mod_prpslGPCC_JJA, 
                             lambda_sel,
                             startyear,
                             multi = T,
                             CR20_psl_5d00_JJA_anom_XAX)


#========================================
obs_trend_frame_recon <- rbind(data.table(GHCN_JJA_Recon[[1]], seas = "JJA", mask = "Pr", obs = "GHCN"),
                               data.table(GPCC_JJA_Recon[[1]], seas = "JJA", mask = "Pr", obs = "GPCC"),
                               data.table(GHCN_DJF_Recon[[1]], seas = "DJF", mask = "Pr", obs = "GHCN"),
                               data.table(GPCC_DJF_Recon[[1]], seas = "DJF", mask = "Pr", obs = "GPCC"),
                               data.table(GHCN_DJF_PSLRecon[[1]], seas = "DJF", mask = "PrPsl", obs = "GHCN"),
                               data.table(GPCC_DJF_PSLRecon[[1]], seas = "DJF", mask = "PrPsl", obs = "GPCC"),
                               data.table(GHCN_JJA_PSLRecon[[1]], seas = "JJA", mask = "PrPsl", obs = "GHCN"),
                               data.table(GPCC_JJA_PSLRecon[[1]], seas = "JJA", mask = "PrPsl", obs = "GPCC")
                               )

LENS_trend_frame_recon <- rbind(data.table(GHCN_JJA_Recon[[2]], seas = "JJA", mask = "Pr", obs = "GHCN"),
                                data.table(GPCC_JJA_Recon[[2]], seas = "JJA", mask = "Pr", obs = "GPCC"),
                                data.table(GHCN_DJF_Recon[[2]], seas = "DJF", mask = "Pr", obs = "GHCN"),
                                data.table(GPCC_DJF_Recon[[2]], seas = "DJF", mask = "Pr", obs = "GPCC"),
                                data.table(GHCN_DJF_PSLRecon[[2]], seas = "DJF", mask = "PrPsl", obs = "GHCN"),
                                data.table(GPCC_DJF_PSLRecon[[2]], seas = "DJF", mask = "PrPsl", obs = "GPCC"),
                                data.table(GHCN_JJA_PSLRecon[[2]], seas = "JJA", mask = "PrPsl", obs = "GHCN"),
                                data.table(GPCC_JJA_PSLRecon[[2]], seas = "JJA", mask = "PrPsl", obs = "GPCC")
                                )

PIC_trend_frame_recon <- rbind(data.table(GHCN_JJA_Recon[[3]], seas = "JJA", mask = "Pr", obs = "GHCN"),
                               data.table(GPCC_JJA_Recon[[3]], seas = "JJA", mask = "Pr", obs = "GPCC"),
                               data.table(GHCN_DJF_Recon[[3]], seas = "DJF", mask = "Pr", obs = "GHCN"),
                               data.table(GPCC_DJF_Recon[[3]], seas = "DJF", mask = "Pr", obs = "GPCC"),
                               data.table(GHCN_DJF_PSLRecon[[3]], seas = "DJF", mask = "PrPsl", obs = "GHCN"),
                               data.table(GPCC_DJF_PSLRecon[[3]], seas = "DJF", mask = "PrPsl", obs = "GPCC"),
                               data.table(GHCN_JJA_PSLRecon[[3]], seas = "JJA", mask = "PrPsl", obs = "GHCN"),
                               data.table(GPCC_JJA_PSLRecon[[3]], seas = "JJA", mask = "PrPsl", obs = "GPCC")
                               )


LENS_trend_quant_recon <- LENS_trend_frame_recon[
  ,.(q975  = quantile(trend, 0.975), 
     q025 = quantile(trend, 0.025),
     q050 = quantile(trend, 0.5)),
  by = .(lat, seas, obs, mask, mod)][
    ,.(q975 = mean(q975), q025 = mean(q025), q050 = mean(q050)), 
    by = .(lat, seas, obs, mask)]

PIC_trend_quant_recon <- PIC_trend_frame_recon[
  , .(q975  = quantile(trend, 0.975), q025 = quantile(trend, 0.025)),
  by = .(lat, seas, obs, mask, mod)][
    ,.(q975 = mean(q975), q025 = mean(q025)), 
    by = .(lat, seas, obs, mask) 
  ]




ggplot() +
  geom_ribbon(data = PIC_trend_quant_recon[mask == "PrPsl",], aes(x = lat, ymax = q975, ymin = q025, fill = "PiC"), alpha = 0.3)+
  geom_ribbon(data = LENS_trend_quant_recon[mask == "PrPsl",], aes(x = lat, ymax = q975, ymin = q025, fill = "LENS"), alpha = 0.3)+
  geom_line(data = obs_trend_frame_recon[mask == "PrPsl",], aes(lat, trend, color = obs))+
  facet_grid(forcats::fct_rev(seas) ~ obs)  +
  scale_color_manual(values = c("forestgreen", "#bf9a3b")) +
  scale_fill_manual(values = c("#3986fa", "grey50")) +
  ggtitle("Masked Pr + SLP") +
  ylab(bquote('ZMP Trend [mm'~month^-1~decade^-1~']')) +  
  scale_x_continuous(breaks = c(-45, 0, 45))  +
  xlab("Latitude [°N]") +
  theme(strip.text.y.right = element_text(angle = 0)) +
  pub_theme

ggsave(file = "./figures/recon_trends_prpsl_1se.png",
       height = 5,
       width = 8)



ggplot() +
  geom_ribbon(data = PIC_trend_quant_recon[mask == "Pr",], aes(x = lat, ymax = q975, ymin = q025, fill = "PiC"), alpha = 0.3)+
  geom_ribbon(data = LENS_trend_quant_recon[mask == "Pr",], aes(x = lat, ymax = q975, ymin = q025, fill = "LENS"), alpha = 0.3)+
  geom_line(data = obs_trend_frame_recon[mask == "Pr",], aes(lat, trend, color = obs))+
  facet_grid(forcats::fct_rev(seas) ~ obs)  +
  pub_theme +
  scale_color_manual(values = c("forestgreen", "#bf9a3b")) +
  scale_fill_manual(values = c("#3986fa", "grey50")) +
  ggtitle("Masked Pr") +
  ylab(bquote('ZMP Trend [mm'~month^-1~decade^-1~']')) +  
  scale_x_continuous(breaks = c(-45, 0, 45))  +
  xlab("Latitude [°N]") +
  theme(strip.text.y.right = element_text(angle = 0))

ggsave(file = "./figures/recon_trends_pr_1se.png",
       height = 5,
       width = 8)


#============Desnsity Plot==========================
fraw_mat <- LENS_pr_DJF_5d00_XAX_norm$Y$fraw %>% unique()
fraw_mat_JJA <- LENS_pr_JJA_5d00_XAX_norm$Y$fraw %>% unique()


mm_fr <- data.table(fraw_mat, 
                    year = LENS_pr_DJF_5d00_XAX_norm$M$year[LENS_pr_DJF_5d00_XAX_norm$M$ens.mem == "r10i1p1"]) %>% 
  .[year %in% 1950:2080, lapply(.SD, mean), by = year] %>% 
  .[,year := NULL] %>% as.matrix()
mm_fr_trend_DJF <- as.numeric(apply(mm_fr,2, function(x) lm(x ~ seq_along(x))$coefficient[2])  *10)



mm_fr_JJA <- data.table(fraw_mat_JJA, 
                        year = LENS_pr_JJA_5d00_XAX_norm$M$year[LENS_pr_JJA_5d00_XAX_norm$M$ens.mem == "r10i1p1"]) %>% 
  .[year %in% 1950:2080, lapply(.SD, mean), by = year] %>% 
  .[,year := NULL] %>% as.matrix()
mm_fr_trend_JJA <- as.numeric(apply(mm_fr_JJA,2, function(x) lm(x ~ seq_along(x))$coefficient[2])  *10)


mm_fr_trend_DJF <- apply(mm_fr,2, function(x) lm(x ~ seq_along(x))$coefficient[2])  *10
mm_fr_trend_JJA <- apply(mm_fr_JJA,2, function(x) lm(x ~ seq_along(x))$coefficient[2]) *10
mm_fr_frame <- data.table(DJF = mm_fr_trend_DJF,
                          JJA = mm_fr_trend_JJA,
                          lat = latitudes) %>% 
  melt(id.vars = "lat", value.name = "FR", variable.name = "seas")


ggplot(mm_fr_frame) +
  geom_line(aes(lat, FR, group = seas, color = seas))

CovFR_obs_DJF <- obs_trend_frame_recon[seas == "DJF", .(obs_cov = cov(trend, mm_fr_trend_DJF)), by = .(seas, mask, obs)]
CovFR_LENS_DJF <- LENS_trend_frame_recon[seas == "DJF", .(LENS_cov = cov(trend, mm_fr_trend_DJF)), by = .(mod, ens, seas, mask, obs)]
CovFR_PIC_DJF <- PIC_trend_frame_recon[seas == "DJF", .(PIC_cov = cov(trend, mm_fr_trend_DJF)), by = .(mod, seas, mask, obs, group)]


CovFR_obs_JJA <- obs_trend_frame_recon[seas == "JJA", .(obs_cov = cov(trend, mm_fr_trend_JJA)), by = .(seas, mask, obs)]
CovFR_LENS_JJA <- LENS_trend_frame_recon[seas == "JJA", .(LENS_cov = cov(trend, mm_fr_trend_JJA)), by = .(mod, ens, seas, mask, obs)]
CovFR_PIC_JJA <- PIC_trend_frame_recon[seas == "JJA", .(PIC_cov = cov(trend, mm_fr_trend_JJA)), by = .(mod, seas, mask, obs, group)]

CovFR_obs <- rbind(CovFR_obs_JJA, CovFR_obs_DJF)
CovFR_LENS <- rbind(CovFR_LENS_JJA, CovFR_LENS_DJF)
CovFR_PIC <- rbind(CovFR_PIC_JJA, CovFR_PIC_DJF)


ecdf_fun <- function(x,perc) ecdf(x)(perc)
CovFR_obsPiC <-  CovFR_PIC[CovFR_obs, on = .(mask, seas, obs)]
CovQuant_PIC <- CovFR_obsPiC[, ecdf_fun(PIC_cov, obs_cov), by= .(seas,mask,obs)] %>% unique() %>% .[,V1 := round(V1, 2)]
CovQuant_PIC[,x:= 0.5]
CovQuant_PIC[,y:= 4.5]



CovFR_obsLENS <-  CovFR_LENS[CovFR_obs, on = .(mask, seas, obs)]
CovQuant_LENS <- CovFR_obsLENS[, ecdf_fun(LENS_cov, obs_cov), by= .(seas,mask,obs)] %>% unique() %>% .[,V1 := round(V1, 2)]
CovQuant_LENS[,x:= 0.5]
CovQuant_LENS[,y:= 5.5]

maks_sel <- "PrPsl"   #"Pr" or "PrPsl"

ggplot() +
  geom_density(data = CovFR_LENS[mask  == maks_sel], aes(LENS_cov,  fill = "LENS"),alpha = 0.5, color = NA) +
  geom_density(data = CovFR_PIC[mask == maks_sel], aes(PIC_cov, fill= "PiC"), alpha = 0.5, color = NA)+
  geom_vline(data = CovFR_obs[mask == maks_sel], aes(xintercept = obs_cov, color = obs))+
  geom_text(data = CovQuant_LENS[mask == maks_sel],aes(x,y, label = V1), color = "#3986fa")+
  geom_text(data = CovQuant_PIC[mask == maks_sel],aes(x,y, label = V1), color = "grey60")+
  facet_grid(seas ~ obs) +
  ggtitle("Pr + SLP") +
  pub_theme +
  scale_color_manual(values = c("forestgreen", "#bf9a3b")) +
  scale_fill_manual(values = c("#3986fa", "grey50")) +
  ylab("Density") +
  xlab("Covariance of Reconstructed Zonal Trends with Forced Response")

ggsave(file = "./figures/cov_trend_distr_prpsl.png",
       height = 5,
       width = 8)





