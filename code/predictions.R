# This script is for 1) phenology predictions - current (1981-2010)
#                                             - future (2041-2060, 2081-2100) SSP2-4.5 SSP5-8.5
#                    2) temporal synchrony calculation
#                    3) visualization

library(tidyverse)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(sf)

cv.df = read_csv(file="data_processed/Claytonia_virginica_iNat_cleaned_data_all.csv")
ae.df = read_csv(file="data_processed/Andrena_erigeniae_occurence_cleaned_data_all_anm_20220606.csv")

#--------------------- prepare data
# reduce duplicated data at the same location and on the same day
cv.df1 <- cv.df %>% select(-gbifID) %>% distinct() 
ae.df1 <- ae.df %>% select(-gbifID, -datasetName) %>% distinct() 

# remove NAs
cv.df1a = cv.df1 %>% filter(!is.na(MAT))
ae.dfa = ae.df1 %>% filter(!is.na(MAT))

# remove outliers
cv.df1a = filter(cv.df1a, doy<=160, year>=2008)    # 17494
ae.dfa = filter(ae.dfa, doy<=180)    # 1168

# select data by ecoregions
t3=as.data.frame(table(cv.df1a[, c("NA_L3CODE", "color2") ]))
t3=t3[t3$Freq>30,]
t3 = mutate(t3, eco3_color2=paste(NA_L3CODE, color2, sep="_"))
t3

cv.df1c=cv.df1a %>% mutate(eco3_color2=paste(NA_L3CODE, color2, sep="_")) %>% 
  filter(eco3_color2 %in% t3$eco3_color2)

ta3=as.data.frame(table(ae.dfa$NA_L3CODE))
ae.dfc=ae.dfa %>% filter(NA_L3CODE %in% ta3$Var1[ta3$Freq>=20])


# load models
load("data_results/springbeauty_lmer_models_rv.Rdata")


# current year estimation (1981-2010) -------------------------------------------

# calculate 30y normal ecoregion level winter temp
cv.sy=read_csv("data_raw/springbeauty_locations_all_1901-2021SY.csv") 
ae.sy=read_csv("data_raw/minerbee_locations_20220510_1901-2021SY.csv") 
ae.syca=read_csv("data_raw/minerbee_locations_CA_20220606_1901-2021SY.csv") 



colnames(cv.sy)
cv.sy1 = cv.sy[, c(1:6, 11, 19)] %>% 
  distinct() %>%
  rename(year=Year, sc=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) %>%
  filter(year>=1981, year<=2010) 
cv.LT =cv.sy %>% group_by(ID1) %>% summarise(LT_Tmin_wt=mean(Tmin_wt),
                                             LT_PPT_wt=mean(PPT_wt))


unique(ae.sy$ID1)
unique(ae.syca$ID1)
ae.sy$ID1a=ae.sy$ID1
ae.syca$ID1a=ae.syca$ID1+832
unique(ae.syca$ID1a)

ae.sy1 =rbind.data.frame(ae.sy[, c(1:6, 11, 19, 92)], ae.syca[, c(1:6, 11, 19, 92)]) %>% 
  distinct() 

ae.LT =ae.sy1 %>% group_by(ID1a) %>% summarise(LT_Tmin_wt=mean(Tmin_wt),
                                               LT_PPT_wt=mean(PPT_wt))

ae.sy2 = ae.sy1 %>% 
  rename(year=Year, sc=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) %>%
  filter(year>=1981, year<=2010) 


# prepare new data
# spatial and random effect groups
cv.sp=cv.df1c %>% select(lat, lon, elevation, color2, NA_L3CODE) %>% distinct()
ae.sp=ae.dfc %>% select(lat, lon, elevation, NA_L3CODE) %>% distinct()

cv.clm = cv.sy1 %>% group_by(ID1) %>% summarise(Tm=mean(Tmin_wt), Tsd=sd(Tmin_wt), 
                                                Pm=mean(PPT_wt), Psd=sd(PPT_wt))
cv.clm = left_join(cv.clm, distinct(cv.sy1[, c(2:6)]))
cv.clm$Tm.w=cv.clm$Tm+cv.clm$Tsd
cv.clm$Tm.c=cv.clm$Tm-cv.clm$Tsd
cv.clm$Pm.w=cv.clm$Pm+cv.clm$Psd
cv.clm$Pm.d=cv.clm$Pm-cv.clm$Psd
cv.clm = left_join(cv.clm, cv.LT)

cv.clm = cv.clm %>% mutate(Tm_anm = Tm-LT_Tmin_wt, 
                           Tm.w_anm = Tm.w-LT_Tmin_wt, 
                           Tm.c_anm = Tm.c-LT_Tmin_wt,
                           Pm_anm = (Pm-LT_PPT_wt)/LT_PPT_wt,
                           Pm.w_anm = (Pm.w-LT_PPT_wt)/LT_PPT_wt,
                           Pm.d_anm = (Pm.d-LT_PPT_wt)/LT_PPT_wt)

cv.nd = left_join(cv.sp, cv.clm[, c(7:9, 16:21)])


ae.clm = ae.sy2 %>% group_by(ID1a) %>% summarise(Tm=mean(Tmin_wt), Tsd=sd(Tmin_wt), 
                                                 Pm=mean(PPT_wt), Psd=sd(PPT_wt))
ae.clm = left_join(ae.clm, distinct(ae.sy2[, c(4:6,9)]))
ae.clm$Tm.w=ae.clm$Tm+ae.clm$Tsd
ae.clm$Tm.c=ae.clm$Tm-ae.clm$Tsd
ae.clm$Pm.w=ae.clm$Pm+ae.clm$Psd
ae.clm$Pm.d=ae.clm$Pm-ae.clm$Psd
ae.clm = left_join(ae.clm, ae.LT)

ae.clm = ae.clm %>% mutate(Tm_anm = Tm-LT_Tmin_wt, 
                           Tm.w_anm = Tm.w-LT_Tmin_wt, 
                           Tm.c_anm = Tm.c-LT_Tmin_wt,
                           Pm_anm = (Pm-LT_PPT_wt)/LT_PPT_wt,
                           Pm.w_anm = (Pm.w-LT_PPT_wt)/LT_PPT_wt,
                           Pm.d_anm = (Pm.d-LT_PPT_wt)/LT_PPT_wt)

ae.nd = left_join(ae.sp, ae.clm[, c(6:8, 15:20)])



# 30-y average
cv.pred1 = cv.nd %>% rename(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm_anm) 
cv.nd$fit.norm=predict(cm.wmin, newdata=cv.pred1)

ae.pred1 = ae.nd %>% rename(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm_anm)
ae.nd$fit.norm=predict(am.wmin, newdata=ae.pred1)

# warm
cv.pred2 = cv.nd %>% mutate(Tmin_wt_anm=Tm.w_anm, PPT_wt_anm=Pm_anm)
cv.nd$fit.warm=predict(cm.wmin, newdata=cv.pred2)

ae.pred2 = ae.nd %>% mutate(Tmin_wt_anm=Tm.w_anm, PPT_wt_anm=Pm_anm)
ae.nd$fit.warm=predict(am.wmin, newdata=ae.pred2)

# cold
cv.pred3 = cv.nd %>% mutate(Tmin_wt_anm=Tm.c_anm, PPT_wt_anm=Pm_anm)
cv.nd$fit.cold=predict(cm.wmin, newdata=cv.pred3)

ae.pred3 = ae.nd %>% mutate(Tmin_wt_anm=Tm.c_anm, PPT_wt_anm=Pm_anm)
ae.nd$fit.cold=predict(am.wmin, newdata=ae.pred3)

# wet
cv.pred4 = cv.nd %>% mutate(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm.w_anm)
cv.nd$fit.wet=predict(cm.wmin, newdata=cv.pred4)

ae.pred4 = ae.nd %>% mutate(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm.w_anm)
ae.nd$fit.wet=predict(am.wmin, newdata=ae.pred4)

# dry
cv.pred5 = cv.nd %>% mutate(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm.d_anm)
cv.nd$fit.dry=predict(cm.wmin, newdata=cv.pred5)

ae.pred5 = ae.nd %>% mutate(Tmin_wt_anm=Tm_anm, PPT_wt_anm=Pm.d_anm)
ae.nd$fit.dry=predict(am.wmin, newdata=ae.pred5)

summary(cv.nd)  
summary(ae.nd)

# future prediction -------------------------------------------------------------
# spring beauty -----------------------------------------------------------------

cv.gcm4=read_csv("data_raw/springbeauty_locations_all_4 GCMsSY.csv") 

cv.gcm4= cv.gcm4 %>% 
  rename(year=Year, sc=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  select(year, lat , lon, elevation, ID1, Tmin_wt, PPT_wt) 

unique(cv.gcm4$year)

cv.gcm4=left_join(cv.gcm4, cv.LT) %>% 
  mutate(Tmin_wt_anm = Tmin_wt-LT_Tmin_wt, 
         PPT_wt_anm = (PPT_wt-LT_PPT_wt)/LT_PPT_wt)


# 13GCM-ensemble-ssp245-2041-2060
cv.fpd1 = left_join(cv.nd[, c(1:5)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp245_2041-2060.gcm")) %>% distinct()
cv.nd$ssp245_4160=predict(cm.wmin, newdata=cv.fpd1)

# 13GCM-ensemble-ssp585-2041-2060
cv.fpd2 = left_join(cv.nd[, c(1:5)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp585_2041-2060.gcm")) %>% distinct()
cv.nd$ssp585_4160=predict(cm.wmin, newdata=cv.fpd2)

# 13GCM-ensemble-ssp245-2081-2100
cv.fpd3 = left_join(cv.nd[, c(1:5)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp245_2081-2100.gcm")) %>% distinct()
cv.nd$ssp245_8100=predict(cm.wmin, newdata=cv.fpd3)

# 13GCM-ensemble-ssp585-2081-2100
cv.fpd4 = left_join(cv.nd[, c(1:5)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp585_2081-2100.gcm")) %>% distinct()
cv.nd$ssp585_8100=predict(cm.wmin, newdata=cv.fpd4)

summary(cv.nd)


# spring beauty miner bee -----------------------------------------------------------
ae.gcm4a=read_csv("data_raw/minerbee_locations_20220510_4 GCMsSY.csv") 
ae.gcm4b=read_csv("data_raw/minerbee_locations_CA_20220606_4 GCMsSY.csv") 


ae.gcm4a$ID1a = ae.gcm4a$ID1
ae.gcm4b$ID1a = ae.gcm4b$ID1 + 832
ae.gcm4 = rbind.data.frame(ae.gcm4a, ae.gcm4b) %>% 
  rename(year=Year, sc=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  dplyr::select(year, lat , lon, elevation, ID1a, Tmin_wt, PPT_wt) 

unique(ae.gcm4$year)

ae.gcm4=left_join(ae.gcm4, ae.LT) %>% 
  mutate(Tmin_wt_anm = Tmin_wt-LT_Tmin_wt, 
         PPT_wt_anm = (PPT_wt-LT_PPT_wt)/LT_PPT_wt)


# 13GCM-ensemble-ssp245-2041-2060
ae.fpd1 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp245_2041-2060.gcm")) %>% distinct()
ae.nd$ssp245_4160=predict(am.wmin, newdata=ae.fpd1)

# 13GCM-ensemble-ssp585-2041-2060
ae.fpd2 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp585_2041-2060.gcm")) %>% distinct()
ae.nd$ssp585_4160=predict(am.wmin, newdata=ae.fpd2)

# 13GCM-ensemble-ssp245-2081-2100
ae.fpd3 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp245_2081-2100.gcm")) %>% distinct()
ae.nd$ssp245_8100=predict(am.wmin, newdata=ae.fpd3)

# 13GCM-ensemble-ssp585-2081-2100
ae.fpd4 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp585_2081-2100.gcm")) %>% distinct()
ae.nd$ssp585_8100=predict(am.wmin, newdata=ae.fpd4)

summary(ae.nd)

save(cv.nd, ae.nd, file="data_results/predictions_current30y_future4gcms_rv.Rdata")
load("data_results/predictions_current30y_future4gcms_rv.Rdata")


# plot predictions -----------------------------------------------------------------
# long format predictions
cv.pred=cv.nd[, c(1:5, 12:20)] %>% pivot_longer(cols=c(6:14), names_to = "scenario", values_to = "pred")
ae.pred=ae.nd[, c(1:4, 11:19)] %>% pivot_longer(cols=c(5:13), names_to = "scenario", values_to = "pred")
unique(cv.pred$scenario)
unique(ae.pred$scenario)

ae.pred$color2="bee"

pred=rbind.data.frame(cv.pred, ae.pred)
unique(pred$color2)

summary(pred)

ae.eco=unique(ae.pred$NA_L3CODE)
scn5=c("fit.norm","ssp245_4160","ssp245_8100","ssp585_4160","ssp585_8100")

ggplot(filter(pred, NA_L3CODE %in% ae.eco, scenario %in% scn5), aes(x=NA_L3CODE, y=pred, fill=color2))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  scale_fill_manual(name="", values=c("orange","coral1","deepskyblue"))+
  facet_wrap(~scenario, ncol=1)+
  theme_cowplot()


# maps for predictions at L3 ecoregion scale
cv.peco=cv.pred %>% group_by(scenario, color2, NA_L3CODE) %>% 
  summarise(n=n(), pred.m=mean(pred), pred.md=median(pred), pred.sd=sd(pred))

ae.peco=ae.pred %>% group_by(scenario, NA_L3CODE) %>% 
  summarise(n=n(), pred.m=mean(pred), pred.md=median(pred), pred.sd=sd(pred))


eco3=st_read(dsn="ecoregions/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")
eco3=st_transform(eco3, crs=4326)
eco3=st_make_valid(eco3)

eco3_cv=left_join(eco3, cv.peco[, c(1:3, 5)] %>% 
                    pivot_wider(names_from = c(scenario, color2), values_from = pred.m)) 


eco3_ae=left_join(eco3, ae.peco[, c(1:2,4)] %>% 
                    pivot_wider(names_from = scenario, values_from = pred.m)) 




# current phenology - Figure S1
pw1=ggplot(eco3_cv)+
  geom_sf(aes(fill=fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

pp1=ggplot(eco3_cv)+
  geom_sf(aes(fill=fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

pb1=ggplot(eco3_ae)+
  geom_sf(aes(fill=fit.norm))+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80", name="day of year")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

maps=cowplot::plot_grid(pw1, pp1, pb1, labels=c("a","b","c"), label_size=16, ncol=3, nrow=1, 
                        rel_widths=c(1, 1, 1))

ggsave(filename = "figures_rv/maps_current_predictions.tif", maps, device="tiff",
       dpi=300, height= 5, width=17, compression="lzw")


# future predictions - Figure S2
pw2=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_4160_white), show.legend = F)+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp2=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_4160_pink), show.legend = F)+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb2=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp245_4160))+
  scale_fill_viridis_c(limits=c(45, 160), na.value = "grey80", name="day of year")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

pw3=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_8100_white), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp3=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_8100_pink), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb3=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp245_8100), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

pw4=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_4160_white), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp4=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_4160_pink), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb4=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp585_4160), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

pw5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_white), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

pp5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_pink), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

pb5=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp585_8100), show.legend = F)+
  scale_fill_viridis_c( limits=c(45, 160), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

maps_future=cowplot::plot_grid(pw2, pp2, pb2, pw3, pp3, pb3, pw4, pp4, pb4, pw5, pp5, pb5, 
                               labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"), 
                               label_size=16, ncol=3, nrow=4, 
                               rel_widths=c(1, 1, 1), rel_heights = c(1, 1, 1, 1))

ggsave(filename = "figures_rv/maps_future_predictions.tif", maps_future, device="tiff",
       dpi=300, height= 20, width=18, compression="lzw")


# maps of changes in future predictions - Figure S3
summary(eco3_cv$ssp245_4160_white-eco3_cv$fit.norm_white)
summary(eco3_cv$ssp585_8100_white-eco3_cv$fit.norm_white)
summary(eco3_cv$ssp585_8100_pink-eco3_cv$fit.norm_pink)
summary(eco3_ae$ssp585_8100-eco3_ae$fit.norm)

pw2=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_4160_white-fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp2=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_4160_pink-fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb2=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp245_4160-fit.norm))+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80", name="days")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

pw3=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_8100_white-fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp3=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp245_8100_pink-fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb3=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp245_8100-fit.norm), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pw4=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_4160_white-fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pp4=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_4160_pink-fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pb4=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp585_4160-fit.norm), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

pw5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_white-fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

pp5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_pink-fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

pb5=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp585_8100-fit.norm), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

diff_future=cowplot::plot_grid(pw2, pp2, pb2, pw3, pp3, pb3, pw4, pp4, pb4, pw5, pp5, pb5, 
                               labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"), 
                               label_size=16, ncol=3, nrow=4, 
                               rel_widths=c(1, 1, 1), rel_heights = c(1, 1, 1, 1))

ggsave(filename = "figures_rv/maps_changes_future_predictions.tif", diff_future, device="tiff",
       dpi=300, height= 20, width=18, compression="lzw")


# ssp585_8100 only - Figure 5
pw5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_white-fit.norm_white), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

pp5=ggplot(eco3_cv)+
  geom_sf(aes(fill=ssp585_8100_pink-fit.norm_pink), show.legend = F)+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

pb5=ggplot(eco3_ae)+
  geom_sf(aes(fill=ssp585_8100-fit.norm))+
  scale_fill_viridis_c(option="A", limits=c(-27,3), na.value = "grey80", name="days")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

diff_ssp585_8100=cowplot::plot_grid(pw5, pp5, pb5, labels=c("a","b","c"), 
                                    label_size=16, ncol=3, nrow=1, rel_widths=c(1, 1, 1))

ggsave(filename = "figures_rv/maps_changes_future_predictions_ssp585_8100.tif", diff_ssp585_8100, device="tiff",
       dpi=300, height= 5, width=17, compression="lzw")


# calculate changes in temporal mismatch in future predictions

cv.peco.w = cv.peco[, c(1:3, 5)] %>% 
  pivot_wider(names_from = c(scenario, color2), values_from = pred.m)
cv.peco.w$ssp585_8100_pink-cv.peco.w$fit.norm_pink
cv.peco.w$ssp585_8100_white-cv.peco.w$fit.norm_white

ae.peco.w = ae.peco[, c(1:2,4)] %>% 
  pivot_wider(names_from = scenario, values_from = pred.m)
ae.peco.w$ssp585_8100-ae.peco.w$fit.norm


pred1=left_join(cv.peco.w, ae.peco.w) %>% 
  filter(NA_L3CODE %in% unique(ae.peco.w$NA_L3CODE)) %>% 
  na.omit()

# predicted temporal gaps between bee and flower
pred1$ssp245_4160_wdiff=pred1$ssp245_4160-pred1$ssp245_4160_white
pred1$ssp245_8100_wdiff=pred1$ssp245_8100-pred1$ssp245_8100_white
pred1$ssp585_4160_wdiff=pred1$ssp585_4160-pred1$ssp585_4160_white
pred1$ssp585_8100_wdiff=pred1$ssp585_8100-pred1$ssp585_8100_white
pred1$ssp245_4160_pdiff=pred1$ssp245_4160-pred1$ssp245_4160_pink
pred1$ssp245_8100_pdiff=pred1$ssp245_8100-pred1$ssp245_8100_pink
pred1$ssp585_4160_pdiff=pred1$ssp585_4160-pred1$ssp585_4160_pink
pred1$ssp585_8100_pdiff=pred1$ssp585_8100-pred1$ssp585_8100_pink
pred1$fit.norm_wdiff=pred1$fit.norm-pred1$fit.norm_white
pred1$fit.norm_pdiff=pred1$fit.norm-pred1$fit.norm_pink

# change in temporal gaps
pred1$ssp245_4160_wd=abs(pred1$ssp245_4160_wdiff)-abs(pred1$fit.norm_wdiff)
pred1$ssp245_8100_wd=abs(pred1$ssp245_8100_wdiff)-abs(pred1$fit.norm_wdiff)
pred1$ssp585_4160_wd=abs(pred1$ssp585_4160_wdiff)-abs(pred1$fit.norm_wdiff)
pred1$ssp585_8100_wd=abs(pred1$ssp585_8100_wdiff)-abs(pred1$fit.norm_wdiff)
pred1$ssp245_4160_pd=abs(pred1$ssp245_4160_pdiff)-abs(pred1$fit.norm_pdiff)
pred1$ssp245_8100_pd=abs(pred1$ssp245_8100_pdiff)-abs(pred1$fit.norm_pdiff)
pred1$ssp585_4160_pd=abs(pred1$ssp585_4160_pdiff)-abs(pred1$fit.norm_pdiff)
pred1$ssp585_8100_pd=abs(pred1$ssp585_8100_pdiff)-abs(pred1$fit.norm_pdiff)

# proportional changes in temporal gaps
pred1$ssp245_4160_wdp=pred1$ssp245_4160_wd/pred1$fit.norm_wdiff
pred1$ssp245_8100_wdp=pred1$ssp245_8100_wd/pred1$fit.norm_wdiff
pred1$ssp585_4160_wdp=pred1$ssp585_4160_wd/pred1$fit.norm_wdiff
pred1$ssp585_8100_wdp=pred1$ssp585_8100_wd/pred1$fit.norm_wdiff
pred1$ssp245_4160_pdp=pred1$ssp245_4160_pd/pred1$fit.norm_pdiff
pred1$ssp245_8100_pdp=pred1$ssp245_8100_pd/pred1$fit.norm_pdiff
pred1$ssp585_4160_pdp=pred1$ssp585_4160_pd/pred1$fit.norm_pdiff
pred1$ssp585_8100_pdp=pred1$ssp585_8100_pd/pred1$fit.norm_pdiff

write_csv(pred1, file="data_results/springbeauty_predictions_gaps_rv.csv")
pred1=read_csv(file="data_results/springbeauty_predictions_gaps_rv.csv")

colnames(pred1)
pred2=pred1[, c(1, 39:46)] 
colnames(pred2)[2:9]=substr(colnames(pred2)[2:9], 1, 11)
pred2=rbind.data.frame(pred2[,c(1:5)], pred2[, c(1, 6:9)])
pred2$color=c(rep("white",12), rep("pink",12))
pred2= pred2 %>% pivot_longer(cols=2:5, names_to = "scenario", values_to = "diff")


# Figure 6
labs=c("2041-2060 SSP2-4.5", "2081-2100 SSP2-4.5", "2041-2060 SSP5-8.5", "2081-2100 SSP5-8.5")
names(labs) = c("ssp245_4160", "ssp245_8100", "ssp585_4160", "ssp585_8100" )

tiff(filename="figures_rv/change_gaps_predictions.tif", width=2000, height=2000, res=300, compression="lzw")
ggplot(pred2, aes(x=NA_L3CODE, y=diff, fill=color))+
  geom_bar(stat = "identity", position="dodge")+
  scale_fill_manual(values=c("coral1","deepskyblue"),name="color morph")+
  facet_wrap(~scenario, ncol=1, labeller = labeller(scenario = labs))+
  labs(x="ecoregion", y="days")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.spacing.y = unit(0.3, 'cm'),
        legend.box.margin=margin(10,10,10,20))+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))
dev.off()

range(pred2$diff)
pred2 %>% filter(NA_L3CODE != "8.5.1") %>% dplyr::select(diff) %>% range()
pred2 %>% filter(NA_L3CODE != "8.5.1") %>% dplyr::select(diff) %>% view()


#
pred3=pred1[, c(1, 6:7, 12:19, 22,25:28)] 
colnames(pred3)[2:3]="fit.norm"
colnames(pred3)[4:11]=substr(colnames(pred3)[4:11], 1, 11)
pred3=rbind.data.frame(pred3[, c(1, 3, 5, 7, 9, 11)], pred3[, c(1, 2, 4, 6, 8, 10)],
                       pred3[, c(1, 12:16)])
pred3$color=c(rep("white",12), rep("pink",12), rep("bee",12))
pred3= pred3 %>% pivot_longer(cols=2:6, names_to = "scenario", values_to = "doy")



tiff(filename="figures_rv/gaps_predictions1.tif", width=2700, height=2000, res=300, compression="lzw")
ggplot(pred3, aes(x=scenario, y=doy, color=color))+
  geom_point()+
  facet_wrap(~NA_L3CODE, scale="free_y")+
  labs(x="scenario", y="day of year")+
  scale_color_manual(values=c("orange","coral1","deepskyblue"), name="")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.spacing.y = unit(0.3, 'cm'),
        legend.box.margin=margin(10,10,10,20))+
  ## important additional element
  scale_x_discrete(labels=c("fit.norm"="1981-2010", "ssp245_4160"="2041-2060 SSP2-4.5",
                            "ssp245_8100"="2081-2100 SSP2-4.5", 
                            "ssp585_4160"="2041-2060 SSP5-8.5", 
                            "ssp585_8100"="2081-2100 SSP5-8.5"))+
  
  guides(fill = guide_legend(byrow = TRUE))
dev.off()



pred4 = pred3 %>% filter(scenario %in% c("fit.norm","ssp585_8100"))

ggplot(pred4, aes(x=scenario, y=doy, color=color, shape=scenario))+
  geom_point(position=position_dodge(width=0.3))+
  facet_wrap(~NA_L3CODE, nrow=2, scale="free_y")+
  labs(x="scenario", y="day of year")+
  scale_color_manual(values=c("orange","coral1","deepskyblue"), name="")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.spacing.y = unit(0.3, 'cm'),
        legend.box.margin=margin(10,10,10,20))+
  ## important additional element
  scale_x_discrete(labels=c("fit.norm"="1981-2010", "ssp585_8100"="2081-2100 SSP5-8.5"))+
  guides(fill = guide_legend(byrow = TRUE))


df_seg=pivot_wider(pred4, names_from = scenario, values_from = doy)

tiff(filename="figures_rv/gaps_predictions2b.tif", width=2700, height=1500, res=300, compression="lzw")
ggplot(pred4, aes(x=color, y=doy, color=color))+
  geom_point(aes(shape=scenario), size=2)+
  geom_segment(data=df_seg, aes(x = color, y = fit.norm, xend = color, yend = ssp585_8100), arrow = arrow(length = unit(0.2,"cm")), linetype=2, color="grey50")+
  facet_wrap(~NA_L3CODE, nrow=2)+
  labs(x="taxa/color morph", y="day of year")+
  scale_color_manual(values=c("orange","coral1","deepskyblue"), name="taxa/color")+
  scale_shape_manual(values=c(19, 17), name="scenario", 
                     labels=c("fit.norm"="1981-2010", "ssp585_8100"="2081-2100 SSP5-8.5"))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.spacing.y = unit(0.3, 'cm'),
        legend.box.margin=margin(10,10,10,20))+
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE))
dev.off()



# map changes in gaps
eco3_gaps=left_join(eco3, pred1[, c(1, 39:46)]) 

wd1=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp245_4160_wd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  #scale_fill_gradient2(high="red", low="blue", limits=c(-2, 16))+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pd1=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp245_4160_pd))+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80", 
                       name="days")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

wd2=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp245_8100_wd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  #scale_fill_gradient2(high="red", low="blue", limits=c(-2, 16))+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pd2=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp245_8100_pd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

wd3=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp585_4160_wd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  #scale_fill_gradient2(high="red", low="blue", limits=c(-2, 16))+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="latitude")+
  theme_cowplot()

pd3=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp585_4160_pd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="", y="")+
  theme_cowplot()

wd4=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp585_8100_wd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  #scale_fill_gradient2(high="red", low="blue", limits=c(-2, 16))+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

pd4=ggplot(eco3_gaps)+
  geom_sf(aes(fill=ssp585_8100_pd), show.legend = F)+
  scale_fill_viridis_c(option="B", limits=c(-26, 9), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="")+
  theme_cowplot()

diff_gaps=cowplot::plot_grid(wd1, pd1, wd2, pd2, wd3, pd3, wd4, pd4,  
                             labels=c("a","b","c","d","e","f","g","h"), 
                             label_size=16, ncol=2, nrow=4, 
                             rel_widths=c(1, 1), rel_heights = c(1, 1, 1, 1))

ggsave(filename = "figures_rv/maps_gaps_changes_future_predictions.tif", diff_gaps, device="tiff",
       dpi=300, height= 18, width=10, compression="lzw")



# temporal gaps between color morphs
colnames(cv.peco.w)
cv.peco.w$fit.norm_gap=cv.peco.w$fit.norm_pink-cv.peco.w$fit.norm_white
cv.peco.w$ssp245_4160_gap=cv.peco.w$ssp245_4160_pink-cv.peco.w$ssp245_4160_white
cv.peco.w$ssp245_8100_gap=cv.peco.w$ssp245_8100_pink-cv.peco.w$ssp245_8100_white
cv.peco.w$ssp585_4160_gap=cv.peco.w$ssp585_4160_pink-cv.peco.w$ssp585_4160_white
cv.peco.w$ssp585_8100_gap=cv.peco.w$ssp585_8100_pink-cv.peco.w$ssp585_8100_white

cv.peco.w$ssp245_4160_gd=cv.peco.w$ssp245_4160_gap-cv.peco.w$fit.norm_gap
cv.peco.w$ssp245_8100_gd=cv.peco.w$ssp245_8100_gap-cv.peco.w$fit.norm_gap
cv.peco.w$ssp585_4160_gd=cv.peco.w$ssp585_4160_gap-cv.peco.w$fit.norm_gap
cv.peco.w$ssp585_8100_gd=cv.peco.w$ssp585_8100_gap-cv.peco.w$fit.norm_gap

mean(abs(cv.peco.w$fit.norm_gap), na.rm=T) # 2.46

write.csv(cv.peco.w, file="data_results/springbeauty_predictions_2colors_gaps_rv.csv")

colnames(cv.peco.w)
cv.peco.w2=cv.peco.w[, c(1, 25:28)] 
colnames(cv.peco.w2)[2:5]=substr(colnames(cv.peco.w2)[2:5], 1, 11)
cv.peco.w2=na.omit(cv.peco.w2)
cv.peco.w2= cv.peco.w2 %>% pivot_longer(cols=2:5, names_to = "scenario", values_to = "gd")


# supplementary Figure S4
labs=c("2041-2060 SSP2-4.5", "2081-2100 SSP2-4.5", "2041-2060 SSP5-8.5", "2081-2100 SSP5-8.5")

tiff(filename="figures_rv/change_gaps_2colors_predictions.tif", width=2700, height=2400, res=300, compression="lzw")
ggplot(cv.peco.w2, aes(x=scenario, y=gd, fill=scenario))+
  geom_bar(stat = "identity", position="dodge")+
  scale_fill_manual(values=c("chocolate1","chocolate4","brown1","brown4"),
                    labels=labs, name="period & scenario")+
  facet_wrap(~NA_L3CODE, ncol=5)+
  labs(x="scenario & time period in predictions", y="change in temporal gap (days)")+
  theme_cowplot(15)+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
dev.off()


