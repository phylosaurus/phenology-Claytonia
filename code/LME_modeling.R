# This script is for linear mixed effect modeling and plot model outputs


library(tidyverse)
library(sf)
library(lme4)
library(lmerTest)
library(MuMIn)

# read data
cv.df = read_csv(file="data_processed/Claytonia_virginica_iNat_cleaned_data_all.csv")
ae.df = read_csv(file="data_processed/Andrena_erigeniae_occurence_cleaned_data_all_anm_20220606.csv")

# reduce duplicated data at the same location and on the same day
cv.df1 <- cv.df %>% select(-gbifID) %>% distinct() 

ae.df1 <- ae.df %>% select(-gbifID, -datasetName) %>% distinct() 

# remove NAs
cv.df1a = cv.df1 %>% filter(!is.na(MAT))
ae.dfa = ae.df1 %>% filter(!is.na(MAT))

# remove outliers
cv.df1a = filter(cv.df1a, doy<=180, year>=2008)    # 17494
ae.dfa = filter(ae.dfa, doy<=180)    # 1168


# spring beauty -----------------------------------------------------------------
# modeling with L3 ecoregions ------------------------------------------------
t3=as.data.frame(table(cv.df1a[, c("NA_L3CODE", "color2") ]))

t3=t3[t3$Freq>30,]

t3 = mutate(t3, eco3_color2=paste(NA_L3CODE, color2, sep="_"))

t3

cv.df1c=cv.df1a %>% mutate(eco3_color2=paste(NA_L3CODE, color2, sep="_")) %>% 
  filter(eco3_color2 %in% t3$eco3_color2)

unique(cv.df1c$eco3_color2)

hist(cv.df1c$doy)

# check correlations
colnames(cv.df1c)
cor(cv.df1c[, c(13, 1:2, 11, 74:96)])  
cor(cv.df1c[, c(1:2, 11, 76, 80)])  # all corr coef <0.5

# preliminary modeling
cm.wave <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE:color2), 
                data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wave)

cm.wmax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE:color2), 
                data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wmax)

cm.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE:color2), 
                data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wmin)

AICc(cm.wave, cm.wmax, cm.wmin)


# model selection
cm.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                  (1+scale(Tmin_wt_anm)|NA_L3CODE:color2), 
                data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wmin)

cm.wmin1 <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_wt_anm)+scale(PPT_wt_anm)+
                   (1+scale(Tmin_wt_anm)|NA_L3CODE:color2), 
                 data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wmin1)

cm.wmin2 <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_wt_anm)+
                   (1+scale(Tmin_wt_anm)|NA_L3CODE:color2), 
                 data = cv.df1c, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cm.wmin2)

AICc(cm.wmin, cm.wmin1, cm.wmin2)

# cm.wmin is the best model


# spring beuaty miner bee --------------------------------------------------------
# by ecoregion
ta3=as.data.frame(table(ae.dfa$NA_L3CODE))

ae.dfc=ae.dfa %>% filter(NA_L3CODE %in% ta3$Var1[ta3$Freq>=20])

table(ae.dfc$NA_L3CODE)
hist(ae.dfc$doy, breaks=50)

# check correlations
colnames(ae.dfc)
cor(ae.dfc[, c(1:2, 10, 75, 79)])  # all corr coef <0.3
cor(ae.dfc[, c(63:64, 76, 80)])    # all corr coef <0.4


# preliminary modeling with L3 ecoregions
am.wave <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfc, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wave)

am.wmax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfc, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wmax)

am.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfc, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wmin)

AICc(am.wave, am.wmax, am.wmin)


# model selection
am.wmin1 <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                   (1+scale(Tmin_wt_anm)|NA_L3CODE), 
                 data = ae.dfc, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(am.wmin1)

am.wmin2 <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_wt_anm)+scale(PPT_wt_anm)+
                   (1+scale(Tmin_wt_anm)|NA_L3CODE), 
                 data = ae.dfc, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(am.wmin2)

am.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)+
                  (1+scale(Tmin_wt_anm)|NA_L3CODE), 
                data = ae.dfc, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(am.wmin)

AICc(am.wmin1, am.wmin2, am.wmin)

# am.wmin is the best model

#------------------------------------------------------------------------------
# fixed coefs from models
cmsm=summary(cm.wmin)
cmsm$coefficients[2, 1]/sd(cv.df1c$lat) # 5.1
cmsm$coefficients[2, 2]/sd(cv.df1c$lat) # 0.089
cmsm$coefficients[3, 1]/sd(cv.df1c$elevation) # 0.017
cmsm$coefficients[3, 2]/sd(cv.df1c$elevation) # 0.001
cmsm$coefficients[4, 1]/sd(cv.df1c$Tmin_wt_anm) # -1.28
cmsm$coefficients[4, 2]/sd(cv.df1c$Tmin_wt_anm) # 0.17
cmsm$coefficients[5, 1]/sd(cv.df1c$PPT_wt_anm) # 4.68
cmsm$coefficients[5, 2]/sd(cv.df1c$PPT_wt_anm) # 0.47
cmsm$coefficients[6, 1]/sd(cv.df1c$Tmin_wt_anm)/sd(cv.df1c$PPT_wt_anm) # -2.42
cmsm$coefficients[6, 2]/sd(cv.df1c$Tmin_wt_anm)/sd(cv.df1c$PPT_wt_anm) # 0.39

amsm=summary(am.wmin)
amsm$coefficients[2, 1]/sd(ae.dfc$lat) # 4.39
amsm$coefficients[2, 2]/sd(ae.dfc$lat) # 0.54
amsm$coefficients[3, 1]/sd(ae.dfc$elevation) # 0.017
amsm$coefficients[3, 2]/sd(ae.dfc$elevation) # 0.004
amsm$coefficients[4, 1]/sd(ae.dfc$Tave_wt_anm) # -1.80
amsm$coefficients[4, 2]/sd(ae.dfc$Tave_wt_anm) # 0.53



# plot model output
cfcm=coef(cm.wmin)$`NA_L3CODE:color2`
cfcm$NA_L3CODE=substr(rownames(cfcm),1, 5)
cfcm$color=substr(rownames(cfcm),7, length(rownames(cfcm)))
cfcm$NA_L3CODE[6:7]="8.1.10"
cfcm$color[6:7]=c("pink", "white")
cfcm
cm.tm=cv.df1c %>% group_by(NA_L3CODE) %>% summarise(tm=mean(LT_Tmin_wt))
cm.tm$NA_L3CODE=as.character(cm.tm$NA_L3CODE)
cfcm=left_join(cfcm, cm.tm)
cfcm
cfcm$uns_cf=cfcm$`scale(Tmin_wt_anm)`/sd(cv.df1c$Tmin_wt_anm)

cfam=coef(am.wmin)$`NA_L3CODE`
cfam$NA_L3CODE=rownames(cfam)
am.tm=ae.dfc %>% group_by(NA_L3CODE) %>% summarise(tm=mean(LT_Tmin_wt))
am.tm$NA_L3CODE=as.character(am.tm$NA_L3CODE)
cfam=left_join(cfam, am.tm)
cfam
cfam$uns_cf=cfam$`scale(Tmin_wt_anm)`/sd(ae.dfc$Tmin_wt_anm)


# random slopes coefficients
rdmcf=rbind.data.frame(cfcm[, c("uns_cf","tm","NA_L3CODE")], 
                       cfam[, c("uns_cf","tm","NA_L3CODE")])
rdmcf$color=c(cfcm$color, rep("bee", nrow(cfam)))
rdmcf$taxa=c(rep("spring beauty", nrow(cfcm)), rep("miner bee", nrow(cfam)))
rdmcf

table(rdmcf$color)
# bee  pink white 
# 13    38    22 


cor.test(rdmcf$tm[which(rdmcf$color=="white")], 
         rdmcf$uns_cf[which(rdmcf$color=="white")])
# r= -0.44, p=0.04

cor.test(rdmcf$tm[which(rdmcf$color=="pink")], 
         rdmcf$uns_cf[which(rdmcf$color=="pink")])
# r= -0.53, p<0.001

cor.test(rdmcf$tm[which(rdmcf$color=="bee")], 
         rdmcf$uns_cf[which(rdmcf$color=="bee")])
# r=0.76, p=0.003


# Figure 4
ylab=expression("phenological" ~ "sensitivity" ~ "(day" ~ paste("\u00b0C" ^ "-1", ")"))

tiff(filename="figures_rv/pheno_sensitivity_coefs_Tm_L3ecoregions_r1a.tif", width=2000, height=1500, res=300, compression="lzw")
ggplot(rdmcf, aes(y=uns_cf, x=tm, color=color, group=color))+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  scale_color_manual(name="taxa", values=c("orange","coral1","deepskyblue"))+
  labs(x="long-term average winter minimum temperature (\u00b0C)", y= ylab)+
  theme_cowplot()+
  theme(legend.box.margin=margin(10,10,10,10))+
  guides(color=guide_legend(
    keywidth=0.3,
    keyheight=0.3,
    default.unit="inch"))
dev.off()




#---- map coef 
eco3=st_read(dsn="C:/Users/xie412/Documents/Data/shapefiles/ecoregions/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")
eco3=st_transform(eco3, crs=4326)
eco3=st_make_valid(eco3)


eco3_p=left_join(eco3, rdmcf[rdmcf$color=="pink", c(1, 3)])
colnames(eco3_p)[12]="ps_pink"
eco3_p=left_join(eco3_p, rdmcf[rdmcf$color=="white", c(1, 3)])
colnames(eco3_p)[13]="ps_white"
eco3_p=left_join(eco3_p, rdmcf[rdmcf$color=="bee", c(1, 3)])
colnames(eco3_p)[14]="ps_bee"


# Figure 3
p1=ggplot(eco3_p)+
  geom_sf(aes(fill=ps_white), show.legend = F)+
  scale_fill_distiller(palette = "OrRd", limits=c(-4.5, 0.15), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude", y="latitude")+
  theme_cowplot()

p2=ggplot(eco3_p)+
  geom_sf(aes(fill=ps_pink), show.legend = F)+
  scale_fill_distiller(palette = "OrRd", limits=c(-4.5, 0.15), na.value = "grey80")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude")+
  theme_cowplot()

p3=ggplot(eco3_p)+
  geom_sf(aes(fill=ps_bee))+
  scale_fill_distiller(palette = "OrRd", limits=c(-4.5, 0.15), na.value = "grey80", 
                       name="phenological \nsensitivity")+
  coord_sf(xlim=c(-100, -66), ylim=c(25,49))+
  labs(x="longitude")+
  theme_cowplot()+
  theme(legend.position = c(0.75, 0.25))

cfmaps=cowplot::plot_grid(p1, p2, p3, labels=c("a","b","c"), label_size=16, ncol=3, nrow=1, 
                          rel_widths=c(1, 1, 1))

ggsave(filename = "figures_rv/maps_phenological_sensitivity_L3ecoregion.tif", cfmaps, device="tiff",
       dpi=300, height= 5, width=18, compression="lzw")
