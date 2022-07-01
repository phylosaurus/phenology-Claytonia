# This script is for comparing phenology among Claytonia color morphs and Andrena occurrence by 5*5km grids


library(tidyverse)
library(sf)
library(rstatix)
library(cowplot)

# read data
cv.df = read_csv(file="data_processed/Claytonia_virginica_iNat_cleaned_data_all.csv")
ae.df = read_csv(file="data_processed/Andrena_erigeniae_occurence_cleaned_data_all_anm_20220606.csv")



# reduce duplicated data at the same location and on the same day and remove outliers
ae <- ae.df %>% select(-gbifID, -datasetName) %>% distinct() %>% 
  select(1:12)

ae <- ae %>% filter(doy<=180)
summary(ae)

cv <- cv.df %>% select(-gbifID) %>% distinct() %>% 
  select(1:13, 102)

cv <- cv %>% filter(doy<=180, year>=2008)


# create 5*5km grids in the geographic range and compare phenology at grid scale
# convert WGS84 coordinates to equal area Albers crs 

# cv %>% select(lat, lon) %>% distinct()

cv_sf = st_as_sf(cv, coords = c("lon","lat"), crs=4326)

cv_sf = st_transform(cv_sf, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")

cv_sf$geometry

st_crs(cv_sf)

st_bbox(cv_sf)
# xmin      ymin      xmax      ymax 
# -299128.9 -881537.4 2067041.5 1172938.7


library(raster)

(2070000-(-300000))/5000
(1175000-(-885000))/5000
r=raster(ncol=474, nrow=412, xmn=-300000, xmx=2070000, ymn=-885000, ymx=1175000)
res(r)
crs(r)= "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
r[]=1:ncell(r)

plot(r)

r_sp = as(r, "SpatialPolygonsDataFrame")
r_sf = st_as_sf(r_sp)

ggplot(r_sf)+
  geom_sf(color="black")+
  geom_sf(data=cv.df1_sf,color="blue")

cv_joined = st_join(r_sf, cv_sf)
cv_jdf=st_set_geometry(cv_joined, NULL)
summary(cv_jdf)

# calculate grid-average date for all data
cv_jdf_sm=na.omit(cv_jdf) %>% group_by(layer, color2) %>% 
  summarize(doy.mean= mean(doy), doy.sd=sd(doy), 
            doy.min=min(doy), doy.max=max(doy), n=n()) %>% 
  filter(n>=3) # select grids with at least 3 observations
summary(cv_jdf_sm)


ggplot(cv_jdf_sm, aes(y=doy.mean, x=factor(color2), group=factor(color2)))+
  geom_boxplot()  

cv_jdf_sm1 = cv_jdf_sm %>% pivot_wider(-n, names_from = color2, 
                                       values_from = c(doy.mean, doy.sd, doy.min, doy.max))
summary(cv_jdf_sm1)

# quantify the difference between white and pink flowering dates
lm(I(doy.mean_pink-1*doy.mean_white)~1, data=cv_jdf_sm1)
# Coefficients:
#   (Intercept)  
# 3.143


cv_jdf_sm1$diff_pw=cv_jdf_sm1$doy.mean_pink-cv_jdf_sm1$doy.mean_white
summary(cv_jdf_sm1$diff_pw)
cv_jdf_sm1 %>% filter(! is.na(diff_pw))

# aggregate miner bees to 5*5km grids
ae_sf = st_as_sf(ae, coords = c("lon","lat"), crs=4326)

ae_sf = st_transform(ae_sf, crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m no_defs")

ae_joined = st_join(r_sf, ae_sf)
ae_jdf=st_set_geometry(ae_joined, NULL)
ae_jdf_sm=na.omit(ae_jdf) %>% group_by(layer) %>% 
  summarize(doy.mean= mean(doy), doy.sd=sd(doy), 
            doy.min=min(doy), doy.max=max(doy), n=n()) %>% 
  filter(n>=3)

summary(ae_jdf_sm)

cv_jdf_sm = cv_jdf_sm %>% 
  select(layer, color2, doy.mean, doy.sd) %>% 
  rename(mdoy.cv=doy.mean, sdoy.cv=doy.sd)

df_joined=left_join(cv_jdf_sm, ae_jdf_sm)


cv_jdf_sm1=left_join(cv_jdf_sm1, ae_jdf_sm)

# calculate the average difference between flowering mean date and bee occurrence mean date
summary(cv_jdf_sm1$doy.mean-cv_jdf_sm1$doy.mean_pink) # mean 5.772
summary(cv_jdf_sm1$doy.mean-cv_jdf_sm1$doy.mean_white) # mean 9.319


# Figure 2 -----------------------------------------------------------------------
p1=ggplot(cv_jdf_sm1, aes(x=doy.mean_white, y=doy.mean_pink))+
  geom_point(aes(x=doy.mean_white, y=doy.mean_pink)) +
  geom_abline(intercept=0, slope=1, linetype=2)+
  geom_abline(intercept=3.1, slope=1, col="red", size=1.2)+
  xlim(50,140) + ylim(50, 140)+
  labs(x="white flowering date (day of year)", y="pink flowering date (day of year)")+
  theme_cowplot(16)

colors=c("pink flower"="coral1", "white flower"="deepskyblue")

p2=ggplot(filter(cv_jdf_sm1, n>=3))+
  geom_point(aes(x=doy.mean, y=doy.mean_pink, color="pink flower")) +
  geom_point(aes(x=doy.mean, y=doy.mean_white, color="white flower")) +
  scale_color_manual(name="", values=colors) +
  geom_abline(intercept=0, slope=1, linetype=2)+
  geom_abline(intercept=-9.3, slope=1, col="deepskyblue", size=1.2)+
  geom_abline(intercept=-5.8, slope=1, col="coral1", size=1.2)+
  xlim(70,130) + ylim(70, 130)+
  labs(x="bee occurrence date (day of year)", y="flowering date (day of year)")+
  theme_cowplot(16)+
  theme(legend.position = c(0.1, 0.9))

plots=cowplot::plot_grid(p1, p2, labels=c("a","b"), label_size=16, ncol=2, nrow=1, 
                         rel_widths=c(1, 1))

ggsave(filename = "figures/Figure2.tif", plots, device="tiff",
       dpi=300, height= 5, width=10, compression="lzw")
