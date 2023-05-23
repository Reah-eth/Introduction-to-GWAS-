#########################################################################

library(ggplot2)
library(scatterpie)
library(dplyr)

list.files()
df <- read.table("Phenotypes_ALL_passport_data.txt", sep="\t",header=T)

df_wm <- data.frame(region=unique(df$Region), lat=unique(df$lat),
                    long=unique(df$lon), Virulent=0, Intermediate=0, Avirulent=0)

ph_col <- 7

for (i in 1:dim(df)[1]) {
  idx <- which(df_wm$region==df[i,2])
  if (!is.na(df[i,ph_col])) {
    if (df[i,ph_col]==0) {
      df_wm[idx,"Avirulent"] = df_wm[idx,"Avirulent"] +1
    }
    if (df[i,ph_col]>0 & df[i,ph_col]<1) {
      df_wm[idx,"Intermediate"] = df_wm[idx,"Intermediate"] +1
    }
    if (df[i,ph_col]==1) {
      df_wm[idx,"Virulent"] = df_wm[idx,"Virulent"] +1
    }
  }
}

world <- map_data('world')

df_wm <- mutate(df_wm,radius=rowSums(df_wm[,4:6]))

ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="grey80", color="black") +
  geom_scatterpie(aes(x=long, y=lat, group=region, r=2*log(radius)), data=df_wm,
                  cols=colnames(df_wm)[4:6], alpha=0.7)  +
  coord_cartesian(xlim = c(-120,150), ylim=c(-50,70)) +
  scale_fill_manual(values = c("#d7191c","#ffffbf","#2c7bb6")) +
  theme_void() +
  theme(legend.title = element_blank(),
        text = element_text(size=15))
