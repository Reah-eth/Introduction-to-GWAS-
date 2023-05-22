### read in the libraries
library(ggplot2)
library(ggrepel)
library(patchwork)

# set correct working directory
setwd("~/data/....")

# load GWAS resutls into indivdual vectors
GWAS_pm1 <- read.table("Pm1a/GAPIT.Association.GWAS_Results.MLM.Pm1a.csv", sep = ",", header = T)
GWAS_pm2 <- read.table("Pm2/GAPIT.Association.GWAS_Results.MLM.Pm2.csv", sep = ",", header = T)
GWAS_pm4 <- read.table("Pm4a/GAPIT.Association.GWAS_Results.MLM.Pm4a.csv", sep = ",", header = T)


# calculate the bonferroni correction to set the sig. threshold
alpha = 0.05 # this is the normal significance threshold
test_number = dim(GWAS_pm1)[1] # should be the same for all GWAS in your case
bonf = alpha/test_number


############## global plots ###################################################

manhattan_pm1 <- ggplot(GWAS_pm1, aes(x=Pos/1000000, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm1a-Tester line") +
  geom_point(data=GWAS_pm1[GWAS_pm1$P.value<bonf,], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )

manhattan_pm2 <- ggplot(GWAS_pm2, aes(x=Pos/1000000, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm2-Tester line") +
  geom_point(data=GWAS_pm2[GWAS_pm2$P.value<bonf,], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )

manhattan_pm4 <- ggplot(GWAS_pm4, aes(x=Pos/1000000, y=-log10(P.value), label=SNP)) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm4a-Tester line") +
  geom_point(data=GWAS_pm4[GWAS_pm4$P.value<bonf,], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )



########################### plots for only one chromosome ######################
# make plots for the chromosome that has the best peak

manhattan_pm1_chr8 <- ggplot(GWAS_pm1[GWAS_pm1$Chr=="BGT_CHR-08",], aes(x=Pos/1000000, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm1a-Tester line") +
  geom_point(data=GWAS_pm1[GWAS_pm1$P.value<bonf&GWAS_pm1$Chr=="BGT_CHR-08",], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )

manhattan_pm2_chr7 <- ggplot(GWAS_pm2[GWAS_pm1$Chr=="BGT_CHR-07",], aes(x=Pos/1000000, y=-log10(P.value))) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm2-Tester line") +
  geom_point(data=GWAS_pm2[GWAS_pm2$P.value<bonf&GWAS_pm1$Chr=="BGT_CHR-07",], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )

manhattan_pm4_chr8 <- ggplot(GWAS_pm4[GWAS_pm1$Chr=="BGT_CHR-08",], aes(x=Pos/1000000, y=-log10(P.value), label=SNP)) +
  geom_point(size=0.3, color="Grey40") +
  ggtitle("GWAS using Pm4a-Tester line") +
  geom_point(data=GWAS_pm4[GWAS_pm4$P.value<bonf&GWAS_pm1$Chr=="BGT_CHR-08",], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]") + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15)
  )


########################### zoomed in plots ####################################



manhattan_pm1_zoom <- ggplot(GWAS_pm1[GWAS_pm1$Chr=="BGT_CHR-08",], aes(x=Pos/1000000, y=-log10(P.value), label=SNP)) +
  geom_point(size=0.3, color="Grey40") +
  geom_point(data=GWAS_pm1[GWAS_pm1$P.value<bonf& GWAS_pm1$Chr=="BGT_CHR-08",], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  #geom_text_repel(data=GWAS_pm1[GWAS_pm1$P.value<bonf& GWAS_pm1$Chr=="BGT_CHR-08",], point.padding = 10) +
  scale_x_continuous(name="Position [Mb]", limits = c(10.6,10.9), expand = c(0,0)) + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15),
    axis.title.x = element_blank()
  )

manhattan_pm2_zoom <- ggplot(GWAS_pm2[GWAS_pm2$Chr=="BGT_CHR-07",], aes(x=Pos/1000000, y=-log10(P.value), label=SNP)) +
  geom_point(size=0.3, color="Grey40") +
  geom_point(data=GWAS_pm2[GWAS_pm2$P.value<bonf& GWAS_pm2$Chr=="BGT_CHR-07",], size = 0.6, color="red") +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]", limits = c(0.8,1.2), expand = c(0,0)) + #put your limits here
  #geom_text_repel(data=GWAS_pm2[GWAS_pm2$P.value<bonf& GWAS_pm2$Chr=="BGT_CHR-07",], point.padding = 10) +
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15),
    axis.title.x = element_blank()
  )

manhattan_pm4_zoom <- ggplot(GWAS_pm4[GWAS_pm4$Chr=="BGT_CHR-08",], aes(x=Pos/1000000, y=-log10(P.value), label=SNP)) +
  geom_point(size=0.3, color="Grey40") +
  geom_point(data=GWAS_pm4[GWAS_pm4$P.value<bonf& GWAS_pm4$Chr=="BGT_CHR-08",], size = 0.6, color="red") +
  #geom_text_repel(data=GWAS_pm4[GWAS_pm4$P.value<bonf& GWAS_pm4$Chr=="BGT_CHR-08",], point.padding = 10) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed", size = 0.3) +
  scale_x_continuous(name="Position [Mb]", limits = c(3.5,6), expand = c(0,0)) + #put your limits here
  facet_wrap(~Chr) +
  theme_classic() +
  theme(
    text = element_text(size=15),
    axis.title.x = element_blank()
  )


######################### add effector positions ###############################

eff <- read.table("~/data/admin/teaching/2023_GWAS/dir_BlockCourse_2023/Bgt_effector_genes.gff", sep = ",")

eff_Pm1 <- ggplot(eff[eff$V1=="Bgt_chr-08",], aes(x=V4/1000000,y=1, label=V9)) +
  geom_rect(data=eff[eff$V1=="Bgt_chr-08",], mapping=aes(xmin=V4/1000000,
                                                         xmax=V5/1000000, 
                                                         ymin=1, ymax=2,
                                                         fill=T ), 
            color="black",fill="forestgreen") +
  
  scale_x_continuous(name="Position [Mb]", limits = c(10.6,10.9), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1,2)) +
  geom_text_repel(nudge_y = -5, nudge_x = -0.01) +
  facet_wrap(~V1) +
  theme_void() +
  theme(strip.text = element_blank(),
)

eff_Pm2 <- ggplot(eff[eff$V1=="Bgt_chr-07",], aes(x=V4/1000000,y=1, label=V9)) +
  geom_rect(data=eff[eff$V1=="Bgt_chr-07",], mapping=aes(xmin=V4/1000000,
                                                         xmax=V5/1000000, 
                                                         ymin=1, ymax=2,
                                                         fill=T ), 
            color="black",fill="forestgreen") +
  
  scale_x_continuous(name="Position [Mb]", limits = c(0.8,1.2), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1,2)) +
  #geom_text_repel(nudge_y = -5, nudge_x = -0.01) +
  geom_text_repel() +
  facet_wrap(~V1) +
  theme_void() +
  theme(strip.text = element_blank())



eff_Pm4 <- ggplot(eff[eff$V1=="Bgt_chr-08",], aes(x=V4/1000000,y=1, label=V9)) +
  geom_rect(data=eff[eff$V1=="Bgt_chr-08",], mapping=aes(xmin=V4/1000000,
                                                         xmax=V5/1000000, 
                                                         ymin=1, ymax=2,
                                                         fill=T ), 
            color="black",fill="forestgreen") +
  
  scale_x_continuous(name="Position [Mb]", limits = c(3.5,6), expand = c(0,0)) +
  scale_y_continuous(limits = c(-0.1,2)) +
  geom_text_repel(nudge_y = -5, nudge_x = -0.01) +
  facet_wrap(~V1) +
  theme_void() +
  theme(strip.text = element_blank())
  

###########################  make panels for zoom in plots #####################


layout <- "
A
A
A
A
A
A
B

"

final_zoom_in_Pm1 <- manhattan_pm1_zoom + eff_Pm1 + plot_layout(design = layout)
final_zoom_in_Pm2 <- manhattan_pm2_zoom + eff_Pm2 + plot_layout(design = layout)
final_zoom_in_Pm4 <- manhattan_pm4_zoom + eff_Pm4 + plot_layout(design = layout)




