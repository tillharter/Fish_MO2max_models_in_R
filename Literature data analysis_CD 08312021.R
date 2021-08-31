require(lubridate)
require(readxl)
require(ggplot2)
require(ape)
require(phytools)
require(nlme)
require(geiger)
require(ggsignif)
require(scales)
require(cowplot)
require(scales)



#########################
### GILL SURFACE AREA ###
#########################

setwd("/Users/au231308/Dropbox/Projects/SWFW/")

# Import data
df<-read_xlsx("./data.xlsx",sheet = "GSA")
df<-as.data.frame(df)

df<-subset(x = df,select = c("Species","Salinity","Mass (g)","SA (cm^2)","SA/g (cm2/g)"))
colnames(df)<-c("sp","sal","bm","sa","sag") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore



# Plot raw mass-specific metabolic rate without correcting for allometric scaling and phylogenetic non-independence
df.0<-
  cbind(
    aggregate(sag~sp,df,mean),
    aggregate(sal~sp,df,unique)[,2]
  )

colnames(df.0)[3]<-"sal"

p.value<-t.test(sag~sal,df.0)$p.value



p0<-
  ggplot(
    data = df.0,
    mapping = aes(x=sal,y=sag,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("Mass-specific gill surface area (mm"^"2"*"g"^"-1"*")"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*df.0$sag), vjust=0,annotation=ifelse(test = p.value<0.001, no = paste("P = ",signif(p.value,3),sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black");p0




df<-df[!is.na(df$bm),] # remove studies that do not report body mass

# in species that do not report total gill surface area, calculate total surface area from body mass and mass specific surface area
for (i in which(is.na(df$sa))){
  df$sa[i]<-df$sag[i]*df$bm[i]
}

# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Ambloplites_rupestri")]<-"Ambloplites_rupestris"
df$sp[which(df$sp=="Amerurus_nebulosis")]<-"Ameiurus_nebulosus"
df$sp[which(df$sp=="Anguilla_vulgaris")]<-"Anguilla_anguilla"
df$sp[which(df$sp=="Botia_ohachata")]<-"Botia_lohachata"
df$sp[which(df$sp=="Catostoumus_commersonii")]<-"Catostomus_commersonii"
df$sp[which(df$sp=="Channa_punctatus")]<-"Channa_punctata"
df$sp[which(df$sp=="Clarias_mossambicus")]<-"Clarias_gariepinus"
df$sp[which(df$sp=="Comephorus_baicalensis")]<-"Comephorus_baikalensis"
df$sp[which(df$sp=="Gnathonemus_victoriae")]<-"Marcusenius_victoriae"
df$sp[which(df$sp=="Hypostomus_plecostomus.")]<-"Hypostomus_plecostomus"
df$sp[which(df$sp=="Ictalurus_nebulosus")]<-"Ameiurus_nebulosus"
df$sp[which(df$sp=="Noemacheilus_barbatulus")]<-"Barbatula_barbatula"
df$sp[which(df$sp=="Notropis_cornutus)")]<-"Luxilus_cornutus"
df$sp[which(df$sp=="Oncorynchus_mykiss")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Oncorhynchus_mykis")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Oreochromis_alcalicus_grahami")]<-"Alcolapia_grahami"
df$sp[which(df$sp=="Prognathochromis_venator")]<-"Haplochromis_venator"
df$sp[which(df$sp=="Stizostedion_vitreum")]<-"Sander_vitreus"
df$sp[which(df$sp=="Acanthocottus_scorpius")]<-"Myoxocephalus_scorpius"
df$sp[which(df$sp=="Anoplomoma_fimbria")]<-"Anoplopoma_fimbria"
df$sp[which(df$sp=="Archosargus_probaiocephalus")]<-"Archosargus_probatocephalus"
df$sp[which(df$sp=="Centropristis_striatus")]<-"Centropristis_striata"
df$sp[which(df$sp=="Chilomycterus_schoepfi")]<-"Chilomycterus_schoepfii"
df$sp[which(df$sp=="Clinottus_globiceps")]<-"Clinocottus_globiceps"
df$sp[which(df$sp=="Dicentrarchuslabrax")]<-"Dicentrarchus_labrax"
df$sp[which(df$sp=="Errex_zachirus")]<-"Glyptocephalus_zachirus"
df$sp[which(df$sp=="Gadus_virens")]<-"Pollachius_virens"
df$sp[which(df$sp=="Gymnosarda_alleterata")]<-"Euthynnus_alletteratus"
df$sp[which(df$sp=="Leptocephalus_conger")]<-"Conger_conger"
df$sp[which(df$sp=="Lophopsetta_maculata")]<-"Mancopsetta_maculata"
df$sp[which(df$sp=="Palinurichthyes_perciformis")]<-"Hyperoglyphe_perciformis"
df$sp[which(df$sp=="Peprilus_alepidatus")]<-"Peprilus_paru"
df$sp[which(df$sp=="Poecilialatipinna")]<-"Poecilia_latipinna"
df$sp[which(df$sp=="Pomatomus_saliatrix")]<-"Pomatomus_saltatrix"
df$sp[which(df$sp=="Poronotus_triacanthus")]<-"Peprilus_triacanthus"
df$sp[which(df$sp=="Prionotus_strigatus")]<-"Prionotus_evolans"
df$sp[which(df$sp=="Roccus_lineatus")]<-"Morone_saxatilis"
df$sp[which(df$sp=="Sarda_chiliensis")]<-"Sarda_chiliensis_chiliensis"
df$sp[which(df$sp=="Spheroides_maculatus")]<-"Sphoeroides_maculatus"
df$sp[which(df$sp=="Tautoga_onitus")]<-"Tautoga_onitis"

# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

# Prune tree to only represent species in the data set
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

# Set up glmm model
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)

glmm.sa.0<-
  MCMCglmm::MCMCglmm(
    log10(sa)~log10(bm),
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)
summary(glmm.sa.0)

# extract intercept and slope
int<-summary(glmm.sa.0)$solutions[1,1];int
slo<-summary(glmm.sa.0)$solutions[2,1];slo


# Calculate residuals
df$resid<-log10(df$sa)-(log10(df$bm)*slo+int)

# Genereate named vector for the residuals
resid<-aggregate(resid~sp,df,median)
resid<-setNames(resid$resid,resid$sp)
resid

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal

#phyloaov
paov<-phylANOVA(tree = tree,sal,resid,nsim = 50000)
phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]
col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]


# Set up axis text format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Plot gill surface area vs body mass scatter plot
p1<-ggplot(data = df, mapping = aes(x = bm, y = sa, col = sal))+
  scale_x_continuous(trans = "log10",labels = comma,limits = c(1,200000))+
  scale_y_continuous(trans = "log10",labels = comma)+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = c(0.35, 0.90),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"),
    legend.background = element_rect(fill=alpha('blue', 0))
    )+
  geom_abline(slope = slo,intercept = int)+
  geom_point(shape = 16, size = 1)+
  scale_color_manual(values = c("#006D2C","#08519C"),breaks = c("FW","SW"),labels = c("Freshwater","Seawater"))+
  
  
  labs(col = "Salinity",
       x = expression("Body mass (g)"), 
       y = expression("Gill surface area (mm"^"2"*")"));p1


## Jitter plot of residuals
p2<-
  ggplot(
    data = data.frame(
      resid = resid,
      sal = sal
    ),
    mapping = aes(x=sal,y=resid,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("residual log"[10]*"[gill surface area (mm"^"2"*")]"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*resid), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black");p2

#pdf("./Figures/GSA_1.pdf",width = 2.25,height = 4.5,useDingbats = F)
#plot_grid(p1,p2,nrow = 2,ncol =1,align = "v",labels = "AUTO",label_size = 12)
#dev.off()

#pdf("./Figures/GSA_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#plotTree.wBars(
#  tree = tree,
#  x = resid,
#  tip.labels = T,method = "plotTree",
#  col = col,
# args.plotTree=list(font = 1,edge.width = 0.5), 
# args.barplot=list(),
#  scale = 90)
#dev.off()

#dev.off()


## Model estimate of gill surface area
# specify body mass in g
bodymass <- 100

# absolute gill surface area
10^(log10(bodymass)*slo+int)

# mass-specific gill surface area
10^(log10(bodymass)*slo+int)/bodymass




##################################
### HEMOGLOBIN OXYGEN AFFINITY ###
##################################

# Import data
df<-read_xlsx("./data.xlsx",sheet = "P50")
df<-as.data.frame(df)

df<-subset(x = df,select = c("Species","Salinity","P50 (mmHg)","Temp"))
colnames(df)<-c("sp","sal","p50","temp") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore



# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Aguilla_rostrata")]<-"Anguilla_rostrata"
df$sp[which(df$sp=="Ambloplites_rupestri")]<-"Ambloplites_rupestris"
df$sp[which(df$sp=="Salmo_gairdneri")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Arius_leptaspis")]<-"Neoarius_leptaspis"
df$sp[which(df$sp=="Fundulus_heteroclitus")]<-"Fundulus_heteroclitus_heteroclitus"
df$sp[which(df$sp=="Pangasius")]<-"Pangasianodon_hypophthalmus"
df$sp[which(df$sp=="Salmo_trutta_fario")]<-"Salmo_trutta"
df$sp[which(df$sp=="Euthynnus_affinish")]<-"Euthynnus_affinis"
df$sp[which(df$sp=="F._malcomi")]<-"Forsterygion_malcolmi"
df$sp[which(df$sp=="F._varium")]<-"Forsterygion_varium"
df$sp[which(df$sp=="G._capito")]<-"Forsterygion_capito"
df$sp[which(df$sp=="Gadus_morhua;")]<-"Gadus_morhua"
df$sp[which(df$sp=="Pleuronectes_platessa;")]<-"Pleuronectes_platessa"
df$sp[which(df$sp=="Scomber_japonicus)")]<-"Scomber_japonicus"
df$sp[which(df$sp=="Leucoraja_ocellata)")]<-"Scomber_japonicus"
df$sp[which(df$sp=="Helcogramma_medium")]<-"Bellapiscis_medius"


# I have removed these species, as they are not teleosts
df <- df[which(df$sp!="Potamotrygon_motoro"),] 
df <- df[which(df$sp!="Leucoraja_ocellata"),] 
df <- df[which(df$sp!="Myliobatis_californica"),] 

# I chose a random species from these genera
df$sp[which(df$sp=="Hypostomus_sp")]<-"Hypostomus_agna" 

# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

df<-df[!is.na(df$temp),] # remove studies that do not report temperature
df$inv.temp<-1/(273+df$temp) # calculate 1/(abs Temp)

# Prune tree to only represent species in the data set
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)


## Does P50 depend on temperature?


# Set up glmm model
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)

glmm<-
  MCMCglmm::MCMCglmm(
    log10(p50)~inv.temp,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)
summary(glmm) 

# Extract solutions
int<-summary(glmm)$solutions[1,1];int
slo<-summary(glmm)$solutions[2,1];slo

# Calculate residuals
df$resid<-log10(df$p50)-(df$inv.temp*slo+int)


# Plot data
p3<-ggplot(data = df, mapping = aes(x = inv.temp, y = p50, col = sal))+
    scale_x_continuous(breaks = c(0.0033,0.0035,0.0037),
    sec.axis = sec_axis(
      trans = ~ .^-1 - 273.15,
      name = "Temperature (°C)"))+
  scale_y_continuous(trans = "log10",labels = comma)+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = "none",
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 8),
    #legend.position = c(0.7, 0.90),
    #legend.spacing.y = unit(0.01,"cm"),
    #legend.spacing.x = unit(0.001,"cm"),
    #legend.key.height = unit(0.01,"cm")
    )+
  geom_abline(slope = slo, intercept = int)+
  geom_point(shape = 16, size = 1)+
  scale_color_manual(values = c("#006D2C","#08519C"),breaks = c("FW","SW"),labels = c("Freshwater","Seawater"))+
  labs(col = "Salinity",
       x = expression("1/T"), 
       y = expression("P"[50]*" (mmHg)"));p3

# Genereate named vector for p50-residuals
resid<-aggregate(resid~sp,df,median)
resid<-setNames(resid$resid,resid$sp)
resid

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal

#phyloaov
paov<-phylANOVA(tree = tree,sal,resid,nsim = 10000)
paov
phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]
col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]


## Jitter plot of P50-residuals

p4<-
  ggplot(
    data = data.frame(
      resid = resid,
      sal = sal
    ),
    mapping = aes(x=sal,y=resid,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("residual log"[10]*"[P"[50]*" (mmHg)]"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*resid), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") ;p4


#pdf("./Figures/P50_1.pdf",width = 2.25,height = 4.5,useDingbats = F)
#plot_grid(p3,p4,nrow = 2,ncol =1,align = "v",labels = "AUTO",label_size = 12)
#dev.off()

#pdf("./Figures/P50_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#plotTree.wBars(
#  tree = tree,
#  x = resid,
#  tip.labels = T,method = "plotTree",
#  col = col,
#  args.plotTree=list(font = 1,edge.width = 0.5), 
#  args.barplot=list(),
#  scale = 50)
#dev.off()

#dev.off()





# Model estimates for p50
FWSW.logP50.effects<-
  aggregate(
    resid~sal,
    data.frame(resid = resid,sal = sal),
    mean)

temperature = 10

# logP50 and P50 without correcting for FW/SW
1/(273.15+temperature)*slo+int
10^(1/(273.15+temperature)*slo+int)


# logP50 and P50 in FW
1/(273.15+temperature)*slo+int+FWSW.logP50.effects[1,2]
10^(1/(273.15+temperature)*slo+int+FWSW.logP50.effects[1,2])

# logP50 and P50 in SW
1/(273.15+temperature)*slo+int+FWSW.logP50.effects[2,2]
10^(1/(273.15+temperature)*slo+int+FWSW.logP50.effects[2,2])


###############
### MO2 MAX ###
###############

# Import data
df<-read_xlsx("./data.xlsx",sheet = "MO2max")
df<-as.data.frame(df)

df<-subset(x = df,select = c("Species","Salinity","Body mass (g)","Uncorrected MO2max (mg O2/h)","MO2max (umol O2/h*g)","MO2 temp"))
colnames(df)<-c("sp","sal","bm","mo2","mo2g","temp") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore

#convert from mg to µmol
df$mo2<-df$mo2/32*1000

df<-df[!is.na(df$bm),] # remove studies that do not report body mass
df

# in species that do not report total gill surface area, calculate total surface area from body mass and mass specific surface area
for (i in which(is.na(df$mo2))){
  df$mo2[i]<-df$mo2g[i]*df$bm[i]
}

# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Ambloplites_rupestri")]<-"Ambloplites_rupestris"
df$sp[which(df$sp=="Ictalurus_nebulosus")]<-"Ameiurus_nebulosus"
df$sp[which(df$sp=="Aristichthys_nobilis")]<-"Hypophthalmichthys_nobilis"
df$sp[which(df$sp=="Gasterosteus_aculeautus")]<-"Gasterosteus_aculeatus"
df$sp[which(df$sp=="Hypophthalmichthys_sp")]<-"Hypophthalmichthys_molitrix"
df$sp[which(df$sp=="Maccullochella_peelii_peelii")]<-"Maccullochella_peelii"
df$sp[which(df$sp=="Onychostoma_sp")]<-"Onychostoma_barbatulum"
df$sp[which(df$sp=="Pelteobagrus_vachelli")]<-"Pseudobagrus_vachellii"
df$sp[which(df$sp=="Salmo_gairdneri")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Salvelinus_alpinus")]<-"Salvelinus_alpinus_alpinus"
df$sp[which(df$sp=="Chromis_veridis")]<-"Chromis_viridis"
df$sp[which(df$sp=="Leiostomus_xanthrus")]<-"Leiostomus_xanthurus"
df$sp[which(df$sp=="Macrozoraces_americanus")]<-"Zoarces_americanus"
df$sp[which(df$sp=="Neopomacentrus_benkieri")]<-"Neopomacentrus_bankieri"
df$sp[which(df$sp=="Pagrus_auratus")]<-"Sparus_aurata"
df$sp[which(df$sp=="Pomacentrus_lipedogenys")]<-"Pomacentrus_lepidogenys"
df$sp[which(df$sp=="Fundulus_heteroclitus")]<-"Fundulus_heteroclitus_heteroclitus"

# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

# Prune tree to only represent species in the data set
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

# Set up glmm model
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)


# Model selection

# Step 1: most complex model
# 1.1: BM and TEMP
glmm1.1<-
  MCMCglmm::MCMCglmm(
    log10(mo2)~log10(bm)+temp,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)


# Step 2: remove one variable
# Remove temp
glmm2.1<-
  MCMCglmm::MCMCglmm(
    log10(mo2)~log10(bm),
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)


# Remove temp
glmm2.2<-
  MCMCglmm::MCMCglmm(
    log10(mo2)~temp,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)

summary(glmm1.1) # DIC = -95
summary(glmm2.1) # DIC = -77 - model 1.1 is better
summary(glmm2.2) # DIC = 230 - model 1.1 is better

asdf<-summary(glmm1.1)
asdf$solutions
# Conclusion on model selection: go ahead with model 1.1


# extract intercept and slope
summary(glmm1.1)$solutions

int<-summary(glmm1.1)$solutions[1,1];int
slo1<-summary(glmm1.1)$solutions[2,1];slo1
slo2<-summary(glmm1.1)$solutions[3,1];slo2



# Calculate residuals
df$resid<-log10(df$mo2)-(log10(df$bm)*slo1+df$temp*slo2+int)

# Genereate named vector for the residuals
resid<-aggregate(resid~sp,df,median)
resid<-setNames(resid$resid,resid$sp)
resid

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal

#phyloaov
paov<-phylANOVA(tree = tree,sal,resid,nsim = 50000)
paov
phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]
col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]

# Set up axis text format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Plot bm and mo2max scatter plot
p5<-ggplot(data = df, mapping = aes(x = bm, y = mo2, col = sal))+
  #scale_y_continuous(trans = "log10",label=scientific_10)+
  #scale_x_continuous(trans = "log10",label=scientific_10)+
  scale_x_continuous(trans = "log10",labels = comma,limits = c(1,200000))+
  scale_y_continuous(trans = "log10",labels = comma)+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = c(0.3, 0.90),
    legend.spacing.y = unit(0.01,"cm"),
    legend.spacing.x = unit(0.001,"cm"),
    legend.key.height = unit(0.01,"cm"),
    legend.background = element_rect(fill=alpha('blue', 0)))+
  
#  geom_abline(slope = slo1,intercept = int)+
  geom_point(shape = 16, size = 1)+
  scale_color_manual(values = c("#006D2C","#08519C"),breaks = c("FW","SW"),labels = c("Freshwater","Seawater"))+
  
  
  labs(col = "Salinity",
       x = expression("Body mass (g)"), 
       y = expression("MO"["2,max"]*" (µmol h"^"-1"*")"));p5



## Jitter plot of residuals
p6<-
  ggplot(
    data = data.frame(
      resid = resid,
      sal = sal
    ),
    mapping = aes(x=sal,y=resid,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("residual log"[10]*"[MO"["2,max"]*" (µmol h"^"-1"*")]"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*resid), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") 
p6

#pdf("./Figures/MO2_1.pdf",width = 2.25,height = 4.5,useDingbats = F)
#plot_grid(p5,p6,nrow = 2,ncol =1,align = "v",labels = "AUTO",label_size = 12)
#dev.off()

#pdf("./Figures/MO2_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#plotTree.wBars(
#  tree = tree,
#  x = resid,
#  tip.labels = T,method = "plotTree",
#  col = col,
#  args.plotTree=list(font = 1,edge.width = 0.5), 
#  args.barplot=list(),
#  scale = 90)
#dev.off()

#dev.off()



# Model estimates for MO2,max
FWSW.MO2.effects<-
  aggregate(
    resid~sal,
    data.frame(resid = resid,sal = sal),
    mean)
FWSW.MO2.effects
temperature = 10
bodymass = 100

# absolute and mass-specific MO2,max without correcting for FW/SW
10^(int+log10(bodymass)*slo1+temperature*slo2)
10^(int+log10(bodymass)*slo1+temperature*slo2)/bodymass


# absolute and mass-specific MO2,max in FW
10^(int+log10(bodymass)*slo1+temperature*slo2+FWSW.MO2.effects[1,2])
10^(int+log10(bodymass)*slo1+temperature*slo2+FWSW.MO2.effects[1,2])/bodymass

# absolute and mass-specific MO2,max in SW
10^(int+log10(bodymass)*slo1+temperature*slo2+FWSW.MO2.effects[2,2])
10^(int+log10(bodymass)*slo1+temperature*slo2+FWSW.MO2.effects[2,2])/bodymass





###############
### PCRIT ###
###############

# Import data
df<-read_xlsx("./data.xlsx",sheet = "Pcrit")
df<-as.data.frame(df)

df<-subset(x = df,select = c("Species","Salinity","Pcrit (mmHg)","Temp C","Body mass (g)"))
colnames(df)<-c("sp","sal","pcrit","temp","bm") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore


# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Ambloplites_rupestri")]<-"Ambloplites_rupestris"
df$sp[which(df$sp=="Aristichthys_nobilis")]<-"Hypophthalmichthys_nobilis"
df$sp[which(df$sp=="Carassius_auratus_grandoculis")]<-"Carassius_auratus"
df$sp[which(df$sp=="Ctenopharyngodon_piceus")]<-"Ctenopharyngodon_idella"
df$sp[which(df$sp=="Cyprinodon_ariegatus")]<-"Cyprinodon_variegatus_variegatus"
df$sp[which(df$sp=="Cyprinodon_variegatus")]<-"Cyprinodon_variegatus_variegatus"
df$sp[which(df$sp=="Apogon_compressus")]<-"Ostorhinchus_compressus"
df$sp[which(df$sp=="Apogon_cyanosoma")]<-"Ostorhinchus_cyanosoma"
df$sp[which(df$sp=="Apogon_exostigma")]<-"Pristiapogon_exostigma"
df$sp[which(df$sp=="Apogon_fragilis")]<-"Zoramia_fragilis"
df$sp[which(df$sp=="Apogon_leptacanthus")]<-"Zoramia_leptacantha"
df$sp[which(df$sp=="Asterropteryx_semipunctatus")]<-"Asterropteryx_semipunctata"
df$sp[which(df$sp=="Chaetopsylla_globiceps\r\n")]<-"Clinocottus_globiceps"
df$sp[which(df$sp=="Fundulus_heteroclitus")]<-"Fundulus_heteroclitus_heteroclitus"
df$sp[which(df$sp=="Ostorhinchus_doederleini")]<-"Apogon_doederleini"
df$sp[which(df$sp=="Paragobiodon_xanthosomus")]<-"Paragobiodon_xanthosoma"
df$sp[which(df$sp=="Ctenopharyngodon_idellus")]<-"Ctenopharyngodon_idella"
df$sp[which(df$sp=="Amblygobius_rainfordi")]<-"Koumansetta_rainfordi"
df$sp[which(df$sp=="Helcogramma_medium")]<-"Bellapiscis_medius"
df$sp[which(df$sp=="Lumpeninae_lumpenus")]<-"Lumpenus_lampretaeformis"
df$sp[which(df$sp=="Pagrus_auratus")]<-"Sparus_aurata"

# I chose a random species from these genera
df$sp[which(df$sp=="Pagothenia")]<-"Pagothenia_borchgrevinki" 
df$sp[which(df$sp=="Forsterygion_Sp.")]<-"Forsterygion_capito" 

# I have removed this species, as I cannot find it in the phylogeny
df <- df[which(df$sp!="Gobiodon_erythrospilus"),] 

# I have removed these species, as they are reported both as SW and FW in the data set
df <- df[which(df$sp!="Cottus_asper"),] 
df <- df[which(df$sp!="Cyprinodon_variegatus_variegatus"),] 
df <- df[which(df$sp!="Leptocottus_armatus"),] 
df <- df[which(df$sp!="Oncorhynchus_mykiss"),] 

# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

# Step 1: Evaluate if pcrit depends on body mass
df.na<-na.omit(df)
tree<-keep.tip(tree,df.na$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)

glmm.bm<-
  MCMCglmm::MCMCglmm(
    pcrit~bm,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df.na,
    nitt=50000,
    burnin=1000,
    thin=500)
summary(glmm.bm) 
# conclusion: pcrit does not depend on body mass. P = 0.918. Include studies that do not report body mass


# Step 2: Evaluate if pcrit depends on temperature
tree<-read.tree(file = "./mcc.nexus")
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)
glmm.temp<-
  MCMCglmm::MCMCglmm(
    pcrit~temp,
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)
summary(glmm.temp) 
# conclusion: no effect of temperature on pcrit (P = 0.102). Analyze all pcrit without taking temperature into account. 


# Step 3: Evaluate if pcrit depend on salinity

# Genereate named vector for the residuals
pcrit<-aggregate(pcrit~sp,df,median)
pcrit<-setNames(pcrit$pcrit,pcrit$sp)
pcrit

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal

#phyloaov
paov<-phylANOVA(tree = tree,sal,pcrit,nsim = 50000)
paov
phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]

col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]


## Jitter plot of Pcrit

p7<-
  ggplot(
    data = data.frame(
      pcrit = pcrit,
      sal = sal
    ),
    mapping = aes(x=sal,y=pcrit,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("P"[crit]*" (mmHg)"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*pcrit), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") 
p7

#pdf("./Figures/Pcrit_1.pdf",width = 2.25,height = 2.25,useDingbats = F)
#p7
#dev.off()

#pdf("./Figures/Pcrit_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#plotTree.wBars(
#  tree = tree,
#  x = pcrit,
#  tip.labels = T,method = "plotTree",
#  col = col,
#  args.plotTree=list(font = 1,edge.width = 0.5), 
#  args.barplot=list(),
#  scale = 5)
#dev.off()

#dev.off()


##########
### Gd ###
##########

# Import data
df<-read_xlsx("./data.xlsx",sheet = "Thickness")
df<-as.data.frame(df)
head(df)
df<-subset(x = df,select = c("Species","Salinity","Body mass (g) for GSA","Gd (umol mmHg-1   min-1  kg-1)"))
colnames(df)<-c("sp","sal","bm","gd") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore


df<-df[!is.na(df$bm),] # remove studies that do not report body mass
df

# Calculate absolute Gd
df$gd<-df$gd*df$bm
head(df)

# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Salmo_gairdneri")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Ictalurus_nebulosus")]<-"Ameiurus_nebulosus"
df$sp[which(df$sp=="Scornber_scombrus")]<-"Scomber_scombrus"
df$sp[which(df$sp=="Katsuuionis_pelamis")]<-"Katsuwonus_pelamis"
df$sp[which(df$sp=="Notropis_cornutus")]<-"Luxilus_cornutus"


# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

dim(df)

# Prune tree to only represent species in the data set
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

# Set up glmm model
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
tree$node.label<-seq(1,length(tree$node.label),1)
inv.phylo<-MCMCglmm::inverseA(tree,nodes="TIPS",scale=TRUE)


glmm1.1<-
  MCMCglmm::MCMCglmm(
    log10(gd)~log10(bm),
    random=~sp,
    family="gaussian",
    ginverse=list(sp=inv.phylo$Ainv),
    prior=prior,
    data=df,
    nitt=50000,
    burnin=1000,
    thin=500)
summary(glmm1.1) # There is a significant positive effect of body mass on Gd


# extract intercept and slope
summary(glmm1.1)$solutions
int<-summary(glmm1.1)$solutions[1,1];int
slo<-summary(glmm1.1)$solutions[2,1];slo


# Calculate residuals
df$resid<-log10(df$gd)-(log10(df$bm)*slo+int)

# Genereate named vector for the residuals
resid<-aggregate(resid~sp,df,median)
resid<-setNames(resid$resid,resid$sp)
resid

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal

#phyloaov
paov<-phylANOVA(tree = tree,sal,resid,nsim = 50000)
paov
phylosig(tree = tree,x = resid,method = "lambda",test = T)


# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]
col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]

# Set up axis text format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Plot bm and Gd scatter plot
p8<-ggplot(data = df, mapping = aes(x = bm, y = gd, col = sal))+
  #scale_y_continuous(trans = "log10",label=scientific_10)+
  #scale_x_continuous(trans = "log10",label=scientific_10)+
  scale_x_continuous(trans = "log10",labels = comma)+
  scale_y_continuous(trans = "log10",labels = comma)+
  theme_classic()+
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.position = "none",
    #legend.title = element_text(size = 8),
    #legend.position = c(0.3, 0.90),
    #legend.spacing.y = unit(0.01,"cm"),
    #legend.spacing.x = unit(0.001,"cm"),
    #legend.key.height = unit(0.01,"cm")
    )+
  
  #  geom_abline(slope = slo1,intercept = int)+
  geom_point(shape = 16, size = 1)+
  scale_color_manual(values = c("#006D2C","#08519C"),breaks = c("FW","SW"),labels = c("Freshwater","Seawater"))+
  geom_abline(slope = slo,intercept = int)+
  
  labs(col = "Salinity",
       x = expression("Body mass (g)"), 
       y = expression("G"[d]*" (µmol mmHg"^"-1"*" min"^"-1"*")"));p8



## Jitter plot of residuals
p9<-
  ggplot(
    data = data.frame(
      resid = resid,
      sal = sal
    ),
    mapping = aes(x=sal,y=resid,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        #legend.text = element_text(size = 8),
        #legend.title = element_text(size = 8),
        )+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("residual log"[10]*"[G"[d]*" (µmol mmHg"^"-1"*" min"^"-1"*")]"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*resid), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") 
p9

#pdf("./Figures/Gd_1.pdf",width = 2.25,height = 4.5,useDingbats = F)
#plot_grid(p8,p9,nrow = 2,ncol =1,align = "v",labels = "AUTO",label_size = 12)
#dev.off()

#pdf("./Figures/Gd_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#lotTree.wBars(
# tree = tree,
#  x = resid,
# tip.labels = T,method = "plotTree",
# col = col,
#  args.plotTree=list(font = 1,edge.width = 0.5), 
#  args.barplot=list(),
#  scale = 90)
#dev.off()

#dev.off()



# Model estimates for Gd
FWSW.Gd.effects<-
  aggregate(
    resid~sal,
    data.frame(resid = resid,sal = sal),
    mean)
FWSW.Gd.effects

bodymass = 100

# absolute and mass-specific Gd without correcting for FW/SW
10^(int+log10(bodymass)*slo)
10^(int+log10(bodymass)*slo)/bodymass


# absolute and mass-specific Gd in FW
10^(int+log10(bodymass)*slo+FWSW.Gd.effects[1,2])
10^(int+log10(bodymass)*slo+FWSW.Gd.effects[1,2])/bodymass

# absolute and mass-specific Gd in SW
10^(int+log10(bodymass)*slo+FWSW.Gd.effects[2,2])
10^(int+log10(bodymass)*slo+FWSW.Gd.effects[2,2])/bodymass








#################
### Thickness ###
#################

# Import data
df<-read_xlsx("./data.xlsx",sheet = "Thickness")
df<-as.data.frame(df)

df<-subset(x = df,select = c("Species","Salinity","Water-blood thickness (um)"))
colnames(df)<-c("sp","sal","thick") # abbreviated column names
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate genus and species by underscore
df$sp             <-sub("\\ ", "_", df$sp)                      # Separate species and subspecies by underscore

# THere was not effect of log10(bm) (p=0.08) or bm (p=0.102) on thickness, so all thickness data from species without body mass reported was included in the analysis

# Import tree
tree<-read.tree(file = "./mcc.nexus")

# Correct recent changes in species names
df$sp[which(df$sp=="Salmo_gairdneri")]<-"Oncorhynchus_mykiss"
df$sp[which(df$sp=="Ictalurus_nebulosus")]<-"Ameiurus_nebulosus"
df$sp[which(df$sp=="Scornber_scombrus")]<-"Scomber_scombrus"
df$sp[which(df$sp=="Channa_punctatus")]<-"Channa_punctata"
df$sp[which(df$sp=="Clarias_mossambicus")]<-"Clarias_gariepinus"
df$sp[which(df$sp=="Euthynnus_afinis")]<-"Euthynnus_affinis"
df$sp[which(df$sp=="Microstornus_kitt")]<-"Microstomus_kitt"
df$sp[which(df$sp=="Solea_variegata")]<-"Microchirus_variegatus"
df$sp[which(df$sp=="Tetrapturus_audax")]<-"Kajikia_audax"
df$sp[which(df$sp=="Katsuuionis_pelamis")]<-"Katsuwonus_pelamis"
df$sp[which(df$sp=="Notropis_cornutus")]<-"Luxilus_cornutus"

# Remove Chaenocephalus_aceratus
df <- df[which(df$sp!="Chaenocephalus_aceratus"),]

# Has any species changed name?
df$sp[is.na(match(df$sp,tree$tip.label))]

# Prune tree to only represent species in the data set
tree<-keep.tip(tree,df$sp)
tree<-ladderize(tree)
tree<-force.ultrametric(tree)

# Genereate named vector for the residuals
thick<-aggregate(thick~sp,df,median)
thick<-setNames(thick$thick,thick$sp)
thick

# Genereate named vector for salinity
sal<-aggregate(sal~sp,df,unique)
sal<-setNames(sal$sal,sal$sp)
sal


#phyloaov
paov<-phylANOVA(tree = tree,sal,thick,nsim = 50000)
paov
phylosig(tree = tree,x = thick,method = "lambda",test = T)

# Set colors
col = to.matrix(sal[tree$tip.label],seq=sort(unique(sal)))[,1]
col[which(col == 0)]<-"#08519C"
col[which(col == 1)]<-"#006D2C"
col<-col[tree$tip.label]




## Jitter plot of residuals
p10<-
  ggplot(
    data = data.frame(
      thick = thick,
      sal = sal
    ),
    mapping = aes(x=sal,y=thick,color=sal)
  )+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  geom_jitter(shape =16, width = 0.3)+
  scale_color_manual(values = c("#006D2C","#08519C"))+
  scale_x_discrete(labels = c("Freshwater","Seawater"))+
  labs(x = expression("Salinity"),
       y = expression("Thickness (µm)"))+
  geom_signif(comparisons=list(c("FW", "SW")),y_position = max(1.08*thick), vjust=0,annotation=ifelse(test = paov$Pf<0.001, no = paste("P = ",paov$Pf,sep=""), yes = paste("P < 0.001")),size = 0.2,textsize = 6/3,color="black") 
p10

#pdf("./Figures/Thick_1.pdf",width = 2.25,height = 2.25,useDingbats = F)
#p10
#dev.off()

#pdf("./Figures/Thick_2.pdf",width = 10,height = 10,pointsize = 8,useDingbats = F)
#plotTree.wBars(
#  tree = tree,
#  x = thick,
#  tip.labels = T,method = "plotTree",
#  col = col,
#  args.plotTree=list(font = 1,edge.width = 0.5), 
#  args.barplot=list(),
#  scale = 90)
#dev.off()

#dev.off()



pdf("./Figures/Gill_morphology.pdf",width = 18/2.54,height = 18/2.54*2/3,useDingbats = F)
plot_grid(p10,p1,p8,p0,p2,p9,nrow = 2,ncol = 3,align = "hv",labels = c("A","B","C","D","E","F"),label_size = 12,hjust = -3)
dev.off()

pdf("./Figures/Metabolism_P50_new data.pdf",width = 18/2.54,height = 18/2.54*2/3,useDingbats = F)
plot_grid(p5,p3,NULL,p6,p4,p7,nrow = 2,ncol = 3,align = "hv",labels = c("A","B","","C","D","E"),label_size = 12,hjust = -3)
dev.off()

