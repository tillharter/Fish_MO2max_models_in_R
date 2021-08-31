library(ggplot2)
library(bvpSolve)

#-----------------------------------------------------------------------------------------
#Pcrit model BVP Gdt ExpGrid
#-----------------------------------------------------------------------------------------
# Set working directory
setwd("E:/Uni/Postdoc/Data/GSA vs Hb P50_Matt 2020/Fig. S1 Model validation")
#-----------------------------------------------------------------------------------------

#gill characteristics
A=2.4                                          #gill surface area, cm2 g-1: Hughes 1972
t_ept=4.92                                     #thickness of the epithelium, um, rainbow trout 3 um Perry 1998; 6 um Hughes 1972; Altman and Dittmer 1971
beta_ept=0.00177                               #solubility coefficient of O2 in epithelium; umol cm-3 mmHg-1; Boutilier et al., 1984; from Jensen and Malte Thessis 2017
#D_ept=0.00002                                 #diffusion coefficient of O2 in epithelium; cm2 s-1; Dejours, 1975; from Jensen and Malte Thessis 2017
#D_ept=0.00000108                              #cm2 sec-1; 300 g rainbow trout; Hills and Hughes 1970
#D_ept=0.000011                                #cm2 sec-1; 2 kg anaesthetised dogfish; Hills and Hughes 1970;Piiper and Baumgarten-Schuman 1968
D_ept=0.0000044                                #determined to match the PaO2, PvO2 and MO2 from Kiceniuk and Jones 1977

#gill diffusive conductance
Gd=(D_ept*beta_ept*A)/(t_ept/10000)*60*1000    #umol mmHg-1   min-1  kg-1; max in RT 1.6 umol mmHg-1   min-1  kg-1 Malte and Wang
Gd

# tissue diffusive conductance
Gdt=6.4                                       #umol min -1 kg-1 mmHg-1 Malte and Wang
Pt=0.01                                        #mitochondrial PO2 Malte and Wang

#solubilities @ 10C
alpha_b = 1.99                              #solubility of O2 in human plasma at 10C (Boutillier 1984), umol L-1 mmHg-1  
alpha_w = 2.24                              #solubility of O2 in water at 0 ppt, 10C (Boutillier 1984), umol L-1 mmHg-1
#alpha_w = 1.7943                             #solubility of O2 in water at 35ppt, 10C (Boutillier 1984), umol L-1 mmHg-1

#solubilities @ 20C
#alpha_b = 1.6276                              #solubility of O2 in human plasma at 20C (Boutillier 1984), umol L-1 mmHg-1  
#alpha_w = 1.8230                              #solubility of O2 in water at 0 ppt, 20C (Boutillier 1984), umol L-1 mmHg-1
#alpha_w = 1.4842                             #solubility of O2 in water at 35ppt, 20C (Boutillier 1984), umol L-1 mmHg-1


#cardiovascular characteristics
#Q=Gd/24                                      #umol L-1 mmHg-1 according to Malte and Weber 1987         
Q=53/1000                                     #Cardiac output, to L kg-1 min-1
cHb= 1050                                     #Hemoglobin content, umol L-1
P50=22.2                                      #Hemoglobin P50, mmHg; 24.1 mmHg Tetens and Lykkeboe 1981; 24.8 Rummer and Brauner 2015
n=2.09                                        #Hill coefficient; 2.05 Tetens and Lykkeboe 1981; 1.3 Rummer and Brauner 2015
B=-0.82                                       #Bohr coefficient; Rummer and Brauner 2015
dpH=-0.12                                     #A-V pH shift; Kiceniuk and Jones 1977
P50t=P50*10^(B*dpH)                           #Hb P50 at the tissues


#ventilatory conductance
Vw=Q*42.051-0.52869                          #water flow set to 10*Q; Davis and Cameron 1971: max in RT 1730 mL kg-1 min-1 Malte and Wang 
Gv=Vw*alpha_w                                #umol kg-1 min-1 mmHg-1


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

# guess initial conditions for analysis
iniPvO2=30                                   #guess of initial venous PO2, mmHg
inspPwO2=150                              #target inspired PO2 that will be reached through interation

#-----------------------------------------------------------------------------------------
#analysis with bvp approach

repeat{                                 #this is the outer loop that adjusts venous PO2 until it changes by less than 0.1 mmHg from the previous run

x <- seq(0, 1, by= 0.01)                #length sequence along the gill

parameters <- c(Gd=Gd,	                #umol mmHg-1   min-1  kg-1
                Gv=Gv,                  #umol mmHg-1 min-1 kg-1  
                Q= Q,                   #L kg-1 min-1
                cHb= cHb,               #mmol L-1
                alpha_b = alpha_b,      #mmol L-1 mmHg-1
                alpha_w = alpha_w,      #mmol L-1 mmHg-1
                P50=P50,                #mmHg
                n=n,                    #
                Gdt=Gdt,                # umol mmHg-1   min-1  kg-1
                Pt=Pt,                 # mmHg
                P50t)                  # mmHg
   
gas_x_pre <- function(x, state, parameters){             #functions describing the differential equations to be solved simultaneously (see Malte and Weber 1985)
   with(as.list(c(state, parameters)), {
      dPw=Gd/Gv*(Pw-Pb)
      dPb=Gd/(Q*(alpha_b+cHb*4*(Pb^(n - 1) * n/(Pb^n + 1^n) - Pb^n * (Pb^(n - 1) * n)/(Pb^n +1^n)^2)))*(Pw-Pb)
      return(list(c(dPw, dPb)))
   })
}

#print(system.time(                                    #prints the time required for the analysis

out_pre <- bvptwp(x = x, func = gas_x_pre,
                  yini = c(Pw = NA, Pb = iniPvO2),     #sets the boundary values, i.e. venous PO2 and inspired water PO2
                  yend = c(Pw=inspPwO2, Pb=NA))

xguess <- out_pre[,1]
yguess <- t(out_pre[,c(2,3)])

gas_x <- function(x, state, parameters){             #functions describing the differential equations to be solved simultaneously (see Malte and Weber 1985)
   with(as.list(c(state, parameters)), {
      dPw=Gd/Gv*(Pw-Pb)
      dPb=Gd/(Q*(alpha_b+cHb*4*(Pb^(n - 1) * n/(Pb^n + P50^n) - Pb^n * (Pb^(n - 1) * n)/(Pb^n + P50^n)^2)))*(Pw-Pb)
      return(list(c(dPw, dPb)))
   })
}

out <- bvptwp(x = x, func = gas_x, xguess = xguess, yguess = yguess, 
              yini = c(Pw = NA, Pb = iniPvO2),     #sets the bounday values, i.e. venous PO2 and inspired water PO2
              yend = c(Pw=inspPwO2, Pb=NA), nmax=10000, atol=1e-8)
#))

out <- as.data.frame(out)                             #converts output to data frame
names(out)[names(out)=="x"] <- "l"                    #use length instead of x (preset for bvp function)

#Calculate arterial blood characteristics based on the outcome of gill gas exchange
PaO2<-out$Pb[out$l==1]
SaO2<-(PaO2^n)/(PaO2^n+P50^n)
CaO2<-(alpha_b*PaO2+cHb*4*SaO2)/1000 #mmol L-1


#diffusion of O2 at the tissues

state <- c(Pb = PaO2)                                 #runs the differential equations with PaO2
x <- seq(1, 0, by= -0.01)                              #length sequence along the capillary

gas_xt <- function(x, state, parameters){             #functions describing the differential equations to be solved simultaneously (see Malte and Weber 1985)
   with(as.list(c(state, parameters)), {
      dPb=Gdt/(Q*(alpha_b+cHb*4*(Pb^(n - 1) * n/(Pb^n + P50t^n) - Pb^n * (Pb^(n - 1) * n)/(Pb^n + P50t^n)^2)))*(Pb-Pt)
      return(list(c(dPb)))
   })
}

## Evaluate the system of differential equations with ODE solver (Ordinary differential equations)
out_t <- ode(y = state, times = x, func = gas_xt, parms = parameters, method = rkMethod("rk4"))   #solves the first guess for water and blood PO2 after gas exchange
out_t <- as.data.frame(out_t)                                            #converts output to data frame
names(out_t)[names(out_t)=="time"] <- "l"                                #use length instead of time (preset for ODE function)


#Calculate venous blood characteristics based on the outcome of gill gas exchange

PvO2<-out_t$Pb[out_t$l==0]
SvO2<-(PvO2^n)/(PvO2^n+P50t^n)
CvO2<-(alpha_b*PvO2+cHb*4*SvO2)/1000 #mmol L-1
O2_extr<-(SaO2-SvO2)

#switch P50 back during venous transit to the gill
PvO2_g<-(P50^n*((PaO2^n-(SaO2-SvO2)*(PaO2^n+P50^n))/(P50^n+(SaO2-SvO2)*(PaO2^n+P50^n))))^(1/n)  # From Wilford et al. 1982 


if (abs(PvO2_g-iniPvO2)<0.1){     #if venous PvO2 hasn't changed since the last iteration it breaks the repeat function

   break
}

iniPvO2<-PvO2_g                    #otherwise it uses the new PvO2 from the gas exchange simulation to run the next iteration from "repeat" above
}

# calculate the efficiency of gas exchange according to Malte and Weber Encyclopedia of Fish Physiology chapter
CIO2<-alpha_w*out$Pw[out$l==1]/1000  #mmol l
CEO2<-alpha_w*out$Pw[out$l==0]/1000  #mmol l
EE<-(CIO2-CEO2)/CIO2                        #extraction efficiency
EE

Ew<-(out$Pw[out$l==1]-out$Pw[out$l==0])/(out$Pw[out$l==1]-out$Pb[out$l==0])      #efficiency of gas exchange on the water side
Ew

Eb<-(out$Pb[out$l==1]-out$Pb[out$l==0])/(out$Pw[out$l==1]-out$Pb[out$l==0])      #efficiency of gas exchange on the blood side
Eb

Et<-1-Ew-Eb                                                                      #transfer efficiency across the gill
Et

#Calculates the MO2 that is supported by gas exchange at the gills
 
MO2_w<-Vw*(CIO2-CEO2)*60                        #based on O2 extraction from the water
MO2_w # umol g-1 h-1
MO2_b<-Q*(CaO2-CvO2)*60                         #based on oxygenation of the blood
MO2_b # umol g-1 h-1

VR<-Vw/MO2_b

#MO2 label
MO2_lab<-paste("MO2 = ",signif(MO2_b,3), " umol g-1 h-1")

#axis labels
yaxis_lab<-expression(paste("PO"[2], " (mmHg)"))
xaxis_lab<-expression(paste("Relative length along the gill"))

#plot
g1<-ggplot(data = out, aes(x = l, colour))+
   scale_colour_manual(name="Medium", labels=c("Water", "Blood"), values = c("blue","red"))+
   geom_point(aes(x=l, y=Pw, colour="blue"), size=1)+
   geom_point(aes(x=l, y=Pb, colour= "red"), size=1)+
   scale_y_continuous(name=yaxis_lab, limits=c(0, 160))+
   scale_x_continuous(name = xaxis_lab, limits=c(0, 1))+
   geom_text(label=MO2_lab, x=0.6, y=160, check_overlap = TRUE)+
   theme_bw()+
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.title.y = element_text(family="sans", face="bold", colour="black", size = 10),
         axis.title.x = element_text(family="sans", face="bold", colour="black", size = 10),
         axis.text.x = element_text(family="sans", face="plain", colour="black", size = 8),
         axis.text.y = element_text(family="sans", face="plain", colour="black", size = 8),
         strip.text.x = element_blank(),
         legend.text = element_text(family="sans", face="plain", colour="black", size = 8),
         legend.title = element_text(family="sans", face="bold", colour="black", size = 8),
         legend.key.height=unit(1,"line"),
         legend.key.width=unit(1,"line"),
         legend.position = c(0.1, 0.9),
         legend.margin = NULL,
         axis.line = element_line(colour = "black", size=0.5),
         panel.border = element_blank(),
         panel.background = element_blank())
g1

#export figure
ggsave(filename="Fig. S1_Model validation.tiff", device="tiff", unit="cm", width = 10, height = 8, plot=g1)
ggsave(filename="Fig. S1_Model validation.pdf", device="pdf", unit="cm", width = 10, height = 8, plot=g1)

#summary tables
input<-data.frame(parameters=c("Gill surface area", "Gill thickness", "beta epithelium", "Diffusion coefficient", "Gill diffusive conductance (Gd)", "Tissue diffusive conductance (Gdt)", "mitochondrial PO2","alpha water", "alpha plasma", "Water PO2", "Haemoglobin P50", "Hill coefficient", "Bohr coefficient","A-V pH shift","Haemoglobin", "Cardiac output", "Ventilation volume"),
                  input=c(A, t_ept, beta_ept, D_ept, Gd, Gdt, Pt, alpha_w, alpha_b, inspPwO2, P50, n, B, dpH, cHb/1000, Q*1000, Vw*1000),
                  units=c("cm2 g-1", "um", "umol cm-3 mmHg-1", "cm2 s-1","umol mmHg-1   min-1  kg-1", "umol mmHg-1   min-1  kg-1","mmHg", "umol L-1 mmHg-1",  "umol L-1 mmHg-1", "mmHg",  "mmHg", "","","", "mM", "mL kg-1 min-1", "mL kg-1 min-1"))
output<-data.frame(parameters=c("MO2", "PaO2", "SaO2","CaO2", "PvO2", "SvO2", "CvO2", "expired water PO2", "water O2 extraction efficiency", "blood oxygenation efficiency", "O2 transfer efficiency",  "Tissue O2 extraction"),
                    output=c(signif(MO2_b,4), signif(PaO2,4), signif(SaO2*100,4), signif(CaO2,4), signif(PvO2,4), signif(SvO2*100,4), signif(CvO2,4), signif(out$Pw[out$l==0],4), signif(Ew*100,4), signif(Eb*100,4), signif(EE,4), O2_extr*100),
                    units=c("umol kg-1 min -1", "mmHg", "%", "mmol L-1", "mmHg", "%", "mmol L-1", "mmHg", "%", "%", "", "%"))

# Export tables
write.table(input, sep = ",", quote = FALSE, "Fig. S1_input.txt", row.names = F)
write.table(output, sep = ",", quote = FALSE, "Fig. S1_output.txt", row.names = F)

