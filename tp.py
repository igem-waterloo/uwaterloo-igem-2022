import math

#METABOLITE NOMENCLATURE

#[GAP]=x(1)          D-Glyceraldehyde 3-Phosphate
#[Pyr]=x(2)          Pyruvate
#[DOXP]=x(3)         1-Deoxy-D-xylulose 5-phosphate
#[ME4P]=x(4)         2-C-Methyl-D-erythritol-4-phosphate
#[CDPME]=x(5)        4-(Cytidine 5'-diphospho)-2-C-methyl-D-erythritol
#[CDPME2P]=x(6)      2-Phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol
#[MEcPP]=x(7)        2-C-Methyl-D-erythritol-2,4-cyclodiphosphate
#[HMBPP]=x(8)        1-Hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate
#[DMAPP]=x(9)        Dimethylallyl-pyrophosphate
#[IPP]=x(10)         Isopentenyl diphosphate
#[GPP]=x(11)         Geranyl diphosphate
#[LIM]=x(12)         (-)-Limonene
#[IPPol]=x(13)       (-)-trans-Isopiperitenol
#[IPPone]=x(14)      (-)-Isopiperitenone
#[CIPUL]=x(15)       (+)-cis-Isopulegone
#[PUL]=x(16)         (+)-Pulegone
#[MF]=x(17)          (+)-Menthofuran
#[IMone]=x(18)       (+)-Isomenthone
#[Mone]=x(19)        (-)-Menthone 
#[NMol]=x(20)        (+)-Neomenthol
#[Mol]=x(21)         (-)-Menthol
#[IMol]=x(22)        (+)-Isomenthol
#[NIMol]=x(23)       (+)-Neoisomenthol


#KINETIC PARAMETERS

#kc units: [1/s]    (kc = Kcat)
#KM units: [uM]
#Ki units: [uM]

KM1a = 68           #1-Deoxy-D-xylulose-5-phosphate synthase (DXS) for GAP 
kc1a = 1.9 
KM1b = 440          #1-Deoxy-D-xylulose-5-phosphate synthase (DXS) for Pyr
kc1b = 1.9
Kia = 16            #Dissociation constant for Pyr
KM2f = 132          #1-Deoxy-D-xylulose-5-phosphate reductoisomerase (DXR; forward reaction) 
kc2f = 4.4
KM2r = 972          #1-Deoxy-D-xylulose-5-phosphate reductoisomerase (DXR; reverse reaction) 
kc2r = 1.6
KM3 = 500           #2-C-Methyl-D-erythritol 4-phosphate cytidylyltransferase (MCT)
kc3 = 26 
KM4 = 100           #4-(Cytidine 5'-diphospho)-2-C-methyl-D-erythritol kinase (CMK)    
kc4 = 1
KM5 = 252           #2-C-Methyl-D-erythritol 2,4-cyclodiphosphate synthase (MECPS)
kc5 = 3.4 
KM6 = 420           #4-Hydroxy-3-methylbut-2-en-1-yl diphosphate synthase (HDS)
kc6 = 0.4
KM7 = 30            #4-Hydroxy-3-methylbut-2-en-1-yl diphosphate reductase (HDR)
kc7 = 3.7
KM8f = 5.1          #Isopentenyl-diphosphate delta-isomerase for IPP (IPPI; forward reaction)
kc8f =0.018
KM8r = 17           #Isopentenyl-diphosphate delta-isomerase for DMAPP (IPPI; rev reaction)
kc8r = 0.89
KM9a = 54           #Geranyl diphosphate synthase (GPPS; DMAPP as substrate)
kc9a = 48
KM9b = 26           #Geranyl diphosphate synthase (GPPS; IPP as substrate)
kc9b = 48
KM10 = 20           #(-)-Limonene synthase (LS)
kc10 = 0.3
KM11 = 18           #(-)-Limonene 3-hyroxylase (L3H)
kc11 = 1.8 
KM12 = 72           #(-)-trans-Isopiperitenol dehydrogenase (IsoDH)
kc12 = 0.002
KM13 = 1            #(-)-Isopiperitenone reductase (IsoR)
kc13 = 1.3
KM14 = 270          #(+)-cis-Isopulegone isomerase (IsoI)
kc14 = 2.5
KM15 = 30           #(+)-Menthofuran synthase (MFS)
kc15 = 2.0 
KM16a = 2.3         #(+)-Pulegone reductase (PR; product: (-)-menthone)
kc16a = 1.8 
KM16b = 2.3         #(+)-Pulegone reductase (PR; product: (+)-isomenthone)
kc16b = 1.8 
KM17a = 3           #(-)-Menthone:(-)-menthol reductase (MMR; substrate: (-)-menthone) 
kc17a = 0.6
KM17b = 41          #(-)-Menthone:(-)-menthol reductase (MMR; substrate: (+)-isomenthone)
kc17b = 0.6
KM18af  = 674       #(-)-Menthone:(+)-neomenthol reductase (MNR; substrate: (-)-menthone) forward reaction)
kc18af = 0.06
KM18ar = 1200       #(-)-Menthone:(+)-neomenthol reductase (MNR; substrate: (-)-menthone) backward reaction)
kc18ar  = 0.06      #estimated
KM18b = 1000        #(-)-Menthone:(+)-neomenthol reductase (MNR; substrate: (+)-isomenthone)
kc18b = 0.06
Kic1=96             #Product inhibition constant (Geranyl diphosphate acting on IPPI)
Kic2=300            #Product inhibition constant ((+)-menthofuran acting on PR)
                    #Competitive inhibition mechanism
Kis=112             #Substrate Inhibition constant ((+)-pulegone acting on PR)
                    #Uncompetitive inhibition mechanism
z=100               #Factor to account for the actual concentration of (+)-menthofuran in secretory cells of glandular trichomes
w=0.05              #Factor to account for the actual concentration of (+)-pulegone in secretory cells of glandular trichomes

#FIRST PEAK OF ENZYME ACTIVITY
#f(x) = Comp * a * exp((-(t-b).^2)/(2*(c)^2))

#where   Comp = Factor to adjust for the volume density of the compartment in which a 
                #particular enzyme is active [Dimensionless]
#        a    = Factor defining the height of the Gaussian peak for enzyme activity
                #[µM]
#       t    = Time [s]
#       b    = Factor defining the position of the center of the Gaussian peak for
                #enzyme activity [s]
#        c    = Factor defining the width of the Gaussian peak for enzyme activity at
                #half maximum [s]


b1=1296000;   #Defines the position of the center of the Gaussian peak for enzyme
              # activity.
              #% Relevant to the following enzyme activities: LS, L3H, IsoDH, IsoR, IsoI, 
              #MFS, PR

c1=800000;    # Defines the width of the Gaussian peak for enzyme activity at half 
            # maximum.    
              #Relevant to the following enzyme activities: LS, L3H, IsoDH, IsoR, IsoI, 
               # MFS, PR
  
b5=1800000;   #Defines the position of the center of the Gaussian peak for enzyme 
            # activity.
              # Relevant to the following enzyme activities: MMR, MNR

c5=900000;    # Defines the width of the Gaussian peak for enzyme activity at half
            # maximum.     
              # Relevant to the following enzyme activities: MMR, MNR     
              # 
E1=(0.139)*0.03*math.exp((-(t-b1)^2)/(2*(c1)^2))       # DXS
E2=(0.139)*0.0225*math.exp((-(t-b1)^2)/(2*(c1)^2))    # DXR
E3=(0.139)*0.5*math.exp((-(t-b1)^2)/(2*(c1)^2))     # MCT 
E4=(0.139)*0.0225*math.exp((-(t-b1)^2)/(2*(c1)^2))     # CMK 
E5=(0.139)*0.5*math.exp((-(t-b1)^2)/(2*(c1)^2))        # MECPS 
E6=(0.139)*0.5*math.exp((-(t-b1)^2)/(2*(c1)^2))       # HDS  
E7a=(0.139)*0.2*math.exp((-(t-b1)^2)/(2*(c1)^2))      # HDR (product: DMAPP)
E7b=(0.139)*0.04*math.exp((-(t-b1)^2)/(2*(c1)^2))     # HDR (product: IPP)
E8=(0.139)*0.3*math.exp((-(t-b1)^2)/(2*(c1)^2))       # IPPI 
E9=(0.139)*0.1*math.exp((-(t-b1)^2)/(2*(c1)^2))       # GPPS
E10= (0.139)*0.017*math.exp((-(t-b1)^2)/(2*(c1)^2))   # LS    
E11= (0.365)*0.003*math.exp((-(t-b1)^2)/(2*(c1)^2))   # L3H 
E12= (0.044)*10*math.exp((-(t-b1)^2)/(2*(c1)^2))     # IsoDH  
E13= (0.204)*0.34*math.exp((-(t-b1)^2)/((2*c1)^2))     # IsoR  
E14= (0.204)*0.34*math.exp((-(t-b1)^2)/((2*c1)^2))    # IsoI  
E15= (0.365)*0.00007*math.exp((-(t-b1)^2)/(2*(c1)^2))  # MFS  
E16a=(0.204)*0.0015*math.exp((-(t-b1)^2)/(2*(c1)^2))  # PR  (product: (-)-menthone)
E16b=(0.204)*0.00015*math.exp((-(t-b1)^2)/(2*(c1)^2))  # PR  (product: (+)-isomenthone)
E17a=(0.204)*0.0011*math.exp((-(t-b5)^2)/(2*(c5)^2))  # MMR (product: (-)-menthol)
E17b=(0.204)*0.0011*math.exp((-(t-b5)^2)/(2*(c5)^2))  # MMR (product: (+)-neoisomenthol)
E18a=(0.204)*0.00001*math.exp((-(t-b5)^2)/(2*(c5)^2))  # MNR (product: (+)-neomenthol)
E18b=(0.204)*0.00001*math.exp((-(t-b5)^2)/(2*(c5)^2)) # MNR (product: (+)-isomenthol)
                  
# The model also takes into account that the glandular trichome density (GN) changes over time.  This behavior is approximated using a logistic function:

c=5*10^5;    #  parameter approximating slope of exponential phase of sigmoid curve
k=1/8*10^4;  #  parameter approximating shape of sigmoid curve

GN = 1+ 1/(1+c*math.exp(-k*t)); # at day 15, gland number is 86.7 % of total gland number at
                           # day 30



#SPECIES EQUATION

#SECOND PEAK OF ENZYME ACTIVTY