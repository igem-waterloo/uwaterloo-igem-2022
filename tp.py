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

#SPECIES EQUATION

#SECOND PEAK OF ENZYME ACTIVTY