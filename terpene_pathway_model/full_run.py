from mevalonate_pathway.full_dynamics import full_dynamics as MVA
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from terpene_pathway_model.terpene_precursors.tp import terpene_precursors


def run_sans_acetyl_coa(acetyl_coa):
    y0 = np.ones(8, dtype=float) * 1e-6
    y0[0] = 1e-3
    y0[1] = acetyl_coa
    T, Y, k, ep, pyruvate, Tmax = MVA(y0)

    # Y is concentrations post MVA of substrates in order:
    # ['Pyruvate', 'Acetyl CoA', 'Acetoacetyl CoA',  'HMG CoA', 'Mevalonate', 'Mevalonate-phosphate','Mevalonate diphosphate','Isopentyl diphosphate']

    x = np.ones(22, dtype=float) * 1e-6
    #[GAP]=x[0]
    # D-Glyceraldehyde 3-Phosphate (G3P) 🚨
    #[Pyr]=x[1]          Pyruvate 🚨
    #[DOXP]=x[2]         1-Deoxy-D-xylulose 5-phosphate (DXS)🚨
    #[ME4P]=x[3]         2-C-Methyl-D-erythritol-4-phosphate
    #[CDPME]=x[4]        4-(Cytidine 5'-diphospho)-2-C-methyl-D-erythritol 🚨 (CDP-ME is an enzyme in MEP)
    #[CDPME2P]=x[5]      2-Phospho-4-(cytidine 5'-diphospho)-2-C-methyl-D-erythritol
    #[MEcPP]=x[6]        2-C-Methyl-D-erythritol-2,4-cyclodiphosphate (MEcPP) 🚨
    #[HMBPP]=x[7]        1-Hydroxy-2-methyl-2-(E)-butenyl 4-diphosphate(HMB-PP)🚨
    #[DMAPP]=x[8]        Dimethylallyl-pyrophosphate (DMAPP) 🚨
    #[IPP]=x[9]         Isopentenyl diphosphate (IPP) 🚨
    #[GPP]=x[10]         Geranyl diphosphate (Geranyl diphosphate) 🚨
    #[LIM]=x[11]         (-)-Limonene
    #[IPPol]=x[12]       (-)-trans-Isopiperitenol
    #[IPPone]=x[13]      (-)-Isopiperitenone
    #[CIPUL]=x[14]       (+)-cis-Isopulegone
    #[PUL]=x[15]         (+)-Pulegone
    #[MF]=x[16]          (+)-Menthofuran
    #[IMone]=x[17]       (+)-Isomenthone
    #[Mone]=x[18]        (-)-Menthone
    #[NMol]=x[19]        (+)-Neomenthol
    #[Mol]=x[20]         (-)-Menthol
    #[IMol]=x[21]        (+)-Isomenthol
    #[NIMol]=x[22]       (+)-Neoisomenthol

    x[1] = Y[0][-1]
    x[9] = Y[7][-1]

    sol = solve_ivp(terpene_precursors,[0, 10], x, method='Radau')
    plt.plot(sol.t, np.transpose(sol.y))
    plt.show()


run_sans_acetyl_coa(0)
