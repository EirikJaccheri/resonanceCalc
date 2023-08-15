import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import pickle

SIZE = 20

N= 500


#miunsCol = "#ff7f0e"
#plusCol = "#1f77b4"
minusCol = "Blue"
plusCol = "red"

def omega(s,pm,omegaq,h,c):
    return omegaq / 2 - s * h / 2  + pm * s / 2 * np.sqrt((s* omegaq + h)**2 - 4*c*h*s)

def Edown(eps,DeltaDown):
    return np.sqrt(eps**2 + DeltaDown**2)

def Eup(eps,DeltaUp,h):
    return np.sqrt((eps - h)**2 + DeltaUp**2)

#for now i take the low temperature limit
def nf(eps):
    return (1 - np.heaviside(eps, 0.5))

#legg til s
def arg(eps,omega, s, h, DeltaDown, DeltaUp):
    return (
    -2 * (-Edown(eps,DeltaDown) - eps) * (-s * omega - Edown(eps,DeltaDown) - eps + h) * nf(-Edown(eps,DeltaDown)) 
    / (Edown(eps,DeltaDown) * (-s * omega - Edown(eps,DeltaDown) - Eup(eps,DeltaUp,h)) * (-s * omega - Edown(eps,DeltaDown) + Eup(eps,DeltaUp,h)))  

    #er + h over brøkstrek riktig
    -2 * (s * omega - Eup(eps,DeltaUp,h) - eps) * (-Eup(eps,DeltaUp,h) - eps + h) * nf(-Eup(eps,DeltaUp,h))
    / (Eup(eps,DeltaUp,h) * (s * omega - Eup(eps,DeltaUp,h) - Edown(eps,DeltaDown)) * (s * omega - Eup(eps,DeltaUp,h) + Edown(eps,DeltaDown)))
    )

#skal det være - mu i denne defnisjinen
def ReAnCont(eps,omega, s, h, delta, DeltaDown, DeltaUp):
    return np.real(arg(eps,omega + delta, s, h, DeltaDown, DeltaUp))

#skal det være - mu i denne defnisjinen
def ImAnCont(eps,omega, s, h, delta, DeltaDown, DeltaUp):
    return np.imag(arg(eps,omega + delta, s, h, DeltaDown, DeltaUp))

def det(omegaq,omega,s, h, c, mu, delta, DeltaDown, DeltaUp):
    Re = quad(lambda eps : ReAnCont(eps,omega, s, h, delta, DeltaDown, DeltaUp),-mu,np.inf)
    Im = quad(lambda eps : ImAnCont(eps,omega, s, h, delta, DeltaDown, DeltaUp),-mu,np.inf)
    return omegaq - omega - delta + c / 4 * (Re[0] + 1j * Im[0])

#regner ut resonanse numerisk of lagrer listen som en pickle i pickleDir
def calculate_resonance(pickleDir,omegaq,h,c,mu, delta, DeltaDown, DeltaUp, DelOmega):
    omegaL = np.linspace((1 - DelOmega)*omegaq, (1 + DelOmega)*omegaq, N)
    detLplus = np.zeros(N) * 1j
    detLminus = np.zeros(N) * 1j
    for i in range(len(omegaL)):
        detLplus[i] = np.abs(delta / det(omegaq,omegaL[i],1, h, c, mu, delta, DeltaDown, DeltaUp))
        detLminus[i] = np.abs(delta / det(omegaq,omegaL[i],-1, h, c, mu, delta, DeltaDown, DeltaUp))
    pickleL = [omegaL,detLplus,detLminus]
    with open(pickleDir,"wb") as f:
        pickle.dump(pickleL,f)

#lager resonans plot for data som den finner i pickleDir
def plot_resonance(plotName, pickleDir, omegaq, h, c):
    with open(pickleDir,'rb') as p:
        [omegaL,detLplus,detLminus] = pickle.load(p)
    fig, ax = plt.subplots()
    ax.plot(omegaL, detLplus,label="L",color=plusCol)
    ax.plot(omegaL, detLminus,label="R",color=minusCol,linestyle="dashed")
    
    #plotting analytic solution
    for S in [-1,1]:
        for PM in [-1,1]:
            if S == -1:
                Label = "analytic -"
                Col = minusCol
            else:
                Label = "analytic +"
                Col = plusCol
            if (omega(S,PM,omegaq,h,c) > 0):
                ax.axvline(x=omega(S,PM,omegaq,h,c), color=Col, linestyle="--",alpha=0.6,dashes=(3,5))
            print("S:",S,"PM:",PM,"omega:",omega(S,PM,omegaq,h,c))

    ax.legend(fontsize = SIZE)
    ax.set_xlabel("$\hbar\omega [\hbar \omega_0]$",size=SIZE)
    ax.set_ylabel("$\delta \cdot |\mathcal{G}(\omega)|$",size=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.xaxis.set_major_formatter('{:.3f}'.format)

    plt.tight_layout()
    plt.savefig(f"plots/{plotName}.pdf")
    plt.show()

def plot_bands(h, DeltaUp,DeltaDown,plotName):
    k_list = np.linspace(-2,2,100)
    
    fig, ax = plt.subplots()
    ax.plot(k_list,  Edown(k_list**2,DeltaDown))
    ax.plot(k_list,  Eup(k_list**2,DeltaUp,h))
    ax.plot(k_list,  -Edown(k_list**2,DeltaDown))
    ax.plot(k_list,  -Eup(k_list**2,DeltaUp,h))

    ax.set_xlabel("$k$",size=SIZE)
    ax.set_ylabel("energy",size=SIZE)
    plt.savefig(f"plots/gap{plotName}.pdf")
    plt.show()

DELOMEGA = 0.001 #how far from the resonance to plot (percentage of resonance frequency)

DELTAUP = 0.0
DELTADOWN = 0.0
OMEGAQ = 1
H = 1
C = 5.7e-7
DELTA = 5e-6 * 1j
MU = 0.3
PICKLE_DIR = "pickles/pickle"
PLOT_NAME = f"resonanceFerSCDelDown{str(DELTADOWN)}DelUp{str(DELTAUP)}OmegaQ{str(OMEGAQ)}DelOmega{str(DELOMEGA)}"
#calculate_resonance(PICKLE_DIR,OMEGAQ,H,C,MU,DELTA,DELTAUP,DELTADOWN,DELOMEGA)
plot_resonance(PLOT_NAME, PICKLE_DIR,OMEGAQ,H,C)
#plot_bands(1, DELTAUP, DELTADOWN, PLOT_NAME)

DELTAUP = 0.
DELTADOWN = 0.05 
#OMEGAQ = np.sqrt(1 + DELTADOWN**2)
OMEGAQ = 1
H = 1
C = 5.7e-7
DELTA = 5e-6 * 1j
MU = 0.3
PICKLE_DIR = "pickles/pickleGapedDown"
PLOT_NAME = f"resonanceFerSCDelDown{str(DELTADOWN)}DelUp{str(DELTAUP)}OmegaQ{str(OMEGAQ)}DelOmega{str(DELOMEGA)}"
#calculate_resonance(PICKLE_DIR,OMEGAQ,H,C,MU,DELTA,DELTAUP,DELTADOWN,DELOMEGA)
plot_resonance(PLOT_NAME, PICKLE_DIR,OMEGAQ,H,C)
#plot_bands(1, DELTAUP, DELTADOWN, PLOT_NAME)

DELTAUP = 0.05
DELTADOWN = 0.
#OMEGAQ = np.sqrt(1 + DELTAUP**2)
OMEGAQ = 1
H = 1
C = 5.7e-7
DELTA = 5e-6 * 1j
MU = 0.3
PICKLE_DIR = "pickles/pickleGapedUp"
PLOT_NAME = f"resonanceFerSCDelDown{str(DELTADOWN)}DelUp{str(DELTAUP)}OmegaQ{str(OMEGAQ)}DelOmega{str(DELOMEGA)}"
#calculate_resonance(PICKLE_DIR,OMEGAQ,H,C,MU,DELTA,DELTAUP,DELTADOWN,DELOMEGA)
plot_resonance(PLOT_NAME, PICKLE_DIR,OMEGAQ,H,C)
#plot_bands(1, DELTAUP, DELTADOWN, PLOT_NAME)

DELTAUP = 0.05
DELTADOWN = 0.05
#OMEGAQ = np.sqrt(1 + DELTAUP**2)
OMEGAQ = 1
H = 1
C = 5.7e-7
DELTA = 5e-6 * 1j
MU = 0.3
PICKLE_DIR = "pickles/pickleGapedUp"
PLOT_NAME = f"resonanceFerSCDelDown{str(DELTADOWN)}DelUp{str(DELTAUP)}OmegaQ{str(OMEGAQ)}DelOmega{str(DELOMEGA)}"
#calculate_resonance(PICKLE_DIR,OMEGAQ,H,C,MU,DELTA,DELTAUP,DELTADOWN,DELOMEGA)
plot_resonance(PLOT_NAME, PICKLE_DIR,OMEGAQ,H,C)
#plot_bands(1, DELTAUP, DELTADOWN, PLOT_NAME)
