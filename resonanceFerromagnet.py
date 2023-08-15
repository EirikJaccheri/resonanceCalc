import numpy as np
import matplotlib.pyplot as plt

SIZE = 20

omegaq = 1
h = 1.
c = 5.7e-7
delta = 5e-6 * 1j
omegaLim = 0.001

#miunsCol = "#ff7f0e"
#plusCol = "#1f77b4"
minusCol = "blue"
plusCol = "red"

def det(omega,s):
    return -omega - delta + omegaq - c * h / (s * (omega + delta) + h)

def omega(s,pm):
    return omegaq / 2 - s * h / 2  + pm * s / 2 * np.sqrt((s* omegaq + h)**2 - 4*c*h*s)

def print_analytic():
    for S in [-1,1]:
        for PM in [-1,1]:
            print("S:",S,"PM:",PM,"omega:",omega(S,PM))


def plot_resonance():
    omegaL = np.linspace((1 - omegaLim)*omegaq,(1 + omegaLim)*omegaq,500)
    fig, ax = plt.subplots()
    for S in [-1,1]:
        for PM in [-1,1]:
            if S == -1:
                Label = "analytic -"
                Col = minusCol
            else:
                Label = "analytic +"
                Col = plusCol
            if (omega(S,PM) > 0):
                ax.axvline(x=omega(S,PM),color=Col,alpha=0.6, dashes=(3, 5))

    ax.plot(omegaL, np.abs(delta / det(omegaL,1)),label="L",color=plusCol)
    ax.plot(omegaL, np.abs(delta / det(omegaL,-1)),label="R",color=minusCol,linestyle="dashed")
    ax.legend(fontsize = SIZE)
    ax.set_xlabel("$\hbar\omega [\hbar \omega_0] $",size=SIZE)
    ax.set_ylabel("$\delta \cdot |\mathcal{G}(\omega)|$",size=SIZE)
    plt.xticks(fontsize=SIZE)
    plt.yticks(fontsize=SIZE)
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.xaxis.set_major_formatter('{:.3f}'.format)
    
    
    plt.tight_layout()
    plt.savefig(f"plots/resonanceFerromagnetC{c}.pdf")
    plt.show()

print_analytic()
plot_resonance()
