'''Determinando a função de Hubble em diferentes cenários'''

from numcosmo_py import Nc, Ncm
import matplotlib.pyplot as plt
import numpy as np
Ncm.cfg_init()

def Hubble() -> None:

    cosmo1 = Nc.HICosmoDECpl(massnu_length=0)
    cosmo2 = Nc.HICosmoDEXcdm(massnu_length=0)
    cosmo3 = Nc.HICosmoLCDM()

    '''Modelo w_0,w_a'''

    cosmo1.props.H0 = 74.5
    cosmo1.props.Omegab = 0.05
    cosmo1.props.Omegac = 0.26
    cosmo1.props.Omegax = 0.68
    cosmo1.props.Tgamma0 = 2.72
    cosmo1.props.w0 = -0.672
    cosmo1.props.w1 = -1.2

    '''Modelo DEX'''

    cosmo2.props.H0 = 67.5
    cosmo2.props.Omegab = 0.05
    cosmo2.props.Omegac = 0.26
    cosmo2.props.Omegax = 0.68
    cosmo2.props.Tgamma0 = 2.72
    cosmo2.props.w = -1.2

    '''Modelo LCDM'''

    cosmo3.props.H0 = 70
    cosmo3.props.Omegab = 0.05
    cosmo3.props.Omegac = 0.26
    cosmo3.props.Omegax = 0.68
    cosmo3.props.Tgamma0 = 2.72

    #Definindo variáveis de nome curto:

    H0_1       = cosmo1.props.H0
    H0_2       = cosmo2.props.H0
    H0_3       = cosmo3.props.H0
    
    Omegab_3   = cosmo3.props.Omegab
    Omegac_3   = cosmo3.props.Omegac
    Omegax_3   = cosmo3.props.Omegax
    Omega_m_3 = Omegab_3 + Omegac_3

    N = 500

    z_list = []
    H1_list = []
    H2_list = []
    H3_list = []
    H3t_list = []
    dH1_list = []
    dH2_list = []
    dH3_list = []

    '''Loop de Redshift z'''

    for i in range(0,N):
        z = 6.0/(N - 1.0) * i
        H1 = cosmo1.H(z)
        H2 = cosmo2.H(z)
        H3t = cosmo3.props.H0 * np.sqrt(Omega_m_3*(1+z)**3 + Omegax_3)

        H3 = cosmo3.H(z)

        z_list.append(z)
        H1_list.append(H1)
        H2_list.append(H2)
        H3_list.append(H3)
        H3t_list.append(H3t)

        dH1 = cosmo1.dH_dz(z)
        dH2 = cosmo2.dH_dz(z)
        dH3 = cosmo3.dH_dz(z)

        dH1_list.append(dH1)
        dH2_list.append(dH2)
        dH3_list.append(dH3)
        

    return z_list, H1_list, H2_list, H3_list, H3t_list, H0_1, H0_2, H0_3, dH1_list, dH2_list, dH3_list

z_list, H1_list, H2_list, H3_list, H3t_list, H0_1, H0_2, H0_3, dH1_list, dH2_list, dH3_list = Hubble() 

'''Plotando'''

plt.rcdefaults()
plt.rcParams['figure.dpi'] = 150

plt.rcParams.update({
    'font.size': 16,
})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

colors = ["#960707", "#0450ac", "#15af31", "#000000"]

# ====== PLOT H(z) ======
ax1.plot(z_list, H1_list,
         label=rf'$H_1(z)$ ($H_0={H0_1:.1f}\,\mathrm{{km\,s^{{-1}}\,Mpc^{{-1}}}}$)',
         color=colors[0], linewidth=2.0)

ax1.plot(z_list, H2_list,
         label=rf'$H_2(z)$ ($H_0={H0_2:.1f}\,\mathrm{{km\,s^{{-1}}\,Mpc^{{-1}}}}$)',
         color=colors[1], linewidth=2.0)

ax1.plot(z_list, H3_list,
         label=rf'$H_3(z)$ ($H_0={H0_3:.1f}\,\mathrm{{km\,s^{{-1}}\,Mpc^{{-1}}}}$)',
         color=colors[2], linewidth=2.0)

ax1.plot(z_list, H3t_list,
         label=rf'$H_3(z)$ Teste Numérico ($H_0={H0_3:.0f}\,\mathrm{{km\,s^{{-1}}\,Mpc^{{-1}}}}$)',
         color=colors[3], linewidth=1.5, linestyle='--')

ax1.set_ylabel(r'$H(z)\,[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$')
ax1.set_title('Evolução de $H(z)$')

ax1.set_xlim(0, 1)
ax1.set_ylim(0.001, 150)

ax1.legend(
    frameon=True,
    facecolor='white',
    edgecolor='0.7',
    framealpha=1.0,
    loc='upper left',
    fontsize=12
)

ax1.grid(True, linestyle='--', alpha=0.7)

# ====== PLOT dH/dz ======
ax2.plot(z_list, dH1_list,
         label=r'$dH_1/dz$',
         color=colors[0], linewidth=2.0)

ax2.plot(z_list, dH2_list,
         label=r'$dH_2/dz$',
         color=colors[1], linewidth=2.0)

ax2.plot(z_list, dH3_list,
         label=r'$dH_3/dz$',
         color=colors[2], linewidth=2.0)

ax2.set_xlabel(r'Redshift $z$')
ax2.set_ylabel(r'$dH/dz\,[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$')

ax2.set_ylim(0, 150)
ax2.set_xlim(0, 1)

ax2.legend(
    frameon=True,
    facecolor='white',
    edgecolor='0.7',
    framealpha=1.0,
    loc='upper left',
    fontsize=14
)

ax2.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()