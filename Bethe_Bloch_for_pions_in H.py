import math
import matplotlib.pyplot as plt
import numpy as np

sqr = lambda x: x**2

""" Constants and formulas to be used for pion """
K = 0.307075 # constant K in MeV cm mol^-1
z = 1 # charge in e
Z = 1 # atomic number Z
A = 1 # atomic mass (g * mol-1)
Z_over_A = np.float(Z/A)
M = 139.6 # mass of pion (MeV)
m_e = 0.511 # mass of electron (MeV)
po = 8.99 * 10**(-5) # density of gas under normal conditions (g*cm^-3)

gamma = lambda p: math.sqrt(1 + sqr(p / M))
beta = lambda p: math.sqrt(1 - 1 / sqr(gamma(p)))
beta_gamma = lambda p: p / M
W_max = lambda p: 2 * m_e * sqr(beta_gamma(p)) / (1 + 2 * gamma(p) * m_e / M + sqr(m_e / M)) # maximum energy (MeV)

T_kin_pion = lambda p: math.sqrt( sqr(p) + sqr(M) ) - M

I = 0.01359844 # ionisation energy (MeV)

""" Thickness of absorber """
thickness_um = 5000
x = thickness_um/10000.0  # (cm)

epsilon = (K * po * Z * x) / (2 * A) # (MeV)
print('epsilon = ', epsilon)
E_mean = lambda p:  epsilon * (math.log(2 * m_e * W_max(p) * sqr(beta_gamma(p)) / sqr(I)) - 2 * sqr(beta(p)) ) / sqr(beta(p)) 

""" Calculation of all formulas in momentum range ps """

ps = np.logspace(1.22, 5, 1000) # (MeV)

dE_over_dX_pions = np.zeros(len(ps))

beta_gammas = np.zeros(len(ps))

betas = np.zeros(len(ps))
ps_GeV = np.zeros(len(ps))

E_means = np.zeros(len(ps))
T_kin = np.zeros(len(ps))

for pi in range(0, len(ps)):

	beta_gammas[pi] = beta_gamma( ps[pi] )
	betas[pi] = beta( ps[pi] )
	E_means[pi] = E_mean( ps[pi] ) * 1000 
	ps_GeV[pi] = ps[pi]/1000.0
	T_kin[pi] = T_kin_pion( ps[pi] )
	
print (' betas = ', betas)

""" Build Bethe Bloch plot versus momentum """
my_font = 14

fig, ax1 = plt.subplots()
fig.set_size_inches(8,6)

plt.tick_params(axis='both', which='major', labelsize=my_font)
ax1.plot(ps_GeV, E_means, label='Bethe Bloch', marker='', color='Red', linestyle='-', linewidth=2)

plt.grid(True, which="both", ls="-", color='0.45')

ax1.set_xscale('log')
plt.ylabel('Energy loss in '+str(thickness_um)+' um H in keV', fontsize = my_font)
plt.xlabel('Momentum p [GeV]', fontsize = my_font)
plt.legend(loc = 'upper center')

plt.tick_params(axis = 'both', which = 'major', labelsize = my_font)
plt.tight_layout()

plt.savefig('Bethe-for-pions_H.png')
plt.show()

plt.clf()
plt.close()