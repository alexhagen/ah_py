import math
import sys
import numpy as np
sys.path.append("/Users/ahagen/code")
from pym import func as ahm
from pyg import twod as ahp
from ah_py.simulation import fluids as ahf
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy import interpolate

dfp = ahf.fluid('dfp');
T = np.linspace(275., 325., 90)
P = np.linspace(50000., 150000., 100)

T_mid = T[:-1] + (T[1:] - T[:-1]) / 2.
P_mid = P[:-1] + (P[1:] - P[:-1]) / 2.

#Ta,Pa = np.meshgrid(T, P)
#T_rho_d = (Ta[:-1, 1:] - Ta[:-1, :-1])/2.0 + Ta[1:, :-1]
#P_rho_d = (Pa[1:, :-1] - Pa[:-1, :-1])/2.0 + Pa[:-1, 1:]
#rho_d = np.zeros_like(T_rho_d)
#for i in range(len(P) - 1):
#    for j in range(len(T) - 1):
#        rho_d[i, j] = dfp.rho(T_rho_d[i, j], P_rho_d[i, j])
#
#mask = dfp.T_b_curve.at(P_rho_d) <= T_rho_d
#rho_d = np.ma.masked_where(mask, rho_d)

#boiling = dfp.T_b_curve
#boiling_plot = boiling.plot()
#boiling_plot.lines_on()
#boiling_plot.markers_off()
#residual = rho_d
#r = interpolate.interp2d(T_mid, P_mid, residual)
#plt.pcolormesh(P_rho_d, T_rho_d, residual,
#               cmap=\
#                    mpl.colors.LinearSegmentedColormap.from_list('PU',["#ffffff","#E3AE24"],\
#    N=1024),edgecolors='face'
#    )
#cbar = plt.colorbar()
#cbar.set_label(r'Density ($\rho_{tait}$) [$\frac{kg}{m^2}$]')
#boiling_plot.ylim(275., 325.)
#boiling_plot.xlim(50000., 150000.)
#boiling_plot.xlabel('Pressure ($p$) [$Pa$]')
#boiling_plot.ylabel('Temperature ($T$) [$^{o} C$]')
#rho_stp = r(298.73, 101325.)
#stp_str = '$\\frac{|\\rho_{tait}-\\rho_{nist}|}{\\rho_{tait}} \left( STP \\right) = %5.3f \%% $' % (rho_stp);
#boiling_plot.export('rho_dfp',
#    formats=['pdf', 'pgf'],sizes=['2']);
#boiling_plot.show()

P_atm = 101.325E3
T_atm = np.linspace(-20, 55., 1000) + 273.15
rho_atm = []
for T in T_atm:
    rho_atm = np.append(rho_atm, dfp.rho(T, P_atm))

atm_curve = ahm.curve(T_atm - 273.15, rho_atm, name=r'$\rho_{dfp}$')
plot = atm_curve.plot(linestyle='-', linecolor='#7299C6')

T = np.arange(-20., 130., 10.)
rho = np.array([1.70, 1.680, 1.660, 1.620, 1.60, 1.57, 1.550, 1.51, 1.49,
                1.460, 1.430, 1.380, 1.34, 1.32, 1.30]) * 1.0E3

b_curve = ahm.curve(T, rho, name=r'Dupont data')
plot = b_curve.plot(linestyle='-', linecolor='#B63F97', addto=plot)

plot.markers_off()

plot.lines_on()
plot.lines['Dupont data0'].set_markersize(6.0)
plot.lines['Dupont data0'].set_linewidth(0.0)

np.savetxt('/Users/ahagen/Desktop/dfp_density.csv', [T_atm.T, rho_atm.T], delimiter=',')

plot.legend()
plot.xlim(-20., 55.)
plot.xlabel(r'Temperature ($T$) [$^{o} C$]')
plot.ylabel(r'Densisty ($\rho$) [$\frac{kg}{m^2}$]')
plot.export('rho_atm_dfp', sizes=['cs'], formats=['pdf', 'pgf'],
            customsize=(4.375, 2.25))
plot.show()
