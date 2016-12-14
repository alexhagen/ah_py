import math
import sys
import numpy as np
sys.path.append("/Users/ahagen/code")
from pym import func as ahm
from pyg import twod as ahp
from ah_py.simulation import fluids as ahf

dfp = ahf.fluid('dfp')
T_b = dfp.T_b
P_b = dfp.P_b
dfp_pos_curve = ahm.curve(P_b, T_b - 273.15, name='dfp')

ace = ahf.fluid('ace')
T_b = ace.T_b
P_b = ace.P_b
ace_pos_curve = ahm.curve(P_b, T_b - 273.15, name='ace')

plot = dfp_pos_curve.plot(linestyle='-', linecolor='#E3AE24')

'''
A = 6.43876
B = 1242.510
C = 46.568
P = np.linspace(1., -700000., 500) / 1.0E3 # in kPa
T_b = - B / (np.log10(P) - A) + C
dfp_neg_curve = ahm.curve(P * 1.0E3, T_b, name='dfp-neg')
plot = dfp_neg_curve.plot(linestyle='--', linecolor='#E3AE24', addto=plot)
'''

dfp_exp_curve = ahm.curve([-3.2 * 1.0E5], [24.], name='dfp expt')
plot = dfp_exp_curve.plot(linestyle=None, linecolor='#E3AE24', addto=plot)

plot = ace_pos_curve.plot(linestyle='-', linecolor='#2EAFA4', addto=plot)
plot.add_data_pointer(101325., curve=dfp_pos_curve, string='atmospheric',
                      place=(-200000., 45.))

plot.legend(loc=2)
plot.xlabel(r"Pressure ($p$) [$Pa$]")
plot.ylabel(r"Temperature ($T$) [$^{o}C$]")
plot.markers_off()
plot.lines_on()
plot.lines["dfp expt0"].set_alpha(1.0)
plot.lines["dfp expt0"].set_markersize(6)
plot.ylim(15., 60.)
plot.xlim(-400000, 200000)
def barfunc(p):
    return p/1.0E5

plot.add_xx(barfunc)
plot.xlabel(r"Pressure ($p$) [$bar$]", axes=plot.ax2)



plot.export('boiling_curve', formats=['pdf'], sizes=['cs'], customsize=(6, 4))
plot.show()
