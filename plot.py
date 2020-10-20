import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 20})

p_a = 1.01325 #[bar] = 1[atm]
l_nozzle = 120 #[mm]
l_duct = 1.500 #m
l_throat = 30 #[mm]
#a) homework a
x_an = np.loadtxt("results\homework_a_nozzle.txt")[0]
Mx_an = np.loadtxt("results\homework_a_nozzle.txt")[1]
px_an = np.loadtxt("results\homework_a_nozzle.txt")[2]

x_ad = np.loadtxt("results\homework_a_duct.txt")[0]
Mx_ad = np.loadtxt("results\homework_a_duct.txt")[1]
px_ad = np.loadtxt("results\homework_a_duct.txt")[2]
pox_ad = np.loadtxt("results\homework_a_duct.txt")[3]

#b) homework b

x_bn = np.loadtxt("results\homework_b_nozzle.txt")[0]
Mx_bn = np.loadtxt("results\homework_b_nozzle.txt")[1]
px_bn = np.loadtxt("results\homework_b_nozzle.txt")[2]

x_bd = np.loadtxt("results\homework_b_duct.txt")[0]
Mx_bd = np.loadtxt("results\homework_b_duct.txt")[1]
px_bd = np.loadtxt("results\homework_b_duct.txt")[2]
pox_bd = np.loadtxt("results\homework_b_duct.txt")[3]

#c) homework c

x_cn = np.loadtxt("results\homework_c_nozzle.txt")[0]
Mx_cn = np.loadtxt("results\homework_c_nozzle.txt")[1]
px_cn = np.loadtxt("results\homework_c_nozzle.txt")[2]

x_cd = np.loadtxt("results\homework_c_duct.txt")[0]
Mx_cd = np.loadtxt("results\homework_c_duct.txt")[1]
px_cd = np.loadtxt("results\homework_c_duct.txt")[2]
pox_cd = np.loadtxt("results\homework_c_duct.txt")[3]

"""
Do nice plots
"""
#a) graphe du nozzle comparer les Me
fig, ax = plt.subplots()
#ax.set_title('Mach number in the nozzle')
ax.plot(x_an/l_nozzle,Mx_an,'-b',label='a')
ax.plot(x_bn/l_nozzle,Mx_bn,'-r',label='b')
ax.plot(x_cn/l_nozzle,Mx_cn,'-g',label='c')
ax.set_xlabel("x/l [-]")
ax.set_ylabel("M[-]")
ax.legend()
ax.grid()


#b) graphe du nozzle comparer les p
fig2, ax2 = plt.subplots()
#ax2.set_title('Pressure in the nozzle')
ax2.plot(x_an/l_nozzle,px_an/p_a,'-b',label='a')
ax2.plot(x_bn/l_nozzle,px_bn/p_a,'-r',label='b')
ax2.plot(x_cn/l_nozzle,px_cn/p_a,'-g',label='c')
ax2.set_xlabel("x/l [-]")
ax2.set_ylabel("$p/p_{a}$ [-]")
ax2.legend()
ax2.grid()

#c) graphe du duct comparer les Me
fig, ax = plt.subplots()
#ax.set_title('Mach number in the duct')
ax.plot(x_ad/l_duct,Mx_ad,'-b',label='a')
ax.plot(x_bd/l_duct,Mx_bd,'-r',label='b')
ax.plot(x_cd/l_duct,Mx_cd,'-g',label='c')
ax.set_xlabel("x/l [-]")
ax.set_ylabel("M[-]")
ax.legend()
ax.grid()


#d) graphe du duct comparer les p
fig, ax = plt.subplots()
#ax.set_title('Pressure in the duct')
ax.plot(x_ad/l_duct,px_ad,'-b',label='a')
ax.plot(x_bd/l_duct,px_bd,'-r',label='b')
ax.plot(x_cd/l_duct,px_cd,'-g',label='c')
ax.set_xlabel("x/l [-]")
ax.set_ylabel("$p/p_{a}$ [-]")
ax.legend()
ax.grid()

#e) graphe du duct comparer les p0
fig, ax = plt.subplots()
#ax.set_title('Total pressure in the duct')
ax.plot(x_ad/l_duct,pox_ad,'-b',label='a')
ax.plot(x_bd/l_duct,pox_bd,'-r',label='b')
ax.plot(x_cd/l_duct,pox_cd,'-g',label='c')
ax.set_xlabel("x/l [-]")
ax.set_ylabel("$p_{0}/p_{a}$ [-]")
ax.legend()
ax.grid()
plt.show()
