import numpy as np
from nozzle import nozzle_height
import useful
from matplotlib import pyplot as plt

#données de base

h_throat = 0.003 #3mm
h_duct = 0.0075 #7.5 mm
length_nozzle = 0.12
length_duct = 1.5#58.29722
d = 2*h_duct
b = 0.05 #50mm
Dh_duct = 4*h_duct*b/(b+2*h_duct) #4*A/perimeter
print('Dh_duct',Dh_duct)

p_a = 1.0 #atm
gamma = 1.4
R=287.1 #J/kg/K
T0 = 300 #K
"""
#a) find po in order to have sonic throat and subsonic after and before

"""

#a.1) Determine the  sonic conditions A_t = A*
T_star = T0/(gamma+1)*2
c_star = np.sqrt(gamma*R*T_star)
#a.2) Determine the conditions at the end of the nozzle
air_ratio = h_throat/h_duct#ratio between throat and end_nozzle with iteration

Mi = useful.inv_air_ratio(air_ratio,np.array([0.5])) # calculates Mi using the air ratio with iteration
print(Mi)
Ti = T0/(1+0.2*Mi**2)
ui = Mi * np.sqrt(gamma*R*Ti)

#a.3) we need to find pi but we have no information about the fanno flow after so lets find the conditions at the exit of the duct and then return to the nozzle
pe = 1.01325 #atm safe assumption pa=pe

#a.3.1)  trouver M_e avec critical length et les perte de charge lambda* (x2-x1)/D = f(M1)-f(M2), x1=0 , x2 = 1.5 M1=Mi et donc f(M1)-lamda*1.5/D = f(M2)

# a3.2) pour trouver les conditions à l outlet on va iterer en prenant comme variable le flux massique qui doit etre conservée
# massflow_i =1000*ui #first estimation rho_i = 1000 kg/m**3 air = eau
# iter2=1
# error2=1
# while error2>0.0001 and iter2 <100 :
#
#     Re = useful.Re(massflow_i/ui,ui,d,Ti)
#     print(Re)
#     lambda_duct = useful.perte_lambda(Re,0.05)
#     print(lambda_duct)
#     a = useful.critical_lenghth(Mi)-0.5*lambda_duct*length_duct/Dh_duct
#     print(a)
#     Me = useful.inv_critical_length(a,np.array([0.5]))
#     Te = T0/(1+0.2*Me**2) #throughout d-the duct T0 is conserved as well as po*
#     rho_e = pe*10**(5)/(287.1*Te)
#     ue = Me*np.sqrt(gamma*R*Te)
#     massflow_new = rho_e*ue
#     error2 = abs(massflow_i-massflow_new)
#     iter2 +=1
#     massflow_i= massflow_new
Me = 0.7 #first estimation rho_i = 1000 kg/m**3 air = eau
iter2=1
error2=1
while error2>0.0001 and iter2 <50 :
    Te = T0/(1+0.2*Me**2) #throughout d-the duct T0 is conserved as well as po*
    rho_e = pe*10**(5)/(287.1*Te)
    ue = Me*np.sqrt(gamma*R*Te)
    #Re = useful.Re(massflow_i/ui,ui,d,Ti)
    Re = useful.Re(rho_e,ue,d,Te)
    lambda_duct = useful.perte_lambda(Re,0.05)
    print('lambda_e',lambda_duct)
    print('Re',Re)
    a = useful.critical_lenghth(Mi)-0.5*lambda_duct*length_duct/Dh_duct
    M_new = useful.inv_critical_length(a,np.array([0.5]))

    error2 = abs(M_new-Me)
    iter2 +=1
    Me = M_new


#a.4) resulats a la sortie
print('Me',Me)
# print('massflow chequ',ue*rho_e,massflow_i)
poe = pe*(1+0.2*Me**2)**(1.4/0.4)
poe_star = poe/(1/Me*((1+0.2*Me**2)/1.2)**(3))
p_star = pe*Me/(1.2/(1+0.2*Me**2))**(1.2/0.4)
print('conditions outlet duct: ',Me,ue,Te,pe,poe)

#a.5) resultats au inlet duct
#rho_i = massflow_i/ui

#pi = rho_i*R*Ti*10**(-5)
poi = poe_star*(1/Mi*((1+0.2*Mi**2)/1.2)**(3))
pi=poi*(Ti/T0)**(1.4/0.4)
p_stari =  pi*Mi/(1.2/(1+0.2*Mi**2))**(1.2/0.4)
#poi = pi*(T0/Ti)**(1.4/0.4)
poi_star = poi/(1/Mi*((1+0.2*Mi**2)/1.2)**(3))
#print('chequ',poi,poi2)
print('p_star chequ',p_star,p_stari)
print('po_star chequ', poe_star,poi_star)
print('conditions inlet duct', Mi,ui,Ti,pi,poi)
print('reponse exercice a', poi)

#a.6) conditions sonic au throat

p_star = poi*(T_star/T0)**(1.4/0.4)
rho_star = p_star*10**(5)/(R*T_star)
print('conditions at throat', 1, c_star, T_star, p_star, poi)

#a.7) massflow chequ
Qm = ue*rho_e*b*2*h_duct#kg/s
print('massflow', Qm)

"""
3.a) graphe  the behaviour of p(x),po(x),M(x) in the nozzle and duct

"""
# askip po* est conservé dans un fanno flow et il peut etre calculé avec le inlet et le outlet ==> verifier qu a on a la meme chose
po_star = poe / (1/Me*((1+0.2*Me**2)/1.2)**(1.2/0.4)) #poi / (1/Mi*((1+0.2*Mi**2)/1.2)**(1.2/0.4))
#3.a_nozzle
x=np.linspace(0,120,1000)
ones=np.ones(len(x))
h= np.zeros(len(x))
for i in range(len(h)):
     h[i]= nozzle_height(x[i])/1000
#h[:]= nozzle_height(x[:])/1000
air_ratios = h_throat/h
Mx = useful.inv_air_ratio(air_ratios,0.5*np.ones(len(x)))
px = poi/((ones+0.2*Mx**2)**(1.4/0.4))

fig, ax = plt.subplots()
ax.set_title('a nozzle')
ax.plot(x,Mx,'-b',label='M')
ax.plot(x,px,'-r',label='p[bar]')
ax.plot(x,poi*ones,'-c',label='p0[bar]')
ax.legend()
ax.grid()


#3.a_duct
x2=np.linspace(0,length_duct,1000)
ones=np.ones(len(x2))
a = -0.5*lambda_duct*x2/Dh_duct+ones*useful.critical_lenghth(Mi)
Mx_d = useful.inv_critical_length(a,0.2*ones)
pox = po_star*(1/Mx_d*((ones+0.2*Mx_d**2)/1.2)**(1.2/0.4))
px = pox/(ones+0.2*Mx_d**2)**(1.4/0.4)

fig2=plt.figure()
ax2 =plt.subplot(311)
ax2.set_title('a duct')
ax2.plot(x2,Mx_d,'-b',label='M')
ax22=plt.subplot(312)
ax22.plot(x2,px,'-r',label='p[bar]')
ax23=plt.subplot(313)
ax23.plot(x2,pox,'-c',label='p0[bar]')
ax2.legend()
ax2.grid()
ax22.legend()
ax22.grid()
ax23.legend()
ax23.grid()

plt.show()
