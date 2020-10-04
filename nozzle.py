import numpy as np
from matplotlib import pyplot as plt
pi= np.pi


yadd = 3 #mm
r1 = 254.3 #mm
r3 = 153.7 #mm
alpha = 3.29 # deg
x_throat = 30 # mm

x1_fin = x_throat + np.sqrt(np.tan(alpha*pi/180)**2*r1**2/(1+np.tan(alpha*pi/180)**2))

x1= np.linspace(0,x1_fin,1000)
y1 = yadd + r1 - np.sqrt(r1**2-(x1-x_throat)**2)

Oy3 = 4.5 - r3
Ox3 = 120
x3_deb = Ox3 - np.sqrt(np.tan(alpha*pi/180)**2*r3**2/(1+np.tan(alpha*pi/180)**2))

x2 = np.linspace(x1_fin,x3_deb,1000)
y2 = y1[-1]+np.tan(alpha*pi/180)*(x2-x1_fin)

x3= np.linspace(x3_deb,Ox3,1000)
y3 = yadd + Oy3 + np.sqrt(r3**2-(x3-Ox3)**2)

def nozzle_height(x):
    if x<=x1_fin:
        y= yadd + r1 - np.sqrt(r1**2-(x-x_throat)**2)
    elif(x>=x3_deb):
        y=yadd + Oy3 + np.sqrt(r3**2-(x-Ox3)**2)
    else:
        y=y1[-1]+np.tan(alpha*pi/180)*(x-x1_fin)
    return y
x=np.linspace(0,120,1000)
y= np.zeros(len(x))
for i in range(len(y)):
    y[i]= nozzle_height(x[i])

# plt.plot(x1,y1,'-b')
# plt.plot(x1,-y1,'-b')
# plt.plot(x2,y2,'-r')
# plt.plot(x2,-y2,'-r')
# plt.plot(x3,y3,'-g')
# plt.plot(x3,-y3,'-g')
plt.plot(x,y,'-c')
plt.plot(x,-y,'-c')
plt.grid()
plt.xlim((0,120))
plt.ylim((-8,8))
plt.show()
