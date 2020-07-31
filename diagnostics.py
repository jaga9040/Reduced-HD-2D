import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf

mode_scan = ['1','2','3','3','4','5','6','7','8'] 

def linear_fit(x, a, b):
    return a * x + b
L = 2.0*np.pi
i = 0
sz=np.size(mode_scan)
kz = np.zeros(sz)
gama = np.zeros(sz)
fsz = 10

for vv in mode_scan:
	data = np.loadtxt('KHI_'+vv+'/fort.1')
	t = data[:,0]
	KEy  = data[:,2]
	log_KEy = np.log(KEy)	
	if int(vv)==1:
	  param, param1 = cf(linear_fit,t[1500:2000],log_KEy[1500:2000])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==2:
	  param, param1 = cf(linear_fit,t[1000:1800],log_KEy[1000:1800])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==3:
	  param, param1 = cf(linear_fit,t[810:1490],log_KEy[810:1490])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==4:
	  param, param1 = cf(linear_fit,t[680:1350],log_KEy[680:1350])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==5:
	  param, param1 = cf(linear_fit,t[570:1120],log_KEy[570:1120])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==6: 
	  param, param1 = cf(linear_fit,t[500:1100],log_KEy[500:1100])
	  ans = param[0] * np.array(t) + param[1]
	elif int(vv)==7: 
	  param, param1 = cf(linear_fit,t[480:1000],log_KEy[480:1000])
	  ans = param[0] * np.array(t) + param[1]
	else: 
	  param, param1 = cf(linear_fit,t[370:880],log_KEy[370:880])
	  ans = param[0] * np.array(t) + param[1]

	kz[i] = 2.0*np.pi*float(vv)/L
	gama[i] = param[0]/2
	i += 1
	
	plt.figure(1)
	plt.plot(t,log_KEy,label='MODE ='+vv+'; growth_rate = '+str(param[0]), linewidth=1.5)
	plt.plot(t,ans, linestyle='dashed',color='gray',linewidth=1.0)
	plt.xticks(fontsize=fsz)
	plt.yticks(fontsize=fsz)
	plt.ylabel(r'$\log(\langle KE_x(t) \rangle ) $', fontsize=fsz)
	plt.title('growth rate of KHI calculation for Linearised NS-equation', fontsize=fsz)
	plt.legend(fontsize=fsz)
	plt.xlabel(r'Time', fontsize=fsz)


plt.figure()
U0 = 1.840
nu = 1.5e-3
d = 0.1 
Re = U0 * d/nu
sz = int(7.1/0.01)
m = np.zeros(sz)
gamma = np.zeros(sz)
k_z = np.zeros(sz)
L = 1.0*np.pi
for i in range(0, sz, 1):
    m[i] = 1.0 + i*0.01
    k_z[i] = 2.0 * np.pi * m[i]/L
    ans = param[0] * np.array(t) + param[1] 
    gamma[i] = (k_z[i]*U0/3.0)*(np.sqrt(3.0) - 2.0*k_z[i]/Re - 2.0*((k_z[i]/Re)**2+2.0*np.sqrt(3.0)*k_z[i]/Re)**(0.5))
#np.savetxt('k_v_gamma.txt',(kz,gama))
print(m)
plt.plot(m, gamma, label='Analytical growth rate for d=0.03',linewidth=2)
plt.plot(kz, gama,label='numerical growth rate for d=0.015',linewidth=2, marker = 'd')
plt.title('Growth rate comparison', fontsize=fsz)
plt.ylabel(r'$Growth~ rate~(\gamma)$', fontsize=fsz)
plt.xlabel(r'$Mode ~number~(K_z)$', fontsize=fsz)
plt.legend(fontsize=fsz)
plt.xticks(fontsize=fsz)
plt.yticks(fontsize=fsz)
plt.grid()
#plt.savefig('gamma_kz.png', dpi=300)
plt.show()
