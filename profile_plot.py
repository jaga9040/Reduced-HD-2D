import matplotlib.pyplot as plt
import numpy as np
import f90nml

# reading data from fortran namelist 'input.inp'
nml = f90nml.read('input.inp')
t_param = nml['time_param']
mesh = nml['mesh']
#phy = nml['phy']

# reading grid size from namelist 
nx=mesh['Nx']
ny=mesh['Ny']

# calculating time_unit (time difference between each snap file) and 
# file_max (max number of file present + 1 + 10, because first file is fort.10)
time_unit = (t_param['dt']*t_param['ifile_write'])
file_max = 11 + int(t_param['time_max']/time_unit)

for i in range(10, file_max, 1):
 f=np.loadtxt('KHI_3/fort.'+ str(i))#, unpack = True)
 xlist=f[:,0]
 ylist=f[:,1]
 xi=np.linspace(np.min(xlist), np.max(xlist), nx)
 yi=np.linspace(np.min(ylist), np.max(ylist), ny)
 Y, X = np.meshgrid(yi, xi)
 plt.figure()
# plt.xlabel('x-axis')
# plt.ylabel('y-axis')
# z=f[:,6]
# zi=np.reshape(z,(Nx,Ny))
# plt.ylim([-5,5])
 #plt.xlim([5.28,7.28])
# levels=200
# cp =plt.contour(X,Y,zi,levels,linewidths=0.7)
# plt.colorbar(cp)
# plt.title('profile of Az at t='+str((i-10)*2))
# plt.savefig('Az_' + str((i-9)) + '.png', format='png', dpi=300)
# plt.clf()
# plt.figure()
 plt.xlabel('x-axis')
 plt.ylabel('y-axis')
 z=f[:,5]
 zi=np.reshape(z,(nx,ny))
# plt.ylim([-5,5])
 #plt.xlim([5.28,7.28])
 levels=200
 cmap1=plt.cm.get_cmap("jet")
 cp =plt.contourf(X,Y,zi,levels,cmap=cmap1)
 plt.colorbar(cp)
 plt.title('color plot of Wz at t='+str((i-10)*time_unit))
 plt.savefig('wz_' + str((i-9)) + '.png', format='png', dpi=300)
 plt.clf()
 plt.close('all')
