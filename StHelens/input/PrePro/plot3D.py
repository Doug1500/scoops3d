from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Read the file. 
f2 = open('sthel_res100mDEM.asc', 'r')

# read the whole file into a single variable, which is a list of every row of the file.
lines = f2.readlines()
f2.close()

xF = []
yF = []
zF = []

colour= 'grey'

# grid size
gdsx = 91
gdsy = 96

#cell size
csix = 100
csiy = 100

# initialize some variable to be lists:
x = np.linspace(557974.631687, 557974.631687+gdsx*csix, gdsx)
y = np.linspace(5111446.1291019, 5111446.1291019+gdsy*csiy, gdsy)
# x = np.linspace(0.0, gdsx*csix, gdsx)
# y = np.linspace(0.0, gdsy*csiy, gdsy)

xm = []
ym = []
pm = []


for i in range(len(x)):
	# counter checking
	# print i
	for j in range(len(y)):
		p = lines[i].split()
		# if float(p[j]) == -9999.0:
			# print x[i],",", y[j],",", 900.0
		# else: 
			# print x[i],",", y[j],",", p[j]
		xm.append(x[i])
		ym.append(y[j])
		if float(p[j]) == -9999.0:
			pm.append(1000)
		if float(p[j]) > 0.0:
			pm.append(p[j])

axm=np.array(xm, dtype=float)
aym=np.array(ym, dtype=float)
apm=np.array(pm, dtype=float)

ax.set_title('Mount St.Helens', fontsize='20')
plt.xlabel('Wdith (m)')
plt.ylabel('Distance (m)')
ax.set_zlabel('Height (m)')

ax.scatter(axm, aym, apm, c='brown', s=100.0)
# ax.plot([564126.0, 564126.0],[5.11825e+06, 5.11825e+06],zs=[3837.82, 1866.98])

# plt.savefig('/home/doug/shared/'+NAME+'.png',dpi=150)
for ii in xrange(0, 360, 10):
	ax.view_init(elev=10.0, azim=ii)
	fig.set_size_inches(20.0, 10.0)
	# plt.savefig("/home/doug/mygo/src/SCOOP/movies/movie"+str(ii)+".png", dpi=150)

ax.axis('equal')

plt.show()

# for i in range(len(x)):
# 	# counter checking
# 	print i
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		q = lines1[i].split()
# 		ax.plot([float(x[i])],[float(y[j])],[float(p[j])], markerfacecolor= colour, marker='o')
# 		ax.plot([float(x[i])],[float(y[j])],[float(q[j])], markerfacecolor= colour1, marker='o')


# X, Y = np.meshgrid(x, y)
# appm = X + Y
# print apm
# k = 0
# for i in range(len(X)):
# 	for j in range(len(X[0])):
# 		print appm[i][j]
# 		appm[i][j] = apm[k]
# 		print k
# 		k =  k + 1


# print X
# print Y
# print appm
	
# ax.plot_surface(X, Y, appm, rstride = 3, cstride = 3, alpha = 0.3, cmap = cm.BuPu)

# plt.show()