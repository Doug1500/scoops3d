from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

# Read the file. 
# f1 = open('lattice.asc', 'r')
f2 = open('cone1000DEM.asc', 'r')
# f3 = open('emb20layer_2plot.asc', 'r')
# f4 = open('emb20layer_1plot.asc', 'r')
# read the whole file into a single variable, which is a list of every row of the file.

lines = f2.readlines()
f2.close()

# lines1 = f3.readlines()
# f3.close()

# lines2 = f4.readlines()
# f4.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


xF = []
yF = []
zF = []
# FF = []


# for i in range(len(lines0)):
# 	F = lines0[i].split()
# 	xF.append(F[0])
# 	yF.append(F[1])
# 	zF.append(F[2])
# 	if float(F[3]) > 10.0:
# 		FF.append(0.0)
# 	if float(F[3]) < 10.0:
# 		FF.append(F[3])

axF = np.array(xF, dtype=float)
ayF = np.array(yF, dtype=float)
azF = np.array(zF, dtype=float)
# aFF = np.array(FF, dtype=float)

cmap = plt.matplotlib.cm.jet

ax.scatter(axF, ayF, azF, cmap=cmap)
ax.set_title('Embankment Problem', fontsize='20')

plt.xlabel('Wdith (m)')
plt.ylabel('Distance (m)')
ax.set_zlabel('Height (m)')

colour = 'green'
# colour1= 'red'
# colour2= 'grey'

# grid size
gdsx = 97
gdsy = 97

#cell size
csix = 43.3014
csiy = 43.3014 

# initialize some variable to be lists:
x = np.linspace(0, gdsx*csix, gdsx)
y = np.linspace(0, gdsy*csiy, gdsy)

xm = []
ym = []
pm = []
qm = []
rm = []

for i in range(len(x)):
	# counter checking
	print i
	for j in range(len(y)):
		p = lines[i].split()
		# q = lines1[i].split()
		# r = lines2[i].split()
		xm.append(x[i])
		ym.append(y[j])
		pm.append(p[j])
		print x[i], y[j], p[j]
		# qm.append(q[j])
		# if float(r[j]) > 0.0:
		# 	rm.append(r[j])
		# if float(r[j]) < 0.0:
		# 	rm.append(p[j])# 	rm.append(0.0)

axm=np.array(xm, dtype=float)
aym=np.array(ym, dtype=float)
apm=np.array(pm, dtype=float)
# aqm=np.array(qm, dtype=float)
# arm=np.array(rm, dtype=float)

# ax.scatter(axm, aym, aqm, color= colour1)
# ax.scatter(axm, aym, arm, color= colour2)
ax.scatter(axm, aym, apm, color= colour)
# ax.text(50.0, 30.0, 40.0, "HAY", color="red")

plt.show()

# for i in range(len(x)):
# 	# counter checking
# 	print i
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		q = lines1[i].split()
# 		ax.plot([float(x[i])],[float(y[j])],[float(p[j])], markerfacecolor= colour, marker='o')
# 		ax.plot([float(x[i])],[float(y[j])],[float(q[j])], markerfacecolor= colour1, marker='o')