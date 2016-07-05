# ############################################################################################3
# ############################################################################################3
# ############################################################################################3
# ############################################################################################3

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')


# f2 = open('emb20DEM.asc', 'r')
# f3 = open('emb20DEM.asc', 'r')
# # f2 = open('A_emb20 - Copy_newDEM_out.asc', 'r')
# # f3 = open('A_emb20 - Copy_newDEM_out.asc', 'r')

# lines = f2.readlines()
# f2.close()

# liness = f3.readlines()
# f3.close()

# xF = []
# yF = []
# zF = []

# # grid size
# gdsx = 86
# gdsy = 152

# #cell size
# csix = 20.0
# csiy = 20.0

# # initialize some variable to be lists:
# x = np.linspace(0, gdsx*csix, gdsx)
# y = np.linspace(0, gdsy*csiy, gdsy)

# xm = []
# ym = []
# qm = []
# qqm = []

# for i in range(len(x)):
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		q = liness[i].split()
# 		if abs(float(p[j]) - float(q[j])) > 0.0:
# 			xm.append(x[i])
# 			ym.append(y[j])
# 			qm.append(q[j])
# 			qqm.append(p[j])
# 			# print x[i],",", y[j],",", q[j]
# 			# print x[i],",", y[j],",", p[j]
# 		# elif float(p[j]) - float(q[j]) < 0.0:
# 		# 	pm.append(0)
# 		# else:
# 		# 	pm.append(0)


# axm=np.array(xm, dtype=float)
# aym=np.array(ym, dtype=float)
# aqm=np.array(qm, dtype=float)
# aqqm=np.array(qqm, dtype=float)

# ax.scatter(axm, aym, aqm, c='red', s=500.0)
# ax.scatter(axm, aym, aqqm, c='yellow', s=500.0)


# bbxm = []
# bbym = []
# bbpm = []


# for i in range(len(x)):
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		bbxm.append(x[i])
# 		bbym.append(y[j])
# 		print y[j],",", x[i],",", p[j]
# 		if float(p[j]) <= -9999.0:
# 			bbpm.append(2000)
# 		else:
# 			bbpm.append(p[j])

# bxm=np.array(bbxm, dtype=float)
# bym=np.array(bbym, dtype=float)
# bpm=np.array(bbpm, dtype=float)

# ax.scatter(bxm, bym, bpm, c='yellow', s=2)

# ccxm = []
# ccym = []
# ccpm = []


# for i in range(len(x)):
# 	for j in range(len(y)):
# 		p = liness[i].split()
# 		ccxm.append(x[i])
# 		ccym.append(y[j])
# 		if float(p[j]) == -9999.0:
# 			ccpm.append(1000)
# 		if float(p[j]) > 0.0:
# 			ccpm.append(p[j])
			
			

# cxm=np.array(ccxm, dtype=float)
# cym=np.array(ccym, dtype=float)
# cpm=np.array(ccpm, dtype=float)

# # ax.scatter(cxm, cym, cpm, c='brown', s=8)

# # ax.axis('equal')
# # ax.auto_scale_xyz([min(bxm), max(bxm)], [min(bym), max(bym)], [min(bpm), max(bpm)])

# plt.show()

############################################################################################3
############################################################################################3
############################################################################################3
############################################################################################3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


f2 = open('emb20DEM.asc', 'r')
f3 = open('emb20DEM.asc', 'r')
# f2 = open('A_emb20 - Copy_newDEM_out.asc', 'r')
# f3 = open('A_emb20 - Copy_newDEM_out.asc', 'r')

lines = f2.readlines()
f2.close()

liness = f3.readlines()
f3.close()

xF = []
yF = []
zF = []

# grid size
gdsx = 162
gdsy = 142

#cell size
csix = 20.0
csiy = 20.0

# initialize some variable to be lists:
x = np.linspace(0, gdsx*csix, gdsx)
y = np.linspace(0, gdsy*csiy, gdsy)

xm = []
ym = []
qm = []
qqm = []

for i in range(len(x)):
	for j in range(len(y)):
		p = lines[i].split()
		q = liness[i].split()
		if abs(float(p[j]) - float(q[j])) > 0.0:
			xm.append(x[i])
			ym.append(y[j])
			qm.append(q[j])
			qqm.append(p[j])
			# print x[i],",", y[j],",", q[j]
			# print x[i],",", y[j],",", p[j]
		# elif float(p[j]) - float(q[j]) < 0.0:
		# 	pm.append(0)
		# else:
		# 	pm.append(0)


axm=np.array(xm, dtype=float)
aym=np.array(ym, dtype=float)
aqm=np.array(qm, dtype=float)
aqqm=np.array(qqm, dtype=float)

ax.scatter(axm, aym, aqm, c='red', s=500.0)
ax.scatter(axm, aym, aqqm, c='yellow', s=500.0)


bbxm = []
bbym = []
bbpm = []


for i in range(len(x)):
	for j in range(len(y)):
		p = lines[i].split()
		bbxm.append(x[i])
		bbym.append(y[j])
		print y[j],",", x[i],",", p[j]
		if float(p[j]) <= -9999.0:
			bbpm.append(2000)
		else:
			bbpm.append(p[j])

bxm=np.array(bbxm, dtype=float)
bym=np.array(bbym, dtype=float)
bpm=np.array(bbpm, dtype=float)

ax.scatter(bxm, bym, bpm, c='yellow', s=2)

ccxm = []
ccym = []
ccpm = []


for i in range(len(x)):
	for j in range(len(y)):
		p = liness[i].split()
		ccxm.append(x[i])
		ccym.append(y[j])
		if float(p[j]) == -9999.0:
			ccpm.append(1000)
		if float(p[j]) > 0.0:
			ccpm.append(p[j])
			
			

cxm=np.array(ccxm, dtype=float)
cym=np.array(ccym, dtype=float)
cpm=np.array(ccpm, dtype=float)

# ax.scatter(cxm, cym, cpm, c='brown', s=8)

# ax.axis('equal')
# ax.auto_scale_xyz([min(bxm), max(bxm)], [min(bym), max(bym)], [min(bpm), max(bpm)])

plt.show()

