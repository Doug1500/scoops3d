# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Read the file. 
# f2 = open('sthel_res100mDEM.asc', 'r')
# # f2 = open('test0_newDEM_out.asc', 'r')

# # read the whole file into a single variable, which is a list of every row of the file.
# lines = f2.readlines()
# f2.close()

# xF = []
# yF = []
# zF = []

# colour= 'grey'

# # grid size
# gdsx = 91
# gdsy = 96

# #cell size
# csix = 100
# csiy = 100

# # initialize some variable to be lists:
# x = np.linspace(557974.631687, 557974.631687+gdsx*csix, gdsx)
# y = np.linspace(5111446.1291019, 5111446.1291019+gdsy*csiy, gdsy)
# # x = np.linspace(0.0, gdsx*csix, gdsx)
# # y = np.linspace(0.0, gdsy*csiy, gdsy)

# xm = []
# ym = []
# pm = []


# for i in range(len(x)):
# 	# counter checking
# 	# print i
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		# if float(p[j]) == -9999.0:
# 			# print x[i],",", y[j],",", 900.0
# 		# else: 
# 			# print x[i],",", y[j],",", p[j]
# 		xm.append(x[i])
# 		ym.append(y[j])
# 		if float(p[j]) == -9999.0:
# 			pm.append(0)
# 		if float(p[j]) > 0.0:
# 			pm.append(p[j])

# axm=np.array(xm, dtype=float)
# aym=np.array(ym, dtype=float)
# apm=np.array(pm, dtype=float)

# ax.scatter(axm, aym, apm, c='green', s=500.0)

# # f3 = open('sthel_res100mDEM.asc', 'r')
# f3 = open('test0_newDEM_out.asc', 'r')

# # read the whole file into a single variable, which is a list of every row of the file.
# liness = f3.readlines()
# f3.close()

# xF = []
# yF = []
# zF = []

# colour= 'grey'

# # grid size
# gdsx = 91
# gdsy = 96

# #cell size
# csix = 100
# csiy = 100

# # initialize some variable to be lists:
# x = np.linspace(557974.631687, 557974.631687+gdsx*csix, gdsx)
# y = np.linspace(5111446.1291019, 5111446.1291019+gdsy*csiy, gdsy)
# # x = np.linspace(0.0, gdsx*csix, gdsx)
# # y = np.linspace(0.0, gdsy*csiy, gdsy)

# xm = []
# ym = []
# pm = []
# qm = []


# for i in range(len(x)):
# 	# counter checking
# 	# print i
# 	for j in range(len(y)):
# 		p = liness[i].split()
# 		# if float(p[j]) == -9999.0:
# 			# print x[i],",", y[j],",", 900.0
# 		# else: 
# 			# print x[i],",", y[j],",", p[j]
# 		xm.append(x[i])
# 		ym.append(y[j])
# 		if float(p[j]) == -9999.0:
# 			pm.append(0)
# 		if float(p[j]) > 0.0:
# 			pm.append(p[j])

# axm=np.array(xm, dtype=float)
# aym=np.array(ym, dtype=float)
# apm=np.array(pm, dtype=float)

# ax.set_title('Mount St.Helens', fontsize='20')
# plt.xlabel('Wdith (m)')
# plt.ylabel('Distance (m)')
# ax.set_zlabel('Height (m)')

# # ax.scatter(axm, aym, apm, c='yellow', s=500.0)

# # ax.plot([564126.0, 564126.0],[5.11825e+06, 5.11825e+06],zs=[3837.82, 1866.98])

# # plt.savefig('/home/doug/shared/'+NAME+'.png',dpi=150)
# for ii in xrange(0, 360, 10):
# 	ax.view_init(elev=10.0, azim=ii)
# 	fig.set_size_inches(20.0, 10.0)
# 	# plt.savefig("/home/doug/mygo/src/SCOOP/movies/movie"+str(ii)+".png", dpi=150)

# # TT = 0
# for i in range(len(x)):
# 	# counter checking
# 	# print i
# 	for j in range(len(y)):
# 		p = lines[i].split()
# 		q = liness[i].split()
# 		# if float(p[j]) == -9999.0:
# 		# 	print x[i],",", y[j],",", 900.0
# 		# else: 
# 		# print x[i],",", y[j],",", q[j]
# 		xm.append(x[i])
# 		ym.append(y[j])
# 		if abs(float(p[j]) - float(q[j])) > 0.0:
# 			pm.append(q[j])
# 			print x[i],",", y[j],",", q[j]
# 		elif float(p[j]) - float(q[j]) < 0.0:
# 			pm.append(0)
# 		else:
# 			pm.append(0)


# axm=np.array(xm, dtype=float)
# aym=np.array(ym, dtype=float)
# apm=np.array(pm, dtype=float)

# # ax.scatter(axm, aym, apm, c='yellow', s=8.0)

# ax.axis('equal')

# plt.show()

# # for i in range(len(x)):
# # 	# counter checking
# # 	print i
# # 	for j in range(len(y)):
# # 		p = lines[i].split()
# # 		q = lines1[i].split()
# # 		ax.plot([float(x[i])],[float(y[j])],[float(p[j])], markerfacecolor= colour, marker='o')
# # 		ax.plot([float(x[i])],[float(y[j])],[float(q[j])], markerfacecolor= colour1, marker='o')


# # X, Y = np.meshgrid(x, y)
# # appm = X + Y
# # print apm
# # k = 0
# # for i in range(len(X)):
# # 	for j in range(len(X[0])):
# # 		print appm[i][j]
# # 		appm[i][j] = apm[k]
# # 		print k
# # 		k =  k + 1


# # print X
# # print Y
# # print appm
	
# # ax.plot_surface(X, Y, appm, rstride = 3, cstride = 3, alpha = 0.3, cmap = cm.BuPu)
# # plt.show()

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

# Read the file. 
f2 = open('sthel_res100mDEM.asc', 'r')
lines = f2.readlines()
f2.close()

f3 = open('PaperFig.asc', 'r')
# f3 = open('test0_newDEM_out.asc', 'r')

liness = f3.readlines()
f3.close()

xF = []
yF = []
zF = []

# grid size
gdsx = 91
gdsy = 96

#cell size
csix = 100
csiy = 100

x = np.linspace(557974.631687, 557974.631687+gdsx*csix, gdsx)
y = np.linspace(5111446.1291019, 5111446.1291019+gdsy*csiy, gdsy)


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
			print x[i],",", y[j],",", q[j]
			print x[i],",", y[j],",", p[j]
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
		if float(p[j]) == -9999.0:
			bbpm.append(1000)
		if float(p[j]) > 0.0:
			bbpm.append(p[j])

bxm=np.array(bbxm, dtype=float)
bym=np.array(bbym, dtype=float)
bpm=np.array(bbpm, dtype=float)

# ax.scatter(bxm, bym, bpm, c='yellow', s=20)

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

##################################MOVIE##############################33
##################################MOVIE##############################33
# plt.savefig('/home/doug/shared/'+NAME+'.png',dpi=150)
# for ii in xrange(0, 360, 10):
# 	ax.view_init(elev=10.0, azim=ii)
# 	fig.set_size_inches(20.0, 10.0)
	# plt.savefig("/home/doug/mygo/src/SCOOP/movies/movie"+str(ii)+".png", dpi=150)
##################################MOVIE##############################33
##################################MOVIE##############################33
ax.axis('equal')
# ax.auto_scale_xyz([min(bxm), max(bxm)], [min(bym), max(bym)], [min(bpm), max(bpm)])

plt.show()

