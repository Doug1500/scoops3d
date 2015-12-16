# Read the file. 
f2 = open('A_emb20.scp', 'r')
f3 = open('A_emb20_out.txt', 'r')

# read the whole file into a single variable, which is a list of every row of the file.
lines = f2.readlines()
outlis = f3.readlines()
f2.close()
f3.close()

FlSu = []
FOS = []

for i in range(len(lines)):
	p = lines[i].split()
	if p[0] == "xcen":
		ap = lines[i+1].split()
		FlSu = ap 
		

for ouli in outlis:
	q = ouli.split()
	# print q
	if len(q)>0:
		if q[0] == "Bishop's":
			FOS.append(q[5])

print "Failure surface geometry", FlSu
print "Simplified Bishop FOS is:", FOS[0]