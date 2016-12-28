#!/usr/bin/python
import sys
import numpy as np

def printf(format, *args):
    sys.stdout.write(format % args)

# The position of the atoms
atoms_x=[]
atoms_y=[]
atoms_z=[]

stri = ''.join(file(sys.argv[1]).readlines()[2])
data = np.fromstring(stri, sep=' ')
num_atoms= int(data[0])  # first number in 3-rd row
start_point=[data[1], data[2], data[3]]

stri = ''.join(file(sys.argv[1]).readlines()[3])
data = np.fromstring(stri, sep=' ')
NX = int(data[0])
dx = [data[1], data[2], data[3]]
stri = ''.join(file(sys.argv[1]).readlines()[4])
data = np.fromstring(stri, sep=' ')
NY =  int(data[0])
dy =  [data[1], data[2], data[3]]
stri = ''.join(file(sys.argv[1]).readlines()[5])
data = np.fromstring(stri, sep=' ')
NZ = int(data[0])
dz = [data[1], data[2], data[3]]
for i in range(num_atoms):
   stri = ''.join(file(sys.argv[1]).readlines()[6+i])
   data = np.fromstring(stri, sep=' ')
   atoms_x.append( data[2] ) #line 7+i, 2. row
   atoms_y.append( data[3] ) #line 7+i, 3. row
   atoms_z.append( data[4] ) #line 7+i, 4. row

linenum=0
pt_num=0

with open(sys.argv[1]) as inp:
    for line in inp:
       if (linenum<6+num_atoms):
         linenum+=1
       else:
          data = np.fromstring(line, sep=' ')
          for k in range(len(data)):
             printf("%.5f   "%(start_point[0]+dx[0]*(pt_num//(NY*NZ))))
             printf("%.5f   "%(start_point[1]+dy[1]*((pt_num//NZ)%NY)))
             printf("%.5f   %.5g\n"%(start_point[2]+dz[2]*(pt_num%NZ),data[k]))
             pt_num+=1


#stri = ''.join(file(sys.argv[1]).readlines()[6+num_atoms:])
#data = np.fromstring(stri, sep=' ')
##data.shape = (NX, NY, NZ)
#x,y,z=np.mgrid[start_point[0]:start_point[0]+dx[0]*NX+dy[0]*NY+dz[0]*NZ:NX*1j,
#               start_point[1]:start_point[1]+dx[1]*NX+dy[1]*NY+dz[1]*NZ:NY*1j,
#               start_point[2]:start_point[2]+dx[2]*NX+dy[2]*NY+dz[2]*NZ:NZ*1j]
#for i in range(x.size):
#   for j in range(y.size):
#      for k in range(z.size):
#         print(x[i],y[j],z[k],data[(i*x.size+j)*z.size+k])
