import numpy

len=220
N=len*len*len

mesh=(len,len,len)
data=numpy.ones(mesh,dtype=numpy.int8)

c = 55
rad = 90
d = 110.01
rad2 = 27
e = 68


## assign a bubble in the middle
for x in range(0,len):
  for y in range(0,len):
    for z in range(0,len):
        data[x,y,z]=2


## assign a bubble in the middle
#for x in range(10,30):
#  for y in range(10,30):
#    for z in range(10,30):
#        data[x,y,z]=1


## assign a bubble in the middle
for x in range(0,len):
  for y in range(0,len):
    for z in range(0,len):
        if ( (x-d)*(x-d) + (y-d)*(y-d) + (z-d)*(z-d) < rad*rad):
           data[x,y,z]=-1

#for x in range(0,len):
#  for y in range(0,len):
#    for z in range(0,len):
#        if ( (x-d)*(x-d) + (y-d)*(y-d) + (z-d)*(z-d) < rad2*rad2):
#           data[x,y,z]=-1


#for x in range(0,len):
#  for y in range(0,len):
#    for z in range(0,len):
#        if ( (x-e)*(x-e) + (y-e)*(y-e) + (z-e)*(z-e) < rad2*rad2):
#           data[x,y,z]=-1

#for x in range(0,len):
#  for y in range(0,len):
#    for z in range(0,len):
#        if ( (x-d)*(x-d) + (y-d)*(y-d) + (z-d)*(z-d) < rad2*rad2):
#           data[x,y,z]=2



#a = int(len/4)
#b = int(3*len/4)
#
## assign a bubble in the middle
#for x in range(0,a):
#  for y in range(0,a):
#    for z in range(0,a):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(b,len):
#  for y in range(0,a):
#    for z in range(0,a):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(0,a):
#  for y in range(b,len):
#    for z in range(0,a):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(0,a):
#  for y in range(0,a):
#    for z in range(b,len):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(b,len):
#  for y in range(0,a):
#    for z in range(b,len):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(b,len):
#  for y in range(b,len):
#    for z in range(0,a):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(b,len):
#  for y in range(b,len):
#    for z in range(b,len):
#         data[x,y,z]=2
#
#
## assign a bubble in the middle
#for x in range(0,a):
#  for y in range(b,len):
#    for z in range(b,len):
#         data[x,y,z]=2
#
#
#
#
#
#
data.tofile("Bubble.raw")
