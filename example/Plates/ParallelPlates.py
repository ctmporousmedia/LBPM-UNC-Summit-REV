import numpy

nx=120
ny=30
nz=30
N=nx*ny*nz

mesh=(nx,ny,nz)
data=numpy.ones(mesh,dtype=numpy.int8)

LabelTop=-1
LabelBottom=-1

#print(data)
print("Writing parllel plates")
print("Mesh size: "+repr(mesh))
print("Top label: "+repr(LabelTop))
print("Bottom label: "+repr(LabelBottom))

# assign a bubble in the middle
for x in range(nx-40,nx):
  for y in range(0,ny):
    for z in range(0,nz):
       data[x,y,z]=2

# solid walls with different component labels
for x in range(0,nx):
  for y in range(0,ny):
    data[x,y,0]=LabelBottom
    data[x,y,1]=LabelBottom
    data[x,y,2]=LabelBottom
    data[x,y,nz-1]=LabelTop
    data[x,y,nz-2]=LabelTop
    data[x,y,nz-3]=LabelTop


for x in range(0,nx):
  for z in range(0,nz):
    data[x,0,z]=LabelBottom
    data[x,1,z]=LabelBottom
    data[x,2,z]=LabelBottom
    data[x,ny-1,z]=LabelTop
    data[x,ny-2,z]=LabelTop
    data[x,ny-3,z]=LabelTop


data.tofile("ParallelPlates.raw")
