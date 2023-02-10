import vtk
from vtk.util import numpy_support as VN
import numpy as np
import matplotlib.pyplot as plt

filename = "../output/out.vtk.0"
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName(filename)
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()
data = reader.GetOutput()

dim = data.GetDimensions()
dim = list(dim)
dim.append(3)

velocity = VN.vtk_to_numpy(data.GetPointData().GetArray("Velocity"))
velocity = velocity.reshape(dim, order='F')

x3D = np.zeros(data.GetNumberOfPoints())
y3D = np.zeros(data.GetNumberOfPoints())
z3D = np.zeros(data.GetNumberOfPoints())
for i in range(data.GetNumberOfPoints()):
        x3D[i],y3D[i],z3D[i] = data.GetPoint(i)
x3D = x3D.reshape(dim[:-1], order='F')
y3D = y3D.reshape(dim[:-1], order='F')
z3D = z3D.reshape(dim[:-1], order='F')

u = velocity[:,:,1,0]
v = velocity[:,:,1,1]
w = velocity[:,:,1,2]
x = x3D[:,:,1]
y = y3D[:,:,1]

uFlip = np.fliplr(u)
vFlip = np.fliplr(v)

plt.figure(1)
plt.contourf(x,y,u-uFlip)
plt.colorbar()
plt.title("u")
plt.figure(2)
plt.contourf(x,y,v+vFlip)
plt.colorbar()
plt.title("v")
plt.figure(3)
plt.contourf(x,y,w)
plt.colorbar()
plt.title("w")
#plt.show()

print("norm, u: {norm:.4f}".format(norm=np.linalg.norm(u-uFlip,ord=2)))
print("norm, v: {norm:.4f}".format(norm=np.linalg.norm(v+vFlip,ord=2)))
print("norm, w: {norm:.4f}".format(norm=np.linalg.norm(w,      ord=2)))
