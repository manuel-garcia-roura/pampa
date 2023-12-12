import random

def main():
   
   dx = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
   dy = [0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5]
   dz = [0.25, 0.25, 0.25, 0.25, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25]
   
   dh = 0.05
   
   nx = len(dx)
   ny = len(dy)
   nz = len(dz)
   
   with open("cartesian_mesh.dat", "w") as f:
      
      f.write("# x-discretization:\n")
      f.write("dx %d\n" % nx)
      for i, d in enumerate(dx):
         f.write("%.3f" % d)
         if i < nx-1:
            f.write(" ")
         else:
            f.write("\n")
      f.write("\n")
      
      f.write("# y-discretization:\n")
      f.write("dy %d\n" % ny)
      for j, d in enumerate(dy):
         f.write("%.3f" % d)
         if j < ny-1:
            f.write(" ")
         else:
            f.write("\n")
      f.write("\n")
      
      f.write("# z-discretization:\n")
      f.write("dz %d\n" % nz)
      for k, d in enumerate(dz):
         f.write("%.3f" % d)
         if k < nz-1:
            f.write(" ")
         else:
            f.write("\n")
   
   with open("unstructured_extruded_mesh.dat", "w") as f:
      
      f.write("# xy-points:\n")
      f.write("points %d\n" % ((nx+1)*(ny+1)))
      x = [0.0] * (nx+1)
      for i in range(nx):
         x[i+1] = x[i] + dx[i]
      y = [0.0] * (ny+1)
      for j in range(ny):
         y[j+1] = y[j] + dy[j]
      for j in range(ny+1):
         for i in range(nx+1):
            f.write("%.3f %.3f\n" % (x[i]+random.uniform(-dh, dh), y[j]+random.uniform(-dh, dh)))
      f.write("\n")
      
      f.write("# xy-cells:\n")
      f.write("cells %d\n" % (nx*ny))
      for j in range(ny):
         for i in range(nx):
            p1 = i + j*(nx+1)
            p2 = (i+1) + j*(nx+1)
            p3 = (i+1) + (j+1)*(nx+1)
            p4 = i + (j+1)*(nx+1)
            f.write("4 %d %d %d %d\n" % (p1, p2, p3, p4))
      f.write("\n")
      
      f.write("# z-discretization:\n")
      f.write("dz %d\n" % len(dz))
      for i, d in enumerate(dz):
         f.write("%.3f" % d)
         if i < len(dz)-1:
            f.write(" ")
         else:
            f.write("\n")

if __name__ == '__main__': main()
