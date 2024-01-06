import random

def main():
   
   dims = 2
   bc_robin = False
   full_core = False
   
   dx = [10.0] * 17
   dy = [10.0] * 17
   if dims == 3:
      dz = [12.0] * 17
      n = 2
   else:
      dz = [1.0]
      n = 8
   
   dh = 1.0 / n
   
   layout_xy = [[3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 1, 1, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4, 4, 4], 
                [3, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 4, 4, 4, 4], 
                [3, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 4, 4, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4], 
                [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4], 
                [2, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
                [2, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4], 
                [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]]
   
   if dims == 3:
      layout_z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4]
   else:
      layout_z = [0]
   
   if full_core:
      dx = dx[::-1] + dx
      dy = dy[::-1] + dy
      if dims == 3:
         dz = dz[::-1] + dz
      layout_xy = [l[::-1] + l for l in layout_xy]
      layout_xy = layout_xy[::-1] + layout_xy
      if dims == 3:
         layout_z = layout_z[::-1] + layout_z
   
   if n > 1:
      dx = [dxi/n for dxi in dx for _ in range(n)]
      dy = [dyi/n for dyi in dy for _ in range(n)]
      if dims == 3:
         dz = [dzi/n for dzi in dz for _ in range(n)]
      for i in range(len(layout_xy)):
         layout_xy[i] = [idx for idx in layout_xy[i] for _ in range(n)]
      layout_xy = [ids for ids in layout_xy for _ in range(n)]
      if dims == 3:
         layout_z = [ids for ids in layout_z for _ in range(n)]
   
   if bc_robin:
      bc_ext = "3 -0.4692"
   else:
      bc_ext = "1"
   
   nx = len(dx)
   ny = len(dy)
   nz = len(dz)
   
   with open("cartesian/mesh.pmp", "w") as f:
      
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
      f.write("\n")
      
      f.write("# boundary conditions:\n")
      if full_core:
         f.write("bc x %s %s\n" % (bc_ext, bc_ext))
         f.write("bc y %s %s\n" % (bc_ext, bc_ext))
         if dims == 3:
            f.write("bc z %s %s\n" % (bc_ext, bc_ext))
         else:
            f.write("bc z 2 2\n")
      else:
         f.write("bc x 2 %s\n" % bc_ext)
         f.write("bc y 2 %s\n" % bc_ext)
         if dims == 3:
            f.write("bc z 2 %s\n" % bc_ext)
         else:
            f.write("bc z 2 2\n")
      f.write("\n")
      
      f.write("# material distribution:\n")
      f.write("materials %d\n" % (nx*ny*nz))
      for k in range(nz):
         f.write("\n")
         for j in range(ny):
            for i in range(nx):
               mat = layout_xy[j][i] if (layout_z[k] == 0) else layout_z[k]
               f.write("%d" % mat)
               if i < nx-1:
                  f.write(" ")
               else:
                  f.write("\n")
   
   with open("unstructured-extruded/mesh.pmp", "w") as f:
      
      f.write("# xy-points:\n")
      f.write("points %d\n" % ((nx+1)*(ny+1)))
      x = [0.0] * (nx+1)
      for i in range(nx):
         x[i+1] = x[i] + dx[i]
      y = [0.0] * (ny+1)
      for j in range(ny):
         y[j+1] = y[j] + dy[j]
      bc_1_points = []
      bc_2_points = []
      for j in range(ny+1):
         for i in range(nx+1):
            f.write("%.3f %.3f\n" % (x[i]+random.uniform(-dh, dh), y[j]+random.uniform(-dh, dh)))
            if (i == nx or j == ny):
               bc_1_points.append(j*(nx+1)+i)
            if (i == 0 or j == 0):
               bc_2_points.append(j*(nx+1)+i)
      f.write("\n")
      
      f.write("# xy-cells:\n")
      f.write("cells %d\n" % (nx*ny))
      for j in range(ny):
         for i in range(nx):
            p1 = i + j*(nx+1)
            p2 = (i+1) + j*(nx+1)
            p3 = (i+1) + (j+1)*(nx+1)
            p4 = i + (j+1)*(nx+1)
            f.write("%d %d %d %d\n" % (p1, p2, p3, p4))
      f.write("\n")
      
      f.write("# z-discretization:\n")
      f.write("dz %d\n" % nz)
      for k, d in enumerate(dz):
         f.write("%.3f" % d)
         if k < nz-1:
            f.write(" ")
         else:
            f.write("\n")
      f.write("\n")
      
      f.write("# zero-flux boundary:\n")
      f.write("boundary %d\n" % (nx+ny+1))
      for i in bc_1_points:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("# reflective boundary:\n")
      f.write("boundary %d\n" % (nx+ny+1))
      for i in bc_2_points:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("# boundary conditions:\n")
      if full_core:
         f.write("bc 1 %s\n" % bc_ext)
         f.write("bc 2 %s\n" % bc_ext)
         if dims == 3:
            f.write("bc z %s %s\n" % (bc_ext, bc_ext))
         else:
            f.write("bc z 2 2\n")
      else:
         f.write("bc 1 %s\n" % bc_ext)
         f.write("bc 2 2\n")
         if dims == 3:
            f.write("bc z 2 %s\n" % bc_ext)
         else:
            f.write("bc z 2 2\n")
      f.write("\n")
      
      f.write("# material distribution:\n")
      f.write("materials %d\n" % (nx*ny*nz))
      for k in range(nz):
         f.write("\n")
         for j in range(ny):
            for i in range(nx):
               mat = layout_xy[j][i] if (layout_z[k] == 0) else layout_z[k]
               f.write("%d" % mat)
               if i < nx-1:
                  f.write(" ")
               else:
                  f.write("\n")

if __name__ == '__main__': main()
