import random
import shutil

def main():
   
   dims = 2
   bc_robin = False
   full_core = True
   
   dx = [10.0] * 17
   dy = [10.0] * 17
   if dims == 3:
      dz = [12.0] * 17
      if full_core:
         n = 2
      else:
         n = 4
   else:
      n = 4
   
   dh = 0.0 / n
   
   layout_xy = [None] * 2
   
   layout_xy[0] = [[3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 1, 1, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4, 4, 4], 
                   [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 4, 4, 4, 4], 
                   [3, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 4, 4, 0, 0], 
                   [3, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 1, 4, 4, 0, 0], 
                   [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0], 
                   [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0], 
                   [2, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0, 0, 0], 
                   [2, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 0, 0, 0, 0], 
                   [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0], 
                   [1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   
   layout_xy[1] = [[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                   [4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   
   if dims == 3:
      layout_z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
   else:
      layout_z = [0]
   
   if full_core:
      dx = dx[::-1] + dx
      dy = dy[::-1] + dy
      if dims == 3:
         dz = dz[::-1] + dz
      for k in range(len(layout_xy)):
         layout_xy[k] = [l[::-1] + l for l in layout_xy[k]]
         layout_xy[k] = layout_xy[k][::-1] + layout_xy[k]
      if dims == 3:
         layout_z = layout_z[::-1] + layout_z
   
   if n > 1:
      dx = [dxi/n for dxi in dx for _ in range(n)]
      dy = [dyi/n for dyi in dy for _ in range(n)]
      if dims == 3:
         dz = [dzi/n for dzi in dz for _ in range(n)]
      for k in range(len(layout_xy)):
         for i in range(len(layout_xy[k])):
            layout_xy[k][i] = [idx for idx in layout_xy[k][i] for _ in range(n)]
         layout_xy[k] = [ids for ids in layout_xy[k] for _ in range(n)]
      if dims == 3:
         layout_z = [ids for ids in layout_z for _ in range(n)]
   
   if bc_robin:
      bc_ext = "3 -0.4692"
   else:
      bc_ext = "1"
   
   nx = len(dx)
   ny = len(dy)
   if dims == 3:
      nz = len(dz)
   else:
      nz = 1
   
   with open("cartesian-diffusion/mesh.pmp", "w") as f:
      
      f.write("# x-discretization:\n")
      f.write("dx %d\n" % nx)
      for i, d in enumerate(dx):
         f.write("%.3f" % d)
         if i < nx-1: f.write(" ")
      f.write("\n\n")
      
      f.write("# y-discretization:\n")
      f.write("dy %d\n" % ny)
      for j, d in enumerate(dy):
         f.write("%.3f" % d)
         if j < ny-1: f.write(" ")
      f.write("\n\n")
      
      if dims == 3:
         f.write("# z-discretization:\n")
         f.write("dz %d\n" % nz)
         for k, d in enumerate(dz):
            f.write("%.3f" % d)
            if k < nz-1: f.write(" ")
         f.write("\n\n")
      
      f.write("# boundary conditions:\n")
      if full_core:
         f.write("bc x %s %s\n" % (bc_ext, bc_ext))
         f.write("bc y %s %s\n" % (bc_ext, bc_ext))
         if dims == 3:
            f.write("bc z %s %s\n" % (bc_ext, bc_ext))
      else:
         f.write("bc x 2 %s\n" % bc_ext)
         f.write("bc y 2 %s\n" % bc_ext)
         if dims == 3:
            f.write("bc z 2 %s\n" % bc_ext)
      f.write("\n")
      
      f.write("# material distribution:\n")
      f.write("materials %d\n" % (nx*ny*nz))
      for k in range(nz):
         f.write("\n")
         for j in range(ny):
            for i in range(nx):
               mat = layout_xy[layout_z[k]][j][i]
               f.write("%d" % mat)
               if i < nx-1: f.write(" ")
            f.write("\n")
   
   with open("unstructured-extruded-diffusion/mesh.pmp", "w") as f:
      
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
            if (i == 0 or j == 0):
               bc_2_points.append(j*(nx+1)+i)
            mats = []
            if (j > 0 and i > 0): mats.append(layout_xy[layout_z[0]][j-1][i-1])
            if (j > 0 and i < nx): mats.append(layout_xy[layout_z[0]][j-1][i])
            if (j < ny and i < nx): mats.append(layout_xy[layout_z[0]][j][i])
            if (j < ny and i > 0):   mats.append(layout_xy[layout_z[0]][j][i-1])
            if (not all(m == 0 for m in mats)):
               if (any(m == 0 for m in mats)):
                  bc_1_points.append(j*(nx+1)+i)
               elif (i == nx or j == ny):
                  bc_1_points.append(j*(nx+1)+i)
      f.write("\n")
      
      f.write("# xy-cells:\n")
      num_xy_cells = 0
      for j in range(ny):
         for i in range(nx):
            if (layout_xy[layout_z[0]][j][i] != 0): num_xy_cells += 1
      f.write("cells %d\n" % num_xy_cells)
      for j in range(ny):
         for i in range(nx):
            if (layout_xy[layout_z[0]][j][i] != 0):
               p1 = i + j*(nx+1)
               p2 = (i+1) + j*(nx+1)
               p3 = (i+1) + (j+1)*(nx+1)
               p4 = i + (j+1)*(nx+1)
               f.write("%d %d %d %d\n" % (p1, p2, p3, p4))
      f.write("\n")
      
      if dims == 3:
         f.write("# z-discretization:\n")
         f.write("dz %d\n" % nz)
         for k, d in enumerate(dz):
            f.write("%.3f" % d)
            if k < nz-1: f.write(" ")
         f.write("\n\n")
      
      f.write("# zero-flux boundary:\n")
      f.write("boundary %d\n" % len(bc_1_points))
      for i in bc_1_points:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("# reflective boundary:\n")
      f.write("boundary %d\n" % len(bc_2_points))
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
         f.write("bc 1 %s\n" % bc_ext)
         f.write("bc 2 2\n")
         if dims == 3:
            f.write("bc z 2 %s\n" % bc_ext)
      f.write("\n")
      
      f.write("# material distribution:\n")
      f.write("materials %d\n" % (num_xy_cells*nz))
      for k in range(nz):
         f.write("\n")
         for j in range(ny):
            for i in range(nx):
               mat = layout_xy[layout_z[k]][j][i]
               if (layout_xy[layout_z[0]][j][i] != 0):
                  f.write("%d" % mat)
                  if i < nx-1: f.write(" ")
            f.write("\n")
   
   shutil.copyfile("cartesian-diffusion/mesh.pmp", "cartesian-sn/mesh.pmp")

if __name__ == '__main__': main()
