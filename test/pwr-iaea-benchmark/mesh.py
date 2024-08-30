import random
import shutil

def main():
   
   methods = ["diffusion", "diffusion", "sn", "precursors", "conduction"]
   dimensions = [2, 3, 2, 3, 3]
   full_core = False
   
   for method, dims in zip(methods, dimensions):
      
      dx = [10.0] * 17
      dy = [10.0] * 17
      if dims == 3:
         dz = [12.0] * 17
         n = 2
      elif dims == 2:
         n = 4
      else:
         raise RuntimeError("Wrong number of dimensions!")
      dh = 0.0 / n
      
      if method == "diffusion":
         bc_ext = "robin 1 -0.4692"
      elif method == "sn":
         bc_ext = "vacuum"
      elif method == "conduction":
         bc_ext = "dirichlet 1 300.0"
      elif method == "precursors":
         bc_ext = None
      else:
         raise RuntimeError("Wrong method!")
      
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
      
      nx = len(dx)
      ny = len(dy)
      if dims == 3:
         nz = len(dz)
      else:
         nz = 1
      
      if dims == 3:
         filename = "cartesian-" + method + "-3d/mesh.pmp"
      else:
         filename = "cartesian-" + method + "/mesh.pmp"
      
      with open(filename, "w") as f:
         
         f.write("dx %d\n" % nx)
         for i, d in enumerate(dx):
            if i > 0: f.write(" ")
            f.write("%.3f" % d)
         f.write("\n\n")
         
         f.write("dy %d\n" % ny)
         for j, d in enumerate(dy):
            if j > 0: f.write(" ")
            f.write("%.3f" % d)
         f.write("\n\n")
         
         if dims == 3:
            f.write("dz %d\n" % nz)
            for k, d in enumerate(dz):
               if k > 0: f.write(" ")
               f.write("%.3f" % d)
            f.write("\n\n")
         
         if bc_ext != None:
            if full_core:
               f.write("bc -x %s\n" % bc_ext)
               f.write("bc +x %s\n" % bc_ext)
               f.write("bc -y %s\n" % bc_ext)
               f.write("bc +y %s\n" % bc_ext)
               if dims == 3:
                  f.write("bc -z %s\n" % bc_ext)
                  f.write("bc +z %s\n" % bc_ext)
            else:
               f.write("bc -x reflective\n")
               f.write("bc +x %s\n" % bc_ext)
               f.write("bc -y reflective\n")
               f.write("bc +y %s\n" % bc_ext)
               if dims == 3:
                  f.write("bc -z reflective\n")
                  f.write("bc +z %s\n" % bc_ext)
            f.write("\n")
         
         f.write("materials %d\n" % (nx*ny*nz))
         for k in range(nz):
            f.write("\n")
            for j in range(ny):
               for i in range(nx):
                  if i > 0: f.write(" ")
                  mat = layout_xy[layout_z[k]][j][i]
                  f.write("%d" % mat)
               f.write("\n")
      
      if method == "conduction" or method == "precursors":
         return
      
      if dims == 3:
         filename = "unstructured-" + method + "-3d/mesh.pmp"
      else:
         filename = "unstructured-" + method + "/mesh.pmp"
      
      with open(filename, "w") as f:
         
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
               if (j < ny and i > 0): mats.append(layout_xy[layout_z[0]][j][i-1])
               if (not all(m == 0 for m in mats)):
                  if (any(m == 0 for m in mats)):
                     bc_1_points.append(j*(nx+1)+i)
                  elif (i == nx or j == ny):
                     bc_1_points.append(j*(nx+1)+i)
         f.write("\n")
         
         num_xy_cells = 0
         for j in range(ny):
            for i in range(nx):
               if (layout_xy[layout_z[0]][j][i] != 0): num_xy_cells += 1
         f.write("cells %d %d\n" % (num_xy_cells, 4*num_xy_cells))
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
            f.write("dz %d\n" % nz)
            for k, d in enumerate(dz):
               if k > 0: f.write(" ")
               f.write("%.3f" % d)
            f.write("\n\n")
         
         if full_core:
            f.write("boundary s1 %d\n" % len(bc_1_points))
         else:
            f.write("boundary exterior %d\n" % len(bc_1_points))
         for i in bc_1_points:
            f.write("%d\n" % i)
         f.write("\n")
         
         if full_core:
            f.write("boundary s2 %d\n" % len(bc_1_points))
         else:
            f.write("boundary interior %d\n" % len(bc_1_points))
         for i in bc_2_points:
            f.write("%d\n" % i)
         f.write("\n")
         
         if full_core:
            f.write("bc 1 %s\n" % bc_ext)
            f.write("bc 2 %s\n" % bc_ext)
            if dims == 3:
               f.write("bc -z %s\n" % bc_ext)
               f.write("bc +z %s\n" % bc_ext)
         else:
            f.write("bc 1 %s\n" % bc_ext)
            f.write("bc 2 reflective\n")
            if dims == 3:
               f.write("bc -z reflective\n")
               f.write("bc +z %s\n" % bc_ext)
         f.write("\n")
         
         f.write("materials %d\n" % (num_xy_cells*nz))
         for k in range(nz):
            f.write("\n")
            for j in range(ny):
               for i in range(nx):
                  mat = layout_xy[layout_z[k]][j][i]
                  if (layout_xy[layout_z[0]][j][i] != 0):
                     if i > 0: f.write(" ")
                     f.write("%d" % mat)
               f.write("\n")

if __name__ == "__main__": main()
