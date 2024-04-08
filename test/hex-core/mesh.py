from typing import NamedTuple
import numpy as np
import matplotlib.pyplot as plt

class Mesh(NamedTuple):
   
   nx: int
   ny: int
   x: list
   y: list
   xc: list
   yc: list
   pt_idx: list
   hex_cells: list
   tri_cells: list
   hex_mats: list
   tri_mats: list
   bc_pts: list

def build_y_hex_mesh(p, layout):
   
   nx = 2*len(layout[0]) + 1
   ny = 2*len(layout) + len(layout[0])%2 - 2
   
   y, x = np.meshgrid(np.arange(float(ny), 0.0, -1.0), np.arange(float(nx)))
   y[1::2, :] += 0.5
   dx = 0.5 * p
   dy = (2.0/np.sqrt(3)) * 0.5 * p
   y = dy * y
   x = dx * x
   
   xc = []; yc = []
   ic = []; jc = []
   hex_cells = []
   tri_cells = []
   hex_mats = []
   tri_mats = []
   n = np.zeros((nx, ny))
   for (j, cj) in enumerate(layout):
      for (i, c) in enumerate(cj):
         if c != 0:
            
            i0 = 2*i + j%2 + 1
            j0 = int(1.5*j) + 1
            xc.append(x[i0, j0])
            yc.append(y[i0, j0])
            ic.append(i0)
            jc.append(j0)
            n[i0, j0] += 3
            
            ip = [i0-1, i0, i0+1, i0+1, i0, i0-1]
            jp = [j0+j%2, j0+1, j0+j%2, j0+j%2-1, j0-1, j0+j%2-1]
            hex_pts = []
            for k in range(6):
               hex_pts.append((ip[k], jp[k]))
               tri_pts = []
               tri_pts.append((ip[k], jp[k]))
               tri_pts.append((ip[(k+1)%6], jp[(k+1)%6]))
               tri_pts.append((i0, j0))
               tri_cells.append(tri_pts)
               tri_mats.append(c)
               n[ip[k], jp[k]] += 1
            hex_cells.append(hex_pts)
            hex_mats.append(c)
   
   pt_idx = -1 * np.ones((nx, ny))
   bc_pts = []
   ip = 0
   for j in range(ny):
      for i in range(nx):
         if n[i, j] > 0:
            pt_idx[i, j] = ip
            if n[i, j] < 3:
               bc_pts.append(ip)
            ip += 1
   
   return Mesh(nx, ny, x, y, xc, yc, pt_idx, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

def build_x_hex_mesh(p, layout):
   
   ny = 2*len(layout[0]) - 3
   nx = 2*len(layout) + len(layout[0])%2 - 6
   
   x, y = np.meshgrid(np.arange(float(nx), 0.0, -1.0), np.arange(float(ny)))
   x[1::2, :] += 0.5
   dy = 0.5 * p
   dx = (2.0/np.sqrt(3)) * 0.5 * p
   x = dx * x
   y = dy * y
   
   yc = []; xc = []
   jc = []; ic = []
   hex_cells = []
   tri_cells = []
   hex_mats = []
   tri_mats = []
   n = np.zeros((ny, nx))
   for (i, ci) in enumerate(layout):
      for (j, c) in enumerate(ci):
         if c != 0:
            
            j0 = 2*j + i%2 - 1
            i0 = int(1.5*i) + 1
            yc.append(y[j0, i0])
            xc.append(x[j0, i0])
            jc.append(j0)
            ic.append(i0)
            n[j0, i0] += 3
            
            jp = [j0-1, j0, j0+1, j0+1, j0, j0-1]
            ip = [i0+i%2, i0+1, i0+i%2, i0+i%2-1, i0-1, i0+i%2-1]
            jp.reverse()
            ip.reverse()
            hex_pts = []
            for j in range(6):
               hex_pts.append((jp[j], ip[j]))
               tri_pts = []
               tri_pts.append((jp[j], ip[j]))
               tri_pts.append((jp[(j+1)%6], ip[(j+1)%6]))
               tri_pts.append((j0, i0))
               tri_cells.append(tri_pts)
               tri_mats.append(c)
               n[jp[j], ip[j]] += 1
            hex_cells.append(hex_pts)
            hex_mats.append(c)
   
   pt_idx = -1 * np.ones((ny, nx))
   bc_pts = []
   jp = 0
   for i in range(nx):
      for j in range(ny):
         if n[j, i] > 0:
            pt_idx[j, i] = jp
            if n[j, i] < 3:
               bc_pts.append(jp)
            jp += 1
   
   return Mesh(ny, nx, x, y, xc, yc, pt_idx, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

def write_mesh(filename, mesh, cells, mats, np):
   
   with open(filename, "w") as f:
      
      f.write("points %d\n" % (mesh.pt_idx >= 0).sum())
      for j in range(mesh.ny):
         for i in range(mesh.nx):
            if (mesh.pt_idx[i, j] >= 0):
               f.write("%.3f %.3f\n" % (mesh.x[i, j], mesh.y[i, j]))
      f.write("\n")
      
      f.write("cells %d %d\n" % (len(cells), np*len(cells)))
      for c in cells:
         for i, p in enumerate(c):
            if i > 0: f.write(" ")
            f.write("%d" % mesh.pt_idx[p])
         f.write("\n")
      f.write("\n")
      
      f.write("dz -16\n")
      f.write("12.0\n")
      f.write("\n")
      
      f.write("boundary %d\n" % len(mesh.bc_pts))
      for i in mesh.bc_pts:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("materials %d\n" % (16*len(cells)))
      for k in range(16):
         f.write("\n")
         for i in range(len(mats)):
            if i > 0: f.write(" ")
            f.write("%d" % mats[i])
         f.write("\n")

def plot_mesh(mesh):
   
   fig, ax = plt.subplots()
   
   ax.scatter(mesh.x, mesh.y, s = 8, c = "blue", alpha = 0.25)
   
   # for (i, (xci, yci)) in enumerate(zip(mesh.xc, mesh.yc)):
   #    ax.annotate(i, (xci, yci))
   
   for c in mesh.tri_cells:
      xp = []; yp = []
      for p in c:
         xp.append(mesh.x[p[0], p[1]])
         yp.append(mesh.y[p[0], p[1]])
      plt.fill(xp, yp, facecolor = "none", edgecolor = "purple", linewidth = 1, alpha = 1.0)
      plt.fill(xp, yp, facecolor = "purple", edgecolor = "none", linewidth = 1, alpha = 0.25)
   
   ax.axis("equal")
   plt.show()

def main():
   
   small = False
   
   pc = 12.0
   if small:
      core = [[0, 1, 1, 1, 0], \
                [1, 1, 1, 1], \
               [1, 1, 1, 1, 1], \
                 [1, 1, 1, 1], \
                [0, 1, 1, 1, 0]]
   else:
      core = [[0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0], \
                [0, 0, 0, 4, 4, 1, 1, 1, 1, 1, 4, 4, 0, 0, 0], \
               [0, 0, 0, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0], \
                 [0, 0, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0], \
                [0, 0, 4, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 4, 0, 0], \
                  [0, 4, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 4, 0], \
                 [0, 4, 1, 1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 1, 4, 0], \
                   [0, 4, 1, 1, 2, 2, 3, 3, 3, 2, 2, 1, 1, 4, 0], \
                  [0, 4, 1, 1, 1, 2, 2, 3, 3, 2, 2, 1, 1, 1, 4, 0], \
                    [0, 4, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 4, 0], \
                   [0, 0, 4, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 4, 0, 0], \
                     [0, 0, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0], \
                    [0, 0, 0, 4, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0], \
                      [0, 0, 0, 4, 4, 1, 1, 1, 1, 1, 4, 4, 0, 0, 0], \
                     [0, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 0]]
   
   pf = 2.0
   if small:
      fa = [[0, 1, 1, 1, 0], \
              [1, 1, 1, 1], \
             [1, 1, 1, 1, 1], \
               [1, 1, 1, 1], \
              [0, 1, 1, 1, 0]]
   else:
      fa = [[0, 0, 1, 1, 1, 1, 0, 0], \
              [0, 1, 1, 1, 1, 1, 0], \
             [0, 1, 1, 1, 1, 1, 1, 0], \
               [1, 1, 1, 1, 1, 1, 1], \
              [0, 1, 1, 1, 1, 1, 1, 0], \
                [0, 1, 1, 1, 1, 1, 0], \
               [0, 0, 1, 1, 1, 1, 0, 0]]
   
   core_mesh = build_x_hex_mesh(pc, core)
   
   write_mesh("hex-cells-diffusion-3d/mesh.pmp", core_mesh, core_mesh.hex_cells, core_mesh.hex_mats, 6)
   write_mesh("tri-cells-diffusion-3d/mesh.pmp", core_mesh, core_mesh.tri_cells, core_mesh.tri_mats, 3)
   
   # plot_mesh(core_mesh)

if __name__ == "__main__": main()
