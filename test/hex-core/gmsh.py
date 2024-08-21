from typing import NamedTuple
import math as m
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import sys

class NodalMesh(NamedTuple):
   
   points: list
   hex_centroids: list
   tri_centroids: list
   hex_cells: list
   tri_cells: list
   hex_mats: list
   tri_mats: list
   bc_pts: list

def build_x_hex_mesh(p, layout):
   
   ny = 2*len(layout[0]) + 1
   nx = int(1.5*len(layout)) + 2
   
   x, y = np.meshgrid(np.arange(float(nx), 0.0, -1.0), np.arange(float(ny)))
   x[1::2, :] += 0.5
   dy = 0.5 * p
   dx = (2.0/np.sqrt(3)) * 0.5 * p
   x = dx * x
   y = dy * y
   
   hex_centroids = []
   tri_centroids = []
   hex_cells = []
   tri_cells = []
   hex_mats = []
   tri_mats = []
   n = np.zeros((ny, nx), dtype = int)
   for (i, ci) in enumerate(layout):
      for (j, c) in enumerate(ci):
         if c != 0:
            
            j0 = 2*j + i%2 + 1
            i0 = int(1.5*i) + 1
            hex_centroids.append([x[j0, i0], y[j0, i0]])
            n[j0, i0] += 3
            
            jp = [j0-1, j0, j0+1, j0+1, j0, j0-1]
            ip = [i0+i%2-1, i0-1, i0+i%2-1, i0+i%2, i0+1, i0+i%2]
            hex_pts = []
            for j in range(6):
               hex_pts.append((jp[j], ip[j]))
               tri_pts = []
               tri_pts.append((jp[j], ip[j]))
               tri_pts.append((jp[(j+1)%6], ip[(j+1)%6]))
               tri_pts.append((j0, i0))
               x0 = 0.0; y0 = 0.0
               for l in range(3):
                  x0 += x[tri_pts[l][0], tri_pts[l][1]]
                  y0 += y[tri_pts[l][0], tri_pts[l][1]]
               x0 /= 3.0; y0 /= 3.0
               tri_centroids.append([x0, y0])
               tri_cells.append(tri_pts)
               tri_mats.append(c)
               n[jp[j], ip[j]] += 1
            hex_cells.append(hex_pts)
            hex_mats.append(c)
   
   pt_idx = -1 * np.ones((ny, nx), dtype = int)
   points = []
   bc_pts = []
   jp = 0
   for i in range(nx):
      for j in range(ny):
         if n[j, i] > 0:
            points.append([x[j, i], y[j, i]])
            pt_idx[j, i] = jp
            if n[j, i] < 3:
               bc_pts.append(jp)
            jp += 1
   
   x0 = 0.0
   y0 = 0.0
   for i in range(len(points)):
      x0 += points[i][0]
      y0 += points[i][1]
   x0 /= len(points)
   y0 /= len(points)
   for i in range(len(points)):
      points[i][0] -= x0
      points[i][1] -= y0
   for i in range(len(tri_centroids)):
      tri_centroids[i][0] -= x0
      tri_centroids[i][1] -= y0
   for i in range(len(hex_centroids)):
      hex_centroids[i][0] -= x0
      hex_centroids[i][1] -= y0
   
   for i in range(len(hex_cells)):
      for j in range(len(hex_cells[i])):
         hex_cells[i][j] = pt_idx[hex_cells[i][j][0], hex_cells[i][j][1]]
   
   for i in range(len(tri_cells)):
      for j in range(len(tri_cells[i])):
         tri_cells[i][j] = pt_idx[tri_cells[i][j][0], tri_cells[i][j][1]]
   
   return NodalMesh(points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

def plot_mesh(mesh):
   
   fig, ax = plt.subplots()
   
   x, y = zip(*mesh.points)
   ax.scatter(x, y, s = 8, c = "blue", alpha = 0.25)
   
   colors = [None, "deeppink", "darkviolet", "darkblue", "darkcyan", "darkgreen", "darkorange", "darkred", "black"]
   for (c, m) in zip(mesh.hex_cells, mesh.hex_mats):
      xp = []; yp = []
      for p in c:
         xp.append(mesh.points[p][0])
         yp.append(mesh.points[p][1])
      plt.fill(xp, yp, facecolor = "none", edgecolor = colors[m], linewidth = 1, alpha = 1.0)
      plt.fill(xp, yp, facecolor = colors[m], edgecolor = "none", linewidth = 1, alpha = 0.25)
   
   ax.axis("equal")
   plt.show()

def build_gmsh_mesh(mesh):
   
   gmsh.initialize()
   
   gmsh.model.add("hex-core")
   
   lc = 0.5
   
   for i, p in enumerate(mesh.points):
      gmsh.model.occ.addPoint(p[0], p[1], 0.0, lc, i+1)
   
   l = 1
   for ic, c in enumerate(mesh.hex_cells):
      lines = []
      for ip, p in enumerate(c):
         gmsh.model.occ.addLine(c[ip]+1, c[(ip+1)%len(c)]+1, l)
         lines.append(l)
         l += 1
      gmsh.model.occ.addCurveLoop(lines, ic+1)
      gmsh.model.occ.addPlaneSurface([ic+1], ic+1)
   
   gmsh.model.occ.synchronize()
   
   gmsh.model.addPhysicalGroup(1, range(1, len(mesh.hex_cells)), 1, name = "all")
   
   gmsh.model.mesh.generate(2)
   
   if '-nopopup' not in sys.argv:
      gmsh.fltk.run()
   
   gmsh.finalize()

def main():
   
   # Core geometry:
   # Fuel-assembly types:
   #    - 1 = fuel 1
   #    - 2 = fuel 2
   #    - 3 = fuel 1 + shutdown rod
   #    - 4 = fuel 2 + shutdown rod
   #    - 5 = moderator
   #    - 6 = reflector
   pc = 12.0
   small = False
   if small:
      core = [[0, 0, 0, 6, 6, 6, 6, 0, 0, 0], \
               [0, 6, 6, 2, 2, 2, 6, 6, 0], \
              [0, 6, 2, 4, 2, 2, 4, 2, 6, 0], \
                [6, 2, 2, 1, 3, 1, 2, 2, 6], \
               [6, 2, 2, 3, 1, 1, 3, 2, 2, 6], \
                 [6, 4, 1, 1, 5, 1, 1, 4, 6], \
                [6, 2, 2, 3, 1, 1, 3, 2, 2, 6], \
                  [6, 2, 2, 1, 3, 1, 2, 2, 6], \
                 [0, 6, 2, 4, 2, 2, 4, 2, 6, 0], \
                   [0, 6, 6, 2, 2, 2, 6, 6, 0], \
                  [0, 0, 0, 6, 6, 6, 6, 0, 0, 0]]
   else:
      core = [[0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0], \
                [0, 0, 6, 6, 2, 2, 2, 2, 2, 6, 6, 0, 0], \
               [0, 0, 6, 2, 2, 4, 2, 2, 4, 2, 2, 6, 0, 0], \
                 [0, 6, 2, 4, 2, 2, 3, 2, 2, 4, 2, 6, 0], \
                [0, 6, 2, 2, 2, 3, 1, 1, 3, 2, 2, 2, 6, 0], \
                  [6, 2, 2, 3, 1, 1, 1, 1, 1, 3, 2, 2, 6], \
                 [6, 2, 4, 2, 1, 1, 5, 5, 1, 1, 2, 4, 2, 6], \
                   [6, 2, 2, 3, 1, 5, 5, 5, 1, 3, 2, 2, 6], \
                  [6, 2, 4, 2, 1, 1, 5, 5, 1, 1, 2, 4, 2, 6], \
                    [6, 2, 2, 3, 1, 1, 1, 1, 1, 3, 2, 2, 6], \
                   [0, 6, 2, 2, 2, 3, 1, 1, 3, 2, 2, 2, 6, 0], \
                     [0, 6, 2, 4, 2, 2, 3, 2, 2, 4, 2, 6, 0], \
                    [0, 0, 6, 2, 2, 4, 2, 2, 4, 2, 2, 6, 0, 0], \
                      [0, 0, 6, 6, 2, 2, 2, 2, 2, 6, 6, 0, 0], \
                     [0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0]]
   
   # Fuel-assembly geometry:
   # Pin types:
   #    - 1 = fuel
   #    - 2 = graphite
   #    - 3 = heat pipe
   #    - 4 = shutdown rod
   pf = 1.8
   fas = [None] * 6
   fas[0] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                [0, 2, 1, 1, 3, 1, 1, 2, 0], \
               [0, 0, 1, 3, 1, 1, 3, 1, 0, 0], \
                 [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                [0, 0, 1, 3, 1, 1, 3, 1, 0, 0], \
                  [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                 [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[1] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                [0, 2, 1, 1, 3, 1, 1, 2, 0], \
               [0, 0, 1, 3, 1, 1, 3, 1, 0, 0], \
                 [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                [0, 0, 1, 3, 1, 1, 3, 1, 0, 0], \
                  [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                 [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[2] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                [0, 2, 1, 1, 3, 1, 1, 2, 0], \
               [0, 0, 1, 3, 4, 4, 3, 1, 0, 0], \
                 [0, 2, 1, 4, 4, 4, 1, 2, 0], \
                [0, 0, 1, 3, 4, 4, 3, 1, 0, 0], \
                  [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                 [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[3] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                [0, 2, 1, 1, 3, 1, 1, 2, 0], \
               [0, 0, 1, 3, 4, 4, 3, 1, 0, 0], \
                 [0, 2, 1, 4, 4, 4, 1, 2, 0], \
                [0, 0, 1, 3, 4, 4, 3, 1, 0, 0], \
                  [0, 2, 1, 1, 3, 1, 1, 2, 0], \
                 [0, 0, 0, 2, 1, 1, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[4] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 2, 2, 2, 0, 0, 0], \
                [0, 2, 2, 2, 3, 2, 2, 2, 0], \
               [0, 0, 2, 3, 2, 2, 3, 2, 0, 0], \
                 [0, 2, 2, 2, 3, 2, 2, 2, 0], \
                [0, 0, 2, 3, 2, 2, 3, 2, 0, 0], \
                  [0, 2, 2, 2, 3, 2, 2, 2, 0], \
                 [0, 0, 0, 2, 2, 2, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[5] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 2, 0, 0, 0, 0], \
              [0, 0, 0, 2, 2, 2, 2, 0, 0, 0], \
                [0, 2, 2, 2, 3, 2, 2, 2, 0], \
               [0, 0, 2, 3, 2, 2, 3, 2, 0, 0], \
                 [0, 2, 2, 2, 3, 2, 2, 2, 0], \
                [0, 0, 2, 3, 2, 2, 3, 2, 0, 0], \
                  [0, 2, 2, 2, 3, 2, 2, 2, 0], \
                 [0, 0, 0, 2, 2, 2, 2, 0, 0, 0], \
                   [0, 0, 0, 0, 2, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   
   core_mesh = build_x_hex_mesh(pc, core)
   # plot_mesh(core_mesh)
   
   build_gmsh_mesh(core_mesh)

if __name__ == "__main__": main()
