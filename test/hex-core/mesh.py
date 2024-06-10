from typing import NamedTuple
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

class Mesh(NamedTuple):
   
   points: list
   hex_centroids: list
   tri_centroids: list
   hex_cells: list
   tri_cells: list
   hex_mats: list
   tri_mats: list
   bc_pts: list

class Grid(NamedTuple):
   
   points: list
   centroids: list
   cells: list
   faces: dict
   mats: list
   bc_pts: list

def build_y_hex_mesh(p, layout):
   
   nx = 2*len(layout[0]) + 1
   ny = int(1.5*len(layout)) + 2
   
   y, x = np.meshgrid(np.arange(float(ny), 0.0, -1.0), np.arange(float(nx)))
   y[1::2, :] += 0.5
   dx = 0.5 * p
   dy = (2.0/np.sqrt(3)) * 0.5 * p
   y = dy * y
   x = dx * x
   
   hex_centroids = []
   tri_centroids = []
   hex_cells = []
   tri_cells = []
   hex_mats = []
   tri_mats = []
   tri_mats = []
   n = np.zeros((nx, ny), dtype = int)
   for (j, cj) in enumerate(layout):
      for (i, c) in enumerate(cj):
         if c != 0:
            
            i0 = 2*i + j%2 + 1
            j0 = int(1.5*j) + 1
            hex_centroids.append([x[i0, j0], y[i0, j0]])
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
               x0 = 0.0; y0 = 0.0
               for l in range(3):
                  x0 += x[tri_pts[l][0], tri_pts[l][1]]
                  y0 += y[tri_pts[l][0], tri_pts[l][1]]
               x0 /= 3.0; y0 /= 3.0
               tri_centroids.append([x0, y0])
               tri_cells.append(tri_pts)
               tri_mats.append(c)
               n[ip[k], jp[k]] += 1
            hex_cells.append(hex_pts)
            hex_mats.append(c)
   
   pt_idx = -1 * np.ones((nx, ny), dtype = int)
   points = []
   bc_pts = []
   ip = 0
   for j in range(ny):
      for i in range(nx):
         if n[i, j] > 0:
            points.append([x[i, j], y[i, j]])
            pt_idx[i, j] = ip
            if n[i, j] < 3:
               bc_pts.append(ip)
            ip += 1
   
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
   
   return Mesh(points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

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
   
   return Mesh(points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

def build_x_hex_grid(p, layout, h, dh):
   
   ny = 2*len(layout[0]) + 1
   nx = int(1.5*len(layout)) + 2
   
   x, y = np.meshgrid(np.arange(float(nx), 0.0, -1.0), np.arange(float(ny)))
   x[1::2, :] += 0.5
   dy = 0.5 * p
   dx = (2.0/np.sqrt(3)) * 0.5 * p
   x = dx * x
   y = dy * y
   
   centroids = []
   hex_cells = []
   hex_mats = []
   n = np.zeros((ny, nx), dtype = int)
   for (i, ci) in enumerate(layout):
      for (j, c) in enumerate(ci):
         if c != 0:
            
            j0 = 2*j + i%2 + 1
            i0 = int(1.5*i) + 1
            centroids.append([x[j0, i0], y[j0, i0]])
            n[j0, i0] += 3
            
            jp = [j0-1, j0, j0+1, j0+1, j0, j0-1]
            ip = [i0+i%2-1, i0-1, i0+i%2-1, i0+i%2, i0+1, i0+i%2]
            hex_pts = []
            for j in range(6):
               hex_pts.append((jp[j], ip[j]))
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
   for i in range(len(centroids)):
      centroids[i][0] -= x0
      centroids[i][1] -= y0
   
   for i in range(len(hex_cells)):
      for j in range(len(hex_cells[i])):
         hex_cells[i][j] = pt_idx[hex_cells[i][j][0], hex_cells[i][j][1]]
   
   hex_faces = dict()
   for i in range(len(hex_cells)):
      m = len(hex_cells[i])
      for j in range(m):
         p1 = hex_cells[i][j]
         p2 = hex_cells[i][(j+1)%m]
         if not ((p1, p2) in hex_faces) and not ((p2, p1) in hex_faces):
            face = (p1, p2)
            p1 = points[p1]
            p2 = points[p2]
            d = sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2)
            dx = (p2[0]-p1[0]) / d
            dy = (p2[1]-p1[1]) / d
            d0 = 0.5*(d-(h-1)*dh)
            pts = []
            for k in range(h):
               s = d0 + k*dh
               x = p1[0] + s*dx
               y = p1[1] + s*dy
               pts.append(len(points))
               if (face[0] in bc_pts) and (face[1] in bc_pts):
                  bc_pts.append(len(points))
               points.append([x, y])
            hex_faces[face] = [face[0]] + pts + [face[1]]
         else:
            hex_faces[(p1, p2)] = hex_faces[(p2, p1)][::-1]
   
   return Grid(points, centroids, hex_cells, hex_faces, hex_mats, bc_pts)

def inside(point, polygon):
   
   n = len(polygon)
   for i in range(n):
      d = (polygon[(i+1)%n][0]-polygon[i][0]) * (point[1]-polygon[i][1]) - \
             (point[0]-polygon[i][0]) * (polygon[(i+1)%n][1]-polygon[i][1])
      if i == 0:
         d0 = d
      else:
         if d*d0 < 0.0: return False
   return True

def distance(p1, p2):
   
   return sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2)

def build_x_hex_nested_mesh(grid, meshes):
   
   points = []
   cells = []
   mats = []
   
   for p in grid.points:
      points.append(p)
   
   for grid_centroid, grid_cell, grid_mat in zip(grid.centroids, grid.cells, grid.mats):
      
      grid_points = [grid.points[i] for i in grid_cell]
      
      p0 = len(points)
      c0 = len(cells)
      mesh = meshes[grid_mat-1]
      
      for p in mesh.points:
         points.append([grid_centroid[0]+p[0], grid_centroid[1]+p[1]])
      
      for centroid, cell, mat in zip(mesh.tri_centroids, mesh.tri_cells, mesh.tri_mats):
         add = True
         for p in cell:
            p = mesh.points[p]
            if not inside([grid_centroid[0]+p[0], grid_centroid[1]+p[1]], grid_points):
               add = False
               break
         if add:
            pts = []
            for p in cell:
               pts.append(p0+p)
            cells.append(pts)
            mats.append(mat)
      
      n = len(grid_cell)
      for i in range(n):
         face = grid.faces[(grid_cell[i], grid_cell[(i+1)%n])]
         m = len(face)
         for j in range(m-1):
            pts = [face[j]]
            for k in range(2):
               dmin = 1.0e16
               p = None
               for p2 in range(len(mesh.points)):
                  if inside([grid_centroid[0]+mesh.points[p2][0], grid_centroid[1]+mesh.points[p2][1]], grid_points):
                     d = distance(grid.points[face[j+k]], (grid_centroid[0]+mesh.points[p2][0], grid_centroid[1]+mesh.points[p2][1]))
                     if d < dmin:
                        dmin = d
                        p = p2
               pts.append(p0+p)
            pts.append(face[j+1])
            # for cell, mat in zip(cells[c0:], mats[c0:]):
            #    if (pts[1] in cell) and (pts[2] in cell):
            #        mats.append(mat)
            mats.append(3)
            pts.reverse()
            cells.append(pts)
   
   return Mesh(points, None, None, None, cells, None, mats, grid.bc_pts)

def clean_up(mesh):
   
   n = [0] * len(mesh.points)
   for cell in mesh.tri_cells:
      for p in cell:
         n[p] += 1
   
   i = 0
   ids = [None] * len(mesh.points)
   for p in range(len(mesh.points)):
      if n[p] > 0:
         ids[p] = i
         i += 1
   
   points = [p for (p, i) in zip(mesh.points, ids) if i != None]
   bc_pts = [ids[p] for p in mesh.bc_pts]
   
   cells = []
   for cell in mesh.tri_cells:
      pts = []
      for p in cell:
         pts.append(ids[p])
      cells.append(pts)
   
   return Mesh(points, None, None, None, cells, None, mesh.tri_mats, bc_pts)

def write_mesh(filename, mesh, cells, mats, nzb, nz, nzt, ref_mat):
   
   with open(filename, "w") as f:
      
      f.write("points %d\n" % len(mesh.points))
      for p in mesh.points:
         f.write("%.3f %.3f\n" % (p[0], p[1]))
      f.write("\n")
      
      np = 0
      for c in cells:
         np += len(c)
      f.write("cells %d %d\n" % (len(cells), np))
      for c in cells:
         for i, p in enumerate(c):
            if i > 0: f.write(" ")
            f.write("%d" % p)
         f.write("\n")
      f.write("\n")
      
      nztot = nzb + nz + nzt
      if nztot > 1:
         f.write("dz -%d\n" % nztot)
         f.write("10.0\n")
         f.write("\n")
      
      f.write("boundary %d\n" % len(mesh.bc_pts))
      for i in mesh.bc_pts:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("materials %d\n" % (nztot*len(cells)))
      for k in range(nztot):
         f.write("\n")
         for i in range(len(mats)):
            if i > 0: f.write(" ")
            if k >= nzb and k < nzb+nz:
               f.write("%d" % mats[i])
            else:
               f.write("%d" % ref_mat)
         f.write("\n")

def plot_mesh(mesh):
   
   fig, ax = plt.subplots()
   
   x, y = zip(*mesh.points)
   ax.scatter(x, y, s = 8, c = "blue", alpha = 0.25)
   
   colors = [None, "deeppink", "darkviolet", "darkblue", "darkcyan", "darkgreen", "darkorange", "darkred", "black"]
   for (c, m) in zip(mesh.tri_cells, mesh.tri_mats):
      xp = []; yp = []
      for p in c:
         xp.append(mesh.points[p][0])
         yp.append(mesh.points[p][1])
      plt.fill(xp, yp, facecolor = "none", edgecolor = colors[m], linewidth = 1, alpha = 1.0)
      plt.fill(xp, yp, facecolor = colors[m], edgecolor = "none", linewidth = 1, alpha = 0.25)
   
   ax.axis("equal")
   plt.show()

def plot_grid(grid):
   
   fig, ax = plt.subplots()
   
   x, y = zip(*grid.points)
   ax.scatter(x, y, s = 8, c = "blue", alpha = 0.25)
   
   colors = [None, "deeppink", "darkviolet", "darkblue", "darkcyan", "darkgreen", "darkorange", "darkred", "black"]
   for (c, m) in zip(grid.cells, grid.mats):
      xp = []; yp = []
      for p in c:
         xp.append(grid.points[p][0])
         yp.append(grid.points[p][1])
      plt.fill(xp, yp, facecolor = "none", edgecolor = colors[m], linewidth = 1, alpha = 1.0)
      plt.fill(xp, yp, facecolor = colors[m], edgecolor = "none", linewidth = 1, alpha = 0.25)
   
   ax.axis("equal")
   plt.show()

def main():
   
   build_nodal_mesh = False
   build_pin_mesh = False
   build_multiscale_mesh = True
   
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
   #    - 1 = fuel 1
   #    - 2 = fuel 2
   #    - 3 = fuel 1 fill
   #    - 4 = fuel 2 fill
   #    - 5 = heat pipe
   #    - 6 = shutdown rod
   #    - 7 = moderator
   #    - 8 = reflector
   pf = 1.8
   fas = [None] * 6
   fas[0] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 3, 0, 0, 0, 0], \
              [0, 0, 0, 3, 1, 1, 3, 0, 0, 0], \
                [0, 3, 1, 1, 5, 1, 1, 3, 0], \
               [0, 0, 1, 5, 1, 1, 5, 1, 0, 0], \
                 [0, 3, 1, 1, 5, 1, 1, 3, 0], \
                [0, 0, 1, 5, 1, 1, 5, 1, 0, 0], \
                  [0, 3, 1, 1, 5, 1, 1, 3, 0], \
                 [0, 0, 0, 3, 1, 1, 3, 0, 0, 0], \
                   [0, 0, 0, 0, 3, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[1] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 4, 0, 0, 0, 0], \
              [0, 0, 0, 4, 2, 2, 4, 0, 0, 0], \
                [0, 4, 2, 2, 5, 2, 2, 4, 0], \
               [0, 0, 2, 5, 2, 2, 5, 2, 0, 0], \
                 [0, 4, 2, 2, 5, 2, 2, 4, 0], \
                [0, 0, 2, 5, 2, 2, 5, 2, 0, 0], \
                  [0, 4, 2, 2, 5, 2, 2, 4, 0], \
                 [0, 0, 0, 4, 2, 2, 4, 0, 0, 0], \
                   [0, 0, 0, 0, 4, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[2] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 3, 0, 0, 0, 0], \
              [0, 0, 0, 3, 1, 1, 3, 0, 0, 0], \
                [0, 3, 1, 1, 5, 1, 1, 3, 0], \
               [0, 0, 1, 5, 6, 6, 5, 1, 0, 0], \
                 [0, 3, 1, 6, 6, 6, 1, 3, 0], \
                [0, 0, 1, 5, 6, 6, 5, 1, 0, 0], \
                  [0, 3, 1, 1, 5, 1, 1, 3, 0], \
                 [0, 0, 0, 3, 1, 1, 3, 0, 0, 0], \
                   [0, 0, 0, 0, 3, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[3] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 4, 0, 0, 0, 0], \
              [0, 0, 0, 4, 2, 2, 4, 0, 0, 0], \
                [0, 4, 2, 2, 5, 2, 2, 4, 0], \
               [0, 0, 2, 5, 6, 6, 5, 2, 0, 0], \
                 [0, 4, 2, 6, 6, 6, 2, 4, 0], \
                [0, 0, 2, 5, 6, 6, 5, 2, 0, 0], \
                  [0, 4, 2, 2, 5, 2, 2, 4, 0], \
                 [0, 0, 0, 4, 2, 2, 4, 0, 0, 0], \
                   [0, 0, 0, 0, 4, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[4] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 7, 0, 0, 0, 0], \
              [0, 0, 0, 7, 7, 7, 7, 0, 0, 0], \
                [0, 7, 7, 7, 7, 7, 7, 7, 0], \
               [0, 0, 7, 7, 7, 7, 7, 7, 0, 0], \
                 [0, 7, 7, 7, 7, 7, 7, 7, 0], \
                [0, 0, 7, 7, 7, 7, 7, 7, 0, 0], \
                  [0, 7, 7, 7, 7, 7, 7, 7, 0], \
                 [0, 0, 0, 7, 7, 7, 7, 0, 0, 0], \
                   [0, 0, 0, 0, 7, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   fas[5] = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0], \
               [0, 0, 0, 0, 8, 0, 0, 0, 0], \
              [0, 0, 0, 8, 8, 8, 8, 0, 0, 0], \
                [0, 8, 8, 8, 8, 8, 8, 8, 0], \
               [0, 0, 8, 8, 8, 8, 8, 8, 0, 0], \
                 [0, 8, 8, 8, 8, 8, 8, 8, 0], \
                [0, 0, 8, 8, 8, 8, 8, 8, 0, 0], \
                  [0, 8, 8, 8, 8, 8, 8, 8, 0], \
                 [0, 0, 0, 8, 8, 8, 8, 0, 0, 0], \
                   [0, 0, 0, 0, 8, 0, 0, 0, 0], \
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
   
   if build_nodal_mesh:
      
      core_mesh = build_x_hex_mesh(pc, core)
      # plot_mesh(core_mesh)
      
      write_mesh("hex-cells-diffusion-3d/mesh.pmp", core_mesh, core_mesh.hex_cells, core_mesh.hex_mats, 2, 12, 2, 6)
      write_mesh("tri-cells-diffusion-3d/mesh.pmp", core_mesh, core_mesh.tri_cells, core_mesh.tri_mats, 2, 12, 2, 6)
   
   if build_pin_mesh:
      
      core_grid = build_x_hex_grid(pc, core, 5, (2.0/np.sqrt(3))*0.5*pf)
      # plot_grid(core_grid)
      
      fa_meshes = []
      for fa in fas:
         fa_meshes.append(build_x_hex_mesh(pf, fa))
         # plot_mesh(fa_meshes[-1])
      
      pin_mesh = build_x_hex_nested_mesh(core_grid, fa_meshes)
      pin_mesh = clean_up(pin_mesh)
      # plot_mesh(pin_mesh)
      
      write_mesh("pin-level-diffusion-3d/mesh.pmp", pin_mesh, pin_mesh.tri_cells, pin_mesh.tri_mats, 2, 12, 2, 8)
   
   if build_multiscale_mesh:
      
      core_mesh = build_x_hex_mesh(pc, core)
      # plot_mesh(core_mesh)
      
      write_mesh("multiscale-conduction-3d/mesh.pmp", core_mesh, core_mesh.hex_cells, core_mesh.hex_mats, 2, 12, 2, 6)
      
      for i, fa in enumerate(fas, 1):
         
         fa_mesh = build_x_hex_mesh(pf, fa)
         # plot_mesh(fa_mesh)
         
         fa_grid = build_x_hex_grid(pc, [[1]], 5, (2.0/np.sqrt(3))*0.5*pf)
         
         pin_mesh = build_x_hex_nested_mesh(fa_grid, [fa_mesh])
         pin_mesh = clean_up(pin_mesh)
         # plot_mesh(pin_mesh)
         
         write_mesh("multiscale-conduction-3d/submesh-%d.pmp" % i, pin_mesh, pin_mesh.tri_cells, pin_mesh.tri_mats, 0, 1, 0, 8)

if __name__ == "__main__": main()
