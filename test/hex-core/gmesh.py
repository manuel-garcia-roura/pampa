from typing import NamedTuple
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import sys

class Mesh(NamedTuple):
   
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
   
   return Mesh(points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts)

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

def build_gmsh_mesh(core_mesh, fa_meshes, d, r):
   
   gmsh.initialize()
   
   gmsh.model.add("hex-core")
   
   lc1 = 0.5
   lc2 = 0.5
   
   for p in core_mesh.points:
      gmsh.model.occ.addPoint(p[0], p[1], 0.0, lc1)
   
   hxs = []; fas = []; fps = []; mps = []; hps = []; srs = []
   for c0, c, fa in zip(core_mesh.hex_centroids, core_mesh.hex_cells, core_mesh.hex_mats):
      
      sides = []
      for p in range(len(c)):
         
         il = gmsh.model.occ.addLine(c[p]+1, c[(p+1)%len(c)]+1)
         sides.append(il)
      
      holes = []
      for c00, m in zip(fa_meshes[fa-1].hex_centroids, fa_meshes[fa-1].hex_mats):
         
         gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]+0.5*d[m-1], 0.0, lc2)
         gmsh.model.occ.addPoint(c0[0]+c00[0]-0.5*d[m-1], c0[1]+c00[1], 0.0, lc2)
         gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]-0.5*d[m-1], 0.0, lc2)
         gmsh.model.occ.addPoint(c0[0]+c00[0]+0.5*d[m-1], c0[1]+c00[1], 0.0, lc2)
         
         il = gmsh.model.occ.addCircle(c0[0]+c00[0], c0[1]+c00[1], 0.0, 0.5*d[m-1], angle1 = 0.0, angle2 = 2*np.pi)
         icl = gmsh.model.occ.addCurveLoop([il])
         ips = gmsh.model.occ.addPlaneSurface([icl])
         
         holes.append(ips)
         if m == 1:
            fps.append(ips)
         elif m == 2:
            mps.append(ips)
         elif m == 3:
            hps.append(ips)
         elif m == 4:
            srs.append(ips)
      
      icl = gmsh.model.occ.addCurveLoop(sides)
      ips = gmsh.model.occ.addPlaneSurface([icl] + [-x for x in holes])
      fas.append(ips)
      
      icl = gmsh.model.occ.addCurveLoop(sides)
      ips = gmsh.model.occ.addPlaneSurface([icl])
      hxs.append(ips)
   
   gmsh.model.occ.addPoint(0.0, r, 0.0, lc1)
   gmsh.model.occ.addPoint(-r, 0.0, 0.0, lc1)
   gmsh.model.occ.addPoint(0.0, -r, 0.0, lc1)
   gmsh.model.occ.addPoint(r, 0.0, 0.0, lc1)
   
   boundary = gmsh.model.occ.addCircle(0.0, 0.0, 0.0, r, angle1 = 0.0, angle2 = 2*np.pi)
   icl = gmsh.model.occ.addCurveLoop([boundary])
   reflector = gmsh.model.occ.addPlaneSurface([icl])
   gmsh.model.occ.cut([(2, reflector)], [(2, x) for x in hxs])
   
   gmsh.model.occ.synchronize()
   
   gmsh.model.addPhysicalGroup(2, fas + [reflector], name = "graphite")
   gmsh.model.addPhysicalGroup(2, fps, name = "fuel")
   gmsh.model.addPhysicalGroup(2, mps, name = "moderator")
   gmsh.model.addPhysicalGroup(2, hps, name = "heat-pipe")
   gmsh.model.addPhysicalGroup(2, srs, name = "shutdown-rod")
   gmsh.model.addPhysicalGroup(1, [boundary], name = "boundary")
   
   gmsh.model.mesh.generate(2)
   
   gmsh.write("mesh.msh")
   
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
   small = True
   if small:
      core = [[0, 0, 6, 6, 0, 0], \
               [6, 4, 3, 4, 6], \
              [6, 3, 1, 2, 3, 6], \
                [4, 2, 5, 1, 4], \
               [6, 3, 1, 2, 3, 6], \
                 [6, 4, 3, 4, 6], \
                [0, 0, 6, 6, 0, 0]]
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
   r = 0.5 * pc * len(core)
   
   # Fuel-assembly geometry:
   # Pin types:
   #    - 1 = fuel
   #    - 2 = moderator
   #    - 3 = heat pipe
   #    - 4 = shutdown rod
   pf = 1.8
   d = [1.0, 1.5, 1.2, 4.0]
   fas = [None] * 6
   fas[0] = [[0, 0, 1, 1, 0, 0], \
               [1, 1, 3, 1, 1], \
              [1, 3, 1, 1, 3, 1], \
                [1, 1, 3, 1, 1], \
               [1, 3, 1, 1, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[1] = [[0, 0, 1, 1, 0, 0], \
               [1, 1, 3, 1, 1], \
              [1, 3, 1, 1, 3, 1], \
                [1, 1, 3, 1, 1], \
               [1, 3, 1, 1, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[2] = [[0, 0, 1, 1, 0, 0], \
               [1, 1, 3, 1, 1], \
              [1, 3, 0, 0, 3, 1], \
                [1, 0, 4, 0, 1], \
               [1, 3, 0, 0, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[3] = [[0, 0, 1, 1, 0, 0], \
               [1, 1, 3, 1, 1], \
              [1, 3, 0, 0, 3, 1], \
                [1, 0, 4, 0, 1], \
               [1, 3, 0, 0, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[4] = [[0, 0, 2, 2, 0, 0], \
               [2, 2, 2, 2, 2], \
              [2, 2, 2, 2, 2, 2], \
                [2, 2, 2, 2, 2], \
               [2, 2, 2, 2, 2, 2], \
                 [2, 2, 2, 2, 2], \
                [0, 0, 2, 2, 0, 0]]
   fas[5] = [[0, 0, 0, 0, 0, 0], \
               [0, 0, 3, 0, 0], \
              [0, 3, 0, 0, 3, 0], \
                [0, 0, 3, 0, 0], \
               [0, 3, 0, 0, 3, 0], \
                 [0, 0, 3, 0, 0], \
                [0, 0, 0, 0, 0, 0]]
   
   core_mesh = build_x_hex_mesh(pc, core)
   fa_meshes = [build_x_hex_mesh(pf, fa) for fa in fas]
   
   build_gmsh_mesh(core_mesh, fa_meshes, d, r)

if __name__ == "__main__": main()
