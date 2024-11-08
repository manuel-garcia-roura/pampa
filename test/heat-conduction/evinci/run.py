import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import gmsh
import subprocess
import os
import shutil

class Mesh:
   
   def __init__(self, points, cells, centroids, mats, boundaries, nodes):
      
      self.points = points
      self.cells = cells
      self.centroids = centroids
      self.mats = mats
      self.boundaries = boundaries
      self.nodes = nodes

def parse_vtk_file(file_path):
   
   scalar_fields = {}
   with open(file_path, "r") as file:
      
      lines = file.readlines()
      
      is_reading_scalars = False
      current_scalar_name = None
      current_scalar_data = []
      
      for line in lines:
         
         if line.startswith("SCALARS"):
            
            if current_scalar_name:
               scalar_fields[current_scalar_name] = np.array(current_scalar_data)
            
            parts = line.split()
            current_scalar_name = parts[1]
            current_scalar_data = []
            is_reading_scalars = True
         
         elif is_reading_scalars and line.startswith("LOOKUP_TABLE"):
            continue
         
         elif is_reading_scalars and not line.startswith(("SCALARS", "LOOKUP_TABLE")):
            current_scalar_data.extend(map(float, line.split()))
      
      if current_scalar_name:
         scalar_fields[current_scalar_name] = np.array(current_scalar_data)
   
   return scalar_fields

def build_x_hex_mesh(pitch, layout):
   
   nx = int(1.5*len(layout)) + 2
   ny = 2*len(layout[0]) + 1
   
   x, y = np.meshgrid(np.arange(float(nx), 0.0, -1.0), np.arange(float(ny)))
   x[1::2, :] += 0.5
   dy = 0.5 * pitch
   dx = (2.0/np.sqrt(3.0)) * 0.5 * pitch
   
   x = dx * x
   y = dy * y
   
   cells = []; centroids = []; mats = []
   n = np.zeros((ny, nx), dtype = int)
   for (i, row) in enumerate(layout):
      for (j, c) in enumerate(row):
         if c != 0:
            
            j0 = 2*j + i%2 + 1
            i0 = int(1.5*i) + 1
            
            jp = [j0-1, j0, j0+1, j0+1, j0, j0-1]
            ip = [i0+i%2-1, i0-1, i0+i%2-1, i0+i%2, i0+1, i0+i%2]
            
            cells.append(list(zip(jp, ip)))
            centroids.append([x[j0, i0], y[j0, i0]])
            mats.append(c)
            
            n[j0, i0] += 3
            for j in range(6):
               n[jp[j], ip[j]] += 1
   
   points = []; boundaries = []
   idx = -1 * np.ones((ny, nx), dtype = int)
   jp = 0
   for i in range(nx):
      for j in range(ny):
         if n[j, i] > 0:
            points.append([x[j, i], y[j, i]])
            idx[j, i] = jp
            if n[j, i] < 3:
               boundaries.append(jp)
            jp += 1
   
   x0, y0 = map(sum, zip(*points))
   x0 /= len(points); y0 /= len(points)
   points = [(x-x0, y-y0) for (x, y) in points]
   centroids = [(x-x0, y-y0) for (x, y) in centroids]
   
   for i in range(len(cells)):
      for j in range(6):
         cells[i][j] = idx[cells[i][j][0], cells[i][j][1]]
   
   return Mesh(points, cells, centroids, mats, boundaries, None)

def build_gmsh_mesh(core_mesh, fa_meshes, d, r, lc1, lc2, lc3, pitch):
   
   write_vtk = False
   run_fltk = False
   
   gmsh.initialize()
   
   gmsh.model.add("evinci")
   
   pin_pts = [False] * len(core_mesh.points)
   for c, fa in zip(core_mesh.cells, core_mesh.mats):
      if not fa_meshes[fa-1] is None:
         for p in c:
            pin_pts[p] = True
   
   for ip, p in enumerate(core_mesh.points):
      if ip in core_mesh.boundaries or not pin_pts[ip]:
         lc = lc3
      else:
         lc = lc2
      gmsh.model.occ.addPoint(p[0], p[1], 0.0, lc)
   
   fas = []; hxs = []; pins = [[] for _ in range(4)]; regions = [None]
   for c0, c, fa in zip(core_mesh.centroids, core_mesh.cells, core_mesh.mats):
      
      sides = []
      for p in range(len(c)):
         
         p1 = c[p] + 1
         p2 = c[(p+1)%len(c)] + 1
         il = gmsh.model.occ.addLine(p1, p2)
         sides.append(il)
      
      holes = []
      if not fa_meshes[fa-1] is None:
         
         for c00, m in zip(fa_meshes[fa-1].centroids, fa_meshes[fa-1].mats):
            
            r1 = 0.5 * d[m-1]
            r2 = 0.25 * d[m-1]
            r3 = 0.75 * d[m-1]
            
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]+r1, 0.0, lc1)
            gmsh.model.occ.addPoint(c0[0]+c00[0]-r1, c0[1]+c00[1], 0.0, lc1)
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]-r1, 0.0, lc1)
            gmsh.model.occ.addPoint(c0[0]+c00[0]+r1, c0[1]+c00[1], 0.0, lc1)
            il1 = gmsh.model.occ.addCircle(c0[0]+c00[0], c0[1]+c00[1], 0.0, r1, angle1 = 0.0, angle2 = 2*np.pi)
            icl1 = gmsh.model.occ.addCurveLoop([il1])
            
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]+r2, 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]-r2, c0[1]+c00[1], 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]-r2, 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]+r2, c0[1]+c00[1], 0.0, lc2)
            il2 = gmsh.model.occ.addCircle(c0[0]+c00[0], c0[1]+c00[1], 0.0, r2, angle1 = 0.0, angle2 = 2*np.pi)
            icl2 = gmsh.model.occ.addCurveLoop([il2])
            
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]+r3, 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]-r3, c0[1]+c00[1], 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]-r3, 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]+r3, c0[1]+c00[1], 0.0, lc2)
            il3 = gmsh.model.occ.addCircle(c0[0]+c00[0], c0[1]+c00[1], 0.0, r3, angle1 = 0.0, angle2 = 2*np.pi)
            icl3 = gmsh.model.occ.addCurveLoop([il3])
            
            ips2 = gmsh.model.occ.addPlaneSurface([icl2])
            holes.append(ips2)
            pins[m-1].append(ips2)
            
            ips1 = gmsh.model.occ.addPlaneSurface([icl1, -icl2])
            holes.append(ips1)
            pins[m-1].append(ips1)
            
            ips3 = gmsh.model.occ.addPlaneSurface([icl3, -icl1])
            holes.append(ips3)
            pins[3].append(ips3)
      
      else:
         
         pts = []
         for ci in c:
            p0 = core_mesh.points[ci]
            x = 0.5 * (p0[0]+c0[0])
            y = 0.5 * (p0[1]+c0[1])
            ip = gmsh.model.occ.addPoint(x, y, 0.0, lc3)
            pts.append(ip)
         
         lines = []
         for p1, p2 in zip(pts, pts[1:] + [pts[0]]):
            il = gmsh.model.occ.addLine(p1, p2)
            lines.append(il)
         
         icl = gmsh.model.occ.addCurveLoop(lines)
         ips = gmsh.model.occ.addPlaneSurface([icl])
         
         holes.append(ips)
         pins[3].append(ips)
      
      icl = gmsh.model.occ.addCurveLoop(sides)
      ips = gmsh.model.occ.addPlaneSurface([icl] + [-x for x in holes])
      fas.append(ips)
      regions.append([ips] + [x for x in holes])
      
      icl = gmsh.model.occ.addCurveLoop(sides)
      ips = gmsh.model.occ.addPlaneSurface([icl])
      hxs.append(ips)
   
   reflector = []
   if not r is None:
      
      gmsh.model.occ.addPoint(0.0, r, 0.0, lc3)
      gmsh.model.occ.addPoint(-r, 0.0, 0.0, lc3)
      gmsh.model.occ.addPoint(0.0, -r, 0.0, lc3)
      gmsh.model.occ.addPoint(r, 0.0, 0.0, lc3)
      
      boundary = gmsh.model.occ.addCircle(0.0, 0.0, 0.0, r, angle1 = 0.0, angle2 = 2*np.pi)
      curve = gmsh.model.occ.addCurveLoop([boundary])
      reflector = gmsh.model.occ.addPlaneSurface([curve])
      gmsh.model.occ.cut([(2, reflector)], [(2, x) for x in hxs])
      
      reflector = [reflector] + pins[3]
   
   gmsh.model.occ.synchronize()
   
   materials = []
   materials.append(gmsh.model.addPhysicalGroup(2, fas + pins[3] + reflector, name = "graphite"))
   materials.append(gmsh.model.addPhysicalGroup(2, pins[0], name = "fuel"))
   if pins[1]:
      materials.append(gmsh.model.addPhysicalGroup(2, pins[1], name = "shutdown-rod"))
   materials.append(gmsh.model.addPhysicalGroup(2, pins[2], name = "heat-pipe"))
   if not r is None:
      boundary = gmsh.model.addPhysicalGroup(1, [boundary], name = "boundary")
      regions[0] = gmsh.model.addPhysicalGroup(2, reflector, name = "reflector")
   for ir, reg in enumerate(regions[1:], 1):
      regions[ir] = gmsh.model.addPhysicalGroup(2, reg, name = "node-" + str(ir))
   
   gmsh.model.mesh.generate(2)
   
   gmsh.model.mesh.removeDuplicateNodes()
   gmsh.model.occ.removeAllDuplicates()
   
   if write_vtk:
      gmsh.write("mesh.vtk")
   
   if run_fltk:
      gmsh.fltk.run()
   
   tags, pts, _ = gmsh.model.mesh.getNodes()
   points = [(pts[i], pts[i+1]) for i in range(0, len(pts), 3)]
   points = [p for _, p in sorted(zip(tags, points))]
   
   cells = None
   types, _, tags = gmsh.model.mesh.getElements()
   for t, pts in zip(types, tags):
      n = gmsh.model.mesh.getElementProperties(t)[3]
      if n == 3:
         cells = [tuple(pts[i:i+n]) for i in range(0, len(pts), n)]
         break
   cells = [[int(i-1) for i in c] for c in cells]
   
   n = 0
   for mat in materials:
      entities = gmsh.model.getEntitiesForPhysicalGroup(2, mat)
      for entity_tag in entities:
         element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
         for tags in element_tags:
            if len(tags) > 0:
               n = max(n, max(tags))
   
   mats = [None] * (int(n) + 1)
   for mat in materials:
      entities = gmsh.model.getEntitiesForPhysicalGroup(2, mat)
      for entity_tag in entities:
         element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
         for tags in element_tags:
            if len(tags) > 0:
               for tag in tags:
                  mats[tag] = mat
   mats = [x for x in mats[1:] if x is not None]
   cells = cells[:len(mats)]
   
   boundaries = None
   if not r is None:
      boundaries = [i-1 for i in gmsh.model.mesh.getNodesForPhysicalGroup(1, boundary)[0]]
   
   n = 0
   for reg in regions:
      if not reg is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(2, reg)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  n = max(n, max(tags))
   
   nodes = [None] * (int(n) + 1)
   for reg in regions:
      if not reg is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(2, reg)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  for tag in tags:
                     nodes[tag] = reg
   n = len(materials) + 1
   if not r is None: n += 2
   nodes = [x-n for x in nodes[1:] if x is not None]
   
   gmsh.finalize()
   
   return Mesh(points, cells, None, mats, boundaries, nodes)

def write_mesh(filename, mesh, hb, h, ht, nzb, nz, nzt, rb_mat, rt_mat):
   
   with open(filename, "w") as f:
      
      f.write("points %d\n" % len(mesh.points))
      for p in mesh.points:
         f.write("%.9e %.9e\n" % (p[0], p[1]))
      f.write("\n")
      
      np = 0
      for c in mesh.cells:
         np += len(c)
      f.write("cells %d %d\n" % (len(mesh.cells), np))
      for c in mesh.cells:
         for i, p in enumerate(c):
            if i > 0: f.write(" ")
            f.write("%d" % p)
         f.write("\n")
      f.write("\n")
      
      if not mesh.boundaries is None:
         f.write("boundary exterior %d\n" % len(mesh.boundaries))
         for i in mesh.boundaries:
            f.write("%d\n" % i)
         f.write("\n")
      else:
         f.write("boundary xy 0\n\n")
      
      nztot = nzb + nz + nzt
      if nztot > 1:
         f.write("dz %d\n" % nztot)
         dz = [hb/max(nzb, 1)] * nzb + [h/max(nz, 1)] * nz + [ht/max(nzt, 1)] * nzt
         for i, d in enumerate(dz):
            if i > 0: f.write(" ")
            f.write("%.9e" % d)
         f.write("\n\n")
      
      f.write("materials %d\n" % (nztot*len(mesh.mats)))
      for k in range(nztot):
         f.write("\n")
         for i in range(len(mesh.mats)):
            if i > 0: f.write(" ")
            if k < nzb:
               f.write("%d" % rb_mat[mesh.mats[i]-1])
            elif k < nzb+nz:
               f.write("%d" % mesh.mats[i])
            else:
               f.write("%d" % rt_mat[mesh.mats[i]-1])
         f.write("\n")
      
      if not mesh.nodes is None:
         num_2d_nodes = max(mesh.nodes) + 1
         f.write("\n")
         f.write("nodal-indices %d\n" % (nztot*len(mesh.cells)))
         for k in range(nztot):
            f.write("\n")
            for i in range(len(mesh.nodes)):
               if i > 0: f.write(" ")
               if mesh.nodes[i] < 0:
                  f.write("-1")
               else:
                  f.write("%d" % (mesh.nodes[i]+k*num_2d_nodes))
            f.write("\n")

def main():
   
   # Axial dimensions:
   if dims == 2:
      h = 1.0
      hb = 0.0
      ht = 0.0
   else:
      l = 280.0
      h = 182.0
      hb = 0.5 * (l-h)
      ht = 0.5 * (l-h)
   
   # Case number:
   #    - 1 = 2D convergence for lc1, fuel type 1 (no shutdown rod).
   #    - 2 = 2D convergence for lc2 with fixed lc1, fuel type 1 (no shutdown rod).
   #    - 3 = 2D convergence for lc1, fuel type 3 (with shutdown rod).
   #    - 4 = 2D convergence for lc2 with fixed lc1, fuel type 3 (with shutdown rod).
   #    - 5 = 3D convergence for nz, fuel type 1 (no shutdown rod).
   #    - 6 = 3D convergence for nzb and nzt with fixed nz, fuel type 1 (no shutdown rod).
   #    - 7 = 3D convergence for nz, fuel type 3 (with shutdown rod).
   #    - 8 = 3D convergence for nzb and nzt with fixed nz, fuel type 3 (with shutdown rod).
   #    - 9 = 3D convergence for lc3, full core.
   #    - 10 = 3D convergence for all parameters with fixed ratios, full core.
   case = 5
   if case == 1:
      
      # Case type:
      single = 1
      dims = 2
      
      # Mesh size at the pins (lc1) and the hexagonal grid (lc2):
      l1 = 0.05
      l2 = 0.8
      dl = 0.05
      lc0 = [None, None, None]
      l0 = 0.15
   
   elif case == 2:
      
      # Case type:
      single = 1
      dims = 2
      
      # Mesh size at the pins (lc1) and the hexagonal grid (lc2):
      l1 = 0.15
      l2 = 0.9
      dl = 0.05
      lc0 = [0.15, None, None]
      l0 = 0.3
   
   elif case == 3:
      
      # Case type:
      single = 3
      dims = 2
      
      # Mesh size at the pins (lc1) and the hexagonal grid (lc2):
      l1 = 0.05
      l2 = 0.8
      dl = 0.05
      lc0 = [None, None, None]
      l0 = 0.15
   
   elif case == 4:
      
      # Case type:
      single = 3
      dims = 2
      
      # Mesh size at the pins (lc1) and the hexagonal grid (lc2):
      l1 = 0.15
      l2 = 0.9
      dl = 0.05
      lc0 = [0.15, None, None]
      l0 = 0.3
   
   elif case == 5:
      
      # Case type:
      single = 1
      dims = 3
      
      # Mesh size at the pins (lc1) and the hexagonal grid (lc2):
      lc0 = [0.15, 0.3, None]
      
      # Axial discretization:
      l1 = h / 100
      l2 = h / 10
      dl = h / 100
      l0 = 4 * l1
      nz0 = [None, None, None]
   
   elif case == 6:
      
      # Case type:
      single = 1
      dims = 3
   
   elif case == 7:
      
      # Case type:
      single = 3
      dims = 3
   
   elif case == 8:
      
      # Case type:
      single = 3
      dims = 3
   
   elif case == 9:
      
      # Mesh size at the pins (lc1), the hexagonal grid (lc2) and the outer reflector boundary (lc3):
      l1 = 0.6
      l2 = 1.0
      dl = 0.1
      lc0 = [None, None, None]
      
      # Case type:
      single = 0
      dims = 3
      small = False
   
   elif case == 10:
      
      # Case type:
      single = 0
      dims = 3
      small = False
   
   else:
      raise Exception("Case not implemented!")
   
   # Core geometry:
   # Fuel-assembly types:
   #    - 1 = fuel 1
   #    - 2 = fuel 2
   #    - 3 = fuel 1 + shutdown rod
   #    - 4 = fuel 2 + shutdown rod
   #    - 5 = moderator
   #    - 6 = reflector
   pc = 18.0
   rb_mat = [1, 1, 3, 1]
   rt_mat = [1, 1, 1, 5]
   if single:
      r = None
      core = [[single]]
      if single in [1, 2]:
         case_name = "fa_no_sdr"
         rb_mat = [1, 1, 1]
         rt_mat = [1, 1, 4]
      else:
         case_name = "fa_with_sdr"
   else:
      if small:
         r = 60.0
         core = [[0, 0, 6, 6, 0, 0], \
                   [6, 4, 3, 4, 6], \
                  [6, 3, 1, 2, 3, 6], \
                    [4, 2, 5, 1, 4], \
                   [6, 3, 1, 2, 3, 6], \
                     [6, 4, 3, 4, 6], \
                    [0, 0, 6, 6, 0, 0]]
         case_name = "mini_core"
      else:
         r = 130.0
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
         case_name = "full_core"
   if dims == 2: case_name += "_2d"
   
   # Fuel-assembly geometry:
   # Pin types:
   #    - 1 = fuel
   #    - 2 = shutdown rod
   #    - 3 = heat pipe
   pf = 2.86
   d = [1.7, 4.0, 1.6]
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
                [1, 0, 2, 0, 1], \
               [1, 3, 0, 0, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[3] = [[0, 0, 1, 1, 0, 0], \
               [1, 1, 3, 1, 1], \
              [1, 3, 0, 0, 3, 1], \
                [1, 0, 2, 0, 1], \
               [1, 3, 0, 0, 3, 1], \
                 [1, 1, 3, 1, 1], \
                [0, 0, 1, 1, 0, 0]]
   fas[4] = None
   fas[5] = None
   
   # Axial discretization:
   if dims == 2:
      nzb = 0
      nz = 1
      nzt = 0
   else:
      nzb = 10
      nz = 40
      nzt = 10
   
   # Get the pin power:
   P0 = 15.0e6; n0 = 2592; h0 = 182.0
   n = 0
   core_flat = [x for row in core for x in row]
   for i in range(4):
      fa_flat = [x for row in fas[i] for x in row]
      n += core_flat.count(i+1) * fa_flat.count(1)
   P = (float(n)/n0) * P0
   if dims == 2:
      P /= h0
      print("P = %.3e W/cm" % P)
   else:
      print("P = %.3e W" % P)
   
   # Mesh size:
   l = np.arange(l1, l2+0.1*dl, dl)
   print("Mesh-size points:", len(l))
   
   # Nodal mesh:
   core_ref_mat = [6, 6, 6, 6, 6, 6]
   core_mesh = build_x_hex_mesh(pc, core)
   write_mesh("mesh-nodal.pmp", core_mesh, hb, h, ht, nzb, nz, nzt, core_ref_mat, core_ref_mat)
   
   # Pin mesh for each fuel-assembly type:
   fa_meshes = [build_x_hex_mesh(pf, fa) if not fa is None else None for fa in fas]
   
   case_dir = "case-" + str(case)
   if os.path.exists(case_dir): shutil.rmtree(case_dir)
   os.mkdir(case_dir)
   
   T = [None] * 2; T0 = [None] * 2
   Tmax = [None] * 2; Tmax0 = [None] * 2
   dT = np.full((8, len(l)), np.nan)
   for i, lc in enumerate(l):
      
      lc1 = lc if lc0[0] is None else lc0[0]
      lc2 = lc if lc0[1] is None else lc0[1]
      lc3 = lc if lc0[2] is None else lc0[2]
      mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, lc1, lc2, lc3, pf)
      write_mesh("mesh.pmp", mesh, hb, h, ht, nzb, nz, nzt, rb_mat, rt_mat)
      
      print("lc =", lc)
      subprocess.run(["./run.sh", "input_" + case_name + ".pmp"])
      
      if len(l) == 1:
         return
      
      nodal_fields = parse_vtk_file("output_nodal_0.vtk")
      T[0] = nodal_fields["fuel_temperature"]
      T[1] = nodal_fields["graphite_temperature"]
      Tmax[0] = nodal_fields["fuel_temperature_max"]
      Tmax[1] = nodal_fields["graphite_temperature_max"]
      
      if i == 0:
         T0[0] = np.copy(T[0])
         T0[1] = np.copy(T[1])
         Tmax0[0] = np.copy(Tmax[0])
         Tmax0[1] = np.copy(Tmax[1])
      else:
         error = [None] * 4
         error[0] = np.abs(T[0]-T0[0])
         error[1] = np.abs(T[1]-T0[1])
         error[2] = np.abs(Tmax[0]-Tmax0[0])
         error[3] = np.abs(Tmax[1]-Tmax0[1])
         for j in range(4):
            dT[j][i] = np.max(error[j])
            dT[j+4][i] = np.linalg.norm(error[j]) / np.sqrt(len(error[j]))
      
      if not l0 is None and np.abs(lc-l0) < 1.0e-6:
         os.rename("output_0.vtk", case_dir + "/output.vtk")
      else:
         os.remove("output_0.vtk")
      os.remove("output_nodal_0.vtk")
   
   matplotlib.rcParams.update({'font.size': 12})
   plot_labels = ["Fuel temperature", "Graphite temperature"]
   if len(core_mesh.cells) == 1:
      plot_data = [[0, 1], [2, 3]]
      y_labels = ["Error (K)", "Error (K)"]
      filenames = ["dT.png", "dTmax.png"]
   else:
      plot_data = [[0, 1], [2, 3], [4, 5], [6, 7]]
      y_labels = ["Maximum error (K)", "Maximum error (K)", "L2-norm error (K)", "L2-norm error (K)"]
      filenames = ["dTmax.png", "dTmaxmax.png", "dTmean.png", "dTmaxmean.png"]
   for i, (data, y_label, filename) in enumerate(zip(plot_data, y_labels, filenames)):
      
      fig, ax = plt.subplots()
      
      for j, (k, label) in enumerate(zip(data, plot_labels)):
         ax.plot(l, dT[k], label = label)
      ax.axhline(y = 5.0, linestyle = "dotted")
      
      ax.set_xlabel("Mesh size (cm)")
      ax.set_ylabel(y_label)
      
      ax.legend()
      plt.gca().invert_xaxis()
      plt.savefig(case_dir + "/dT-" + str(i) + ".png")
      plt.clf()

if __name__ == "__main__": main()
