import numpy as np
import matplotlib.pyplot as plt
import gmsh
import subprocess

class Mesh:

   def __init__(self, points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts, nodes):

      self.points = points
      self.hex_centroids = hex_centroids
      self.tri_centroids = tri_centroids
      self.hex_cells = hex_cells
      self.tri_cells = tri_cells
      self.hex_mats = hex_mats
      self.tri_mats = tri_mats
      self.bc_pts = bc_pts
      self.nodes = nodes

def build_x_hex_mesh(p, layout):
   
   ny = 2*len(layout[0]) + 1
   nx = int(1.5*len(layout)) + 2
   
   x, y = np.meshgrid(np.arange(float(nx), 0.0, -1.0), np.arange(float(ny)))
   x[1::2, :] += 0.5
   dy = 0.5 * p
   dx = (2.0/np.sqrt(3.0)) * 0.5 * p
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
   
   return Mesh(points, hex_centroids, tri_centroids, hex_cells, tri_cells, hex_mats, tri_mats, bc_pts, None)

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

def build_gmsh_mesh(core_mesh, fa_meshes, d, r, lc1, lc2, lc3):
   
   write_vtk = True
   run_fltk = False
   
   gmsh.initialize()
   
   gmsh.model.add("hex-core")
   
   pin_pts = [False] * len(core_mesh.points)
   for c, fa in zip(core_mesh.hex_cells, core_mesh.hex_mats):
      if not fa_meshes[fa-1] is None:
         for p in c:
            pin_pts[p] = True
   
   for ip, p in enumerate(core_mesh.points):
      if ip in core_mesh.bc_pts or not pin_pts[ip]:
         lc = lc3
      else:
         lc = lc1
      gmsh.model.occ.addPoint(p[0], p[1], 0.0, lc)
   
   fas = []; hxs = []; pins = [[] for _ in range(4)]; regions = [None]; pin_bcs = [[] for _ in range(4)]
   for c0, c, fa in zip(core_mesh.hex_centroids, core_mesh.hex_cells, core_mesh.hex_mats):
      
      sides = []
      for p in range(len(c)):
         
         p1 = c[p] + 1
         p2 = c[(p+1)%len(c)] + 1
         il = gmsh.model.occ.addLine(p1, p2)
         sides.append(il)
      
      holes = []
      if not fa_meshes[fa-1] is None:
         
         for c00, m in zip(fa_meshes[fa-1].hex_centroids, fa_meshes[fa-1].hex_mats):
            
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]+0.5*d[m-1], 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]-0.5*d[m-1], c0[1]+c00[1], 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0], c0[1]+c00[1]-0.5*d[m-1], 0.0, lc2)
            gmsh.model.occ.addPoint(c0[0]+c00[0]+0.5*d[m-1], c0[1]+c00[1], 0.0, lc2)
            
            il = gmsh.model.occ.addCircle(c0[0]+c00[0], c0[1]+c00[1], 0.0, 0.5*d[m-1], angle1 = 0.0, angle2 = 2*np.pi)
            icl = gmsh.model.occ.addCurveLoop([il])
            ips = gmsh.model.occ.addPlaneSurface([icl])
            
            holes.append(ips)
            pins[m-1].append(ips)
      
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
   
   materials = [None] * 4
   materials[0] = (2, gmsh.model.addPhysicalGroup(2, fas + pins[3] + reflector, name = "graphite"))
   materials[1] = (2, gmsh.model.addPhysicalGroup(2, pins[0], name = "fuel"))
   if pins[1]:
      materials[2] = (2, gmsh.model.addPhysicalGroup(2, pins[1], name = "shutdown-rod"))
   else:
      materials[2] = (2, None)
   materials[3] = (2, gmsh.model.addPhysicalGroup(2, pins[2], name = "heat-pipe"))
   if not r is None:
      boundary = gmsh.model.addPhysicalGroup(1, [boundary], name = "boundary")
      regions[0] = (2, gmsh.model.addPhysicalGroup(2, reflector, name = "reflector"))
   else:
      regions[0] = (2, None)
   for ir, reg in enumerate(regions[1:], 1):
      regions[ir] = (2, gmsh.model.addPhysicalGroup(2, reg, name = "node-" + str(ir)))
   
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
   cells = [[i-1 for i in c] for c in cells]
   
   n = 0
   for dim, mat in materials:
      if not mat is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(dim, mat)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  n = max(n, max(tags))
   
   mats = [None] * (int(n) + 1)
   for dim, mat in materials:
      if not mat is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(dim, mat)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  for tag in tags:
                     mats[tag] = mat
   mats = [x for x in mats[1:] if x is not None]
   
   cells = cells[:len(mats)]
   
   bc_pts = None
   if not r is None:
      bc_pts = [i-1 for i in gmsh.model.mesh.getNodesForPhysicalGroup(1, boundary)[0]]
   
   n = 0
   for dim, reg in regions:
      if not reg is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(dim, reg)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  n = max(n, max(tags))
   
   nodes = [None] * (int(n) + 1)
   for dim, reg in regions:
      if not reg is None:
         entities = gmsh.model.getEntitiesForPhysicalGroup(dim, reg)
         for entity_tag in entities:
            element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
            for tags in element_tags:
               if len(tags) > 0:
                  for tag in tags:
                     nodes[tag] = reg
   if not r is None:
      nodes = [x-7 for x in nodes[1:] if x is not None]
   else:
      nodes = [x-5 for x in nodes[1:] if x is not None]
   
   gmsh.finalize()
   
   return Mesh(points, None, None, None, cells, None, mats, bc_pts, nodes)

def write_mesh(filename, mesh, cells, mats, hb, h, ht, nzb, nz, nzt, rb_mat, rt_mat, nodes):
   
   num_2d_cells = len(cells)
   
   with open(filename, "w") as f:
      
      f.write("points %d\n" % len(mesh.points))
      for p in mesh.points:
         f.write("%.3f %.3f\n" % (p[0], p[1]))
      f.write("\n")
      
      np = 0
      for c in cells:
         np += len(c)
      f.write("cells %d %d\n" % (num_2d_cells, np))
      for c in cells:
         for i, p in enumerate(c):
            if i > 0: f.write(" ")
            f.write("%d" % p)
         f.write("\n")
      f.write("\n")
      
      nztot = nzb + nz + nzt
      if nztot > 0:
         f.write("dz %d\n" % nztot)
         dz = [hb/max(nzb, 1)] * nzb + [h/max(nz, 1)] * nz + [ht/max(nzt, 1)] * nzt
         for i, d in enumerate(dz):
            if i > 0: f.write(" ")
            f.write("%.3f" % d)
         f.write("\n\n")
      
      if not mesh.bc_pts is None:
         f.write("boundary exterior %d\n" % len(mesh.bc_pts))
         for i in mesh.bc_pts:
            f.write("%d\n" % i)
         f.write("\n")
      else:
         f.write("boundary xy 0\n\n")
      
      f.write("materials %d\n" % (nztot*len(mats)))
      for k in range(nztot):
         f.write("\n")
         for i in range(len(mats)):
            if i > 0: f.write(" ")
            if k < nzb:
               f.write("%d" % rb_mat[mats[i]-1])
            elif k < nzb+nz:
               f.write("%d" % mats[i])
            else:
               f.write("%d" % rt_mat[mats[i]-1])
         f.write("\n")
      
      if not nodes is None:
         num_2d_nodes = max(nodes) + 1
         f.write("\n")
         f.write("nodal-indices %d\n" % (nztot*num_2d_cells))
         for k in range(nztot):
            f.write("\n")
            for i in range(len(nodes)):
               if i > 0: f.write(" ")
               if nodes[i] < 0:
                  f.write("-1")
               else:
                  f.write("%d" % (nodes[i]+k*num_2d_nodes))
            f.write("\n")

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
               scalar_fields[current_scalar_name] = current_scalar_data
            
            parts = line.split()
            current_scalar_name = parts[1]
            current_scalar_data = []
            is_reading_scalars = True
         
         elif is_reading_scalars and line.startswith("LOOKUP_TABLE"):
            continue
         
         elif is_reading_scalars and not line.startswith(("SCALARS", "LOOKUP_TABLE")):
            current_scalar_data.extend(map(float, line.split()))
      
      if current_scalar_name:
         scalar_fields[current_scalar_name] = current_scalar_data
   
   return scalar_fields

def main():
   
   # Core geometry:
   # Fuel-assembly types:
   #    - 1 = fuel 1
   #    - 2 = fuel 2
   #    - 3 = fuel 1 + shutdown rod
   #    - 4 = fuel 2 + shutdown rod
   #    - 5 = moderator
   #    - 6 = reflector
   pc = 18.0
   single = True
   small = False
   dims = 2
   if single:
      r = None
      core = [[1]]
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
   rb_mat = [1, 1, 3, 1]
   rt_mat = [1, 1, 1, 5]
   
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
   
   # Axial discretization:
   if dims == 2:
      nzb = 0
      nz = 1
      nzt = 0
   else:
      nzb = 10
      nz = 40
      nzt = 10
   
   # Mesh size at the hexagonal grid (lc1), the pins (lc2) and the outer reflector boundary (lc3):
   lc1 = 0.01
   lc2 = 0.10
   dlc = 0.01
   lc = np.arange(lc1, lc2+0.1*dlc, dlc)
   lc = [0.5]
   
   # Nodal mesh:
   core_ref_mat = [6, 6, 6, 6, 6, 6]
   core_mesh = build_x_hex_mesh(pc, core)
   write_mesh("mesh-nodal.pmp", core_mesh, core_mesh.hex_cells, core_mesh.hex_mats, hb, h, ht, nzb, nz, nzt, core_ref_mat, core_ref_mat, None)
   
   # Pin mesh for each fuel-assembly type:
   fa_meshes = [build_x_hex_mesh(pf, fa) if not fa is None else None for fa in fas]
   
   T = np.empty(shape = (3, len(lc)))
   for i, l in enumerate(lc):
      
      mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, 2.0*l, l, 3.0*l)
      write_mesh("mesh.pmp", mesh, mesh.tri_cells, mesh.tri_mats, hb, h, ht, nzb, nz, nzt, rb_mat, rt_mat, mesh.nodes)
      
      print("lc = ", l)
      subprocess.run(["./run.sh"])
      
      if len(lc) == 1:
         return
      
      nodal_fields = parse_vtk_file("output_nodal_0.vtk")
      T[0][i] = np.max(np.array(nodal_fields["fuel_temperature_max"]))
      T[1][i] = np.max(np.array(nodal_fields["graphite_temperature_max"]))
      T[2][i] = np.max(np.array(nodal_fields["heat-pipe-active_temperature_max"]))
   
   labels = ["Peak fuel temperature (K)", "Peak graphite temperature (K)", "Peak heat-pipe temperature (K)"]
   filenames = ["Tf.png", "Tg.png", "Thp.png"]
   for i in range(3):
      
      dT = np.append(np.abs(T[i][:-1]-T[i][1:]), np.nan)
      
      fig, ax1 = plt.subplots()
      ax2 = ax1.twinx()
      
      ax1.plot(lc, T[i], "b-")
      ax2.plot(lc, dT, "r-")
      ax2.axhline(y = 1.0, color = "r", linestyle = "--")
      
      ax1.set_ylim([np.min(T[i]), np.max(T[i])+1.0])
      ax2.set_ylim([np.min(dT[:-1])-1.0, np.max(dT[:-1])])
      
      ax1.set_xlabel("Mesh size (cm)")
      ax1.set_ylabel(labels[i])
      ax2.set_ylabel("Difference (K)")
      
      plt.gca().invert_xaxis()
      plt.savefig(filenames[i])
      plt.clf()

if __name__ == "__main__": main()
