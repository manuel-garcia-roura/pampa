import numpy as np
import gmsh
import subprocess
import os

class Mesh:
   
   def __init__(self, points, cells, centroids, mats, boundaries):
      
      self.points = points
      self.cells = cells
      self.centroids = centroids
      self.mats = mats
      self.boundaries = boundaries

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
   
   return Mesh(points, cells, centroids, mats, boundaries)

def add_circle(x0, y0, r, lc):
   
   gmsh.model.occ.addPoint(x0, y0+r, 0.0, lc)
   gmsh.model.occ.addPoint(x0-r, y0, 0.0, lc)
   gmsh.model.occ.addPoint(x0, y0-r, 0.0, lc)
   gmsh.model.occ.addPoint(x0+r, y0, 0.0, lc)
   
   il = gmsh.model.occ.addCircle(x0, y0, 0.0, r, angle1 = 0.0, angle2 = 2*np.pi)
   ic = gmsh.model.occ.addCurveLoop([il])
   
   return il, ic

def build_gmsh_mesh(core_mesh, fa_meshes, d, r, r0, lc1, lc2, lc3, quad, write_vtk, run_fltk):
   
   gmsh.initialize()
   
   gmsh.model.add("eVinci")
   
   pins = [[] for _ in range(4)]; holes = []
   for c0, fa in zip(core_mesh.centroids, core_mesh.mats):
      
      if not fa_meshes[fa-1] is None:
         
         for c00, m in zip(fa_meshes[fa-1].centroids, fa_meshes[fa-1].mats):
            
            x0 = c0[0] + c00[0]
            y0 = c0[1] + c00[1]
            
            r1 = 0.25 * d[m-1]
            r2 = 0.5 * d[m-1]
            r3 = 0.75 * d[m-1]
            
            _, ic1 = add_circle(x0, y0, r1, lc2)
            _, ic2 = add_circle(x0, y0, r2, lc1)
            _, ic3 = add_circle(x0, y0, r3, lc2)
            holes.append(ic3)
            
            is1 = gmsh.model.occ.addPlaneSurface([ic1])
            pins[m-1].append(is1)
            
            is2 = gmsh.model.occ.addPlaneSurface([ic2, -ic1])
            pins[m-1].append(is2)
            
            is3 = gmsh.model.occ.addPlaneSurface([ic3, -ic2])
            pins[3].append(is3)
   
   _, ic1 = add_circle(0.0, 0.0, r0, lc3)
   is1 = gmsh.model.occ.addPlaneSurface([ic1])
   
   il2, ic2 = add_circle(0.0, 0.0, r, lc3)
   is2 = gmsh.model.occ.addPlaneSurface([ic2, -ic1] + [-x for x in holes])
   
   gmsh.model.occ.synchronize()
   
   materials = []
   materials.append(gmsh.model.addPhysicalGroup(2, pins[3] + [is1, is2], name = "graphite"))
   materials.append(gmsh.model.addPhysicalGroup(2, pins[0], name = "fuel"))
   if pins[1]:
      materials.append(gmsh.model.addPhysicalGroup(2, pins[1], name = "shutdown-rod"))
   materials.append(gmsh.model.addPhysicalGroup(2, pins[2], name = "heat-pipe"))
   boundary = gmsh.model.addPhysicalGroup(1, [il2], name = "boundary")
   
   if quad:
      gmsh.option.setNumber("Mesh.Algorithm", 8)
      gmsh.option.setNumber("Mesh.RecombineAll", 1)
      gmsh.model.mesh.recombine()
   
   gmsh.model.mesh.generate(2)
   gmsh.model.mesh.removeDuplicateNodes()
   gmsh.model.occ.removeAllDuplicates()
   
   if write_vtk:
      gmsh.write("mesh.vtk")
   
   if run_fltk:
      gmsh.fltk.run()
   
   mesh = get_gmsh_mesh_data(materials, boundary)
   
   gmsh.finalize()
   
   return mesh

def get_gmsh_mesh_data(materials, boundary):
   
   tags, coordinates, _ = gmsh.model.mesh.getNodes()
   points = [p for _, p in sorted(zip(tags, zip(coordinates[0::3], coordinates[1::3])))]
   
   element_types, element_tags, node_tags = gmsh.model.mesh.getElements()
   cell_tags = []
   for etype, etags in zip(element_types, element_tags):
      if etype in [2, 3]:
         cell_tags.extend(etags)
   num_cells = len(cell_tags)
   min_tag = int(min(cell_tags))
   max_tag = int(max(cell_tags))
   cells = [None] * num_cells
   for etype, etags, ntags in zip(element_types, element_tags, node_tags):
      if etype in [2, 3]:
         n = etype + 1
         for i in range(0, len(etags)):
            cells[int(etags[i])-min_tag] = ntags[i*n:(i+1)*n]
   cells = [[int(i-1) for i in c] for c in cells]
   
   mats = [None] * (max_tag+1)
   for mat in materials:
      entities = gmsh.model.getEntitiesForPhysicalGroup(2, mat)
      for entity_tag in entities:
         element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
         for etags in element_tags:
            if len(etags) > 0:
               for i in etags:
                  mats[i] = mat
   mats = [x for x in mats[1:] if x is not None]
   
   boundaries = [i-1 for i in gmsh.model.mesh.getNodesForPhysicalGroup(1, boundary)[0]]
   
   return Mesh(points, cells, None, mats, boundaries)

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

def write_input(filename, core, fas, twodim, power):
   
   num_pins = get_num_pins(core, fas)
   
   with open(filename, "w") as f:
      
      f.write("# mesh definition:\n")
      f.write("mesh unstructured mesh.pmp\n")
      f.write("\n")
      
      f.write("# material definition:\n")
      
      f.write("material graphite {\n")
      f.write("   thermal-properties graphite-h-451\n")
      f.write("   fuel 0\n")
      f.write("}\n")
      
      f.write("material fuel {\n")
      f.write("   thermal-properties graphite-matrix-a3-27\n")
      f.write("   fuel 1\n")
      f.write("}\n")
      
      if num_pins[1] > 0:
         f.write("material shutdown-rod {\n")
         f.write("   split 0\n")
         f.write("   bc 1\n")
         f.write("}\n")
      
      f.write("material heat-pipe-active {\n")
      f.write("   split 1\n")
      f.write("   bc 1\n")
      f.write("}\n")
      
      if not twodim:
         f.write("material heat-pipe-inactive {\n")
         f.write("   split 0\n")
         f.write("   bc 1\n")
         f.write("}\n")
      
      f.write("\n")
      
      f.write("# conduction solver definition:\n")
      f.write("solver conduction {\n")
      
      f.write("   bc exterior convection 5.0e-4 298.0\n")
      
      if not twodim:
         f.write("   bc -z convection 5.0e-4 298.0\n")
         f.write("   bc +z convection 5.0e-4 298.0\n")
      
      if num_pins[1] > 0:
         f.write("   bc shutdown-rod adiabatic\n")
      
      f.write("   bc heat-pipe-active convection 750.0e-4 950.0\n")
      if not twodim:
         f.write("   bc heat-pipe-inactive adiabatic\n")
      f.write("   power %.9e\n" % power)
      f.write("   convergence temperature max 0 1.0e-3\n")
      f.write("   heat-pipe heat-pipe-active 1 1.6 182.0 0.5 1400.0e-4 790.0\n")
      f.write("}\n")

def get_num_pins(core, fas):
   
   num_pins = [0] * 3
   for fa in [fa for row in core for fa in row]:
      if fa > 0 and not fas[fa-1] is None:
         for pin in [pin for row in fas[fa-1] for pin in row]:
            if pin > 0:
               num_pins[pin-1] += 1
   
   return num_pins

def get_pin_power_density(core, fas, d, h, power):
   
   n = get_num_pins(core, fas)
   volume = n[0] * 0.25 * np.pi * d[0]**2 * h
   
   return power / volume

def main():
   
   # Mesh parameters:
   small = True
   twodim = True
   quad = False
   write_vtk = False
   run_fltk = False
   
   # Mesh sizes:
   lc1 = 0.5
   lc2 = 1.0
   lc3 = 5.0
   if twodim:
      nzb = 0
      nz = 1
      nzt = 0
   else:
      nzb = 4
      nz = 16
      nzt = 4
   
   # Axial dimensions:
   l = 280.0
   h = 182.0
   hb = 0.5 * (l-h)
   ht = 0.5 * (l-h)
   
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
   
   # Core geometry:
   # Fuel-assembly types:
   #    - 1 = fuel 1
   #    - 2 = fuel 2
   #    - 3 = fuel 1 + shutdown rod
   #    - 4 = fuel 2 + shutdown rod
   #    - 5 = moderator
   #    - 6 = reflector
   power = 15.0e6
   pc = 18.0
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
   num_pins = get_num_pins(core, fas)
   if small:
      core = [[0, 0, 6, 6, 0, 0], \
                [6, 3, 2, 4, 6], \
               [6, 2, 1, 1, 2, 6], \
                 [4, 1, 5, 1, 3], \
                [6, 2, 1, 1, 2, 6], \
                  [6, 3, 2, 4, 6], \
                 [0, 0, 6, 6, 0, 0]]
      r = 0.5 * len(core) * pc
      r0 = 0.25 * pc
      case = "mini_core"
      num_pins_small = get_num_pins(core, fas)
      power *= float(num_pins_small[0]) / num_pins[0]
   else:
      r = 0.5 * len(core) * pc
      r0 = 0.75 * pc
      case = "full_core"
   bottom_ref_mat = [1, 1, 3, 1]
   top_ref_mat = [1, 1, 1, 5]
   if twodim:
      case += "_2d"
      power /= h
   
   # Build and export the mesh:
   fa_meshes = [build_x_hex_mesh(pf, fa) if not fa is None else None for fa in fas]
   core_mesh = build_x_hex_mesh(pc, core)
   mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, r0, lc1, lc2, lc3, quad, write_vtk, run_fltk)
   write_mesh("mesh.pmp", mesh, hb, h, ht, nzb, nz, nzt, bottom_ref_mat, top_ref_mat)
   
   write_input("input.pmp", core, fas, twodim, power)
   
   # Run the heat-conduction solver:
   subprocess.run(["./run.sh", "input.pmp"])

if __name__ == "__main__": main()
