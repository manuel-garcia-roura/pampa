import os
import shutil
import subprocess
import numpy as np
import gmsh
import torch

class Mesh:
   
   def __init__(self, points, cells, centroids, materials, pins, boundaries):
      
      self.points = points
      self.cells = cells
      self.centroids = centroids
      self.materials = materials
      self.pins = pins
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
   
   cells = []; centroids = []; materials = []
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
            materials.append(c)
            
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
   
   return Mesh(points, cells, centroids, materials, None, boundaries)

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
   for c0, fa in zip(core_mesh.centroids, core_mesh.materials):
      
      if fa_meshes[fa-1] is not None:
         
         for c00, m in zip(fa_meshes[fa-1].centroids, fa_meshes[fa-1].materials):
            
            x0 = c0[0] + c00[0]
            y0 = c0[1] + c00[1]
            
            if lc2 is not None:
               
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
            
            else:
               
               r1 = 0.5 * d[m-1]
               
               _, ic1 = add_circle(x0, y0, r1, lc1)
               holes.append(ic1)
               
               is1 = gmsh.model.occ.addPlaneSurface([ic1])
               pins[m-1].append(is1)
   
   _, ic1 = add_circle(0.0, 0.0, r0, lc3)
   is1 = gmsh.model.occ.addPlaneSurface([ic1])
   
   il2, ic2 = add_circle(0.0, 0.0, r, lc3)
   is2 = gmsh.model.occ.addPlaneSurface([ic2, -ic1] + [-x for x in holes])
   
   gmsh.model.occ.synchronize()
   
   regions = []
   num_materials = 4 if pins[1] else 3
   
   regions.append(gmsh.model.addPhysicalGroup(2, pins[3] + [is1, is2], name = "graphite"))
   regions.append(gmsh.model.addPhysicalGroup(2, pins[0], name = "fuel"))
   if pins[1]:
      regions.append(gmsh.model.addPhysicalGroup(2, pins[1], name = "shutdown-rod"))
   regions.append(gmsh.model.addPhysicalGroup(2, pins[2], name = "heat-pipe"))
   
   regions.append(gmsh.model.addPhysicalGroup(2, pins[1] + pins[2] + pins[3] + [is1, is2], name = "non-fuel"))
   if lc2 is not None:
      for i, pin in enumerate(zip(pins[0][::2], pins[0][1::2]), 1):
         regions.append(gmsh.model.addPhysicalGroup(2, pin, name = "fuel-pin-" + str(i)))
   else:
      for i, pin in enumerate(pins[0], 1):
         regions.append(gmsh.model.addPhysicalGroup(2, [pin], name = "fuel-pin-" + str(i)))
   
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
   
   mesh = get_gmsh_mesh_data(regions, boundary, num_materials)
   
   gmsh.finalize()
   
   return mesh

def get_gmsh_mesh_data(regions, boundary, num_materials):
   
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
   
   materials = [None] * (max_tag+1)
   pins = [None] * (max_tag+1)
   for ir, reg in enumerate(regions):
      entities = gmsh.model.getEntitiesForPhysicalGroup(2, reg)
      for entity_tag in entities:
         element_types, element_tags, _ = gmsh.model.mesh.getElements(2, entity_tag)
         for etags in element_tags:
            if len(etags) > 0:
               for i in etags:
                  if ir < num_materials:
                     materials[i] = reg
                  else:
                     pins[i] = reg - num_materials - 2
   materials = [x for x in materials[1:] if x is not None]
   pins = [x if x > -1 else None for x in pins[1:] if x is not None]
   
   boundaries = [i-1 for i in gmsh.model.mesh.getNodesForPhysicalGroup(1, boundary)[0]]
   
   return Mesh(points, cells, None, materials, pins, boundaries)

def write_mesh(filename, mesh, hb, h, ht, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
   
   num_points = len(mesh.points)
   num_cells = len(mesh.cells)
   num_boundaries = len(mesh.boundaries)
   num_vertices = 0
   for c in mesh.cells:
      num_vertices += len(c)
   
   nztot = nzb + nz + nzt
   dz = [hb/max(nzb, 1)] * nzb + [h/max(nz, 1)] * nz + [ht/max(nzt, 1)] * nzt
   
   with open(filename, "w") as f:
      
      f.write("points %d\n" % num_points)
      for p in mesh.points:
         f.write("%.9e %.9e\n" % (p[0], p[1]))
      f.write("\n")
      
      f.write("cells %d %d\n" % (num_cells, num_vertices))
      for c in mesh.cells:
         for i, p in enumerate(c):
            if i > 0: f.write(" ")
            f.write("%d" % p)
         f.write("\n")
      f.write("\n")
      
      if mesh.boundaries is not None:
         f.write("boundary exterior %d\n" % num_boundaries)
         for p in mesh.boundaries:
            f.write("%d\n" % p)
         f.write("\n")
      else:
         f.write("boundary xy 0\n\n")
      
      if nztot > 1:
         f.write("dz %d\n" % nztot)
         for i, d in enumerate(dz):
            f.write("%.9e\n" % d)
         f.write("\n")
      
      f.write("materials %d\n" % (nztot*num_cells))
      for k in range(nztot):
         f.write("\n")
         for i, m in enumerate(mesh.materials):
            if i > 0: f.write(" ")
            if k < nzb:
               f.write("%d" % bottom_ref_mats[m-1])
            elif k < nzb+nz:
               f.write("%d" % m)
            else:
               f.write("%d" % top_ref_mats[m-1])
         f.write("\n")

def write_input(filename, with_sdr, two_dim, power):
   
   with open(filename, "w") as f:
      
      f.write("mesh unstructured mesh.pmp\n")
      f.write("\n")
      
      f.write("material graphite {\n")
      f.write("   thermal-properties graphite-h-451\n")
      f.write("   fuel 0\n")
      f.write("}\n")
      f.write("material fuel {\n")
      f.write("   thermal-properties graphite-matrix-a3-27\n")
      f.write("   fuel 1\n")
      f.write("}\n")
      if with_sdr:
         f.write("material shutdown-rod {\n")
         f.write("   split 0\n")
         f.write("   bc 1\n")
         f.write("}\n")
      f.write("material heat-pipe-active {\n")
      f.write("   split 1\n")
      f.write("   bc 1\n")
      f.write("}\n")
      if not two_dim:
         f.write("material heat-pipe-inactive {\n")
         f.write("   split 0\n")
         f.write("   bc 1\n")
         f.write("}\n")
      f.write("\n")
      
      f.write("solver conduction {\n")
      f.write("   bc exterior convection 5.0e-4 298.0\n")
      if not two_dim:
         f.write("   bc -z convection 5.0e-4 298.0\n")
         f.write("   bc +z convection 5.0e-4 298.0\n")
      if with_sdr:
         f.write("   bc shutdown-rod adiabatic\n")
      f.write("   bc heat-pipe-active convection 750.0e-4 950.0\n")
      if not two_dim:
         f.write("   bc heat-pipe-inactive adiabatic\n")
      f.write("   power %.9e\n" % power)
      f.write("   data data.pmp\n")
      f.write("   convergence temperature max 0 1.0e-3\n")
      f.write("   heat-pipe heat-pipe-active 1 1.6 182.0 0.75 1400.0e-4 790.0\n")
      f.write("}\n")

def write_heat_source(filename, heat_source, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
   
   num_physical_cells = 0
   for m in mesh.materials:
      if bottom_ref_mats[m-1] == 1:
         num_physical_cells += nzb
      if m < 3:
         num_physical_cells += nz
      if top_ref_mats[m-1] == 1:
         num_physical_cells += nzt
   
   nztot = nzb + nz + nzt
   
   with open(filename, "w") as f:
      
      f.write("heat-source %d\n" % num_physical_cells)
      for k in range(nztot):
         for i, m in enumerate(mesh.materials):
            q = None
            if k < nzb:
               if bottom_ref_mats[m-1] == 1:
                  q = 0.0
            elif k < nzb+nz:
               if m < 3:
                  q = heat_source[k-nzb, i]
            else:
               if top_ref_mats[m-1] == 1:
                  q = 0.0
            if q is not None:
               f.write("%.9e\n" % q)

def get_num_pins(core, fas):
   
   num_pins = [0] * 3
   for fa in [fa for row in core for fa in row]:
      if fa > 0 and fas[fa-1] is not None:
         for pin in [pin for row in fas[fa-1] for pin in row]:
            if pin > 0:
               num_pins[pin-1] += 1
   
   return num_pins

def get_heat_source(heat_source_range, mesh, num_pins, nz):
   
   fp = np.random.uniform(1.0-heat_source_range, 1.0+heat_source_range, (nz, num_pins))
   
   num_cells = len(mesh.cells)
   heat_source = np.zeros((nz, num_cells))
   for k in range(nz):
      for i in range(num_cells):
         if mesh.pins[i] is not None:
            heat_source[k, i] = max(fp[k, mesh.pins[i]], 0.0)
   
   return heat_source

def create_new_directories(paths):
   
   for path in paths:
      if os.path.exists(path):
         shutil.rmtree(path)
      os.mkdir(path)

def parse_vtk_file(filename):
   
   scalar_fields = {}
   with open(filename, "r") as file:
      
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

def main():
   
   # Run parameters:
   run = False
   train = True
   
   # Number of training cases:
   num_cases = 100
   
   # Heat-source parameters:
   power_fraction_range = [0.1, 1.1]
   heat_source_range = 0.5
   
   # Mesh parameters:
   small = True
   two_dim = True
   quad = False
   write_vtk = False
   run_fltk = False
   
   # Mesh sizes:
   lc1 = 1.0
   lc2 = None
   lc3 = 5.0
   nzb = 5
   nz = 20
   nzt = 5
   
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
   pin_power = power / num_pins[0]
   r = 0.5 * len(core) * pc
   r0 = 0.75 * pc
   
   # Small geometry:
   if small:
      core = [[0, 0, 6, 6, 0, 0], \
                [6, 4, 2, 3, 6], \
               [6, 2, 1, 1, 2, 6], \
                 [3, 1, 5, 1, 4], \
                [6, 2, 1, 1, 2, 6], \
                  [6, 4, 2, 3, 6], \
                 [0, 0, 6, 6, 0, 0]]
      num_pins = get_num_pins(core, fas)
      power = pin_power * num_pins[0]
      r = 0.5 * len(core) * pc
      r0 = 0.25 * pc
   
   # Reflector materials:
   with_sdr = num_pins[1] > 0
   if with_sdr:
      bottom_ref_mats = [1, 1, 3, 1]
      top_ref_mats = [1, 1, 1, 5]
   else:
      bottom_ref_mats = [1, 1, 1]
      top_ref_mats = [1, 1, 4]
   
   # Two-dimensional geometry:
   if two_dim:
      nzb = 0
      nz = 1
      nzt = 0
      power /= h
   
   # Fields needed for training:
   field_names = ["power", "temperature"]
   
   # Get the training data:
   if run:
      
      # Build and export the mesh:
      fa_meshes = [build_x_hex_mesh(pf, fa) if fa is not None else None for fa in fas]
      core_mesh = build_x_hex_mesh(pc, core)
      mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, r0, lc1, lc2, lc3, quad, write_vtk, run_fltk)
      write_mesh("mesh.pmp", mesh, hb, h, ht, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
      
      # Run all training cases:
      create_new_directories(field_names)
      for i in range(num_cases):
         
         # Calculate and export the relative volumetric heat-source distribution:
         heat_source = get_heat_source(heat_source_range, mesh, num_pins[0], nz)
         write_heat_source("data.pmp", heat_source, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         
         # Export the main input file:
         power_fraction = np.random.uniform(power_fraction_range[0], power_fraction_range[1])
         write_input("input.pmp", with_sdr, two_dim, power_fraction*power)
         
         # Run the heat-conduction solver:
         subprocess.run(["./run.sh", "input.pmp"])
         
         # Get the temperature and power distributions:
         fields = parse_vtk_file("output_0.vtk")
         for name in field_names:
            fields[name].tofile(name + "/" + str(i) + ".np")
   
   # Train the model:
   if train:
      
      # Load the training data:
      fields = []
      for name in field_names:
         field = []
         for i in range(num_cases):
            field.append(np.fromfile(name + "/" + str(i) + ".np"))
         fields.append(field)
      fields = np.array(fields)
      print(fields.shape)
      
      x = torch.rand(5, 3)
      print(x)

if __name__ == "__main__": main()
