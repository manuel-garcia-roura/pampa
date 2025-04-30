import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import torch
from sklearn.model_selection import train_test_split

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
   
   if r0 is not None:
      
      _, ic1 = add_circle(0.0, 0.0, r0, lc3)
      is1 = gmsh.model.occ.addPlaneSurface([ic1])
      
      il2, ic2 = add_circle(0.0, 0.0, r, lc3)
      is2 = gmsh.model.occ.addPlaneSurface([ic2, -ic1] + [-x for x in holes])
      
      cylinders = [is1, is2]
      
   else:
      
      il2, ic2 = add_circle(0.0, 0.0, r, lc3)
      is2 = gmsh.model.occ.addPlaneSurface([ic2] + [-x for x in holes])
      
      cylinders = [is2]
   
   gmsh.model.occ.synchronize()
   
   regions = []
   num_materials = 4 if pins[1] else 3
   
   regions.append(gmsh.model.addPhysicalGroup(2, pins[3] + cylinders, name = "graphite"))
   regions.append(gmsh.model.addPhysicalGroup(2, pins[0], name = "fuel"))
   if pins[1]:
      regions.append(gmsh.model.addPhysicalGroup(2, pins[1], name = "shutdown-rod"))
   regions.append(gmsh.model.addPhysicalGroup(2, pins[2], name = "heat-pipe"))
   
   regions.append(gmsh.model.addPhysicalGroup(2, pins[1] + pins[2] + pins[3] + cylinders, name = "non-fuel"))
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

def create_directory(path):
   
   if os.path.exists(path):
      shutil.rmtree(path)
   os.mkdir(path)

def remove_file(path):
   
   if os.path.exists(path):
      os.remove(path)

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

class HeatConductionDataNormalization:
   
   def __init__(self, data):
      
      self.dmin = np.min(data, axis = (0, 2), keepdims = True)
      self.dmax = np.max(data, axis = (0, 2), keepdims = True)
   
   def apply(self, data, idx = None):
      
      if idx is None:
         return (data-self.dmin) / (self.dmax-self.dmin)
      else:
         return (data-self.dmin[0, idx]) / (self.dmax[0, idx]-self.dmin[0, idx])
   
   def reverse(self, data, idx = None):
      
      if idx is None:
         return data*(self.dmax-self.dmin) + self.dmin
      else:
         return data*(self.dmax[0, idx]-self.dmin[0, idx]) + self.dmin[0, idx]

class HeatConductionDataset(torch.utils.data.Dataset):
   
   def __init__(self, data, probe_locations):
      
      self.heat_sources = torch.tensor(data[:, 0], dtype = torch.float32)
      self.temperatures = torch.tensor(data[:, 1], dtype = torch.float32)
      self.probe_locations = probe_locations
   
   def __len__(self):
      
      return self.heat_sources.shape[0]
   
   def __getitem__(self, idx):
      
      heat_source = self.heat_sources[idx]
      temperature = self.temperatures[idx]
      probe_temperature = temperature[self.probe_locations]
      
      return probe_temperature, temperature, heat_source

class HeatConductionNN(torch.nn.Module):
   
   def __init__(self, input_size, layer_sizes, output_size):
      
      super(HeatConductionNN, self).__init__()
      
      layers = []
      layer_input_size = input_size
      for size in layer_sizes:
         layers.append(torch.nn.Linear(layer_input_size, size))
         layers.append(torch.nn.ReLU())
         layer_input_size = size
      
      self.layers = torch.nn.Sequential(*layers)
      self.output_temperature = torch.nn.Linear(layer_input_size, output_size)
      self.output_heat_source = torch.nn.Linear(layer_input_size, output_size)
   
   def forward(self, probe_temperature):
      
      layer_output = self.layers(probe_temperature)
      
      return self.output_temperature(layer_output), self.output_heat_source(layer_output)

def main():
   
   # Workflow parameters:
   run = False
   train = True
   test = True
   reproducibility = True
   
   # Number of temperature measurements:
   num_probes = 50
   
   # Training parameters:
   num_samples = 10000
   num_epochs = 200
   batch_size = 500
   test_size = 0.1
   learning_rate = 0.001
   
   # Neural-network architecture:
   # Some tips:
   #    - If performance plateaus early, try adding depth (more layers).
   #    - If training is unstable or slow, try reducing width.
   #    - Overfitting? Consider smaller network, dropout, or regularization.
   # Some possible configurations:
   #    - [2, 1]: simple, decent baseline.
   #    - [4, 2, 1]: wider, more expressive.
   #    - [4, 4, 2, 1]: deeper model.
   layer_sizes = [4, 2, 1]
   
   # Heat-source parameters:
   power_fraction_range = [0.9, 1.1]
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
      # core = [[0, 0, 6, 6, 0, 0], \
      #           [6, 4, 2, 3, 6], \
      #          [6, 2, 1, 1, 2, 6], \
      #            [3, 1, 5, 1, 4], \
      #           [6, 2, 1, 1, 2, 6], \
      #             [6, 4, 2, 3, 6], \
      #            [0, 0, 6, 6, 0, 0]]
      # core = [[0, 1, 1, 0], \
      #           [1, 3, 1], \
      #          [0, 1, 1, 0]]
      core = [[1]]
      num_pins = get_num_pins(core, fas)
      power = pin_power * num_pins[0]
      # r = 0.5 * len(core) * pc
      r = 0.6 * len(core) * pc
      # r0 = 0.25 * pc
      r0 = None
   
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
   
   # Get the training data:
   if run:
      
      # Build and export the mesh:
      fa_meshes = [build_x_hex_mesh(pf, fa) if fa is not None else None for fa in fas]
      core_mesh = build_x_hex_mesh(pc, core)
      mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, r0, lc1, lc2, lc3, quad, write_vtk, run_fltk)
      write_mesh("mesh.pmp", mesh, hb, h, ht, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
      
      # Run all training cases:
      print("\nRun training cases...\n")
      create_directory("data")
      for i in range(num_samples):
         
         # Calculate and export the relative volumetric heat-source distribution:
         heat_source = get_heat_source(heat_source_range, mesh, num_pins[0], nz)
         write_heat_source("data.pmp", heat_source, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         
         # Export the main input file:
         power_fraction = np.random.uniform(power_fraction_range[0], power_fraction_range[1])
         write_input("input.pmp", with_sdr, two_dim, power_fraction*power)
         
         # Run the heat-conduction solver:
         print("Run training case %d..." % i)
         subprocess.run(["./run.sh", "input.pmp"])
         
         # Get the temperature and power distributions:
         data = parse_vtk_file("output_0.vtk")
         data = np.array([data[name] for name in ["power", "temperature"]])
         data.tofile("data/" + str(i) + ".np")
      
      # Get the mesh in .vtk format:
      with open("output_0.vtk", "r") as fin, open("data/mesh.vtk", "w") as fout:
         for line in fin:
            if "CELL_DATA" in line: break
            fout.write(line)
      
      # Clean up solver files:
      for path in ["mesh.pmp", "data.pmp", "input.pmp", "output_0.vtk"]:
         remove_file(path)
      
      print("Done running training cases.")
   
   # Prepare the training data:
   if train or test:
      
      # Set the random seed for reproducibility:
      if reproducibility:
         seed = 42
         np.random.seed(seed)
         torch.manual_seed(seed)
         torch.backends.cudnn.deterministic = True
         torch.backends.cudnn.benchmark = False
      else:
         seed = None
      
      # Load the training data:
      data = []
      for i in range(num_samples):
         data.append(np.fromfile("data/" + str(i) + ".np").reshape(2, -1))
      data = np.array(data)
      num_cells = data.shape[2]
      
      # Normalize the training data:
      normalization = HeatConductionDataNormalization(data)
      data = normalization.apply(data)
      
      # Get random probe locations:
      probe_locations = torch.randint(0, num_cells, (num_probes,))
      
      # Split the data for training and testing:
      train_data, test_data = train_test_split(data, test_size = test_size, random_state = seed)
      
      # Wrap the training and testing data:
      train_dataset = HeatConductionDataset(train_data, probe_locations)
      train_loader = torch.utils.data.DataLoader(train_dataset, batch_size = batch_size, shuffle = True)
      test_dataset = HeatConductionDataset(test_data, probe_locations)
      test_loader = torch.utils.data.DataLoader(test_dataset, batch_size = 1, shuffle = False)
      
      # Initialize the neural network:
      layer_sizes = [int(n*num_cells) for n in layer_sizes]
      model = HeatConductionNN(num_probes, layer_sizes, num_cells)
      
      # Define optimizer and loss functions:
      optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate)
      criterion = torch.nn.MSELoss()
   
   # Train the model:
   if train:
      
      # Set training mode:
      model.train()
      
      # Loss weighting:
      w = 0.1
      dynamic_weighting = True
      
      # Run the training loop:
      print("\nTrain the neural network...\n")
      losses = np.zeros((3, num_epochs))
      for i in range(num_epochs):
         for probe_temperature_batch, temperature_batch, heat_source_batch in train_loader:
            
            # Forward pass:
            optimizer.zero_grad()
            temperature, heat_source = model(probe_temperature_batch)
            
            # Compute the loss (MSE for both temperature and heat source):
            loss_temperature = criterion(temperature, temperature_batch)
            loss_heat_source = criterion(heat_source, heat_source_batch)
            loss = w * loss_temperature + (1.0-w) * loss_heat_source
            
            # Backward pass and optimization:
            loss.backward()
            optimizer.step()
            
            # Accumulate the losses:
            losses[0, i] += loss_temperature.item()
            losses[1, i] += loss_heat_source.item()
            losses[2, i] += loss.item()
         
         # Adjust the loss weighting:
         if dynamic_weighting:
            w = losses[1, i] / (losses[0, i] + losses[1, i])
         
         # Print the loss:
         losses[:, i] /= len(train_loader)
         if i % 10 == 9:
            print("Epoch %d/%d loss: %.3e (T), %.3e (q), %.3e (total), %.3e (w)" % (i+1, num_epochs, losses[0, i], losses[1, i], losses[2, i], w))
      
      # Plot the loss:
      x = np.arange(1, num_epochs+1, dtype = int)
      plt.plot(x, losses[0], label = "Temperature")
      plt.plot(x, losses[1], label = "Heat source")
      plt.plot(x, losses[2], label = "Total")
      plt.xlabel("Epoch")
      plt.ylabel("Loss")
      plt.legend()
      plt.grid(True)
      plt.savefig("loss.png")
      
      # Save the neural network:
      torch.save(model.state_dict(), "model.torch")
      
      print("\nDone training the neural network.")
   
   # Test the model:
   if test:
      
      # Load the neural network:
      model.load_state_dict(torch.load("model.torch", weights_only = True))
      
      # Set testing mode:
      model.eval()
      
      # Run the testing loop:
      print("\nTest the neural network...\n")
      loss_temperature_total = 0.0; loss_heat_source_total = 0.0; loss_total = 0.0
      worse_case = [0.0, None, None, None, None, None, None]
      with torch.no_grad():
         for probe_temperature_batch, temperature_batch, heat_source_batch in test_loader:
            
            # Get the predicted values:
            temperature, heat_source = model(probe_temperature_batch)
            
            # Compute the loss (MSE for both temperature and heat source):
            loss_temperature = criterion(temperature, temperature_batch)
            loss_heat_source = criterion(heat_source, heat_source_batch)
            
            # Compute the total loss:
            loss_temperature_total += loss_temperature.item()
            loss_heat_source_total += loss_heat_source.item()
            loss_total += (loss_temperature + loss_heat_source).item()
            
            # Get the worst temperature prediction:
            if loss_temperature.item() > worse_case[0]:
               worse_case[0] = loss_temperature.item()
               worse_case[1] = heat_source_batch[0]
               worse_case[2] = heat_source[0]
               worse_case[4] = temperature_batch[0]
               worse_case[5] = temperature[0]
      
      # Print the mean loss:
      loss_temperature_mean = loss_temperature_total / len(test_loader)
      loss_heat_source_mean = loss_heat_source_total / len(test_loader)
      loss_mean = loss_total / len(test_loader)
      print("Temperature mean loss: %.3e" % loss_temperature_mean)
      print("Heat-source mean loss: %.3e" % loss_heat_source_mean)
      print("Total mean loss: %.3e" % loss_mean)
      
      # Denormalize the temperature and heat-source data:
      worse_case[1] = normalization.reverse(worse_case[1].numpy(), 0)
      worse_case[2] = normalization.reverse(worse_case[2].numpy(), 0)
      worse_case[3] = 100.0 * (worse_case[2]-worse_case[1]) / np.mean(worse_case[1])
      worse_case[4] = normalization.reverse(worse_case[4].numpy(), 1)
      worse_case[5] = normalization.reverse(worse_case[5].numpy(), 1)
      worse_case[6] = worse_case[5] - worse_case[4]
      
      # Export the worst temperature prediction in .vtk format:
      field_names = ["heat-source-reference", "heat-source", "heat-source-error", "temperature-reference", "temperature", "temperature-error"]
      shutil.copyfile("data/mesh.vtk", "worst-temperature.vtk")
      with open("worst-temperature.vtk", "a") as f:
         f.write("CELL_DATA %d\n" % num_cells)
         for data, name in zip(worse_case[1:], field_names):
            f.write("SCALARS %s double 1\n" % name)
            f.write("LOOKUP_TABLE default\n")
            for x in data:
               f.write("%.9e\n" % x)
            f.write("\n")
      
      print("\nDone evaluating the neural network.")

if __name__ == "__main__": main()
