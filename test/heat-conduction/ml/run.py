import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import gmsh
import torch
import torch_geometric
import torch_scatter
from collections import defaultdict
from itertools import combinations
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

def create_directory(path):
   
   if os.path.exists(path):
      shutil.rmtree(path)
   os.mkdir(path)

def remove_file(path):
   
   if os.path.exists(path):
      os.remove(path)

def parse_vtk_file(filename, field = None):
   
   scalar_fields = {}
   
   with open(filename, "r") as f:
      
      is_reading_scalars = False
      current_scalar_name = None
      current_scalar_data = []
      
      for line in f:
         
         if line.startswith("SCALARS"):
            
            if current_scalar_name:
               scalar_fields[current_scalar_name] = np.array(current_scalar_data)
            
            parts = line.split()
            current_scalar_name = parts[1]
            current_scalar_data = []
            is_reading_scalars = field is None or parts[1] == field
         
         elif is_reading_scalars and line.startswith("LOOKUP_TABLE"):
            continue
         
         elif is_reading_scalars and not line.startswith(("SCALARS", "LOOKUP_TABLE")):
            current_scalar_data.extend(map(float, line.split()))
      
      if current_scalar_name:
         scalar_fields[current_scalar_name] = np.array(current_scalar_data)
   
   if field is not None:
      scalar_fields = scalar_fields[field]
   
   return scalar_fields

def get_num_pins(core, fas):
   
   num_pins = [0] * 3
   for fa in [fa for row in core for fa in row]:
      if fa > 0 and fas[fa-1] is not None:
         for pin in [pin for row in fas[fa-1] for pin in row]:
            if pin > 0:
               num_pins[pin-1] += 1
   
   return num_pins

def get_pin_power(pin_power_range, mesh, num_pins, nz):
   
   pin_power = np.random.uniform(*pin_power_range, (nz, num_pins))
   pin_power /= np.sum(pin_power)
   
   return pin_power

def get_heat_source(pin_power, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
   
   num_physical_cells = 0
   for m in mesh.materials:
      if bottom_ref_mats[m-1] == 1:
         num_physical_cells += nzb
      if m < 3:
         num_physical_cells += nz
      if top_ref_mats[m-1] == 1:
         num_physical_cells += nzt
   
   nztot = nzb + nz + nzt
   
   ic = 0
   heat_source = np.empty(num_physical_cells)
   for k in range(nztot):
      for i, m in enumerate(mesh.materials):
         q = None
         if k < nzb:
            if bottom_ref_mats[m-1] == 1:
               q = 0.0
         elif k < nzb + nz:
            if m < 3:
               if mesh.pins[i] is not None:
                  q = pin_power[k-nzb, mesh.pins[i]]
               else:
                  q = 0.0
         else:
            if top_ref_mats[m-1] == 1:
               q = 0.0
         if q is not None:
            heat_source[ic] = q
            ic += 1
   
   return heat_source

def get_cell_to_pin_mapping(mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
   
   num_pins = max([ip for ip in mesh.pins if ip is not None]) + 1
   num_physical_cells = 0
   for m in mesh.materials:
      if bottom_ref_mats[m-1] == 1:
         num_physical_cells += nzb
      if m < 3:
         num_physical_cells += nz
      if top_ref_mats[m-1] == 1:
         num_physical_cells += nzt
   
   nztot = nzb + nz + nzt
   
   ic = 0
   cell_to_pin_mapping = np.empty(num_physical_cells, dtype = int)
   for k in range(nztot):
      for i, m in enumerate(mesh.materials):
         ip = None
         if k < nzb:
            if bottom_ref_mats[m-1] == 1:
               ip = -1
         elif k < nzb + nz:
            if m < 3:
               if mesh.pins[i] is not None:
                  ip = (k-nzb)*num_pins + mesh.pins[i]
               else:
                  ip = -1
         else:
            if top_ref_mats[m-1] == 1:
               ip = -1
         if ip is not None:
            cell_to_pin_mapping[ic] = ip
            ic += 1
   
   return torch.tensor(cell_to_pin_mapping, dtype = torch.long)

def get_cell_adjacency_graph(filename):
   
   with open(filename, "r") as f:
      
      cells = None
      
      for line in f:
         
         if cells is not None:
            
            if not line or line == "\n" or "CELL_TYPES" in line:
               break
            
            line = line.split()
            n = int(line[0])
            for j in range(n):
               cells[i, j] = int(line[j+1])
            i += 1
         
         if "CELLS" in line:
            
            line = line.split()
            n = int(line[1])
            m = int(line[2])
            cells = np.empty((n, int(m/n-1)), dtype = int)
            i = 0
   
   point_to_cells = defaultdict(set)
   for c, points in enumerate(cells):
      for p in points:
         point_to_cells[p].add(c)
   
   edges = set()
   for cells in point_to_cells.values():
      for c1, c2 in combinations(cells, 2):
         edges.add((c1, c2))
         edges.add((c2, c1))
   
   return torch.tensor(list(edges), dtype = torch.long).t().contiguous()

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
            elif k < nzb + nz:
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

def write_heat_source(filename, pin_power, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
   
   heat_source = get_heat_source(pin_power, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
   num_physical_cells = len(heat_source)
   
   with open(filename, "w") as f:
      
      f.write("heat-source %d\n" % num_physical_cells)
      for q in heat_source:
         f.write("%.9e\n" % q)

def load_training_data(path, num_samples):
   
   data = []
   for i in range(num_samples):
      data.append(np.fromfile(path + "/" + str(i) + ".np"))
   data = np.array(data)
   
   normalization = DataNormalization(data)
   data = normalization.apply(data)
   
   return data, normalization

def wrap_training_data(num_samples, batch_size, test_size, seed, pin_power, temperature, probe_temperature):
   
   indices = np.arange(0, num_samples, dtype = int)
   train_indices, test_indices = train_test_split(indices, test_size = test_size, random_state = seed)
   
   train_pin_power = pin_power[train_indices] if pin_power is not None else None
   train_temperature = temperature[train_indices] if temperature is not None else None
   train_dataset = HeatConductionDataset(train_pin_power, train_temperature, probe_temperature[train_indices])
   train_loader = torch.utils.data.DataLoader(train_dataset, batch_size = batch_size, shuffle = True)
   
   test_pin_power = pin_power[test_indices] if pin_power is not None else None
   test_temperature = temperature[test_indices] if temperature is not None else None
   test_dataset = HeatConductionDataset(test_pin_power, test_temperature, probe_temperature[test_indices])
   test_loader = torch.utils.data.DataLoader(test_dataset, batch_size = 1, shuffle = False)
   
   pin_power = None; temperature = None; probe_temperature = None
   
   return train_loader, test_loader

class DataNormalization:
   
   def __init__(self, data):
      
      self.dmin = np.min(data)
      self.dmax = np.max(data)
   
   def apply(self, data):
      
      return (data-self.dmin) / (self.dmax-self.dmin)
   
   def reverse(self, data):
      
      return data*(self.dmax-self.dmin) + self.dmin

class HeatConductionDataset(torch.utils.data.Dataset):
   
   def __init__(self, pin_power, temperature, probe_temperature):
      
      self.pin_power = torch.tensor(pin_power, dtype = torch.float32) if pin_power is not None else None
      self.temperature = torch.tensor(temperature, dtype = torch.float32) if temperature is not None else None
      self.probe_temperature = torch.tensor(probe_temperature, dtype = torch.float32)
   
   def __len__(self):
      
      return self.probe_temperature.shape[0]
   
   def __getitem__(self, idx):
      
      pin_power = self.pin_power[idx] if self.pin_power is not None else torch.empty(0)
      temperature = self.temperature[idx] if self.temperature is not None else torch.empty(0)
      
      return pin_power, temperature, self.probe_temperature[idx]

def build_mlp(input_size, layer_sizes, output_size = None):
   
   if layer_sizes is None:
      return None
   
   layers = []
   in_size = input_size
   for out_size in layer_sizes:
      layers.append(torch.nn.Linear(in_size, out_size))
      layers.append(torch.nn.ReLU())
      in_size = out_size
   
   if output_size is not None:
      layers.append(torch.nn.Linear(in_size, output_size))
   
   return torch.nn.Sequential(*layers)

def build_gnn(input_size, layer_sizes, output_size = None):
   
   if layer_sizes is None:
      return None
   
   layers = torch.nn.ModuleList()
   in_size = input_size
   for out_size in layer_sizes:
      layers.append(torch_geometric.nn.GCNConv(in_size, out_size))
      in_size = out_size
   
   if output_size is not None:
      layers.append(torch_geometric.nn.GCNConv(in_size, output_size))
   
   return layers

# Pure MLP model:
#    - Shared MLP encoder.
#    - MLP heads for the pin-power and temperature tasks.
#    - No geometrical awareness.
class MLPModel(torch.nn.Module):
   
   def __init__(self, num_probes, num_pins, num_cells, encoder_sizes, pin_power_head_sizes, temperature_head_sizes):
      
      super(MLPModel, self).__init__()
      
      # Shared MLP encoder:
      self.encoder = build_mlp(num_probes, encoder_sizes)
      
      # Pin-power MLP head:
      self.pin_power_head = build_mlp(encoder_sizes[-1], pin_power_head_sizes, num_pins)
      
      # Temperature MLP head:
      self.temperature_head = build_mlp(encoder_sizes[-1], temperature_head_sizes, num_cells)
   
   def forward(self, probe_temperature):
      
      # Shared encoder:
      x = self.encoder(probe_temperature)
      
      # Task-specific heads:
      pin_power = self.pin_power_head(x) if self.pin_power_head is not None else torch.empty(0)
      temperature = self.temperature_head(x) if self.temperature_head is not None else torch.empty(0)
      
      return pin_power, temperature

# GNN + MLP model:
#    - Shared GNN encoder.
#    - MLP head for the pin-power task.
#    - GNN head for the temperature task.
#    - Cell-adjacency for the shared encoder and the temperature head.
#    - Cell-to-pin mapping between the GNN and the MLP for the pin-power head (optional).
class GNNModel(torch.nn.Module):
   
   def __init__(self, num_pins, num_cells, probe_locations, cell_adjacency, cell_to_pin_mapping, encoder_sizes, pin_power_head_sizes, temperature_head_sizes):
      
      super(GNNModel, self).__init__()
      
      self.num_pins = num_pins
      self.num_cells = num_cells
      self.register_buffer('probe_locations', probe_locations)
      self.register_buffer('cell_adjacency', cell_adjacency)
      self.register_buffer('cell_to_pin_mapping', cell_to_pin_mapping)
      
      # Fixed binary mask (1 at probe locations, 0 elsewhere):
      probe_mask = torch.zeros(num_cells, dtype = torch.float32)
      probe_mask[probe_locations] = 1.0
      self.register_buffer('probe_mask', probe_mask)
      
      # Shared GNN encoder:
      #    - Input: temperature + binary mask.
      #    - Output: intermediate representation with F = encoder_sizes[-1] features per cell.
      self.encoder = build_gnn(2, encoder_sizes)
      
      # Pin-power MLP head:
      #    - Input: F features per cell (flattened, optionally collapsed using the cell-to-pin mapping).
      #    - Output: full pin-power distribution.
      input_size = num_pins if cell_to_pin_mapping is not None else num_cells
      self.pin_power_head = build_mlp(input_size*encoder_sizes[-1], pin_power_head_sizes, num_pins)
      
      # Temperature GNN head:
      #    - Input: F features per cell (graph structure preserved).
      #    - Output: full temperature distribution (collapsed using an extra single-feature layer).
      self.temperature_head = build_gnn(encoder_sizes[-1], temperature_head_sizes, 1)
      
      self.activation = torch.nn.ReLU()
   
   def forward(self, probe_temperature):
      
      batch_size = probe_temperature.shape[0]
      
      # Expand the probe values and mask them to the full mesh:
      probe_values = torch.zeros(batch_size, self.num_cells)
      probe_values[:, self.probe_locations] = probe_temperature
      probe_mask = self.probe_mask.unsqueeze(0).repeat(batch_size, 1)
      
      # Stack both input features:
      x = torch.stack([probe_values, probe_mask], dim = 2).reshape(batch_size*self.num_cells, 2)
      
      # Shared encoder (GNN message passing):
      for layer in self.encoder:
         x = layer(x, self.cell_adjacency)
         x = self.activation(x)
      
      # Aggregate full feature vectors per pin:
      if self.pin_power_head is not None and self.cell_to_pin_mapping is not None:
         x = x.view(batch_size, self.num_cells, -1)
         valid_mask = self.cell_to_pin_mapping >= 0
         pin_input = []
         for b in range(batch_size):
            pin_features = torch_scatter.scatter_mean(x[b, valid_mask], self.cell_to_pin_mapping[valid_mask], dim = 0, dim_size = self.num_pins)            
            pin_input.append(pin_features.flatten())
         pin_input = torch.stack(pin_input, dim = 0)
      else:
         pin_input = x.view(batch_size, -1)
      
      # Pin-power head:
      pin_power = self.pin_power_head(pin_input) if self.pin_power_head is not None else torch.empty(0)
      
      # Temperature head:
      if self.temperature_head is not None:
         x = x.view(batch_size*self.num_cells, -1)
         for i, layer in enumerate(self.temperature_head):
            x = layer(x, self.cell_adjacency)
            if i < len(self.temperature_head) - 1:
               x = self.activation(x)
         temperature = x.view(batch_size, self.num_cells)
      else:
         temperature = torch.empty(0)
      
      return pin_power, temperature

class HeatConductionCase:
   
   def __init__(self, loss_pin_power = 0.0, loss_temperature = 0.0, pin_power_ref = None, pin_power = None, temperature_ref = None, temperature = None):
      
      self.loss_pin_power = loss_pin_power
      self.loss_temperature = loss_temperature
      
      self.pin_power_ref = pin_power_ref if pin_power_ref is not None and pin_power_ref.numel() > 0 else None
      self.pin_power = pin_power if pin_power is not None and pin_power.numel() > 0 else None
      self.pin_power_error = None
      self.heat_source_ref = None
      self.heat_source = None
      self.heat_source_error = None
      
      self.temperature_ref = temperature_ref if temperature_ref is not None and temperature_ref.numel() > 0 else None
      self.temperature = temperature if temperature is not None and temperature.numel() > 0 else None
      self.temperature_error = None
   
   def process(self, pin_power_normalization, temperature_normalization, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats):
      
      if self.pin_power_ref is not None:
         self.pin_power_ref = pin_power_normalization.reverse(self.pin_power_ref[0].numpy())
         self.pin_power = pin_power_normalization.reverse(self.pin_power[0].numpy())
         self.pin_power_error = 100.0 * (self.pin_power-self.pin_power_ref) / np.mean(self.pin_power_ref)
         self.heat_source_ref = get_heat_source(self.pin_power_ref.reshape((nz, -1)), mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         self.heat_source = get_heat_source(self.pin_power.reshape((nz, -1)), mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         self.heat_source_error = get_heat_source(self.pin_power_error.reshape((nz, -1)), mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         self.num_cells = self.heat_source_ref.shape[0]
      
      if self.temperature_ref is not None:
         self.temperature_ref = temperature_normalization.reverse(self.temperature_ref[0].numpy())
         self.temperature = temperature_normalization.reverse(self.temperature[0].numpy())
         self.temperature_error = self.temperature - self.temperature_ref
         self.num_cells = self.temperature_ref.shape[0]
   
   def write(self, prefix, filename):
      
      shutil.copyfile(prefix + "/mesh.vtk", filename)
      
      field_data = [self.heat_source_ref, self.heat_source, self.heat_source_error, self.temperature_ref, self.temperature, self.temperature_error]
      field_names = ["pin-power-reference", "pin-power", "pin-power-error", "temperature-reference", "temperature", "temperature-error"]
      
      with open(filename, "a") as f:
         
         f.write("CELL_DATA %d\n\n" % self.num_cells)
         for data, name in zip(field_data, field_names):
            if not data is None:
               f.write("SCALARS %s double 1\n" % name)
               f.write("LOOKUP_TABLE default\n")
               for x in data:
                  f.write("%.9e\n" % x)
               f.write("\n")

def main():
   
   # Workflow parameters:
   run = False
   train = True
   test = True
   reproducibility = True
   reconstruct_pin_power = True
   reconstruct_temperature = True
   
   # Number of temperature measurements:
   num_probes = 50
   
   # Training parameters:
   num_samples = 10000
   num_epochs = 500
   batch_size = 100
   test_size = 0.1
   learning_rate = 0.001
   
   # Neural-network architecture (shared encoder + task-specific heads):
   # Some tips:
   #    - If performance plateaus early, try adding depth (more layers).
   #    - If training is unstable or slow, try reducing width.
   #    - Overfitting? Consider smaller network, dropout, or regularization.
   # Some possible configurations:
   #    - [2, 1]: simple, decent baseline.
   #    - [4, 2, 1]: wider, more expressive.
   #    - [4, 4, 2, 1]: deeper model.
   gnn = True
   if gnn:
      encoder_sizes = [8, 4]
      pin_power_head_sizes = [-2, -1]
      temperature_head_sizes = [2, 1]
      use_cell_to_pin_mapping = True
   else:
      encoder_sizes = [-2, -1]
      pin_power_head_sizes = [-1, -1]
      temperature_head_sizes = [-1, -1]
   
   # Power parameters:
   power_fraction_range = [0.9, 1.1]
   pin_power_range = [0.8, 1.2]
   
   # Mesh parameters:
   case = "single"
   two_dim = True
   quad = False
   write_vtk = False
   run_fltk = False
   
   # Mesh sizes:
   lc1 = 1.0
   lc2 = None
   lc3 = 3.0
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
   power_nominal = 15.0e6
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
   pin_power = power_nominal / num_pins[0]
   r = 0.5 * len(core) * pc
   r0 = 0.75 * pc
   
   # Small geometries:
   if case != "full":
      if case == "mini":
         core = [[0, 0, 6, 6, 0, 0], \
                   [6, 4, 2, 3, 6], \
                  [6, 2, 1, 1, 2, 6], \
                    [3, 1, 5, 1, 4], \
                   [6, 2, 1, 1, 2, 6], \
                     [6, 4, 2, 3, 6], \
                    [0, 0, 6, 6, 0, 0]]
         r = 0.5 * len(core) * pc
         r0 = 0.25 * pc
      elif case == "micro":
         core = [[0, 1, 1, 0], \
                   [1, 3, 1], \
                  [0, 1, 1, 0]]
         r = 0.6 * len(core) * pc
         r0 = None
      elif case == "single":
         core = [[1]]
         r = 0.7 * len(core) * pc
         r0 = None
      else:
         raise ValueError("Wrong case!")
      num_pins = get_num_pins(core, fas)
      power_nominal = pin_power * num_pins[0]
   prefix = "data-" + case
   prefix += "-2d" if two_dim else "-3d"
   
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
      power_nominal /= h
   
   # Build the mesh:
   fa_meshes = [build_x_hex_mesh(pf, fa) if fa is not None else None for fa in fas]
   core_mesh = build_x_hex_mesh(pc, core)
   mesh = build_gmsh_mesh(core_mesh, fa_meshes, d, r, r0, lc1, lc2, lc3, quad, write_vtk, run_fltk)
   
   # Generate the training data:
   if run:
      
      print("\nRun the training cases...\n")
      
      # Export the mesh:
      write_mesh("mesh.pmp", mesh, hb, h, ht, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
      
      # Create the data directories:
      for path in [prefix, prefix + "/temperature", prefix + "/pin-power"]:
         create_directory(path)
      
      # Run all training cases:
      for i in range(num_samples):
         
         # Export the main input file:
         power = np.random.uniform(*power_fraction_range) * power_nominal
         write_input("input.pmp", with_sdr, two_dim, power)
         
         # Calculate the pin power and export the relative volumetric heat-source distribution:
         pin_power = get_pin_power(pin_power_range, mesh, num_pins[0], nz)
         write_heat_source("data.pmp", pin_power, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
         
         # Run the heat-conduction solver:
         print("Run training case %d..." % i)
         subprocess.run(["./run.sh", "input.pmp"])
         
         # Get and export the pin power:
         pin_power *= power
         pin_power.tofile(prefix + "/pin-power/" + str(i) + ".np")
         
         # Get and export the temperature:
         temperature = parse_vtk_file("output_0.vtk", "temperature")
         temperature.tofile(prefix + "/temperature/" + str(i) + ".np")
      
      # Get the mesh in .vtk format:
      with open("output_0.vtk", "r") as fin, open(prefix + "/mesh.vtk", "w") as fout:
         for line in fin:
            if "CELL_DATA" in line: break
            fout.write(line)
      
      # Clean up solver files:
      for path in ["mesh.pmp", "data.pmp", "input.pmp", "output_0.vtk"]:
         remove_file(path)
      
      print("Done running the training cases.")
   
   # Prepare the training data:
   if train or test:
      
      print("\nLoad the training data...")
      
      # Set the random seed for reproducibility:
      if reproducibility:
         seed = 1
         np.random.seed(seed)
         torch.manual_seed(seed)
         torch.backends.cudnn.deterministic = True
         torch.backends.cudnn.benchmark = False
      else:
         seed = None
      
      # Get the mesh connectivity and the cell-to-pin mapping for GNNs:
      if gnn:
         cell_adjacency = get_cell_adjacency_graph(prefix + "/mesh.vtk")
         cell_to_pin_mapping = get_cell_to_pin_mapping(mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats) if use_cell_to_pin_mapping else None
      
      # Load the pin-power data:
      pin_power = None; pin_power_normalization = None; num_pins = 0
      if reconstruct_pin_power:
         pin_power, pin_power_normalization = load_training_data(prefix + "/pin-power", num_samples)
         num_pins = pin_power.shape[1]
      
      # Load the temperature data:
      temperature = None; temperature_normalization = None; num_cells = 0
      if reconstruct_temperature:
         temperature, temperature_normalization = load_training_data(prefix + "/temperature", num_samples)
         num_cells = temperature.shape[1]
      
      # Get the temperature at random probe locations:
      if reconstruct_temperature:
         probe_locations = torch.randint(0, num_cells, (num_probes,))
         probe_temperature = temperature[:, probe_locations]
      else:
         temperature, temperature_normalization = load_training_data(prefix + "/temperature", num_samples)
         num_cells = temperature.shape[1]
         probe_locations = torch.randint(0, num_cells, (num_probes,))
         probe_temperature = temperature[:, probe_locations]
         temperature = None
      
      # Wrap the training and testing data:
      train_loader, test_loader = wrap_training_data(num_samples, batch_size, test_size, seed, pin_power, temperature, probe_temperature)
      
      # Initialize the neural network:
      encoder_sizes = [int(-n*(num_pins+num_cells)) if n < 0 else n for n in encoder_sizes]
      pin_power_head_sizes = [int(-n*num_pins) if n < 0 else n for n in pin_power_head_sizes] if reconstruct_pin_power else None
      temperature_head_sizes = [int(-n*num_cells) if n < 0 else n for n in temperature_head_sizes] if reconstruct_temperature else None
      if gnn:
         model = GNNModel(num_pins, num_cells, probe_locations, cell_adjacency, cell_to_pin_mapping, encoder_sizes, pin_power_head_sizes, temperature_head_sizes)
      else:
         model = MLPModel(num_probes, num_pins, num_cells, encoder_sizes, pin_power_head_sizes, temperature_head_sizes)
      
      # Define optimizer and loss functions:
      optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate)
      criterion = torch.nn.MSELoss()
      
      print("Done loading the training data.")
   
   # Train the model:
   if train:
      
      print("\nTrain the neural network...\n")
      
      # Set training mode:
      model.train()
      
      # Loss weighting:
      w = 0.9 if reconstruct_pin_power and reconstruct_temperature else 1.0
      dynamic_weighting = True and reconstruct_pin_power and reconstruct_temperature
      
      # Run the training loop:
      loss_pin_power_total = np.zeros(num_epochs); loss_temperature_total = np.zeros(num_epochs); loss_total = np.zeros(num_epochs)
      for i in range(num_epochs):
         for pin_power_batch, temperature_batch, probe_temperature_batch in train_loader:
            
            # Forward pass:
            optimizer.zero_grad()
            pin_power, temperature = model(probe_temperature_batch)
            
            # Compute the pin-power and temperature MSE losses:
            loss_pin_power = criterion(pin_power, pin_power_batch) if reconstruct_pin_power else torch.tensor(0.0)
            loss_temperature = criterion(temperature, temperature_batch) if reconstruct_temperature else torch.tensor(0.0)
            
            # Compute the total loss:
            if not reconstruct_pin_power:
               loss = loss_temperature
            elif not reconstruct_temperature:
               loss = loss_pin_power
            else:
               if dynamic_weighting:
                  w = loss_pin_power.item() / (loss_pin_power.item()+loss_temperature.item())
               loss = w * loss_pin_power + (1.0-w) * loss_temperature
            
            # Backward pass and optimization:
            loss.backward()
            optimizer.step()
            
            # Accumulate the losses:
            if reconstruct_pin_power: loss_pin_power_total[i] += loss_pin_power.item()
            if reconstruct_temperature: loss_temperature_total[i] += loss_temperature.item()
            loss_total[i] += loss.item()
         
         # Print the loss:
         loss_pin_power_total[i] /= len(train_loader)
         loss_temperature_total[i] /= len(train_loader)
         loss_total[i] /= len(train_loader)
         if i % 10 == 9:
            print("Epoch %d/%d loss: %.3e (q), %.3e (T), %.3e (total), %.3e (w)" % (i+1, num_epochs, loss_pin_power_total[i], loss_temperature_total[i], loss_total[i], w))
      
      # Plot the loss:
      x = np.arange(1, num_epochs+1, dtype = int)
      if reconstruct_pin_power: plt.plot(x, loss_pin_power_total, label = "Pin power")
      if reconstruct_temperature: plt.plot(x, loss_temperature_total, label = "Temperature")
      if reconstruct_pin_power and reconstruct_temperature: plt.plot(x, loss_total, label = "Total")
      plt.xlabel("Epoch")
      plt.ylabel("Loss")
      plt.yscale("log")
      plt.legend()
      plt.grid(True)
      plt.savefig("loss.png")
      
      # Save the neural network:
      torch.save(model.state_dict(), "model.torch")
      
      print("\nDone training the neural network.")
   
   # Test the model:
   if test:
      
      print("\nTest the neural network...\n")
      
      # Load the neural network:
      model.load_state_dict(torch.load("model.torch", weights_only = True))
      
      # Set testing mode:
      model.eval()
      
      # Run the testing loop:
      loss_pin_power_total = 0.0; loss_temperature_total = 0.0; loss_total = 0.0
      worst_case_pin_power = HeatConductionCase(); worst_case_temperature = HeatConductionCase()
      with torch.no_grad():
         for pin_power_ref, temperature_ref, probe_temperature_ref in test_loader:
            
            # Get the predicted values:
            pin_power, temperature = model(probe_temperature_ref)
            
            # Compute the loss (MSE for both temperature and heat source):
            loss_pin_power = criterion(pin_power, pin_power_ref) if reconstruct_pin_power else torch.tensor(0.0)
            loss_temperature = criterion(temperature, temperature_ref) if reconstruct_temperature else torch.tensor(0.0)
            loss = loss_pin_power + loss_temperature
            
            # Compute the total loss:
            if reconstruct_pin_power: loss_pin_power_total += loss_pin_power.item()
            if reconstruct_temperature: loss_temperature_total += loss_temperature.item()
            loss_total += loss.item()
            
            # Get the worst pin-power prediction:
            if reconstruct_pin_power and loss_pin_power.item() > worst_case_pin_power.loss_pin_power:
               worst_case_pin_power = HeatConductionCase(loss_pin_power.item(), loss_temperature.item(), pin_power_ref, pin_power, temperature_ref, temperature)
            
            # Get the worst temperature prediction:
            if reconstruct_temperature and loss_temperature.item() > worst_case_temperature.loss_temperature:
               worst_case_temperature = HeatConductionCase(loss_pin_power.item(), loss_temperature.item(), pin_power_ref, pin_power, temperature_ref, temperature)
      
      # Print the mean loss:
      loss_pin_power_mean = loss_pin_power_total / len(test_loader)
      loss_temperature_mean = loss_temperature_total / len(test_loader)
      loss_mean = loss_total / len(test_loader)
      print("Mean loss: %.3e (q), %.3e (T), %.3e (total)" % (loss_pin_power_mean, loss_temperature_mean, loss_mean))
      
      # Denormalize and process the pin-power and temperature data:
      if reconstruct_pin_power: worst_case_pin_power.process(pin_power_normalization, temperature_normalization, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
      if reconstruct_temperature: worst_case_temperature.process(pin_power_normalization, temperature_normalization, mesh, nzb, nz, nzt, bottom_ref_mats, top_ref_mats)
      
      # Export the pin-power and temperature data for the worst cases:
      if reconstruct_pin_power: worst_case_pin_power.write(prefix, "worst-pin-power.vtk")
      if reconstruct_temperature: worst_case_temperature.write(prefix, "worst-temperature.vtk")
      
      print("\nDone evaluating the neural network.")

if __name__ == "__main__": main()
