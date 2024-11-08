import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import gmsh
import subprocess
import os

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

def build_gmsh_mesh(r1, r2, lc):
   
   write_vtk = False
   run_fltk = False
   
   gmsh.initialize()
   
   gmsh.model.add("cylinder")
   
   gmsh.model.occ.addPoint(0.0, r1, 0.0, lc)
   gmsh.model.occ.addPoint(-r1, 0.0, 0.0, lc)
   gmsh.model.occ.addPoint(0.0, -r1, 0.0, lc)
   gmsh.model.occ.addPoint(r1, 0.0, 0.0, lc)
   
   inner_circle = gmsh.model.occ.addCircle(0.0, 0.0, 0.0, r1, angle1 = 0.0, angle2 = 2*np.pi)
   inner_curve = gmsh.model.occ.addCurveLoop([inner_circle])
   inner_surface = gmsh.model.occ.addPlaneSurface([inner_curve])
   
   gmsh.model.occ.addPoint(0.0, r2, 0.0, lc)
   gmsh.model.occ.addPoint(-r2, 0.0, 0.0, lc)
   gmsh.model.occ.addPoint(0.0, -r2, 0.0, lc)
   gmsh.model.occ.addPoint(r2, 0.0, 0.0, lc)
   
   outer_circle = gmsh.model.occ.addCircle(0.0, 0.0, 0.0, r2, angle1 = 0.0, angle2 = 2*np.pi)
   outer_curve = gmsh.model.occ.addCurveLoop([outer_circle])
   outer_surface = gmsh.model.occ.addPlaneSurface([outer_curve])
   
   gmsh.model.occ.cut([(2, outer_surface)], [(2, inner_surface)])
   
   gmsh.model.occ.synchronize()
   
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
   
   types, _, tags = gmsh.model.mesh.getElements()
   for t, pts in zip(types, tags):
      n = gmsh.model.mesh.getElementProperties(t)[3]
      if n == 3:
         cells = [tuple(pts[i:i+n]) for i in range(0, len(pts), n)]
         break
   cells = [[int(i-1) for i in c] for c in cells]
   
   boundaries = [[] for i in range(2)]
   for i, p in enumerate(points):
      r = np.sqrt(p[0]*p[0]+p[1]*p[1])
      if abs(r-r1) < 0.1*lc:
         boundaries[0].append(i)
      elif abs(r-r2) < 0.1*lc:
         boundaries[1].append(i)
   
   gmsh.finalize()
   
   return points, cells, boundaries

def write_input(points, cells, boundaries, nz, dz, k, rho, cp, P, T1, h1, T2, h2):
   
   with open("input.pmp", "w") as f:
      
      f.write("mesh unstructured mesh.pmp\n\n")
      
      f.write("material graphite {\n")
      f.write("   thermal-properties constant %.3e %.3e %.3e\n" % (k, rho, cp))
      f.write("   fuel 1\n")
      f.write("}\n\n")
      
      f.write("solver conduction {\n")
      f.write("   power %.3e\n" % P)
      f.write("}")
   
   with open("mesh.pmp", "w") as f:
      
      f.write("points %d\n" % len(points))
      for p in points:
         f.write("%.9e %.9e\n" % (p[0], p[1]))
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
      
      f.write("boundary inner %d\n" % len(boundaries[0]))
      for i in boundaries[0]:
         f.write("%d\n" % i)
      f.write("\n")
      
      f.write("boundary outer %d\n" % len(boundaries[1]))
      for i in boundaries[1]:
         f.write("%d\n" % i)
      f.write("\n")
      
      if nz > 1:
         f.write("dz -%d\n" % nz)
         f.write("%.9e\n\n" % dz)
      
      if T1 is None:
         f.write("bc inner adiabatic\n")
      else:
         f.write("bc inner convection %.3f %.3f\n" % (h1, T1))
      if h2 is None:
         f.write("bc outer dirichlet %.3f\n" % T2)
      else:
         f.write("bc outer convection %.3f %.3f\n" % (h2, T2))
      f.write("\n")
      
      if nz > 1:
         f.write("bc -z reflective\n")
         f.write("bc +z reflective\n")
         f.write("\n")
      
      n = len(cells) * nz
      f.write("materials %d\n" % n)
      for i in range(n):
         if i > 0: f.write(" ")
         f.write("1")

def get_analytical_solution(points, cells, nz, r1, r2, k, P, T1, h1, T2, h2):
   
   q = P / (np.pi*(r2*r2-r1*r1))
   c2 = -0.25*q / k
   if T1 is None and h2 is None:
      
      c0 = T2 + (0.5*q/k)*(0.5*r2*r2-r1*r1*np.log(r2))
      c1 = 0.5*r1*r1*q / k
      
   elif T1 is None and not h2 is None:
      
      c0 = T2 + (0.5*q/k)*(0.5*r2*r2-r1*r1*np.log(r2)) + (0.5*q/h2)*(r2+r1*r1/r2)
      c1 = 0.5*r1*r1*q / k
      
   elif not h1 is None and h2 is None:
      
      c1 = T2 - T1 + (0.25*q/k)*(r2*r2-r1*r1) + 0.5*r1*q/h1
      c1 /= np.log(r2/r1) + k/(r1*h1)
      c0 = T1 + 0.25*r1*r1*q/k - 0.5*r1*q/h1 + (k/(r1*h1)-np.log(r1))*c1
      
   elif not h1 is None and not h2 is None:
      
      c1 = T2 - T1 + (0.25*q/k)*(r2*r2-r1*r1) + 0.5*q*(r1/h1+r2/h2)
      c1 /= np.log(r2/r1) + k*(1.0/(r1*h1)+1.0/(r2*h2))
      c0 = T1 + 0.25*r1*r1*q/k - 0.5*r1*q/h1 + (k/(r1*h1)-np.log(r1))*c1
      
   else:
      raise Exception("Case not implemented!")
   
   nxy = len(cells)
   n = nxy * nz
   T = np.empty(n)
   for i, c in enumerate(cells):
      
      x0 = 0.0; y0 = 0.0
      for p in c:
         x0 += points[p][0]
         y0 += points[p][1]
      x0 /= len(c); y0 /= len(c)
      r = np.sqrt(x0*x0+y0*y0)
      
      for j in range(nz):
         T[j*nxy+i] = c2*r*r + c1*np.log(r) + c0
   
   return T

def main():
   
   # Problem parameters:
   r1 = 5.0
   r2 = 15.0
   k = 0.24
   rho = 0.00225
   cp = 707.7
   q = 1.0
   P = q * np.pi * (r2*r2-r1*r1)
   T1 = 850.0
   h1 = 0.05
   T2 = 950.0
   h2 = 0.075
   dT2 = 250.0
   T = [(None, T2+dT2), (None, T2), (T1, T2+dT2), (T1, T2)]
   h = [(None, None), (None, h2), (h1, None), (h1, h2)]
   case_labels = ["Adiabatic at r$\mathregular{_1}$ + Dirichlet at r$\mathregular{_2}$", 
                  "Adiabatic at r$\mathregular{_1}$ + convection at r$\mathregular{_2}$", 
                  "Convection at r$\mathregular{_1}$ + Dirichlet at r$\mathregular{_2}$", 
                  "Convection at r$\mathregular{_1}$ + convection at r$\mathregular{_2}$"]
   
   # Axial discretization:
   nz = 1
   dz = 1.0 / nz
   
   # Mesh size:
   l1 = 0.025
   l2 = 0.5
   dl = 0.025
   l = np.arange(l1, l2+0.1*dl, dl)
   
   dT = np.empty((2, len(T), len(l)))
   for j, lc in enumerate(l):
      
      points, cells, boundaries = build_gmsh_mesh(r1, r2, lc)
      
      for i, ((T1, T2), (h1, h2)) in enumerate(zip(T, h)):
         
         write_input(points, cells, boundaries, nz, dz, k, rho, cp, P, T1, h1, T2, h2)
         
         print("lc = ", lc)
         subprocess.run(["./run.sh", "input.pmp"])
         
         fields = parse_vtk_file("output_0.vtk")
         T = fields["temperature"]
         T0 = get_analytical_solution(points, cells, nz, r1, r2, k, P, T1, h1, T2, h2)
         error = np.abs(T-T0)
         
         dT[0][i][j] = np.max(error)
         dT[1][i][j] = np.linalg.norm(error) / np.sqrt(len(error))
         
         files = os.listdir(".")
         for f in files:
            if f.endswith(".pmp") or f.endswith(".vtk"):
               os.remove(f)
   
   matplotlib.rcParams.update({'font.size': 12})
   y_labels = ["Maximum error (K)", "L2-norm error (K)"]
   filenames = ["dTmax.png", "dTmean.png"]
   for i, (y_label, filename) in enumerate(zip(y_labels, filenames)):
      
      fig, ax = plt.subplots()
      
      for j, label in enumerate(case_labels):
         ax.plot(l, dT[i][j], label = label)
      
      ax.set_xlabel("Mesh size (cm)")
      ax.set_ylabel(y_label)
      
      ax.legend()
      plt.gca().invert_xaxis()
      plt.savefig(filename)
      plt.clf()

if __name__ == "__main__": main()
