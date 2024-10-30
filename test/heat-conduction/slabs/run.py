import numpy as np
import matplotlib.pyplot as plt
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

def write_input(nx, dx, k, rho, cp, P, T1, h1, T2, h2):
   
   with open("input.pmp", "w") as f:
      
      f.write("mesh cartesian mesh.pmp\n\n")
      
      f.write("material graphite {\n")
      f.write("   thermal-properties constant %.3e %.3e %.3e\n" % (k, rho, cp))
      f.write("   fuel 1\n")
      f.write("}\n\n")
      
      f.write("solver conduction {\n")
      f.write("   power %.3e\n" % P)
      f.write("}")
   
   with open("mesh.pmp", "w") as f:
      
      f.write("dx -%d\n" % nx)
      f.write("%.9e\n\n" % dx)
      
      if T1 is None:
         f.write("bc -x adiabatic\n")
      elif h1 is None:
         f.write("bc -x dirichlet %.3f\n" % T1)
      else:
         f.write("bc -x convection %.3f %.3f\n" % (h1, T1))
      if h2 is None:
         f.write("bc +x dirichlet %.3f\n" % T2)
      else:
         f.write("bc +x convection %.3f %.3f\n" % (h2, T2))
      f.write("\n")
      
      f.write("materials %d\n" % nx)
      for i in range(nx):
         if i > 0: f.write(" ")
         f.write("1")

def get_analytical_solution(nx, dx, L, k, P, T1, h1, T2, h2):
   
   if T1 is None and h2 is None:
      c0 = T2 + 0.5*L*P/k
      c1 = 0.0
   elif T1 is None and not h2 is None:
      c0 = T2 + (0.5*L/k+1.0/h2)*P
      c1 = 0.0
   elif h1 is None and h2 is None:
      c0 = T1
      c1 = 0.5*P/k + (T2-T1)/L
   elif h1 is None and not h2 is None:
      c0 = T1
      c1 = ((1.0+0.5*L*h2/k)*P+h2*(T2-T1)) / (k+L*h2)
   elif not h1 is None and h2 is None:
      c0 = (0.5*L*P/k+L*T1*h1/k+T2) / (1.0+h1*L/k)
      c1 = (h1/k) * (c0-T1)
   elif not h1 is None and not h2 is None:
      c0 = ((0.5*L/k+1.0/h2)*P+(h1/h2+L*h1/k)*T1+T2) / (1.0+h1/h2+L*h1/k)
      c1 = (h1/k) * (c0-T1)
   else:
      raise Exception("Case not implemented!")
   c2 = -0.5*P / (k*L)
   
   T = np.empty(nx)
   for i in range(nx):
      x = (i+0.5) * dx
      T[i] = c2*x*x + c1*x + c0
   
   return T

def main():
   
   L = 10.0
   k = 0.24
   rho = 0.00225
   cp = 707.7
   P = 100.0
   T1 = 600.0
   h1 = 0.1
   T2 = 300.0
   h2 = 0.5
   T = [(None, T2), (None, T2), (T1, T2), (T1, T2)]
   h = [(None, None), (None, h2), (h1, None), (h1, h2)]
   labels = ["adiabatic + Dirichlet", "adiabatic + convection", "convection + Dirichlet", "convection + convection"]
   
   n1 = 500
   n2 = 10000
   dn = 500
   n = np.arange(n1, n2+1, dn, dtype = int)
   d = L / n
   
   dT = np.empty((2, len(T), len(n)))
   for i, ((T1, T2), (h1, h2)) in enumerate(zip(T, h)):
      
      for j, (nx, dx) in enumerate(zip(n, d)):
         
         write_input(nx, dx, k, rho, cp, P, T1, h1, T2, h2)
         subprocess.run(["./run.sh", "input.pmp"])
         
         fields = parse_vtk_file("output_0.vtk")
         temperature = fields["temperature"]
         temperature_analytical = get_analytical_solution(nx, dx, L, k, P, T1, h1, T2, h2)
         error = np.abs(temperature-temperature_analytical)
         
         dT[0][i][j] = np.max(error)
         dT[1][i][j] = np.linalg.norm(error) / nx
         
         files = os.listdir(".")
         for f in files:
            if f.endswith(".pmp") or f.endswith(".vtk"):
               os.remove(f)
   
   filenames = ["dTmax.png", "dTmean.png"]
   for i, filename in enumerate(filenames):
      
      fig, ax = plt.subplots()
      
      for j, label in enumerate(labels):
         ax.plot(d, dT[i][j], label = label)
      
      ax.set_xlabel("Mesh size (cm)")
      ax.set_ylabel("Error (K)")
      
      ax.legend()
      plt.gca().invert_xaxis()
      plt.savefig(filename)
      plt.clf()

if __name__ == "__main__": main()
