import shutil

def main():
   
   lengths = [20.0, 100.0, 15.0]
   materials = [2, 1, 3]
   dx = 0.1
   
   layout = []
   for (l, m) in zip(lengths, materials):
      n = round(l/dx)
      for i in range(n):
         layout.append(m)
   nx = len(layout)
   
   with open("reflected-diffusion/mesh.pmp", "w") as f:
      
      f.write("dx %d\n" % (-nx))
      f.write("%.3f\n\n" % dx)
      
      f.write("bc -x robin -0.4692\n")
      f.write("bc +x robin -0.4692\n")
      f.write("\n")
      
      f.write("materials %d\n" % nx)
      for i in range(nx):
         f.write("%d" % layout[i])
         if i < nx-1:
            f.write(" ")
         else:
            f.write("\n")
   
   with open("reflected-s2/mesh.pmp", "w") as f:
      
      f.write("dx %d\n" % (-nx))
      f.write("%.3f\n\n" % dx)
      
      f.write("bc -x vacuum\n")
      f.write("bc +x vacuum\n")
      f.write("\n")
      
      f.write("materials %d\n" % nx)
      for i in range(nx):
         f.write("%d" % layout[i])
         if i < nx-1:
            f.write(" ")
         else:
            f.write("\n")
   
   shutil.copyfile("reflected-s2/mesh.pmp", "reflected-s4/mesh.pmp")

if __name__ == "__main__": main()
