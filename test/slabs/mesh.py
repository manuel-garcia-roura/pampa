def main():
   
   lengths = [20.0, 100.0, 15.0]
   materials = [2, 1, 3]
   dx = 1.0
   
   layout = []
   for (l, m) in zip(lengths, materials):
      n = round(l/dx)
      for i in range(n):
         layout.append(m)
   nx = len(layout)
   
   with open("reflected/mesh.pmp", "w") as f:
      
      f.write("# x-discretization:\n")
      f.write("dx %d\n" % (-nx))
      f.write("%.3f\n\n" % dx)
      
      f.write("# y-discretization:\n")
      f.write("dy 1\n")
      f.write("1.000\n\n")
      
      f.write("# z-discretization:\n")
      f.write("dz 1\n")
      f.write("1.000\n\n")
      
      f.write("# boundary conditions:\n")
      f.write("bc x 1 1\n")
      f.write("bc y 2 2\n")
      f.write("bc z 2 2\n")
      f.write("\n")
      
      f.write("# material distribution:\n")
      f.write("materials %d\n" % nx)
      for i in range(nx):
         f.write("%d" % layout[i])
         if i < nx-1:
            f.write(" ")
         else:
            f.write("\n")
      f.write("\n")

if __name__ == '__main__': main()
