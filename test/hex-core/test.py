import math as m
import numpy as np
import gmsh
import sys

gmsh.initialize()

gmsh.model.add("hex-core")

p = 1.0

lc = 0.1

dy = 0.5 * p
dx = (2.0/np.sqrt(3)) * 0.5 * p

gmsh.model.occ.addPoint(0.5*dx, -dy, 0.0, lc, 1)
gmsh.model.occ.addPoint(-0.5*dx, -dy, 0.0, lc, 2)
gmsh.model.occ.addPoint(-dx, 0.0, 0.0, lc, 3)
gmsh.model.occ.addPoint(-0.5*dx, dy, 0.0, lc, 4)
gmsh.model.occ.addPoint(0.5*dx, dy, 0.0, lc, 5)
gmsh.model.occ.addPoint(dx, 0.0, 0.0, lc, 6)

gmsh.model.occ.addLine(1, 2, 1)
gmsh.model.occ.addLine(2, 3, 2)
gmsh.model.occ.addLine(3, 4, 3)
gmsh.model.occ.addLine(4, 5, 4)
gmsh.model.occ.addLine(5, 6, 5)
gmsh.model.occ.addLine(6, 1, 6)

gmsh.model.occ.addCircle(0.0, 0.0, 0.0, 0.25, 7, angle1 = 0.0, angle2 = 2*m.pi)
gmsh.model.occ.addCurveLoop([7], 2)
gmsh.model.occ.addPlaneSurface([2], 2)

gmsh.model.occ.addCurveLoop([1, 2, 3, 4, 5, 6], 1)
gmsh.model.occ.addPlaneSurface([1, -2], 1)

gmsh.model.occ.synchronize()

# At this level, Gmsh knows everything to display the rectangular surface 1 and
# to mesh it. An optional step is needed if we want to group elementary
# geometrical entities into more meaningful groups, e.g. to define some
# mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
# material ("steel", "carbon") properties.
#
# Such groups are called "Physical Groups" in Gmsh. By default, if physical
# groups are defined, Gmsh will export in output files only mesh elements that
# belong to at least one physical group. (To force Gmsh to save all elements,
# whether they belong to physical groups or not, set the `Mesh.SaveAll' option
# to 1.) Physical groups are also identified by tags, i.e. stricly positive
# integers, that should be unique per dimension (0D, 1D, 2D or 3D). Physical
# groups can also be given names.
#
# Here we define a physical curve that groups the left, bottom and right curves
# in a single group (with prescribed tag 5); and a physical surface with name
# "My surface" (with an automatic tag) containing the geometrical surface 1:
gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4, 5, 6], 1, name = "out")
# gmsh.model.addPhysicalGroup(2, [2], 2, name = "hp")

# We can then generate a 2D mesh...
gmsh.model.mesh.generate(2)

points = gmsh.model.mesh.getNodes()
point_indices = points[0]
point_coordinates = np.reshape(points[1], (-1, 3))

cells = gmsh.model.mesh.getElements()
print(cells)

# ... and save it to disk
# gmsh.write("mesh.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
