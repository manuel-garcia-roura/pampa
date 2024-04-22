#include "Parser.hxx"

/* Read a plain-text input file: */
int Parser::read(const std::string& filename, Mesh** mesh, Array1D<Material>& materials, 
   Array1D<Solver*>& solvers, Array1D<double>& dt) {
   
   /* Open the input file: */
   std::ifstream file(filename);
   PAMPA_CHECK(!file.is_open(), 1, "unable to open " + filename);
   
   /* Read the file line by line: */
   while (true) {
      
      /* Get the next line: */
      std::vector<std::string> line = utils::get_next_line(file);
      if (line.empty()) break;
      
      /* Get the next keyword: */
      if (line[0] == "mesh") {
         
         /* Create the mesh depending on the mesh type: */
         std::string mesh_type = line[1];
         if (mesh_type == "cartesian")
            *mesh = new CartesianMesh();
         else if (mesh_type == "unstructured")
            *mesh = new UnstructuredExtrudedMesh();
         else if (mesh_type == "partitioned")
            *mesh = new PartitionedMesh();
         else {
            PAMPA_CHECK(true, 2, "wrong mesh type");
         }
         
         /* Read the mesh: */
         std::string mesh_filename = line[2];
         PAMPA_CALL((*mesh)->read(mesh_filename), "unable to read the mesh from " + mesh_filename);
         
         /* Build the mesh: */
         PAMPA_CALL((*mesh)->build(), "unable to build the mesh");
         
         /* Partition the mesh and swap the meshes: */
         if (mpi::size > 1 && !(*mesh)->isPartitioned()) {
            Mesh* submesh = NULL;
            PAMPA_CALL((*mesh)->partition(&submesh), "unable to partition the mesh");
            delete *mesh;
            *mesh = submesh;
         }
         
      }
      else if (line[0] == "material") {
         
         /* Create the material: */
         Material mat;
         
         /* Read the material: */
         std::string mat_filename = line[1];
         PAMPA_CALL(mat.read(mat_filename), "unable to read the material from " + mat_filename);
         
         /* Keep the material definition: */
         materials.pushBack(mat);
         
      }
      else if (line[0] == "solver") {
         
         /* Create the solver depending on the solver type: */
         std::string solver_type = line[1];
         if (solver_type == "diffusion") {
            
            /* Get the number of energy groups: */
            int num_energy_groups;
            PAMPA_CALL(utils::read(num_energy_groups, 1, INT_MAX, line[2]), 
               "wrong number of energy groups");
            
            /* Create the solver: */
            DiffusionSolver* solver = new DiffusionSolver(*mesh, materials, num_energy_groups);
            solvers.pushBack(solver);
            
         }
         else if (solver_type == "sn") {
            
            /* Get the order (N) of the SN method: */
            int order;
            PAMPA_CALL(utils::read(order, 1, INT_MAX, line[2]), "wrong SN order");
            
            /* Get the number of energy groups: */
            int num_energy_groups;
            PAMPA_CALL(utils::read(num_energy_groups, 1, INT_MAX, line[3]), 
               "wrong number of energy groups");
            
            /* Get the weight between upwind and linear schemes for face interpolation: */
            double face_interpolation_delta;
            PAMPA_CALL(utils::read(face_interpolation_delta, 0.0, 1.0, line[4]), 
               "wrong weight between upwind and linear interpolation");
            
            /* Get the switch to use the least-squares gradient for boundary interpolation: */
            bool boundary_interpolation_ls;
            PAMPA_CALL(utils::read(boundary_interpolation_ls, line[5]), 
               "wrong switch for least-squares boundary interpolation");
            
            /* Create the solver: */
            SNSolver* solver = new SNSolver(*mesh, materials, num_energy_groups, order, 
                             face_interpolation_delta, boundary_interpolation_ls);
            solvers.pushBack(solver);
            
         }
         else if (solver_type == "conduction") {
            
            /* Get the convergence tolerance and p-norm for nonlinear problems: */
            double tol = 1.0;
            if (line.size() > 2) {
               PAMPA_CALL(utils::read(tol, 0.0, DBL_MAX, line[2]), 
                  "wrong tolerance for implicit coupling");
            }
            double p = 2.0;
            if (line.size() > 3) {
               PAMPA_CALL(utils::read(p, 0.0, DBL_MAX, line[3]), 
                  "wrong p-norm for implicit coupling");
            }
            
            /* Create the solver: */
            HeatConductionSolver* solver = new HeatConductionSolver(*mesh, materials, tol, p);
            solvers.pushBack(solver);
            
         }
         else if (solver_type == "precursors") {
            
            /* Get the number of delayed-neutron precursor groups: */
            int num_precursor_groups;
            PAMPA_CALL(utils::read(num_precursor_groups, 1, INT_MAX, line[2]), 
               "wrong number of delayed-neutron precursor groups");
            
            /* Create the solver: */
            PrecursorSolver* solver = new PrecursorSolver(*mesh, materials, num_precursor_groups);
            solvers.pushBack(solver);
            
         }
         else if (solver_type == "coupled") {
            
            /* Get the solver name: */
            std::string name = line[2];
            
            /* Get the number of coupled solvers: */
            int num_coupled_solvers;
            PAMPA_CALL(utils::read(num_coupled_solvers, 1, 3, line[3]), 
               "wrong number of coupled solvers");
            
            /* Get the coupled solvers: */
            Array1D<Solver*> coupled_solvers(num_coupled_solvers, NULL);
            int l = 0;
            unsigned int i = 4;
            while (l < num_coupled_solvers) {
               PAMPA_CALL(utils::find(line[i++], solvers, &(coupled_solvers(l++))), 
                  "unable to find coupled solver");
            }
            
            /* Get the switch to use implicit coupling: */
            bool implicit = false;
            if (i < line.size()) {
               PAMPA_CALL(utils::read(implicit, line[i++]), "wrong switch for implicit coupling");
            }
            
            /* Get the convergence tolerance and p-norm for implicit coupling: */
            double tol = 1.0;
            if (i < line.size()) {
               PAMPA_CALL(utils::read(tol, 0.0, DBL_MAX, line[i++]), 
                  "wrong tolerance for implicit coupling");
            }
            double p = 2.0;
            if (i < line.size()) {
               PAMPA_CALL(utils::read(p, 0.0, DBL_MAX, line[i++]), 
                  "wrong p-norm for implicit coupling");
            }
            
            /* Create the solver: */
            CouplingSolver* solver = new CouplingSolver(name, *mesh, coupled_solvers, implicit, 
               tol, p);
            solvers.pushBack(solver);
            
         }
         else {
            PAMPA_CHECK(true, 3, "wrong solver type");
         }
         
      }
      else if (line[0] == "dt") {
         
         /* Get the dt values: */
         int nt;
         PAMPA_CALL(utils::read(nt, -INT_MAX, INT_MAX, line[1]), "wrong number of time steps");
         if (nt > 0) {
            PAMPA_CALL(utils::read(dt, nt, file), "wrong dt data");
         }
         else {
            PAMPA_CALL(utils::read(dt, 1, file), "wrong dt data");
            nt = -nt;
            dt.resize(nt, dt(0));
         }
         
      }
      else if (line[0] == "power") {
         
         /* Get the total power: */
         int np;
         PAMPA_CALL(utils::read(np, 1, INT_MAX, line[1]), "wrong number of power levels");
         Function power;
         if (np == 1) {
            double p;
            PAMPA_CALL(utils::read(p, 0.0, DBL_MAX, line[2]), "wrong power level");
            power = Function(p);
         }
         else {
            Array1D<double> t, p;
            PAMPA_CALL(utils::read(t, np, file), "wrong time data");
            PAMPA_CALL(utils::read(p, np, file), "wrong power data");
            power = Function(t, p);
         }
         
         /* Set the total power in the solvers: */
         for (int i = 0; i < solvers.size(); i++)
            solvers(i)->setPower(power);
         
      }
      else if (line[0] == "bc") {
         
         /* Get the solver name: */
         std::string name = line[1];
         
         /* Get the solver: */
         Solver* solver;
         PAMPA_CALL(utils::find(name, solvers, &solver), "unable to find solver");
         
         /* Get the boundary condition (1-based indexed): */
         int l, i = 2;
         BoundaryCondition bc;
         PAMPA_CALL(utils::read(l, 1, INT_MAX, line[i++]), "wrong boundary condition index");
         PAMPA_CALL(utils::read(bc, line, i), "wrong boundary condition");
         
         /* Add the boundary condition to the solver: */
         PAMPA_CALL(solver->addBoundaryCondition(bc, l), "unable to add the boundary condition");
         
      }
      else if (line[0] == "fixed") {
         
         /* Get the solver name: */
         std::string name = line[1];
         PAMPA_CHECK(name != "conduction", 4, "fixed values only allowed for heat conduction");
         
         /* Get the solver: */
         Solver* solver;
         PAMPA_CALL(utils::find(name, solvers, &solver), "unable to find solver");
         
         /* Get the material and the value: */
         int mat;
         PAMPA_CALL(utils::read(mat, 1, materials.size(), line[2]), "wrong material index");
         double x;
         PAMPA_CALL(utils::read(x, 0.0, DBL_MAX, line[3]), "wrong fixed value");
         
         /* Add the fixed temperature to the solver: */
         PAMPA_CALL(((HeatConductionSolver*)solver)->addFixedTemperature(mat-1, x), 
            "unable to add the fixed value");
         
      }
      else if (line[0] == "include") {
         
         /* Read an included input file: */
         std::string include_filename = line[1];
         PAMPA_CALL(read(include_filename, mesh, materials, solvers, dt), 
            "unable to parse " + include_filename);
         
      }
      else {
         
         /* Wrong keyword: */
         PAMPA_CHECK(true, 5, "unrecognized keyword '" + line[0] + "'");
         
      }
      
   }
   
   return 0;
   
}
