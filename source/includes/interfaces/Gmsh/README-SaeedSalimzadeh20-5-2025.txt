CSMP Gmsh Interface
===================
Saeed Salimzadeh (CSIRO, Clayton, Australia)

TODO: adjust: Gmsh_to_CSMP_elementType_enum to make it work with Open CSMP.
 
This is how I create and initialise the object: CSP_GMSH_Interface to create a mesh automatically from a 3dm input geometry file (Rhino3D file):
 
CSP_GMSH_Interface gmsh(“3dm filename”);
 
              gmsh.SetMeshSizeMax(10);
 
              gmsh.ConvertToQuadraticElements(False/True);
 
              gmsh.SetMinimumElementsPerTwoPi(60);
 
              gmsh.AddMeshPoints(points, points_mesh_sizes);
 
              gmsh.Run();
 
              gmsh.GmshGUI(True/False);
 
 
Here is a simple function to open a Gmsh mesh file (*.msh) directly and convert it to csp format using the OutputMesh function:
 
               
                void CSP_GMSH_Interface::Convert_Gmsh_Msh_to_Csmp_datasc(std::string file_name)
                {
                                // this is just a convenience function that clears everything in Gmsh then
                                // opens a Gmsh msh file then outputs it as csmp asc and dat files
 
                                gmsh::clear();
                                gmsh::open(file_name);
 
                                CSP_GMSH_Interface::OutputMesh(file_name);
                }
 
 
 
There are three ways to impose element sizes in Gmsh:
The maximum mesh size of the model
Number of elements per two PI in a curve
Each point can have a mesh size attached to it
 
You need to close the Gmsh GUI (if true) to continue. 
 
Cheers,
Saeed
 