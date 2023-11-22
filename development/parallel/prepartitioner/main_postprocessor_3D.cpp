#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "mpi.h"
#include "ParallelIO.h"
#include "SuperGroup.h"
#include "ParallelSuperGroup.h"
#include "CSP_VSet.h"
#include "GlobalTimeToInteger.h"


using namespace std;

namespace csp
 {
    csp::ErrorHandler skm_err;
    csp_float         global_time;
    void WriteNode2NodeFile( std::vector<std::pair<int,int> >& node_mapping, unsigned int my_rank );
    void ReadNode2NodeFile( std::vector<std::pair<int,int> >& node_mapping, unsigned int my_rank );
    void CalculateMassAndEnergyFluxes ( ParallelSuperGroup<csp_float,3>& sg );
 }

using namespace csp;

int main()
{
   /*
  Prepartitioner:

  Use this main file to partition a vset file containing a global
  2D triangular mesh into subdomains. 
  vset files, continaining the sub-domains will be stored as well
  as txt files containing the communication data between processors
  These files can be read in using memberfunctions in ParallelIO<fT,3>
  */

  cout<<"\n\nNOTE: THIS METHOD ONLY WORKS FOR NODAL SCALAR VARIABLES FOR NOW!"<< endl;

  MPI::Init();

  char fname[200];
  cout<<"\n\nGive global vset file"<< endl;
  cin >> fname;

  std::string iso;
  bool isoparametric(false);
  cout <<"\nIsoparametric or not (y/n) ?"<< endl;
  cin >> iso;
  if(iso == "y") isoparametric = true;
  else if (iso == "n") isoparametric = false;
  else 
   {
    cout <<"\nWRONG ENTRY...enter 'y' or 'n' please!!" << endl;
    return 0;
   }

  std::string n2n;
  bool do_n2n_mapping(false);
  cout <<"\nAre 'node2node_rank.txt' files available for all partitions (y/n) ? (If not a time-consuming node-to-node mapping has to be performed)"<< endl;
  cin >> n2n;
  if(n2n == "y") do_n2n_mapping = false;
  else if (n2n == "n") do_n2n_mapping = true;
  else 
   {
    cout <<"\nWRONG ENTRY...enter 'y' or 'n' please!!" << endl;
    return 0;
   }

  char phys_var_fname[200];
  cout<<"\n\nGive physical variable file name"<< endl;
  cin >> phys_var_fname;

/*
  int nov;
  cout<<"\n\nHow many variables should be post-processed"<< endl;
  cin >> nov;

  std::vector<char*> phys_vars;
  phys_vars.resize(nov);

  char phys_var[200];
  for(unsigned int i=0; i!=nov; i++)
   {
     cout<<"\n\nGive name variable "<<i+1<< endl;
     cin >> phys_var;
     phys_vars[i] = phys_var;
   }
*/
  std::vector<char*> phys_vars;
  phys_vars.resize(9);
  phys_vars[0]="temperature";    phys_vars[1]="nodal velocity";
  phys_vars[2]="mass flux x";    phys_vars[3]="mass flux y";    phys_vars[4]="mass flux z";
  phys_vars[5]="energy flux x";  phys_vars[6]="energy flux y";  phys_vars[7]="energy flux z";
  phys_vars[8]="fluid pressure";

  int n_procs;
  cout<<"\n\nGive number of partitions"<< endl;
  cin >> n_procs;

  VSet<csp_float,3> local_model, global_model;
  global_model.InputFrom( fname, global_time);
  stl_index         first_outerhalo;

  cout <<"\nPostprocessor 3D: building global SuperGroup"<< endl;
  SuperGroup<csp_float,3> global_sg ( global_model, phys_var_fname, isoparametric );

  // loop over number of processors
  for (unsigned int my_rank=0; my_rank!=n_procs; my_rank++ )
   {
     // read in local vset
     ParallelIO<csp_float,3>  VSetIO(n_procs, my_rank);
     VSetIO.ReadinLocalVSet(fname, local_model, first_outerhalo, global_time );

     cout <<"\nPostprocessor 3D: building local SuperGroup "<<my_rank<< endl;
     ParallelSuperGroup<csp_float,3>  local_sg ( local_model, phys_var_fname, isoparametric );
     std::vector<std::pair<int,int> > node_mapping(0);

     CalculateMassAndEnergyFluxes ( local_sg );

     if(do_n2n_mapping)
      {
        cout <<"\nPostprocessor 3D: creating local-to-global mapping for partition "<<my_rank<< endl;
        // loop over local sg nodes:
        csp_float max_noo = static_cast<csp_float>(local_sg.ReferenceMesh().Nodes() * global_sg.ReferenceMesh().Nodes());

        //cout <<local_sg.ReferenceMesh().Nodes()<< endl;
        //cout <<global_sg.ReferenceMesh().Nodes()<< endl;
        cout <<"\t\t\tMaximum number of operations to do "<<max_noo<< endl;

        if( max_noo < 0 ) max_noo = 1e10;

        int operation(0), operation2(0);
        for (deque<Node<csp_float,3> >::iterator lit = local_sg.ReferenceMesh().NodesBegin(); lit != local_sg.ReferenceMesh().NodesEnd(); lit++ )
         for (deque<Node<csp_float,3> >::iterator git = global_sg.ReferenceMesh().NodesBegin(); git != global_sg.ReferenceMesh().NodesEnd(); git++ )
          {
            operation++;
            operation2++;
            if( operation > static_cast<int> (max_noo/100.) )
             {
               cout <<".";
               operation = 0;
             }
            if( (*lit).x() == (*git).x() && (*lit).y() == (*git).y() && (*lit).z() == (*git).z() )
             {
               node_mapping.push_back( make_pair( (*lit).ID(),(*git).ID() ) );
               break;
             }
          }
        WriteNode2NodeFile( node_mapping, my_rank );
      }
     else // read in file 
      {
        ReadNode2NodeFile( node_mapping, my_rank );
      }

     assert(node_mapping.size()==local_sg.ReferenceMesh().Nodes());

     cout <<"\nPostprocessor 3D: Copying variable data from partition "<<my_rank<<" to the global domain"<< endl;
     // loop over node_mapping:
     ScalarVariable<csp_float> variable;
     VectorVariable<csp_float,3> vv;
     for (unsigned int i=0; i!=node_mapping.size();i++)
      {
        // loop over variables
        for(unsigned int j=0; j!=phys_vars.size(); j++)
         {
           csp::Index local_key  = local_sg.ReferencePropertyDatabase().StorageKey(phys_vars[j]);
           csp::Index global_key = global_sg.ReferencePropertyDatabase().StorageKey(phys_vars[j]);
if(j==1)
 {
           local_sg.ReferencePropertyStorage().Read(   node_mapping[i].first,  NODE, local_key.index,  vv );
//variable() = static_cast<csp_float>(my_rank);
           global_sg.ReferencePropertyStorage().Store( node_mapping[i].second, NODE, global_key.index, vv );
 }
else
 {
           local_sg.ReferencePropertyStorage().Read(   node_mapping[i].first,  NODE, local_key.index,  variable );
//variable() = static_cast<csp_float>(my_rank);
           global_sg.ReferencePropertyStorage().Store( node_mapping[i].second, NODE, global_key.index, variable );
 }
         }
      }

   }

  cout <<"\nPostprocessor 3D: writing output files"<< endl;
  VTK_Interface<csp_float,3> vtk;
  GlobalTimeToInteger<csp_float>    time_converter;
  for(unsigned int j=0; j!=phys_vars.size(); j++)
    vtk.OutputDataToVTK( global_sg, phys_vars[j], phys_vars[j], time_converter.GlobalTimeInYears(global_time));

  MPI::Finalize();

  return 0;


} // end main

namespace csp
 {
void WriteNode2NodeFile( std::vector<std::pair<int,int> >& node_mapping, unsigned int my_rank )
  {
     // write node2node mapping to file
     std::string  filename("node2node");
     char rank[200];
     sprintf(rank, "%ld", static_cast<long>(my_rank));
     filename += "_";
     filename += rank;
     filename += ".txt";
     std::ofstream ofs(filename.c_str());
     // write header
     ofs <<"# This file contains local node to global node mapping for partition "<<my_rank<< endl;
     ofs <<"# Total number of local nodes:";
     ofs <<"\n"<<node_mapping.size();
     // write local node IDs first
     ofs <<"\n# local nodes:";
     for (unsigned int i=0; i!=node_mapping.size(); i++)
       ofs <<"\n"<<node_mapping[i].first;
     // now write global node IDs
     ofs <<"\n# global nodes:";
     for (unsigned int i=0; i!=node_mapping.size(); i++)
       ofs <<"\n"<<node_mapping[i].second;
  }

void ReadNode2NodeFile( std::vector<std::pair<int,int> >& node_mapping, unsigned int my_rank )
  {
     std::string  filename("node2node");
     char rank[200];
     sprintf(rank, "%ld", static_cast<long>(my_rank));
     filename += "_";
     filename += rank;
     filename += ".txt";
     std::ifstream ifs;
     ifs.open(filename.c_str());
     if(!ifs.is_open())
       skm_err.notice(CSP_ERROR, "Postprocessor", "Text file seems to be missing", filename.c_str());
     char* pch;
     char text_line[256];
     int size;
     char* delim ="\n";

     // skip header
     ifs.getline(text_line, 256);// 1
     ifs.getline(text_line, 256);// 2
     // read number of nodes
     ifs.getline(text_line, 256);// 3
     pch = strtok(text_line, delim);
     node_mapping.resize(atoi(pch));
     // read local node IDs
     ifs.getline(text_line, 256);// 4
     for (unsigned int i=0; i!=node_mapping.size(); i++)
      {
        ifs.getline(text_line, 256);
        pch = strtok(text_line, delim);
        node_mapping[i].first = atoi(pch);
      }
     // read global node IDs
     ifs.getline(text_line, 256);
     for (unsigned int i=0; i!=node_mapping.size(); i++)
      {
        ifs.getline(text_line, 256);
        pch = strtok(text_line, delim);
        node_mapping[i].second = atoi(pch);
      }
  }
void CalculateMassAndEnergyFluxes ( ParallelSuperGroup<csp_float,3>& sg )
 {
   csp::Index v_key = sg.ReferencePropertyDatabase().StorageKey("nodal velocity");
   csp::Index h_key = sg.ReferencePropertyDatabase().StorageKey("enthalpy liquid");
   csp::Index r_key = sg.ReferencePropertyDatabase().StorageKey("density liquid");

   csp::Index M_key = sg.ReferencePropertyDatabase().StorageKey("mass flux");
   csp::Index E_key = sg.ReferencePropertyDatabase().StorageKey("energy flux");
   csp::Index Mx_key = sg.ReferencePropertyDatabase().StorageKey("mass flux x");
   csp::Index Ex_key = sg.ReferencePropertyDatabase().StorageKey("energy flux x");
   csp::Index My_key = sg.ReferencePropertyDatabase().StorageKey("mass flux y");
   csp::Index Ey_key = sg.ReferencePropertyDatabase().StorageKey("energy flux y");
   csp::Index Mz_key = sg.ReferencePropertyDatabase().StorageKey("mass flux z");
   csp::Index Ez_key = sg.ReferencePropertyDatabase().StorageKey("energy flux z");

   sg.ExtrapolateElementPropertyToNodeProperty( "velocity", "nodal velocity" );
   csp_float r, h;

   ScalarVariable<csp_float> mass_flux(PLAIN,0.), E_flux(PLAIN,0.);
   VectorVariable<csp_float,3> v;
   for (deque<Node<csp_float,3> >::iterator it = sg.ReferenceMesh().NodesBegin(); it != sg.ReferenceMesh().NodesEnd(); it++ )
    {
     // read variables
     r = (*it).Read( sg.ReferencePropertyStorage(), r_key.index );
     h = (*it).Read( sg.ReferencePropertyStorage(), h_key.index );
     (*it).Read( sg.ReferencePropertyStorage(), v_key.index, v );

     // total fluxes
     mass_flux = r * v.Length();
     E_flux    = r * h * v.Length();
     (*it).Store( sg.ReferencePropertyStorage(), M_key.index, mass_flux );
     (*it).Store( sg.ReferencePropertyStorage(), E_key.index, E_flux );

     // x flux
     mass_flux = r * v[0];
     E_flux    = r * h * v[0];
     (*it).Store( sg.ReferencePropertyStorage(), Mx_key.index, mass_flux );
     (*it).Store( sg.ReferencePropertyStorage(), Ex_key.index, E_flux );

     // y flux
     mass_flux = r * v[1];
     E_flux    = r * h * v[1];
     (*it).Store( sg.ReferencePropertyStorage(), My_key.index, mass_flux );
     (*it).Store( sg.ReferencePropertyStorage(), Ey_key.index, E_flux );

     // z flux
     if(3==3)
      {
       mass_flux = r * v[2];
       E_flux    = r * h * v[2];
       (*it).Store( sg.ReferencePropertyStorage(), Mz_key.index, mass_flux );
       (*it).Store( sg.ReferencePropertyStorage(), Ez_key.index, E_flux );
      }
    }

 } // end function SetBoundaryToMax

 } // csp