#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "mpi.h"
#include "ParallelSuperGroup.h"
#include "DataInterpreter.h"
#include "Triangle.h"
#include "CSP_String.h"
#include "CSP_TextInterface.h"
#include "CSP_VSet.h"
#include "CSP_String.h"
#include "FemToGridVisitor.h"
#include "FemFromGridVisitor.h"
#include "CSP_Standard_IO_Handler.h"

using namespace std;

namespace csp {

extern ErrorHandler  skm_err;



/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  ParallelSuperGroup<fT,dim>::ParallelSuperGroup( VSet& vset, 
                                                  const char* var_file, 
                                                  bool isoparametric_local_mesh )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

Builds 2 and 3D, finite-element local_mesh according to the specifications
given in the supplied VSet. The vset may also contain 
elements of different types, for instance, tetrahedra and prism elements.If,
these are not already stored in the VSet, different material properties 
can be assigned to regions inside the the calculation. <p>

<H4>Input Arguments:</H4><!----------------------------------------------->

A VSet is supplied as constructor argument and should contain a
finite-element local_mesh and associated properties with names which must correspond 
to the properties specified in the Property input file 
(default: CSP_variables.txt).
A bool isoparametric_local_mesh can be supplied optionally (defaults to false) if the
the local_mesh is isoparamteric
A bool flag_irregular_top_boundary can be supplied optionally (defaults to true)
if the top boundary need to be flagged. This switch should be set to false for 
parallel computations.
<p>

<H4>Implementation:</H4><!------------------------------------------------>

The constructor first initializes the MeshManager and then the 
MemoryManager with the data from the VSet. Then property data 
stored in the VSet as FEM_Data are used to initialize some 
of the properties for which memory was allocated for by the 
MemoryManager.<p> 

<H4>Application:</H4><!--------------------------------------------------->

Constructor will build a ParallelSuperGroup for 2D local_meshes from 'Triangle' as well 
as 3D local_meshes from GoCad or NASTRAN input data.<p>

<H4>Messages:</H4><!------------------------------------------------------>

If no properties are supplied, the constructor will issue a Warning message.
In this case, property values have to be assigned before a computation is
carried out since all basic variables are initialized to NAN (not a number)
by default. <p>

<!------------------------------------------------------------------------>
tested: O.K. */
template<typename fT,stl_index dim>
ParallelSuperGroup<fT,dim>::ParallelSuperGroup( VSet<fT,dim>& vset, const char* var_file, bool isoparametric )
  : SuperGroup<fT,dim>( vset, var_file, isoparametric, false ),
    fem_manager(dim,vset.OrderOfFiniteElementInterpolationFunctions(),isoparametric),
    phys_vars(var_file),
    n_procs(MPI::COMM_WORLD.Get_size()), my_rank(MPI::COMM_WORLD.Get_rank())
 {
    vset.CheckFixPData();
    cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::ParallelSuperGroup(VSet): ";
    cout <<"The mesh is treated as adaptively refined and Y points upward." << endl; 

    // 1. building the finite element local_mesh
    local_mesh.BuildConnectivityFrom( fem_manager, vset );
    if ( phys_vars.VariablesFace() != 0U ) local_mesh.AddFaces( fem_manager );
    if ( phys_vars.VariablesConstraintPoint() != 0U ) local_mesh.AddConstraintPoints();

    // 2. building the property storage
    property_collection.SetStorageSpecifications( phys_vars, local_mesh.Nodes(), local_mesh.ConstraintPoints(),
                                                  local_mesh.Faces(), local_mesh.Elements() ); 
    property_collection.BuildPropertyStorage();

    unsigned int boxed_shaped_model(0);
    int rank(100000),send_rank(100000);

    if (local_mesh.BoxShapedModel())
     {
       boxed_shaped_model=1;
       rank = my_rank;
     }
    MPI::COMM_WORLD.Allreduce(&rank,&send_rank,1,MPI::INT,MPI::MIN);
    if(my_rank==send_rank)
      cout <<"\nParallelSuperGroup: Rank "<<my_rank<<": Global mesh is box-shaped: Broadcasting to all processors." << endl;
    if(send_rank!=100000)
      MPI::COMM_WORLD.Bcast(&boxed_shaped_model,1,MPI::UNSIGNED,send_rank);

    if ( boxed_shaped_model==1 )
     {
        cout <<"\nParallelSuperGroup: Rank "<<my_rank<<": Global mesh is box-shaped: Flag corners" << endl;
        FlagCornerElements();
        //if ( dim == 2U && flag_irregular_top_boundary ) FlagIrregularTopBoundary();
     }
    else
        cout <<"\nParallelSuperGroup: Rank "<<my_rank<<": Global mesh is NOT box-shaped: Do not flag corners" << endl;

    // test:
    FlagCornerElements();


    csp::Index       prop_key = phys_vars.StorageKey("permeability");
    ScalarVariable<fT>  perm;

    // 3. assigning properties to local_mesh
    for ( typename map<string,FEM_Data<ScalarVariable<fT> > >::const_iterator
          sit=vset.ScalarPropertiesBegin(); sit!=vset.ScalarPropertiesEnd(); sit++ )
       SuperGroup<fT,dim>::InputVariableFrom( (*sit).first.c_str(), (*sit).second );
    for ( typename map<string,FEM_Data<VectorVariable<fT,dim> > >::const_iterator
          vit=vset.VectorPropertiesBegin(); vit!=vset.VectorPropertiesEnd(); vit++ )
       SuperGroup<fT,dim>::InputVariableFrom( (*vit).first.c_str(), (*vit).second );
    for ( typename map<string,FEM_Data<TensorVariable<fT,dim> > >::const_iterator
          tit=vset.TensorPropertiesBegin(); tit!=vset.TensorPropertiesEnd(); tit++ )
       SuperGroup<fT,dim>::InputVariableFrom( (*tit).first.c_str(), (*tit).second );

    // 4. if no properties are contained in the VSet
    if ( vset.ScalarPropertiesBegin() == vset.ScalarPropertiesEnd() &&
         vset.VectorPropertiesBegin() == vset.VectorPropertiesEnd() &&
         vset.TensorPropertiesBegin() == vset.TensorPropertiesEnd() )
      skm_err.notice( WARNING, "ParallelSuperGroup<csp_float,dim>::ParallelSuperGroup (VSet)", 
                               "Warning no property data were supplied...");
    else
    cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::ParallelSuperGroup (VSet): Mesh has been built successfully..." << endl; 
    
    // local_mesh.Out();
    // property_collection.Out();

 } // end VSet constructor
 



template<typename fT,stl_index dim>
ParallelSuperGroup<fT,dim>::~ParallelSuperGroup()
  {
  } // end 




/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  ParallelSuperGroup<fT,dim>& ParallelSuperGroup::operator=( const ParallelSuperGroup<fT,dim>& sg )
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

The copy-constructor permits to duplicate models which may come in handy
if one wants to compare the results of slightly different computations
at runtime. This has, for instance, the advantage that large datasets 
must not be written to file first before they can be compared.<p>

<!------------------------------------------------------------------------>
tested: fails ? */
template<typename fT,stl_index dim>
ParallelSuperGroup<fT,dim>& ParallelSuperGroup<fT,dim>::operator=( const ParallelSuperGroup<fT,dim>& sg )
 {
    if ( &sg == this ) return *this;
    phys_vars            = sg.phys_vars;
    fem_manager          = sg.fem_manager;
    local_mesh                 = sg.local_mesh;
    property_collection  = sg.property_collection;
    group_list           = sg.group_list;
    
    return *this;
 } // end assignment operator


/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  void ParallelSuperGroup<fT,dim>::FlagCornerElements()
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

In CSP_definitions.h  boundary flags are defined which comprise labels 
for the 4/8 or more corner elements/ and nodes of each triangle / tetrahedra. 
These flags are used to assign boundary conditions etc. The numbering of the 
corners is counter-clockwise starting from the lower left:
<pre>
     
          always                    (in 3D)
  4                  3        8                  7 
   o----------------o          o----------------o
   |                |          |                |
   |     BACK       |          |      FRONT     | 
   |                |          |                |
   o----------------o          o----------------o
  1                  2        5                  6
 
 CNR1 to CNR8
</pre> 
 
 In the sketch above, both, the front and the back are seen from the 
 front. The model corners are always numbered in the displayed
 fashion. <p>

<H4>Implementation:</H4><!------------------------------------------------>

After having determined the dimensions of the model, the method examines 
all the nodes which lie at the model boundary and finds 
those nodes that lie at the spatial extrema and labels them according
to the scheme shown above. Subsequently, the Elements to which the
corner nodes belong are flagged equivalent to the former. <p>

<H4>Application:</H4><!--------------------------------------------------->

Obviously, this method makes sense only for rectangular and/or box-shaped 
models. In irregular shaped models, no corner points will be found and
no assigments will be made. Thus FlagCornerElements() will issue a 
warning. <p>
        
Since constraint points lie along Segments or Faces they never are corner
points. Hence, this method  does not affect their boundary flagging.<p>

<H4>Messages:</H4><!------------------------------------------------------>

The method distinguishes between two- and three-dimensional models and
will report if the number of identified corner nodes is too small.<p>

<!------------------------------------------------------------------------>
tested: O.K. */
template<typename fT,stl_index dim>
void ParallelSuperGroup<fT,dim>::FlagCornerElements()
 {
     extern ErrorHandler  skm_err;
 
//     if ( !local_mesh.BoxShapedModel() )
//       skm_err.notice( FATAL_ERROR, "ParallelSuperGroup<csp_float,dim>::FlagCornerElements", 
//                                    "Corners can only be identified for box-shaped models.");

     typename deque<Node<fT,dim> >::iterator     nit;
     stl_index                          i;
     int counter=0;
     
     // 0. extrema of local_mesh
     // ------------------
     vector<fT> xyz, xyz_loc;
     SuperGroup<fT,dim>::Dimensions( xyz );

     // if parallel dimensions need to be communicated!
     unsigned int entry(0);
     SuperGroup<fT,dim>::Dimensions( xyz_loc );
     xyz.resize(xyz_loc.size()); 
     for (unsigned int i=0;i!=(xyz_loc.size()/2);i++)
      {
        entry = i*2;
        MPI::COMM_WORLD.Allreduce(&xyz_loc[entry],&xyz[entry],1,MPI::DOUBLE,MPI::MIN);
        entry = i*2 + 1;
        MPI::COMM_WORLD.Allreduce(&xyz_loc[entry],&xyz[entry],1,MPI::DOUBLE,MPI::MAX);
      }
     cout <<"\nParallelSuperGroup<fT,dim>::FlagCornerElements: Rank "<<my_rank<<":Global local_mesh dimensions are: "<< endl <<"\t"<<xyz[0]<<" < x < "<<xyz[1]<< endl <<"\t"<< xyz[2]<<" < y < "<<xyz[3]<< endl;
     if(dim==3) cout <<"\t"<< xyz[4]<<" < z < "<<xyz[5]<< endl;

     // 1. Finding the corner nodes loop for 2D and 3D models
     // -----------------------------------------------------
     if ( dim == 2 ) for ( nit=local_mesh.NodesBegin(); nit!=local_mesh.NodesEnd(); nit++ )
       {
          // if the node lies at the model boundary
          if ( (*nit).AtBoundary() != NOT ) 
            {
               // CNR_MIN in the origin of a 2D model
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[2] )
                 {
                    (*nit).AtBoundary( CNR1 );
                    counter++;
                 } 
               // CNR_MIN_MAXX
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[2] )
                 {
                    (*nit).AtBoundary( CNR2 );
                    counter++;
                 } 
               // CNR_MIN_MAXXZ
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[3] )
                 {
                    (*nit).AtBoundary( CNR3 );
                    counter++;
                 } 
               // CNR_MIN_MAXZ
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[3] )
                 {
                    (*nit).AtBoundary( CNR4 );
                    counter++;
                 } 
            }
        } // end node loop

     else if ( dim == 3 ) for ( nit=local_mesh.NodesBegin(); nit!=local_mesh.NodesEnd(); nit++ )
       {
          if ( (*nit).AtBoundary() != NOT ) // if the node lies at the model boundary
            {
               // BACK LOWER LEFT defined as CNR_MIN
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[2] && (*nit).z() == xyz[4] )
                 {
                    (*nit).AtBoundary( CNR1 );
                    counter++;
                 } 
               // BACK LOWER RIGTH defined as CNR_MIN_MAXX
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[2] && (*nit).z() == xyz[4] )
                 {
                    (*nit).AtBoundary( CNR2 );
                    counter++;
                 } 
               // BACK UPPER RIGHT defined as CNR_MAX_MAXX
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[3] && (*nit).z() == xyz[4] )
                 {
                    (*nit).AtBoundary( CNR3 );
                    counter++;
                 } 
               // BACK UPPER LEFT defined as CNR_MAX_MINXZ 
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[3] && (*nit).z() == xyz[4] )
                 {
                    (*nit).AtBoundary( CNR4 );
                    counter++;
                 } 
               // FRONT BOTTOM LEFT defined as CNR_MIN_MAXZ
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[2] && (*nit).z() == xyz[5] )
                 {
                    (*nit).AtBoundary( CNR5 );
                    counter++;
                 } 
               // FRONT BOTTOM RIGHT defined as CNR_MIN_MAXXZ
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[2] && (*nit).z() == xyz[5] )
                 {
                    (*nit).AtBoundary( CNR6 );
                    counter++;
                 } 
               // FRONT TOP RIGHT defined as CNR_MAX
               if ( (*nit).x() == xyz[1] && (*nit).y() == xyz[3] && (*nit).z() == xyz[5] )
                 {
                    (*nit).AtBoundary( CNR7 );
                    counter++;
                 } 
               // FRONT TOP LEFT defined as CNR_MAX_MAXZ
               if ( (*nit).x() == xyz[0] && (*nit).y() == xyz[3] && (*nit).z() == xyz[5] )
                 {
                    (*nit).AtBoundary( CNR8 );
                    counter++;
                 } 
            }
        } // end node loop

    // Summing all corners flagged
    int sum(0);
    MPI::COMM_WORLD.Allreduce(&counter,&sum,1,MPI::INT,MPI::SUM);

    // 2. testing whether method identified all corners
    // ------------------------------------------------   
    if ( dim == 2 ) 
      {
         if ( sum < 4 )
           {
               cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagCornerElements: Warning: Only "<< counter <<" corner points ";
               cout <<"were identified (2D model should have 4)."<< endl;
           }
         else
           cout <<"\nParallelSuperGroup::FlagCornerElements: Rank "<<my_rank<<": Flagged "<<counter<<" corners (total corners = "<<sum<<")"<< endl;
      }  
    else
      {
         if ( sum < 8 )
           {
               cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagCornerElements: Warning: Only "<< counter <<" corner points ";
               cout <<"were identified (3D model should have 8)."<< endl;
           }
         else
           cout <<"\nParallelSuperGroup::FlagCornerElements: Rank "<<my_rank<<": Flagged "<<counter<<" corners (total corners = "<<sum<<")"<< endl;
      }
        
    // 3. flagging elements like their nodes  
    // -------------------------------------
    typename deque<Element<fT,dim> >::iterator  it;
    Node<fT,dim>*                      nd;
    SG_BOUNDARY                        b;
    
    for ( counter=0, it=local_mesh.ElementsBegin(); it!=local_mesh.ElementsEnd(); it++ )    
      for ( i=0; i<(*it).Nodes(); i++ )
        {
           nd = (*it).ConnectedNode(i);
           b  = nd->AtBoundary();
           if ( b == CNR1 || b == CNR2 || b == CNR3 || b == CNR4 ||
                b == CNR5 || b == CNR6 || b == CNR7 || b == CNR8  )
             { 
                (*it).AtBoundary( b );
                counter++;
             }
        }
    cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagCornerElements: "<< counter <<" assignments made."<< endl;  
 } // end FlagCornerElements
 



/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  bool ParallelSuperGroup<fT,dim>::FlagIrregularTopBoundary()
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

For topography-driven flow problems the upper model boundary may be 
irregular with slopes exceeding 45 degrees. In this case lateral straight
model boundaries will still have been identified correctly but elements 
at the top boundary of the model may have been flagged as lateral 
boundaries or as irregular. These elements must be correctly identified
and re-assigned as top-boundary elements. <p>

<H4>Output Arguments &amp; Return Value</H4><!---------------------------->

The method returns the boolean variable 'true' if no assignments had to
be made or if all assignments could be made without difficulties.<p>

<H4>Application:</H4><!--------------------------------------------------->

To model topography-driven flow. <p>

<H4>Messages:</H4><!------------------------------------------------------>

If the method cannot be applied if will issue an INFO message. <p>

<!------------------------------------------------------------------------>
tested:  */
template<typename fT,stl_index dim>
bool ParallelSuperGroup<fT,dim>::FlagIrregularTopBoundary( fT right_left_tolerance )
 {
     // 0. before doing anything check whether boundary elements have been
     //    flagged irregular if not nothing has to be done
     list<stl_index>                  belmts;
     vector<fT>                       minmax_xy;
     typename list<stl_index>::const_iterator  lit;
     bool                             at_top, at_bottom;
     stl_index                        reflagged_e(0);

     // getting the model dimensions for coordinate comparison
     SuperGroup<fT,dim>::Dimensions( minmax_xy );

     if ( !EnlistBoundaryElements( IRREGULAR, belmts ) ) {
          cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagIrregularTopBoundary: "; 
          cout <<"No elements flagged IRREGULAR were found."<< endl;
          return false;
       }
     else
       {  
          // 1. if elements flagged as irregular were found, these must be
          //    tested for which boundary they belong to. If this is not the top
          //    boundary or some internal boundary, an error must be issued 
          //    eventually, but first the top ones are identified by their node
          //    coordinates.
          // Note: If there is topography, only a few of the top facing elements 
          //    are likely to have been identified correctly the others however will
          //    be flagged IRREGULAR 
     
          // looping over the elements, finding their boundary nodes
          for ( lit=belmts.begin(); lit!=belmts.end(); lit++ )
            {
               at_top = false;
       
               for ( stl_index i=0; i<local_mesh.E( (*lit)-1 ).Nodes(); i++ )
                 if ( local_mesh.E( (*lit)-1 ).ConnectedNode(i)->AtBoundary() != NOT )
                   {
                      if ( local_mesh.E( (*lit)-1 ).AtBoundary() == BOTTOM ||
                           local_mesh.E( (*lit)-1 ).AtBoundary() == CNR1   ||
                           local_mesh.E( (*lit)-1 ).AtBoundary() == CNR2   ||
                           local_mesh.E( (*lit)-1 ).ConnectedNode(i)->y() <= (minmax_xy[2]+right_left_tolerance) )
                        {
                           at_bottom = true;
                           break;
                        }
                      else
                       {
                          // test: xmin == x or x == xmax
                          // if yes, Node is a corner node in a model with topography and should
                          // be flagged as such
                          if      ( local_mesh.E( (*lit)-1 ).ConnectedNode(i)->x() <= (minmax_xy[0]+right_left_tolerance) )
                            {
                               local_mesh.E( (*lit)-1 ).ConnectedNode(i)->AtBoundary( CNR4 );
                               local_mesh.E( (*lit)-1 ).AtBoundary( CNR4 );
                            }
                          else if ( local_mesh.E( (*lit)-1 ).ConnectedNode(i)->x() >= (minmax_xy[1]-right_left_tolerance) )
                            { 
                               local_mesh.E( (*lit)-1 ).ConnectedNode(i)->AtBoundary( CNR3 );
                               local_mesh.E( (*lit)-1 ).AtBoundary( CNR3 );
                            }
                          else
                            {
                               // now it is plausible that the irregular boundary node lies at 
                               // the model top
                               local_mesh.E( (*lit)-1 ).ConnectedNode(i)->AtBoundary( TOP );
                               at_top = true;
                            }
                       }
                   }

               if ( at_top ) local_mesh.E( (*lit)-1 ).AtBoundary( TOP );
               reflagged_e++;
            }
       }
       
    // 2. Now, elements of the upper boundary which were accidentially assigned to the 
    // left or right model boundary must be re-assigned to the upper model boundary. 
    // Strategy: Subvertically-bounded wrongly-flagged elements are identified by daughter
    // nodes whose coordinates differ from x-min and x-max.
    // Note: also the upper corner elements must be identified since these 
    //       have most likely not been identified correctly.
    stl_index                    j, bnodes;
    map<fT,stl_index,less<fT> >  side_nodes;

    
    // left model boundary (x=xmin)
    // ----------------------------
    if ( !EnlistBoundaryElements( LEFT, belmts ) )
      skm_err.notice( INFO, "ParallelSuperGroup<csp_float,dim>::FlagIrregularTopBoundary", 
                            "Model contains no LEFT boundary elements ?");
    else
    for ( lit=belmts.begin(); lit!=belmts.end(); lit++ )
      {
         side_nodes.erase( side_nodes.begin(), side_nodes.end() );

         for ( at_bottom=false, bnodes=j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
           if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT )
             {
                // test whether one deals with an element on the lower model boundary 
                if ( local_mesh.E( (*lit)-1 ).AtBoundary() == BOTTOM ||
                     local_mesh.E( (*lit)-1 ).AtBoundary() == CNR1   ||
                     local_mesh.E( (*lit)-1 ).ConnectedNode(j)->y() <= (minmax_xy[2]+right_left_tolerance) )
                  { at_bottom = true; break; }
                  
                // test whether x-coordinate of node is equal to xmin
                if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() <= (minmax_xy[0]+right_left_tolerance) ) 
                  side_nodes[ local_mesh.E( (*lit)-1 ).ConnectedNode(j)->y() ] = j;
                  
                // counting number of boundary nodes
                bnodes++;
             }
             
         // if there are no bottom nodes and if one node at most, lies on the left boundary,
         // the element is interpreted to lie on the top boundary    
         if ( !at_bottom && side_nodes.size() <= 1 && bnodes == 2 ) 
           {
              local_mesh.E( (*lit)-1 ).AtBoundary( TOP );
              for ( j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
                if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT  &&
                     local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() >= (minmax_xy[0]+right_left_tolerance) )
                  local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary( TOP );
              reflagged_e++;              
           }
         // LEFT UPPER CORNER ELEMENT (must share two nodes with left boundary)
         if ( !at_bottom && bnodes == local_mesh.E( (*lit)-1 ).Nodes() && side_nodes.size() == 2 ) 
           {
              // upper left corner
              cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagIrregularTopBoundary: Re-flagged upper left model corner CNR4"<< endl;
              local_mesh.E( (*lit)-1 ).AtBoundary( CNR4 );
                   
              // of the 2 nodes flagged left the one with the greater
              // Y-coordinate must be flagged CNR4 (this node will be the last
              // element in the side_nodes map)
              local_mesh.E( (*lit)-1 ).ConnectedNode( (*side_nodes.rbegin()).second )->AtBoundary( CNR4 );
                  
              // now the remaining node which is not located on the left model 
              // boundary is flagged TOP  
              for ( j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
                 if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT  &&
                      local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() >= (minmax_xy[0]+right_left_tolerance) )
                   local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary( TOP );
              reflagged_e++;              
           }
      }  


    // right model boundary (x=xmax)
    // -----------------------------
    if ( !EnlistBoundaryElements( RIGHT, belmts ) )
      skm_err.notice( INFO, "ParallelSuperGroup<csp_float,dim>::FlagIrregularTopBoundary", 
                            "Model contains no RIGHT boundary elements ?");
    else
    for ( lit=belmts.begin(); lit!=belmts.end(); lit++ )
      {
         side_nodes.erase( side_nodes.begin(), side_nodes.end() );

         for ( at_bottom=false, bnodes=j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
           if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT )
             {
                // test whether one deals with an element on the lower model boundary 
                if ( local_mesh.E( (*lit)-1 ).AtBoundary() == BOTTOM ||
                     local_mesh.E( (*lit)-1 ).AtBoundary() == CNR2   ||
                     local_mesh.E( (*lit)-1 ).ConnectedNode(j)->y() <= (minmax_xy[2]+right_left_tolerance) )
                  { at_bottom = true; break; }
                  
                // test whether x-coordinate of node is equal to xmax
                if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() >= (minmax_xy[1]-right_left_tolerance) ) 
                  side_nodes[ local_mesh.E( (*lit)-1 ).ConnectedNode(j)->y() ] = j;
                  
                // counting number of boundary nodes
                bnodes++;
             }
             
         // if there are no bottom nodes and if one node at most, lies on the left boundary,
         // the element is interpreted to lie on the top boundary    
         if ( !at_bottom && side_nodes.size() <= 1 && bnodes == 2 ) 
           {
              local_mesh.E( (*lit)-1 ).AtBoundary( TOP );
              for ( j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
                if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT  &&
                     local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() <= (minmax_xy[1]-right_left_tolerance) )
                  local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary( TOP );
              reflagged_e++;              
           }
           
         // RIGHT UPPER CORNER ELEMENT (must share two nodes with right boundary)
         if ( !at_bottom && bnodes == local_mesh.E( (*lit)-1 ).Nodes() && side_nodes.size() == 2 ) 
           {
              // upper left corner
              cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagIrregularTopBoundary: Re-flagged upper right model corner CNR3"<< endl;
              local_mesh.E( (*lit)-1 ).AtBoundary( CNR3 );
                   
              // of the 2 nodes flagged right the one with the greater
              // Y-coordinate must be flagged CNR3 (this node will be the last
              // element in the side_nodes map)
              local_mesh.E( (*lit)-1 ).ConnectedNode( (*side_nodes.rbegin()).second )->AtBoundary( CNR3 );
                  
              // now the remaining node which is not located on the right model 
              // boundary is flagged TOP  
              for ( j=0; j<local_mesh.E( (*lit)-1 ).Nodes(); j++ )
                if ( local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary() != NOT &&
                     local_mesh.E( (*lit)-1 ).ConnectedNode(j)->x() <= (minmax_xy[1]-right_left_tolerance) )
                  local_mesh.E( (*lit)-1 ).ConnectedNode(j)->AtBoundary( TOP );
              reflagged_e++;              
           }
      }  
      
   cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::FlagIrregularTopBoundary: Re-assigned "<< reflagged_e;
   cout <<" elements to TOP model boundary."<< endl;
 
   return true;
 
 } // end FlagIrregularTopBoundary


template<typename fT,stl_index dim>
bool ParallelSuperGroup<fT,dim>::EnlistBoundaryElements( SG_BOUNDARY side, list<stl_index>& belmts ) const
 {
     if ( !local_mesh.BoxShapedModel() )
       skm_err.notice( FATAL_ERROR, "ParallelSuperGroup<fT,dim>::EnlistBoundaryElements", 
                                    "this method should only be applied to box-shaped models" );

     belmts.erase( belmts.begin(), belmts.end() );
     
     for ( typename deque<Element<fT,dim> >::const_iterator it=local_mesh.ElementsBegin(); it!=local_mesh.ElementsEnd(); it++ )
       if ( (*it).AtBoundary() == side ) belmts.push_back( ((*it).ID()) );

     AddCornersTo( side, ELEMENT, belmts );  

     if ( belmts.empty() ) return false;
     return true;   
 }


/*M <H4>Method:</H4><CODE>
<!------------------------------------------------------------------------>
  void SuperGroup<fT,dim>::AddCornersTo( SG_BOUNDARY side, PLACEMENT place, 
                                         list<stl_index>& bdata ) const
<!------------------------------------------------------------------------>
</CODE>

<H4>Description:</H4><!--------------------------------------------------->

AddCornersTo() adds the corner enumerations CNR1 to CNR8 to the local_mesh 
hierarchy. Different numbering schemes apply as follows:<P>
   
For 2D regular-gridded local_meshes as created with the Triangulator module, 
the Y-axis points downward (hence BOTTOM is at the smallest Y and TOP 
at the largest y value of the model).<p>
     
For adaptively-refined 2D local_meshes created with Triangle, 
Y points upward as shown below (BOTTOM at smallest Y).<p>
  
The order of the supplied values for BC assignment in 3D is counter 
clockwise from the outside looking in. This is also shown below.<p>
<pre>
           sc2 ^                    y |
   2D:   y     |                3D:   | 
               |                      |
           sc1 +------->              o------->
               sc1     sc2           /       x
                    x            z  / 
</pre>

<H4>Input Arguments:</H4><!----------------------------------------------->

AddCornersTo() needs to know the enumeration value of the target model 
boundary, the variable placement (Node, ConstraintPoint, or Element)
and whether the Y-axis in a two-dimensional model points up or downwards.
The latter is important for models which were created using the 
Triangulator local_meshing tool for regular gridded local_meshes which defines Y
as downward pointing. This definition is useful when model coordinate
values are used to calculate lithostatic overburden etc.<p>

<H4>Application:</H4><!--------------------------------------------------->

AddCornersTo() adds the model corners which were previously identified
by FlagCornerElements() to a list of boundary object ID numbers. This 
list may be for instance a list which was created by EnlistBoundaryNodes().
<p>

<H4>Messages:</H4><!------------------------------------------------------>

A warning will be reported, if no corner objects could be identified.<p>

<!------------------------------------------------------------------------>
tested: O.K. */

template<typename fT,stl_index dim>
void ParallelSuperGroup<fT,dim>::AddCornersTo( SG_BOUNDARY side, PLACEMENT place, 
                                       list<stl_index>& bdata ) const
 {
    extern ErrorHandler  skm_err;
 
    if ( !local_mesh.BoxShapedModel() )
      skm_err.notice( WARNING, "ParallelSuperGroup<csp_float,dim>::AddCornersTo", 
                                   "Corners can only be used in box-shaped models.");

    typename deque<Node<fT,dim> >::const_iterator             nit;
    typename deque<ConstraintPoint<fT,dim> >::const_iterator  cit;
    typename deque<Face<fT,dim> >::const_iterator             fit;
    typename deque<Element<fT,dim> >::const_iterator          eit;
    typename list<stl_index>::iterator                        it1, it2;
    SG_BOUNDARY                                      b1, b2, b3=FRONT, b4=BACK, boundary;
    stl_index                                        points = bdata.size();
    bool                                             at_boundary;
    stl_index                                        i;

      // defining which points are wanted
      // NOTE: Up and down is reverted in 2D !
      if ( dim == 2 ) switch( side )
        {
           case BOTTOM:  b1 = CNR1; b2 = CNR2;
             break;
           case RIGHT:   b1 = CNR2; b2 = CNR3; 
             break;
           case TOP:     b1 = CNR4; b2 = CNR3; 
             break;
           case LEFT:    b1 = CNR1; b2 = CNR4; 
        }

      // numbering counter clockwise looking from the outside
      // starting with the lowest number  
      if ( dim == 3 ) switch( side )
        {
           case BOTTOM:  b1 = CNR1; b2 = CNR5; b3 = CNR6; b4 = CNR2;
             break;
           case RIGHT:   b1 = CNR2; b2 = CNR6; b3 = CNR7; b4 = CNR3;
             break;
           case TOP:     b1 = CNR3; b2 = CNR7; b3 = CNR8; b4 = CNR4;
             break;
           case LEFT:    b1 = CNR1; b2 = CNR4; b3 = CNR8; b4 = CNR5;
             break;
           case FRONT:   b1 = CNR5; b2 = CNR6; b3 = CNR7; b4 = CNR8;
             break;
           case BACK:    b1 = CNR1; b2 = CNR2; b3 = CNR3; b4 = CNR4;
           // duplicate entries will be ignored in the application of boundary conditions
           case EDGE1:   b1 = CNR1; b2 = CNR2; b3 = CNR1; b4 = CNR2;
             break;
           case EDGE2:   b1 = CNR2; b2 = CNR3; b3 = CNR2; b4 = CNR3;
             break;
           case EDGE3:   b1 = CNR3; b2 = CNR4; b3 = CNR3; b4 = CNR4;
             break;
           case EDGE4:   b1 = CNR1; b2 = CNR4; b3 = CNR1; b4 = CNR4;
             break;
           case EDGE5:   b1 = CNR1; b2 = CNR5; b3 = CNR1; b4 = CNR5;
             break;
           case EDGE6:   b1 = CNR2; b2 = CNR6; b3 = CNR2; b4 = CNR6;
             break;
           case EDGE7:   b1 = CNR3; b2 = CNR7; b3 = CNR3; b4 = CNR7;
             break;
           case EDGE8:   b1 = CNR4; b2 = CNR8; b3 = CNR4; b4 = CNR8;
             break;
           case EDGE9:   b1 = CNR5; b2 = CNR6; b3 = CNR5; b4 = CNR6;
             break;
           case EDGE10:  b1 = CNR6; b2 = CNR7; b3 = CNR6; b4 = CNR7;
             break;
           case EDGE11:  b1 = CNR7; b2 = CNR8; b3 = CNR7; b4 = CNR8;
             break;
           case EDGE12:  b1 = CNR5; b2 = CNR8; b3 = CNR5; b4 = CNR8;
        }

    switch( place )
      { 
          case NODE:
               for ( nit=local_mesh.NodesBegin(); nit!=local_mesh.NodesEnd(); nit++ )
                 if ( (*nit).AtBoundary() == b1 || (*nit).AtBoundary() == b2 || 
                      (*nit).AtBoundary() == b3 || (*nit).AtBoundary() == b4 ) bdata.push_back( (*nit).ID() );
            break;
          case CONSTRAINT_POINT:
               for ( cit=local_mesh.ConstraintPointsBegin(); cit!=local_mesh.ConstraintPointsEnd(); cit++ )
                 if ( (*cit).AtBoundary() == b1 || (*cit).AtBoundary() == b2 || 
                      (*cit).AtBoundary() == b3 || (*cit).AtBoundary() == b4 ) bdata.push_back( (*cit).ID() );
            break;
          case FACE:
               for ( fit=local_mesh.FacesBegin(); fit!=local_mesh.FacesEnd(); fit++ )
                 if ( (*fit).AtBoundary() == b1 || (*fit).AtBoundary() == b2 || 
                      (*fit).AtBoundary() == b3 || (*fit).AtBoundary() == b4 ) bdata.push_back( (*fit).ID() );
            break;
          case ELEMENT:
               for ( eit=local_mesh.ElementsBegin(); eit!=local_mesh.ElementsEnd(); eit++ )
                 {
                    for ( at_boundary=false, i=0; i<(*eit).Nodes(); i++ )
                      if ( (*eit).ConnectedNode(i)->AtBoundary() == side ) at_boundary = true;
                    
                    if ( at_boundary )
                      {   
                         boundary = (*eit).AtBoundary();
                         if ( boundary == b1 || boundary == b2 || 
                              boundary == b3 || boundary == b4 ) bdata.push_back( (*eit).ID() );
                      }
                 }
               // SKM addendum 22-6-2000: Delete all those elements from the newly created list which
               // do not have at least one node which lies at the desired boundary  
               for ( it1=bdata.begin(); it1!=bdata.end(); it1++ )
                 {
                    for ( at_boundary=false, i=0; i<local_mesh.E( (*it1)-1 ).Nodes(); i++ )
                      if ( local_mesh.E( (*it1)-1 ).ConnectedNode(i)->AtBoundary() == side ) at_boundary = true;
                    if ( !at_boundary )
                      {
                          it2 = it1++;
                          bdata.erase( it2 );
                      }
                 }
      }
    if ( bdata.size() == points )
      cout <<"\nParallelSuperGroup<"<< typeid(fT).name() <<","<< dim <<">::AddCornersTo: Warning: No corner places could be added. "<< endl;
      
 } // end 

template<typename fT,stl_index dim>
void ParallelSuperGroup<fT,dim>::AssignBoundaryValuesParallel( SG_BOUNDARY b, const char* prop, 
                                                       DATA_STYLE bcond, 
                                                       fT bmin, fT bmax  )
 {
    map<stl_index,fT>                           boundary_list;
    typename map<stl_index,fT>::const_iterator  iter;
    ScalarVariable<fT>                          sc;
    VectorVariable<fT,dim>                       vc;
    TensorVariable<fT,dim>                       ts;
    Point<fT,dim>                                xyz;
    typename deque<Node<fT,dim> >::iterator             nit;
    typename deque<ConstraintPoint<fT,dim> >::iterator  sit;
    typename deque<Face<fT,dim> >::iterator             fit;
    typename deque<Element<fT,dim> >::iterator          eit;
    csp::Index                                         prop_key = phys_vars.StorageKey(prop);
    fT                                                 value,
                                                       step = bmax - bmin;
    stl_index                                          i;

    // Check if boundary exist:
    bool boundary_exist(false);
    for ( nit=local_mesh.NodesBegin(); nit!=local_mesh.NodesEnd(); nit++ )
      if( (*nit).AtBoundary() == b )
        {
          boundary_exist = true;
          break;
        }

    if(boundary_exist)
      cout<<"\nParallelSuperGroup<fT,2U>::AssignBoundaryValues: Assigning boundary values."<< endl;
    else
      cout<<"\nParallelSuperGroup<fT,2U>::AssignBoundaryValues: Boundary does not exist, nothing is done."<< endl;

    if(boundary_exist)
      SuperGroup<fT,dim>::AssignBoundaryValues( b,prop,bcond,bmin,bmax );


} // end AssignBoundaryValues





// explicit instantiations

template class ParallelSuperGroup<csp_float,2>;
template class ParallelSuperGroup<csp_float,3>;

} // end namespace csp
