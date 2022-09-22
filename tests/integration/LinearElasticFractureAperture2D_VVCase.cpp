#include "LinearElasticFractureAperture2D_VVCase.h"
#include "ANSYS_Model2D.h"
#include "VTU_Interface.h"
#include "LinearSolver.h"
#include "PDE_IntegratorExperimental.h"
#include "NumIntegral_BT_C_B_dV.h"
#include "NumIntegral_BT_D_B_dV.h"
#include "ExtractTensorVariableComponent.h"
#include "NumIntegral_PT_op_dS.h"
#include "ModelSubDomain.h"
#include "Fracture.h"
#include <stdio.h>
#include <cmath>


using namespace std;


namespace csmp
{

LinearElasticFractureAperture2D_VVCase::LinearElasticFractureAperture2D_VVCase(const char* prefix)
{
    this->setName("LinearElasticFractureAperture2D_VVCase");
    prefix_=prefix;
}

/**

    Author:  E. Pezzulli & M. Nejati
    Linear Elasticity: Plane stress/strain problem with central crack under uniform tension; testing the aperture profile of the crack
    =================================
    Domain:     Omega: w x h Rectangle (Provided meshes is for a 10x10 square plane)
    Mesh:       Linear/Quadratic Triangles tested OK; Quadrilaterials have to be tested
    Material:   Linear ELastic Isotropic Material: E=Young's Modulus, nu= Poisson's Ratio;
    Loading:    Uniaxial/Biaxial Tension/Compression on RIGHT and TOP Boundaries
    BC:         Zero x and y displacements on LEFT and BOTTOM Edges, rexpectively
    ConfigFiles:Square_20x20_WithCentralCrack_a=1_betha=0_LinearTriangles (.asc, .dat
    =================================

    Assume a rectangular plane of width and height of w, and h, with a horizontal crack at its center. The body is subjected to a uniform tensile stress S1
    in both the x and y directions. Let Uy(x) be the solution for the displacement of the fracture boundaries in the vertical direction y.
    BC: The Neumann conditions SHOULD be Sx=S1, Sy = S1 on the LEFT,TOP,RIGHT,BOTTOM boundaries.
        However what is implemenetd is the displacement boundary conditions Ux=0 on LEFT and Uy=0 on BOTTOM; and the Neumann conditions: Sx=S1 on RIGHT and Sy=S1 on TOP.
        The solution for the vertical displacement Uy(x) for an fracture of half length a is:

    Plane Strain:
    Uy(x)=4*S1/E * (a^2-x^2)^(1/2)

    The following implementation is testing the creation of the Fracture Class, its initialisation and creation of InterFace Objects, constructed from
    opposite face pairs of the master and slave Boundaries of a 2D fracture. The calculation and calibration of a Unit Normal and the creation of midPoint
    Elements within the fracture, which are however not used/tested to work.

    The following implemetation, together with the other elastic tests by MN, uses the capability to assign faces over the boundaries
    of the domain. In this case, all the 2D elements building the regions TOP, BOTTOM, LEFT and RIGHT are deleted, and equivalent
    faces are created instead. Using this functionality, one can apply boundary integrals only on the boundaries with known faces.
    An alternative method is to keep the elements, and apply the boundary interals as rhs terms. One can easily switch between these
    two methodologies below. However, the use of boundary faces has the following advantages:
    1) The lower dimension elements are deleted, and therefore they can not have any contribution to the LHS of the system. In the case of
    the presence of boundary elements, one has to make sure that they do not contribute to the stiffness matrix. This can be done, for example,
    by making Young's modulus zero on these boundary elements.
    2) When using RHS operators to add a boundary integral, an iteration happens over all the elements, including the domain higher dimension
    elements, to accumulate the boundary integrals. Although the contribution of domain elements can be forced to be zero (by making the operand
    zero for all domain elements), this iteration is very expensive compared to the one that just iterates over the boundary faces.

    */



void LinearElasticFractureAperture2D_VVCase::set_min_avg_error(double epsilon){
  min_avg_error_ = epsilon;
}

void LinearElasticFractureAperture2D_VVCase::run()
{
    // Defining constants
    enum{dim=2};
    bool plane_strain = true;
    const double s1(1.0e6) , ym (1.0e9), pr(0.3), Kc(1.0e6);
    bool quarterpoint = false;

    // Model configuration:
    //Input file directory locations
    string input_dir         = "ConfigFiles/InternalCrack/";
    string output_file       = "Results/";
    //Variable inputs
    string vars_file         = input_dir + "SneddonConvergenceTest-variables.txt";
    string regions_file      = input_dir + "SneddonConvergenceTest";
    //string config_file            = input_dir +...;
    //Possible meshes
    string coarse_mesh       = "InternalBoundary_test";
    string coarse_mesh_quad  = "InternalBoundary_Test_quadratic";
    string fine_mesh_x       = "InternalCrack_tri_lin_02";
    string fine_mesh_2x      = "InternalCrack_tri_lin_01";
    string fine_mesh_3x      = "InternalCrack_tri_lin_005";

    //choosing mesh via input
    string mesh_file         = input_dir + prefix_;

    //choosing output file name based on mesh chosen
    output_file              += prefix_;


    ANSYS_Model2D model( mesh_file.c_str() , regions_file.c_str(), vars_file.c_str(),
                         false, true, true, true, false);    // Constractor for empty variables

    model.CreateInternalBoundaryFrom("FRACTURE", true);
    model.CreateSplitBoundaryFrom( model.Boundary("FRACTURE_BOUNDARY0_MATRIX_INTERSECTION")); //relies on split boundary naming convention
    model.RegionsFromSplitBoundaries();

    /// -------------------------------
    /// Setting up Fracture Configuration
    /// -------------------------------
    Fracture<dim> myFracture (model, model.SplitBoundary("FRACTURE_SPLIT_BOUNDARY"), quarterpoint, DC_TIP );            //Creates interface objects in Fracture, and initialises & configures Lubrication region

    // Model Variables:
    model.CreateProperty( "Neumann stress",   "SI", VECTOR, FACE);

    model.InputPropertyValue("Young's modulus", ScalarVariable( PLAIN, ym ) );
    model.InputPropertyValue("Poisson's ratio", ScalarVariable( PLAIN, pr ) );
    model.InputPropertyValue("critical stress intensity factor", ScalarVariable(PLAIN, Kc));

    //initialising
    model.InputPropertyValue("displacement", VectorVariable<dim>(PLAIN, PLAIN, 0.0, 0.0));
    model.InputPropertyValue("Neumann stress", VectorVariable<dim>(PLAIN, PLAIN, 0.0, 0.0));
    model.InputPropertyValue("fluid pressure", ScalarVariable( PLAIN, 0.0 ) );
    model.InputPropertyValue("aperture", ScalarVariable( PLAIN, 0.0 ) );


    // Model Boundary Conditions:
    Node<dim> *cornerNodeLeftBottom, *cornerNodeRightBottom, *n_mid_left, *n_mid_right;
    for (vector<Node<dim>*>::const_iterator node (model.Region("Model").NodesBegin()); node != model.Region("Model").NodesEnd(); ++node )
    {
        if ( (*node)->AtBoundary() == CNR1 ) cornerNodeLeftBottom  = (*node);
        if ( (*node)->AtBoundary() == CNR2 ) cornerNodeRightBottom = (*node);
       }


    model.Boundary("RIGHT").InputPropertyValue("Neumann stress", VectorVariable<dim> (NEUMANN, NEUMANN, 0.0, 0.) );
    model.Boundary("TOP").InputPropertyValue(  "Neumann stress", VectorVariable<dim> (NEUMANN, NEUMANN, 0., s1) );
    model.Boundary("LEFT").InputPropertyValue( "Neumann stress", VectorVariable<dim> (NEUMANN, NEUMANN, -0.0, 0.0));
    model.Boundary("BOTTOM").InputPropertyValue("Neumann stress",VectorVariable<dim> (NEUMANN, NEUMANN, 0.0, -s1));
    //model.Boundary("BOTTOM").InputPropertyValue("displacement",VectorVariable<dim> (PLAIN, DIRICH, 0.0, 0.0));

    //cornerNodeLeftBottom->Store(model.Database().StorageKey("displacement"),VectorVariable<dim> (DIRICH, DIRICH , 0., 0.));
    //cornerNodeRightBottom->Store(model.Database().StorageKey("displacement"), VectorVariable<dim> (PLAIN, DIRICH, 0.0, 0.0));


    printRangeOfVariable(model, "fluid pressure");
    printRangeOfVariable(model, "Young's modulus");
    printRangeOfVariable(model, "Poisson's ratio");
    printRangeOfVariable(model, "displacement");
    printRangeOfVariable(model, "Neumann stress");
    printRangeOfVariable(model, "aperture");
    printRangeOfVariable(model, "critical stress intensity factor");


    ///==================================================================
    /// Elasticity Equations
    ///==================================================================
    // setting up integrator
    SAMG_Settings settings;
    myFracture.SetSolverSettings(settings);
    //settings.Set_ncycle(1000); // Depending on the type of problem, it may take many iterations to converge for anisotropic cases.
    //settings.Set_nxtyp(2);
    SAMG_Solver samgSolver(&settings);
    PDE_IntegratorExperimental<dim,Region> deformation (samgSolver);

    // Adding stiffness to the LHS list:
    NumIntegral_BT_D_B_dV<dim,Element<dim> > stiffnessMatrix( model.Database(), "Young's modulus", "Poisson's ratio", "displacement", "displacement" );
    stiffnessMatrix.PlaneStress(plane_strain);
    deformation.Add( &stiffnessMatrix );

    // Adding surface tractions to the RHS (BoundaryIntegrals) list:
    NumIntegral_PT_op_dS<dim> nodalTractions( model.Database(), "Neumann stress", "displacement" );
    deformation.AddBoundaryIntegral( &nodalTractions );


    ///==================================================================
    /// Displacement Solution & Verification
    ///==================================================================
    //Solve for displacements and move then nodes:
    deformation.IntegrateOver( model, model.Region("Model"), false);

    //Moving Nodes
    myFracture.StoreDisplacementDifferenceAsAperture("displacement" , "aperture");        // assumes no previous existing displacement
    model.MoveNodeCoordinatesBy("displacement");

    //Setting vtu file:
    VTU_Interface<dim> vtu( model );
    vtu.OmitZeroInFileName(true);
    vtu.DeleteConnectivity();

    //Output variables:
    list<string> outputProps;
    outputProps.push_back( "displacement");
    vtu.OutputDataToVTU( output_file, outputProps, "Model", static_cast<int>(1) );



    ///ANALYTICAL SOLUTION
    /// Aperture Profile of crack
    double64 a = 0.5 * myFracture.FractureLength(MIDDLE), max_percent_error(0.), avg_percent_error (0.0) ;        // fracture half length
    auto     F_strain = [ s1, ym, pr, a ](double64 x){return (4.0*s1*(1.0-pr*pr)/(ym))*sqrt(a*a - x*x) ;};
    auto     F_stress = [ s1, ym, a ](double64 x){return (4.0*s1/ym)*sqrt(a*a - x*x) ;};

    std::map<Point<dim>,Node<dim>* > NodeMapAfter = myFracture.NodeMap(MIDDLE);
    size_t ends = 0;
    std::vector<std::vector<double64>> data(4);
    for (typename std::map<Point<dim>,Node<dim>*>::iterator it = NodeMapAfter.begin(); it != NodeMapAfter.end(); it++){
        assert( myFracture.DistanceFromTip( NodeMapAfter, it->second, RIGHT) <= 2.01*a);
        double64 percent_error(0), Sol;
        double64 x        =  a  -  myFracture.DistanceFromTip(NodeMapAfter, it->second, RIGHT) ;
        double64 w_num    = it->second->Read(model.Database().StorageKey("aperture"));
        if (plane_strain)
          Sol = F_strain(x);
        else
          Sol = F_stress(x);

        std::cout << "\nDist from Origin: " << x << " Analytical w: " << Sol << std::endl;
        std::cout << "Dist from Origin: " << x << " Numerical  w: " << w_num << std::endl;


        //Storing solutions in vector object
        data[0].push_back(x);
        data[1].push_back(Sol);
        data[2].push_back(w_num);


        //Calculating percent error
        if (w_num == 0.0){
            ends++;
            percent_error = 0.0;
            data[3].push_back(percent_error);
          }
        else{
            percent_error  = std::fabs( w_num  - Sol ) / Sol * 100.0;
            data[3].push_back(percent_error);
            _equal(percent_error, 0.0, 10.0 );
            avg_percent_error += percent_error;
            if (max_percent_error < percent_error )
              max_percent_error   = percent_error ;
          }

        std::cout << "Percent Error e = " << percent_error << std::endl;
      }

    //safety first
    assert(ends == 2);

    //average percent error
    avg_percent_error = avg_percent_error / ( static_cast<double64>(NodeMapAfter.size()) - static_cast<double64>(ends) );
    std::cout << "Average Percent error found for single node is: " << avg_percent_error << std::endl;
    std::cout << "Maximum Percent error found for single node is: " << max_percent_error << std::endl;

    _test( avg_percent_error < min_avg_error_ );

    //Saving solution to file
    if ( false ){
      WriteSolutionToFile(data, output_file + "_Aperture.txt");
      }



    std::cout << "\nLinearElasticFractureAperture finished running VV_Case" << std::endl;

    return;
}






void LinearElasticFractureAperture2D_VVCase::SetSettings(SAMG_Settings& settings)
{

  cerr << "\nsetSolverSettings";
  settings.Set_isym(2);		/*	1	A is symmetric
                      2	A is not symmetric*/
  settings.Set_itypu(0);		/*	itypu	0 Actual content of u is chosen as first approximation
                      1	First approximation u==0
                      2	First approximation u==1
                      3	First approximation is a random function */
  settings.Set_eps(1.e-11);		/* =0.0	Stopping criteria based in eps is de-activated
                        >0.0	Iteration stops if res <= eps.res0 (res0 = frist residual)
                        <0.0	Iteration stops if res <= |eps| */
  settings.Set_rel_eps(1.e-11); /*Stopping criterion when first guess is set (itypu == 0).
                         As stopping criterion rel_eps * ||rhs|| will be used. Make sure that rel_eps
                         has a negative value, so that the absolute value of rel_eps is chosen as
                         convergence criterion. (See explanation in Set_eps())*/
  settings.Set_napproach(2);	/*	1	Scalar approach (regardless of nsys)
                      2	Unknown-based (if used in scalar system napproach will be reset to 1)
                      3-5	Point32 based approaches - selects type of interpolation to use
                      3: interp. is separate for each unknown
                      4: interp. is same " " "
                      5: interp. is point- (block-) wise*/
  settings.Set_nxtyp(0);      /*	0	Gauss-Seidel relaxation
                      1	ILU(0) (substantial increase in required memory)
                      2	ILUT
                      3	Special box relaxation
                      5	Gauss-Seidel blockwise
                      ILU - best for coupled problems? - but more expensive*/
  settings.Set_igam(1);		/*	igam		1 V-cycle (standard)
                      2	F-cycle
                      3	W-cycle
                      4	WW-cycle (very expensive)*/
  settings.Set_ncgrad(1);		/*	0	default accelerator
                      1	Preconditioner for CG (standard)
                      2	Precon. for BI-CGSTAB
                      3	Precon. for GMRES*/
  settings.Set_nkdim(0);		/*	Subswitch for ncyc.
                      0	Select default dimension
                      1-8	Dimension = nkdim+1
                      9	Dimension = 20*/

  settings.Set_iswit(5);		/*	iswit	Controls re-use of SAMG decompositions during repeated calls.
                      5 	Complete SAMG run. [Upon return, memory is released
                      4	Same as 5 except memory not released
                      3	partial setup: Re-use coarser grids and interp. but update
                      Galerkin operators. Memory not released
                      2	No setup: Re-use coarser grids, interp. and Galerkin from prev. run
                      1	Same as 2 except SAMG assumes matrix A to be the same as prev. run*/
  settings.Set_iextent(1);	/*		iextent	Memory extension switch. Selects beahaviour when limits of initial
                      dimensioning have been reached
                      0	SAMG returns with error code
                      1	SAMG allocates ext. memory and continues (if no core space,
                      writes prev. allocated data to disk
                      2	SAMG allocates ext. memory and continues (if no core space,
                      SAMG terminates)
                      3	SAMG allocates ext. memory and continues (prev. allocated data
                      is written to disk)  */
  settings.Set_norm_typ(0);	/*	norm_typ	Selects type of norm to be used in computing residuals
                      0	L2-norm
                      1	L1-norm
                      2	Maximum norm*/
  settings.Set_idmp(1);		/*  idmp   0 Coarsening history
                      1 Standard print output (coarsening history).*/
  settings.Set_igdp(1);		/* igdp  >1 is only relevant for coupled systems. Otherwise: ignored.
                      0 No particular output.
                      1 Display table on grids (full problem).
                      2 Same for all submatrices (only if nsys>1).*/
  settings.Set_iadp(1);		/*  >1 is only relevant for coupled systems.
                      0 No particular output.
                      1 Display table on coarse-level matrices (full problem).
                      2 In addition: same info for all submatrices.
                      3 In addition: connectivity info between unknowns.*/
  settings.Set_iwdp(1);		/*  iwdp   >1 is only relevant for coupled systems.
                      0 No particular output.
                      1 Display table on interpolation matrices (full problem).
                      2 In addition: same info for all submatrices.
                      3 In addition: connectivity info between unknowns.*/
  settings.Set_iout1(2);		/*	1 Table of input data and work statistics.
                      2 Standard history of cycling process.
                      3 Extended history: including all levels.
                      4 Extended history: including all levels and even partial smoothing*/
  settings.Set_iout2(1);		/*  0 No action.
                      1 Display all of SAMG’s hidden parameters.*/
  settings.Set_ncgtyp(1);		/*	1 Standard process. This is supposed to be used if A has
                      mostly negative off-diagonals. Positive off-diagonal
                      elements (if any) should be small. Variables with only
                      positive couplings will become C-variables.
                      2 Standard process except that variables which have only
                      positive couplings are treated by absolute value.
                      3 Standard process except that, for mixed-sign rows, all
                      "large" positive entries (threshold parameter ewt2) are
                      eliminated before a decision on strong connectivity is
                      made. If, for some variable i, this does not lead to a clear
                      decision, i will become a C-variable.
                      4 Same as 3 except that, if the elimination of positive
                      couplings of variable i fails to give a clear picture, this
                      option temporarily switches to the standard process 1.
                      5 Same as 4 except that variables which have only
                      positive couplings are treated by absolute value.*/
  settings.Set_nred(1);		/*	0 Standard coarsening.
                      1-4 Aggressive coarsening. 1 is most, 4 is least
                      aggressive;
                      2-3 are in between. Recommendation:
                      use 1 for anisotropic and 2 for isotropic problems.
                      5 Cluster coarsening & piecewise constant interpolation.
                      6 Cluster coarsening & multi-pass interpolation.*/
  settings.Set_ncycle(10000);

}










template<typename T>
void LinearElasticFractureAperture2D_VVCase::WriteSolutionToFile(std::vector<std::vector<T>> data, std::string name) {
    std::ofstream SolutionFile (name);
    cout << "Saving file: " << name << " --- must be a .dat file" << endl;
    assert(data.size() >= 2);
    assert(data[0].size() == data[1].size());
    for (size_t ix = 0; ix < data[0].size(); ix++){
        for (size_t iy = 0; iy<data.size(); iy++){
            SolutionFile << data[iy][ix] << "      " ;
          }
        SolutionFile << std::endl;
    }

}





} // csmp
