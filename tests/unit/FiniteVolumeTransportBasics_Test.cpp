
#include "FiniteVolumeTransportBasics_Test.h"

#include "vsetMakers.h"
#include "Boundary.h"
#include "vsetMakers.h"
#include "CSMP_mathUtilities.h"

#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"

#include "finiteVolumeFunctions.h"
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"
#include "NodeCenteredFiniteVolumeTransport.h"

#include "SteadyStateDiffusor.h"
#include "VelocityAndVolumeFlux.h"
#include "ExplicitTransport.h"

#include "ConstantFactor.h"

#include "ImplicitTransport.h"
#include "LinearSystemAccumulator.h"


#include "InputDataManager.h"
#include "Standard_IO_Handler.h"

#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"

#include<cmath>
#include<vector>
#include<iostream>

#define DUMP_VTK_OUTPUT

namespace csmp
{


FiniteVolumeTransportBasics_Test::FiniteVolumeTransportBasics_Test()
{
}
	
FiniteVolumeTransportBasics_Test::~FiniteVolumeTransportBasics_Test()
{
}


void
FiniteVolumeTransportBasics_Test::run()
{
#if 0
  {
    ANSYS_Model3D model( "HeuristicModel1coarse", "HeuristicModel1coarse",  "CSMP-brine-CO2-phase-variables.txt", true, true, true, true );
    model.OutputToBinaryFile("TestThisModel");
    Model<3> m1("TestThisModel");
	  }

#endif

  
#if 0
  {
    ANSYS_Model2D model( "2DCSP", "2DCSP",  "CSMP-brine-CO2-phase-variables.txt", true, true, true, true );
    test_placements(model);
  }

  {
    VSet<3U> vset;
    test_Create_Hexahedra_VSet(vset, true);
    Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
    std::cerr << "Testing hexahedra\n";
    test_placements(model);
  }
  
  {
    VSet<3U> vset;
    test_Create_Prism_VSet(vset, true);
    Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
    std::cerr << "Testing prism\n";
    test_placements(model);
  }
  
  {
    VSet<3U> vset;
    test_Create_Pyramid_VSet(vset, true);
    Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
    std::cerr << "Testing pyramid\n";
    test_placements(model);
  }
#endif

#if 0

    {
        VSet<3U> vset;
        test_Create_Hexahedra_VSet(vset, true);
        Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
        std::cerr << "Testing hexahedra\n";
        test_constant_velocity_field(model);
        model.OutputToBinaryFile("TestThisModel");
        Model<3> m1("TestThisModel");
        test_constant_velocity_field(m1);
    }
    
    {
        VSet<3U> vset;
        test_Create_Prism_VSet(vset, true);
        Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
        std::cerr << "Testing prism\n";
        test_constant_velocity_field(model);
    }

    {
        VSet<3U> vset;
        test_Create_Pyramid_VSet(vset, true);
        Model<3U> model( vset, "CSMP-2phase-variables.txt", true );
        std::cerr << "Testing pyramid\n";
        test_constant_velocity_field(model);
    }
#endif
#if 0
    test_b25();
#endif
}

  
  template<size_t dim>
  void FiniteVolumeTransportBasics_Test::test_placements(Model<dim>& model)
  {
    model.InstantiateFiniteVolumes();
    auto& gref = model.Region("Model");
    VectorVariable<dim> vD;
    switch (dim) {
      case 3:
        vD(2) = 0.3;
      case 2:
        vD(1) = 0.2;
      case 1:
        vD(0) = 0.1;
    }

    for ( auto it=gref.NodesBegin(); it!=gref.NodesEnd(); ++it ) {
      auto& n = **it;
      for (auto fip : n.AllFacetIntegrationPoints()) {
        auto& e = fip.TheElement();
        const size_t iFacet = fip.FacetId();
        const size_t iFacetIp = fip.FacetIp();
        
        const double64 facet_area = e.FacetArea(iFacet);
        const Point<dim> facet_normal = e.FacetNormal(iFacet);
        
        const Point<dim> directed_area_stencil = facet_area * facet_normal;
        const Point<dim> directed_area_parametric = fip.DirectedArea();
        
        for (size_t i = 0; i < dim; ++i) {
          _equal(directed_area_parametric[i], directed_area_stencil[i], 1e-12);
        }
        
#if 0
        const double64 vD_n = vD.DotProduct(facet_area_mapped);
        const double64 facet_flux_prev = vD_n * facet_area;
        
        const double64 facet_flux_new = vD.DotProduct(fip.DirectedArea());
        _equal(facet_flux_prev, facet_flux_new, 1e-15);
#endif
      }
    }
  }


void FiniteVolumeTransportBasics_Test::test_constant_velocity_field(Model<3U>& model)
{
    const double target_velocity = 1.0;
    const double target_k = 0.5;

    model.InstantiateFiniteVolumes();

    const csmp::Index v_key(model.Database().StorageKey("velocity"));
    const csmp::Index p_key(model.Database().StorageKey("fluid pressure"));
    const csmp::Index K_key(model.Database().StorageKey("permeability"));
    Point<3U> xyz_min, xyz_max;
    model.MinMaxCoordinates(xyz_min, xyz_max);
    xyz_max += Point<3U>(0.5,0.5,0.5);
    Point<3U> velocity_vector(0.1,0.2,0.3);
    velocity_vector.NormalizeLengthTo(target_velocity);
    Point<3U> velocity_direction = velocity_vector;
    velocity_direction.NormalizeLengthTo(1.0);

    ScalarVariable k(PLAIN, target_k);
  
    const auto& model_domain = model.Region("Model");

    for (auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it) {
        (*it)->Store(K_key, k);
    }
    for (auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it) {
        auto point = (*it)->Coordinate();
        double pressure = dotProduct(xyz_max - point,velocity_vector) / target_k;
        ScalarVariable p(PLAIN, pressure);
        (*it)->Store(p_key, p);
    }


    size_t nodes = model_domain.Nodes();
    std::vector<double64> cross_section(model_domain.Nodes());
    std::vector<double64> velocity_magnitude(model_domain.Nodes());

     for (auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it) {
         double64 csa = 0.0;
         double64 surface_area = 0.0;
         VectorVariable<3U> vc(PLAIN,PLAIN,PLAIN,0.,0.,0.);

         velocity_magnitude[(*it)->Idx()] = target_velocity;
         for (size_t j = 0; j < (*it)->Parents(); ++j) {
             auto e = (*it)->Parent(j);

             const size_t child = (*it)->ParentNodeNumber(j);
             
             for (size_t k = 0; k < e->FV()->FacetsPerSector(child); ++k) {
                 const auto facet = e->FV()->FacetSurroundingSector(child, k);
                 double64 costheta = std::abs(dotProduct(e->FacetNormal(facet), velocity_direction));
                 csa += costheta * e->FacetArea(facet);
                 surface_area += std::abs(e->FacetArea(facet));
             }

         }
         cross_section[(*it)->Idx()] = csa * 0.5;
     }

#ifdef WRITE_VTK
    VTK_Interface<3U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
    vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
#endif

    std::vector<double64> flux_balance(nodes);
    std::vector<Point<3u>> directed_area_para(nodes);
    std::vector<Point<3u>> directed_area_phys(nodes);

#ifdef WRITE_VTK
    size_t node_to_debug = 3; // change to node of interest
    std::map<Point<3u>,size_t> ptmap;
    std::vector<Point<3u>> pts;
    std::vector<std::vector<size_t>> surfs;
    std::vector<Point<3u>> surf_normal_def;
    std::vector<Point<3u>> surf_normal_outside;
    std::vector<Point<3u>> vd_remapped;
    std::vector<double64> surf_data;
    std::vector<double64> angle_data;
    size_t mesh_size = 0;
#endif

    std::vector<double64> DNR, DNS, DNT;
    for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
    {
        // std::cerr << "Element type: " << parseFiniteElementType((*it)->FE_Type()) << '\n';

        if (!(*it)->IsSurfaceElement()) {
            continue;
        }

        std::cerr << "Element type: " << parseFiniteElementType((*it)->FE_Type()) << '\n';

        const size_t nodes((*it)->Nodes());
        
        // 1. computing facet velocity in parametric space
        // -----------------------------------------------
        // 1.1 getting ipol-functions at barycentre and computing the pressure gradient in parametric space
        const Point<3U> bctr = (*it)->FV()->Barycenter();
        (*it)->FE()->dNr( bctr[0], bctr[1], bctr[2], DNR );
        (*it)->FE()->dNs( bctr[0], bctr[1], bctr[2], DNS );
        (*it)->FE()->dNt( bctr[0], bctr[1], bctr[2], DNT );
        (*it)->CoordinateMatrix();
        (*it)->FE()->Jacobian( DNR, DNS, DNT );
        const double64 jdet_bctr = (*it)->FE()->JacobianInverse();
        const DenseMatrix<DM_MIN> jinv_bctr = (*it)->FE()->JINV;

        // pressure gradients / velocities
        const double64 K((*it)->Read(K_key));
        Point<3U> vD(0.);
        for ( size_t i=0U; i<nodes; ++i ) {
            const double64 p_node((*it)->N(i)->Read(p_key));
            vD[0] += DNR[i] * p_node;
            vD[1] += DNS[i] * p_node;
            vD[2] += DNT[i] * p_node;
        }
        vD *= -K;

        Point<3u> vDremapped(jinv_bctr * vD.Coordinates());
        const double64 vd_scale = vDremapped.Length() / vD.Length();

        const size_t facets = (*it)->Facets();
        for ( size_t i=0U; i<facets; ++i ) {
            size_t inside_node, outside_node;
            (*it)->FV()->FacetEdgeNodes( i, inside_node, outside_node );

            Point<3u> facet_normal_real((*it)->FacetNormal(i));
            Point<3u> v0(0.0);
            Point<3u> v1(0.0);
            for (size_t nn = 0; nn < (*it)->Nodes(); ++nn) {
                auto xform_weights = (*it)->FV()->FacetNormalTransformationNodeWeights(i, nn);
                const Point<3u> n((*it)->N(nn)->Coordinate());
                v0 += xform_weights.first * n;
                v1 += xform_weights.second * n;
            }

            Point<3u> facet_normal_remapped(crossProduct(v1, v0));
            facet_normal_remapped.NormalizeLengthTo(1.0);

            double64 projected_velocity_parametric = dotProduct(facet_normal_remapped, vDremapped);
            double64 projected_velocity_physical = dotProduct((*it)->FacetNormal(i), velocity_vector);

            double64 flux_physical = dotProduct((*it)->FacetNormal(i), velocity_vector) * (*it)->FacetArea(i);
            double64 flux_parametric = projected_velocity_parametric * (*it)->FacetArea(i);
            
            _equal(projected_velocity_parametric, projected_velocity_physical, 1.0e-14);

            _test( (*it)->N(inside_node)->Idx() < model_domain.Nodes() );
            _test( (*it)->N(outside_node)->Idx() < model_domain.Nodes() );
            
#ifdef WRITE_VTK
            if ((*it)->N(inside_node)->Idx() == node_to_debug ||(*it)->N(outside_node)->Idx() == node_to_debug) {
                {
                    double cosAngle = std::min(dotProduct(facet_normal_remapped, (*it)->FacetNormal(i)) / (facet_normal_remapped.Length() *  (*it)->FacetNormal(i).Length()), 1.0);
                    double angle = std::acos(cosAngle) * 180.0 / M_PI;
                    std::cerr << "Angle: " << angle << '\n';
                    angle_data.push_back(angle);
                }

                std::vector<size_t> face;
                face.reserve((*it)->FV()->FacetPoints(i));
                for (size_t k = 0; k < (*it)->FV()->FacetPoints(i); ++k) {
                    const Point<3u> fp = (*it)->RstToXYZ(((*it)->FV()->FacetPoint(i, k)));
                    auto ptit = ptmap.find(fp);
                    if (ptit != ptmap.end()) {
                        face.push_back(ptit->second);
                    }
                    else {
                        ptmap[fp] = pts.size();
                        face.push_back(pts.size());
                        pts.push_back(fp);
                    }
                    std::cerr << "element " << (*it)->Idx() << " facet " << i << " point " << k << " " << fp << "\n";
                }
                surfs.push_back(face);
                mesh_size += (*it)->FV()->FacetPoints(i) + 1;
                
                vd_remapped.push_back(vDremapped);

                if ((*it)->N(inside_node)->Idx() == node_to_debug) {
                    std::cerr << "element " << (*it)->Idx() << " facet " << i << " inside += " << flux_parametric << "\n";
                    surf_data.push_back(flux_parametric);
                    surf_normal_def.push_back((*it)->FacetNormal(i));
                    surf_normal_outside.push_back(facet_normal_remapped);
                }
                else if ((*it)->N(outside_node)->Idx() == node_to_debug) {
                    std::cerr << "element " << (*it)->Idx() << " facet " << i << " outsid -= " << flux_parametric << "\n";
                    surf_data.push_back(-flux_parametric);
                    surf_normal_def.push_back(-1.0 * (*it)->FacetNormal(i));
                    surf_normal_outside.push_back(-1.0 * facet_normal_remapped);
                }
                else {
                    std::cerr << "HUH?!\n";
                }
            }
#endif

            directed_area_para[(*it)->N(inside_node)->Idx()] += facet_normal_remapped * (*it)->FacetArea(i);
            directed_area_para[(*it)->N(outside_node)->Idx()] -= facet_normal_remapped * (*it)->FacetArea(i);
            directed_area_phys[(*it)->N(inside_node)->Idx()] += (*it)->FacetNormal(i) * (*it)->FacetArea(i);
            directed_area_phys[(*it)->N(outside_node)->Idx()] -= (*it)->FacetNormal(i) * (*it)->FacetArea(i);

            flux_balance[(*it)->N(inside_node)->Idx()] += flux_parametric;
            flux_balance[(*it)->N(outside_node)->Idx()] -= flux_parametric;
        }
    }

     for ( auto pit = model_domain.PerimeterNodesBegin(); pit != model_domain.PerimeterNodesEnd(); ++pit )
     {
         _test( (*pit)->Idx() < model_domain.Nodes() );
         flux_balance[(*pit)->Idx()] = 0.;
         directed_area_phys[(*pit)->Idx()] = Point<3u>(0.);
         directed_area_para[(*pit)->Idx()] = Point<3u>(0.);

     }

     double64 fmin = +std::numeric_limits<double64>::max();
     double64 fmax = -std::numeric_limits<double64>::max();
     csmp::Index fb_key = model.CreateProperty("absolute flux balance","m^3 s-1",SCALAR,NODE);
     csmp::Index sfb_key = model.CreateProperty("scaled flux balance","1",SCALAR,NODE);
     csmp::Index da_phys_key = model.CreateProperty("directed area physical","m^2",VECTOR,NODE);
     csmp::Index da_para_key = model.CreateProperty("directed area parametric","m^2",VECTOR,NODE);

     for (size_t i = 0; i < flux_balance.size(); ++i) {
         auto n = model_domain.N(i);

         // 1. Test that the directed area around the node is zero.
         n->Store(da_phys_key, VectorVariable<3u>(directed_area_phys[i]));
         n->Store(da_para_key, VectorVariable<3u>(directed_area_para[i]));
         _test(directed_area_para[i].Length() < 1.0e-15);
         _test(directed_area_phys[i].Length() < 1.0e-15);

         // 2. Test that the flux balance around the node matches.
         double64 numer = flux_balance[i];

         n->Store(fb_key, ScalarVariable(PLAIN, numer));

         double64 denom = cross_section[i] * target_velocity;
         double64 balance = (numer == 0 || denom == 0) ? 0 : numer / denom;
         n->Store(sfb_key, ScalarVariable(PLAIN, balance));

         fmin = std::min(fmin, balance);
         fmax = std::max(fmax, balance);

         _test(std::abs(balance) < 1.0e-13);
     }
    
#ifdef WRITE_VTK
    {
        std::ofstream ofs("debug_data.vtk");
        ofs << "# vtk DataFile Version 2.0\n";
        ofs << "Finite-element dataset (CSMP): variable: facet_normal_remapped, model domain: Model, timestep: 0\n";
        ofs << "ASCII\n\n";
        ofs << "DATASET POLYDATA\n";
        ofs << "POINTS " << pts.size() << " double\n";
        for (auto& p : pts) {
            ofs << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
        ofs << "\nPOLYGONS " << surfs.size() << " " << mesh_size << "\n";
        for (auto& f : surfs) {
            ofs << f.size();
            for (auto& p : f) {
                ofs << " " << p;
            }
            ofs << '\n';
        }
        
        ofs << "\nCELL_DATA " << surf_data.size() << '\n';
        ofs << "SCALARS projected_velocity double\n";
        ofs << "LOOKUP_TABLE default\n";
        for (auto f : surf_data) {
            ofs << f << '\n';
        }
        ofs << "SCALARS angle double\n";
        ofs << "LOOKUP_TABLE default\n";
        for (auto f : angle_data) {
            ofs << f << '\n';
        }
        ofs << "VECTORS surf_normal_def double\n";
        for (auto& f : surf_normal_def) {
            ofs << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
        }
        std::cerr << "VECTORS surf_normal_outside double\n";
        for (auto& f : surf_normal_outside) {
            ofs << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
        }
        ofs << "VECTORS velocity_vector double\n";
        for (auto& f : surf_normal_outside) {
            ofs << velocity_vector[0] << ' ' << velocity_vector[1] << ' ' << velocity_vector[2] << '\n';
        }
    }
    vtk_output.OutputDataToVTK( model, "absolute-flux-balance", "absolute flux balance", 1 );
    vtk_output.OutputDataToVTK( model, "scaled-flux-balance", "scaled flux balance", 1 );
    vtk_output.OutputDataToVTK( model, "directed-area-physical", "directed area physical", 1 );
    vtk_output.OutputDataToVTK( model, "directed-area-parametric", "directed area parametric", 1 );
#endif
    
     std::cerr << "Flux balance: (" << fmin << ", "  << fmax << ")\n";
}
    
    
    void FiniteVolumeTransportBasics_Test::test_b25()
    {
        using namespace std;
#if 0
        string  model_name("cylinder");
        string  variables_name("cylinder-variables.txt");
#endif
#if 1
      string  model_name("cube100");
      string  variables_name("cube100-variables.txt");
#endif
#if 0
      string  model_name("brick");
      string  variables_name("cube100-variables.txt");
#endif


        ANSYS_Model3D  model( model_name.c_str(), model_name.c_str(), variables_name.c_str(), false, true, true, true);

        printModelDimensions( model, true );

        InputDataManager<3U>   model_configuration;
        ComputationalSettings  run_settings;
        model_configuration.ConfigureFromFile( model,
                                              model_name.c_str(),
                                              false,
                                              true,   // 2) default prop.values
                                              true,   // 3) region prop.values
                                              true,   // 4) essential box-boundary conditions
                                              true,   // 5) essential flags
                                              true,   // 6) boundary conditions
                                              run_settings );
        model.InstantiateFiniteVolumes();

        Standard_IO_Handler  stdio;
        printRangeOfVariable( model, stdio, "permeability" );
    
        const double64  fluid_viscosity(1.0e-03);
        variables::Variables_TracerTransfer	 vars(model.Database());
        auto& gref = model.Region("Model");
        {
          csmp::INDEX<TENSOR,ELEMENT> conductivity(model.Database().StorageKey("conductivity"));
          for (auto it = gref.ElementsBegin(); it != gref.ElementsEnd(); ++it) {
            auto e = (*it)->AtBarycenter();
            TensorVariable<3u> k;
            e.Read(vars.key_k, k);
            k *= 1.0 / fluid_viscosity;
            e.Store(conductivity, k);
          }
        } 

        SteadyStateDiffusor<3U,Region> steady_state_pressure( model,
                                                             "conductivity", "fluid pressure",
                                                             "fluid volume source" );
        VelocityAndVolumeFlux<3U,Element<3U> >  postpro0( model, "conductivity", "porosity", "fluid pressure" );
        steady_state_pressure.AddPostProcess( &postpro0 );
        steady_state_pressure.ComputeSteadyState( model );

        
        VTK_Interface<3U>  vtk_output;
        vtk_output.OutputDataToVTK( model, "permeability", "permeability", 0 );
        vtk_output.OutputDataToVTK( model, "fluid-pressure", "fluid pressure", 0 );
        vtk_output.OutputDataToVTK( model, "velocity",       "velocity",       0 );

#if 0
        NodeCenteredFiniteVolumeTransport<3U> advector( "Model", model,
                                                        "porosity", "concentration", "velocity",
                                                        "nodal fluid volume source",false, false );
        PropertyHandle<3U>  fv( model,"finite volume",SCALAR,NODE);
        advector.FiniteVolume( "finite volume" );
        cout <<"\n\nadvectVariableFirstOrderImplicit: Measuring the divergence of fluxes."<< endl;
        advector.Divergence( "velocity", "nodal flux mismatch" );
#endif
      
#if 1
        ExplicitTransport<3U> advector(model, "Model", false);
        advector.StepSizeReductionFactor(0.1);
#endif
#if 0
      CSMP_DEFAULT_LINEAR_SOLVER solver;

      ImplicitTransport<3U> advector(solver, model, "Model", true);
      advector.StepSizeReductionFactor(1.0);
#endif

        // the calculation of fluid pressure
      
#if 1
        vtk_output.OutputDataToVTK( model, "impl_concentration", "concentration", 0 );
        vtk_output.OutputDataToVTK( model, "pressure", "fluid pressure", 0 );

#endif

      gref.RenumberNodes();
        for (unsigned i = 1; i < 200; ++i) {
            advector.AdvectVariable(10.0);
            // advector.AdvectVariable(1000.0);
#if 1
            vtk_output.OutputDataToVTK( model, "velocity", "velocity", i );
            vtk_output.OutputDataToVTK( model, "impl_concentration", "concentration", i );
            vtk_output.OutputDataToVTK( model, "pressure", "fluid pressure", i );
            vtk_output.OutputDataToVTK( model, "fluxbalance", "flux balance", i );
            vtk_output.OutputDataToVTK( model, "ff", "facet flux", i );
            vtk_output.OutputDataToVTK( model, "ffC", "facet flux concentration", i );

            vtk_output.OutputDataToVTK( model, "newconcentration", "new concentration", i );
#endif
        }
    }


} //end namespace csmp
