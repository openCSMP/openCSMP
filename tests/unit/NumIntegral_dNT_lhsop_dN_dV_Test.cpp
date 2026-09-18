// NumIntegral_dNT_lhsop_dN_dV_Test.cpp
#include "NumIntegral_dNT_lhsop_dN_dV_Test.h"

#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "TRIANGLE_Interface.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// =============================================================================
// Constructor / Destructor
// =============================================================================

NumIntegral_dNT_lhsop_dN_dV_Test::NumIntegral_dNT_lhsop_dN_dV_Test( bool verbose )
  : tol_( 0.001 ),
    verbose_( verbose )
 {
    TRIANGLE_Interface triangle_mesh_reader;
    VSet<2> mesh_container;
    triangle_mesh_reader.ReadTriangle2DMesh( "iso.1", mesh_container );
    mesh_container.SingleElementType( ISOPARAMETRIC_LINEAR_TRIANGLE );
    mesh_container.EstablishElementConnectivity2D();

    if ( verbose_ )
      cout << "\nNumIntegral_dNT_lhsop_dN_dV_Test: Building Model..." << endl;

    sg_ = new Model<2U>( mesh_container, "CSMP-1phase-variables.txt" );
    sg_->Region("Model").RenumberNodes();

    sg_->InputPropertyValue( "fluid pressure", makeScalar( PLAIN, 1.0 ) );
    sg_->InputPropertyValue( "diffusivity",    makeScalar( PLAIN, 1.0 ) );
 }


NumIntegral_dNT_lhsop_dN_dV_Test::~NumIntegral_dNT_lhsop_dN_dV_Test()
 {
    delete sg_;
 }


// =============================================================================
// Test runner
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::run()
 {
    scalarTest();
    vectorTest();
    tensorTest();
    symmetryTest();
    rowSumTest();
 }


// =============================================================================
// scalarTest — isotropic scalar diffusivity
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::scalarTest()
 {
    NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                             "diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    if ( verbose_ )
      {
        cout << "\nscalarTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference values
    // -----------------------------------------------------------------------
    // K_jk^(e) = E * (g_j . g_k) / (4A)
    // s = 125000 / (4 * 62500) = 0.5
    //
    // Diagonal:
    //   K_00 = E0*(g0.g0) + E3*(g0.g0) = (1+4)*s*2 ... wait
    //
    // Let me be precise. g_j.g_j for each element:
    //   e0: g0.g0 = 250^2+250^2 = 125000, g1.g1 = 125000, g4.g4 = 250000
    //   e1: g1.g1 = 125000,                g2.g2 = 125000, g4.g4 = 250000
    //   e2: g2.g2 = 125000,                g3.g3 = 125000, g4.g4 = 250000
    //   e3: g3.g3 = 125000,                g0.g0 = 125000, g4.g4 = 250000
    //
    // K_jk^(e) = E * (g_j.g_k) / (4A) = E * (g_j.g_k) / 250000
    //
    // s = 125000/250000 = 0.5
    // s2 = 250000/250000 = 1.0
    //
    // Diagonal:
    //   K_00 = E0*s + E3*s = (1+4)*0.5 = 2.5
    //   K_11 = E0*s + E1*s = (1+2)*0.5 = 1.5
    //   K_22 = E1*s + E2*s = (2+3)*0.5 = 2.5
    //   K_33 = E2*s + E3*s = (3+4)*0.5 = 3.5
    //   K_44 = E0*s2 + E1*s2 + E2*s2 + E3*s2 = (1+2+3+4)*1.0 = 10.0
    //
    // Off-diagonal:
    //   K_01 = E0*(g0.g1)/(4A) = 1*(-250*-250+250*-250)/250000 = 1*(62500-62500)/250000 = 0
    //   K_03 = E3*(g0.g3)/(4A) = 4*(-250*250+250*250)/250000   = 4*(0)/250000           = 0
    //   K_04 = E0*(g0.g4)/(4A) + E3*(g0.g4)/(4A)
    //        = 1*(-250*500+250*0)/250000 + 4*(-250*0+250*-500)/250000
    //        = 1*(-125000)/250000 + 4*(-125000)/250000
    //        = -0.5 + -2.0 = -2.5  ... wait, let me recheck e3
    //
    // In e3: nodes are {3,0,4}, so g0 in e3 = (-250,250), g4 in e3 = (0,-500)
    //   K_04 from e3 = E3*(g0.g4)/(4A) = 4*(-250*0+250*-500)/250000
    //                = 4*(-125000)/250000 = -2.0
    //   K_04 from e0 = E0*(g0.g4)/(4A) = 1*(-250*500+250*0)/250000
    //                = 1*(-125000)/250000 = -0.5
    //   K_04 total = -0.5 + -2.0 = -2.5  ✓ matches -5s = -5*0.5
    //
    //   K_12 = E1*(g1.g2)/(4A) = 2*(-250*250+(-250)*(-250))/250000
    //        = 2*(−62500+62500)/250000 = 0
    //   K_14 = E0*(g1.g4)/(4A) + E1*(g1.g4)/(4A)
    //        e0: g1=(-250,-250), g4=(500,0)  -> g1.g4 = -125000
    //        e1: g1=(-250,-250), g4=(0,500)  -> g1.g4 = -125000
    //        = (1+2)*(-125000)/250000 = 3*(-0.5) = -1.5  ✓ matches -3s
    //   K_23 = E2*(g2.g3)/(4A) = 3*(250*250+(-250)*250)/250000
    //        = 3*(62500-62500)/250000 = 0
    //   K_24 = E1*(g2.g4)/(4A) + E2*(g2.g4)/(4A)
    //        e1: g2=(250,-250), g4=(0,500)   -> g2.g4 = -125000
    //        e2: g2=(250,-250), g4=(-500,0)  -> g2.g4 = -125000
    //        = (2+3)*(-125000)/250000 = 5*(-0.5) = -2.5  ✓ matches -5s
    //   K_34 = E2*(g3.g4)/(4A) + E3*(g3.g4)/(4A)
    //        e2: g3=(250,250),  g4=(-500,0)  -> g3.g4 = -125000
    //        e3: g3=(250,250),  g4=(0,-500)  -> g3.g4 = -125000
    //        = (3+4)*(-125000)/250000 = 7*(-0.5) = -3.5  ✓ matches -7s
    // -----------------------------------------------------------------------

    constexpr double s{ 0.5 };
    const double t{ tol_ * s };

    // diagonal
    _equal( sm.At(0,0),  5.0*s, t );
    _equal( sm.At(1,1),  3.0*s, t );
    _equal( sm.At(2,2),  5.0*s, t );
    _equal( sm.At(3,3),  7.0*s, t );
    _equal( sm.At(4,4), 20.0*s, t );

    // off-diagonal — node 4 couples to all boundary nodes
    _equal( sm.At(0,4), -5.0*s, t );
    _equal( sm.At(1,4), -3.0*s, t );
    _equal( sm.At(2,4), -5.0*s, t );
    _equal( sm.At(3,4), -7.0*s, t );

    _equal( sm.At(4,0), -5.0*s, t );  // symmetric
    _equal( sm.At(4,1), -3.0*s, t );  // symmetric
    _equal( sm.At(4,2), -5.0*s, t );  // symmetric
    _equal( sm.At(4,3), -7.0*s, t );  // symmetric

    // boundary node pairs — zero coupling through stiffness matrix
    // (these entries may not exist in the sparse structure — check first)
    const std::vector<std::pair<size_t,size_t>> zero_pairs{
        {0U,1U}, {0U,2U}, {0U,3U},
        {1U,2U}, {1U,3U},
        {2U,3U}
    };
    for ( const auto& [i,j] : zero_pairs )
      if ( sm.HasEntry(i,j) )
        _equal( sm.At(i,j), 0.0, t );
 }


// =============================================================================
// vectorTest — diagonal anisotropic diffusivity D = diag(Dx, Dy)
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::vectorTest()
 {
    // set diffusivity as a vector variable: D = diag(Dx, Dy)
    sg_->CreateProperty("vector diffusivity", "Dvc", "m2/Pa", VECTOR, ELEMENT );

    NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                             "vector diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    // use Dx=1, Dy=2 for all elements — anisotropic but uniform
    setElementVector( { {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0} },
                      "vector diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    if ( verbose_ )
      {
        cout << "\nvectorTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for D = diag(1,2), uniform across elements
    // -----------------------------------------------------------------------
    // K_jk^(e) = (g_j^T [D] g_k) / (4A)
    //          = (Dx*gx_j*gx_k + Dy*gy_j*gy_k) / (4A)
    //
    // With Dx=1, Dy=2, 4A=250000:
    //
    // e0: g0=(-250,250), g1=(-250,-250), g4=(500,0)
    //   K_00^(0) = (1*250^2 + 2*250^2)/250000 = (62500+125000)/250000 = 0.75
    //   K_11^(0) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_44^(0) = (1*500^2 + 2*0^2)/250000   = 250000/250000 = 1.0
    //   K_01^(0) = (1*(-250)(-250)+2*(250)(-250))/250000 = (62500-125000)/250000 = -0.25
    //   K_04^(0) = (1*(-250)(500)+2*(250)(0))/250000     = -125000/250000 = -0.5
    //   K_14^(0) = (1*(-250)(500)+2*(-250)(0))/250000    = -125000/250000 = -0.5
    //
    // e1: g1=(-250,-250), g2=(250,-250), g4=(0,500)
    //   K_11^(1) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_22^(1) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_44^(1) = (1*0^2   + 2*500^2)/250000 = 500000/250000 = 2.0
    //   K_12^(1) = (1*(-250)(250)+2*(-250)(-250))/250000 = (-62500+125000)/250000 = 0.25
    //   K_14^(1) = (1*(-250)(0)+2*(-250)(500))/250000    = -250000/250000 = -1.0
    //   K_24^(1) = (1*(250)(0)+2*(-250)(500))/250000     = -250000/250000 = -1.0
    //
    // e2: g2=(250,-250), g3=(250,250), g4=(-500,0)
    //   K_22^(2) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_33^(2) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_44^(2) = (1*500^2 + 2*0^2)/250000   = 1.0
    //   K_23^(2) = (1*(250)(250)+2*(-250)(250))/250000 = (62500-125000)/250000 = -0.25
    //   K_24^(2) = (1*(250)(-500)+2*(-250)(0))/250000  = -125000/250000 = -0.5
    //   K_34^(2) = (1*(250)(-500)+2*(250)(0))/250000   = -125000/250000 = -0.5
    //
    // e3: g3=(250,250), g0=(-250,250), g4=(0,-500)
    //   K_33^(3) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_00^(3) = (1*250^2 + 2*250^2)/250000 = 0.75
    //   K_44^(3) = (1*0^2   + 2*500^2)/250000 = 2.0
    //   K_03^(3) = (1*(-250)(250)+2*(250)(250))/250000 = (-62500+125000)/250000 = 0.25
    //   K_04^(3) = (1*(-250)(0)+2*(250)(-500))/250000  = -250000/250000 = -1.0
    //   K_34^(3) = (1*(250)(0)+2*(250)(-500))/250000   = -250000/250000 = -1.0
    //
    // Global assembly:
    //   K_00 = K_00^(0) + K_00^(3) = 0.75 + 0.75 = 1.5
    //   K_11 = K_11^(0) + K_11^(1) = 0.75 + 0.75 = 1.5
    //   K_22 = K_22^(1) + K_22^(2) = 0.75 + 0.75 = 1.5
    //   K_33 = K_33^(2) + K_33^(3) = 0.75 + 0.75 = 1.5
    //   K_44 = K_44^(0) + K_44^(1) + K_44^(2) + K_44^(3) = 1+2+1+2 = 6.0
    //   K_01 = K_01^(0) = -0.25
    //   K_03 = K_03^(3) = 0.25
    //   K_04 = K_04^(0) + K_04^(3) = -0.5 + -1.0 = -1.5
    //   K_12 = K_12^(1) = 0.25
    //   K_14 = K_14^(0) + K_14^(1) = -0.5 + -1.0 = -1.5
    //   K_23 = K_23^(2) = -0.25
    //   K_24 = K_24^(1) + K_24^(2) = -1.0 + -0.5 = -1.5
    //   K_34 = K_34^(2) + K_34^(3) = -0.5 + -1.0 = -1.5
    // -----------------------------------------------------------------------

    const double t{ tol_ * 0.25 };

    // diagonal
    _equal( sm.At(0,0), 1.5, t );
    _equal( sm.At(1,1), 1.5, t );
    _equal( sm.At(2,2), 1.5, t );
    _equal( sm.At(3,3), 1.5, t );
    _equal( sm.At(4,4), 6.0, t );

    // off-diagonal
    _equal( sm.At(0,1), -0.25, t );
    _equal( sm.At(0,3),  0.25, t );
    _equal( sm.At(0,4), -1.5,  t );

    _equal( sm.At(1,0), -0.25, t );  // symmetric
    _equal( sm.At(1,2),  0.25, t );
    _equal( sm.At(1,4), -1.5,  t );

    _equal( sm.At(2,1),  0.25, t );  // symmetric
    _equal( sm.At(2,3), -0.25, t );
    _equal( sm.At(2,4), -1.5,  t );

    _equal( sm.At(3,0),  0.25, t );  // symmetric
    _equal( sm.At(3,2), -0.25, t );  // symmetric
    _equal( sm.At(3,4), -1.5,  t );

    _equal( sm.At(4,0), -1.5,  t );  // symmetric
    _equal( sm.At(4,1), -1.5,  t );  // symmetric
    _equal( sm.At(4,2), -1.5,  t );  // symmetric
    _equal( sm.At(4,3), -1.5,  t );  // symmetric
 }


// =============================================================================
// tensorTest — full tensor anisotropic diffusivity
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::tensorTest()
 {
    // set diffusivity as a vector variable: D = diag(Dx, Dy)
    sg_->CreateProperty("tensor diffusivity", "Dts", "m2/Pa", TENSOR, ELEMENT );

    NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                             "tensor diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    // set diffusivity as a full tensor: D = [[1,0.5],[0.5,2]]
    // uniform across all elements — symmetric positive definite
    setElementTensor( { {1.0, 0.5, 0.5, 2.0},
                        {1.0, 0.5, 0.5, 2.0},
                        {1.0, 0.5, 0.5, 2.0},
                        {1.0, 0.5, 0.5, 2.0} },
                      "tensor diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    if ( verbose_ )
      {
        cout << "\ntensorTest: SparseMatrix after accumulation:\n";
        sm.Out();
      }

    // -----------------------------------------------------------------------
    // Analytical reference values for D = [[1,0.5],[0.5,2]], uniform
    // -----------------------------------------------------------------------
    // K_jk^(e) = (g_j^T [D] g_k) / (4A)
    //
    // [D] g_k = [Dxx*gx_k + Dxy*gy_k, Dyx*gx_k + Dyy*gy_k]
    //         = [gx_k + 0.5*gy_k, 0.5*gx_k + 2*gy_k]
    //
    // g_j^T [D] g_k = gx_j*(gx_k + 0.5*gy_k) + gy_j*(0.5*gx_k + 2*gy_k)
    //               = gx_j*gx_k + 0.5*gx_j*gy_k + 0.5*gy_j*gx_k + 2*gy_j*gy_k
    //
    // For each element, compute g_j^T [D] g_k for all node pairs:
    //
    // e0: g0=(-250,250), g1=(-250,-250), g4=(500,0)
    //   [D]g0 = (-250+125, -125+500)   = (-125, 375)
    //   [D]g1 = (-250-125, -125-500)   = (-375,-625)
    //   [D]g4 = (500+0,    250+0)      = (500,  250)
    //
    //   K_00^(0) = g0.[D]g0 = (-250)(-125)+(250)(375)  = 31250+93750  = 125000
    //   K_11^(0) = g1.[D]g1 = (-250)(-375)+(-250)(-625)= 93750+156250 = 250000
    //   K_44^(0) = g4.[D]g4 = (500)(500)+(0)(250)      = 250000
    //   K_01^(0) = g0.[D]g1 = (-250)(-375)+(250)(-625) = 93750-156250 = -62500
    //   K_04^(0) = g0.[D]g4 = (-250)(500)+(250)(250)   = -125000+62500= -62500
    //   K_14^(0) = g1.[D]g4 = (-250)(500)+(-250)(250)  = -125000-62500= -187500
    //
    // e1: g1=(-250,-250), g2=(250,-250), g4=(0,500)
    //   [D]g1 = (-250-125, -125-500)   = (-375,-625)
    //   [D]g2 = (250-125,  125-500)    = (125, -375)
    //   [D]g4 = (0+250,    0+1000)     = (250, 1000)
    //
    //   K_11^(1) = g1.[D]g1 = (-250)(-375)+(-250)(-625)= 93750+156250 = 250000
    //   K_22^(1) = g2.[D]g2 = (250)(125)+(-250)(-375)  = 31250+93750  = 125000
    //   K_44^(1) = g4.[D]g4 = (0)(250)+(500)(1000)     = 500000
    //   K_12^(1) = g1.[D]g2 = (-250)(125)+(-250)(-375) = -31250+93750 = 62500
    //   K_14^(1) = g1.[D]g4 = (-250)(250)+(-250)(1000) = -62500-250000= -312500
    //   K_24^(1) = g2.[D]g4 = (250)(250)+(-250)(1000)  = 62500-250000 = -187500
    //
    // e2: g2=(250,-250), g3=(250,250), g4=(-500,0)
    //   [D]g2 = (250-125,  125-500)    = (125, -375)
    //   [D]g3 = (250+125,  125+500)    = (375,  625)
    //   [D]g4 = (-500+0,  -250+0)      = (-500,-250)
    //
    //   K_22^(2) = g2.[D]g2 = (250)(125)+(-250)(-375)  = 31250+93750  = 125000
    //   K_33^(2) = g3.[D]g3 = (250)(375)+(250)(625)    = 93750+156250 = 250000
    //   K_44^(2) = g4.[D]g4 = (-500)(-500)+(0)(-250)   = 250000
    //   K_23^(2) = g2.[D]g3 = (250)(375)+(-250)(625)   = 93750-156250 = -62500
    //   K_24^(2) = g2.[D]g4 = (250)(-500)+(-250)(-250) = -125000+62500= -62500
    //   K_34^(2) = g3.[D]g4 = (250)(-500)+(250)(-250)  = -125000-62500= -187500
    //
    // e3: g3=(250,250), g0=(-250,250), g4=(0,-500)
    //   [D]g3 = (250+125,  125+500)    = (375,  625)
    //   [D]g0 = (-250+125,-125+500)    = (-125, 375)
    //   [D]g4 = (0-250,    0-1000)     = (-250,-1000)
    //
    //   K_33^(3) = g3.[D]g3 = (250)(375)+(250)(625)    = 93750+156250 = 250000
    //   K_00^(3) = g0.[D]g0 = (-250)(-125)+(250)(375)  = 31250+93750  = 125000
    //   K_44^(3) = g4.[D]g4 = (0)(-250)+(-500)(-1000)  = 500000
    //   K_03^(3) = g0.[D]g3 = (-250)(375)+(250)(625)   = -93750+156250= 62500
    //   K_04^(3) = g0.[D]g4 = (-250)(-250)+(250)(-1000)= 62500-250000 = -187500
    //   K_34^(3) = g3.[D]g4 = (250)(-250)+(250)(-1000) = -62500-250000= -312500
    //
    // Global assembly (divide all by 4A = 250000):
    //   K_00 = (K_00^(0) + K_00^(3))/250000 = (125000+125000)/250000 = 1.0
    //   K_11 = (K_11^(0) + K_11^(1))/250000 = (250000+250000)/250000 = 2.0
    //   K_22 = (K_22^(1) + K_22^(2))/250000 = (125000+125000)/250000 = 1.0
    //   K_33 = (K_33^(2) + K_33^(3))/250000 = (250000+250000)/250000 = 2.0
    //   K_44 = (250000+500000+250000+500000)/250000                   = 6.0
    //   K_01 = K_01^(0)/250000 = -62500/250000                        = -0.25
    //   K_03 = K_03^(3)/250000 = 62500/250000                         =  0.25
    //   K_04 = (K_04^(0)+K_04^(3))/250000 = (-62500-187500)/250000    = -1.0
    //   K_12 = K_12^(1)/250000 = 62500/250000                         =  0.25
    //   K_14 = (K_14^(0)+K_14^(1))/250000 = (-187500-312500)/250000   = -2.0
    //   K_23 = K_23^(2)/250000 = -62500/250000                        = -0.25
    //   K_24 = (K_24^(1)+K_24^(2))/250000 = (-187500-62500)/250000    = -1.0
    //   K_34 = (K_34^(2)+K_34^(3))/250000 = (-187500-312500)/250000   = -2.0
    // -----------------------------------------------------------------------

    const double t{ tol_ * 0.25 };

    // diagonal
    _equal( sm.At(0,0), 1.0, t );
    _equal( sm.At(1,1), 2.0, t );
    _equal( sm.At(2,2), 1.0, t );
    _equal( sm.At(3,3), 2.0, t );
    _equal( sm.At(4,4), 6.0, t );

    // off-diagonal
    _equal( sm.At(0,1), -0.25, t );
    _equal( sm.At(0,3),  0.25, t );
    _equal( sm.At(0,4), -1.0,  t );

    _equal( sm.At(1,0), -0.25, t );  // symmetric
    _equal( sm.At(1,2),  0.25, t );
    _equal( sm.At(1,4), -2.0,  t );

    _equal( sm.At(2,1),  0.25, t );  // symmetric
    _equal( sm.At(2,3), -0.25, t );
    _equal( sm.At(2,4), -1.0,  t );

    _equal( sm.At(3,0),  0.25, t );  // symmetric
    _equal( sm.At(3,2), -0.25, t );  // symmetric
    _equal( sm.At(3,4), -2.0,  t );

    _equal( sm.At(4,0), -1.0,  t );  // symmetric
    _equal( sm.At(4,1), -2.0,  t );  // symmetric
    _equal( sm.At(4,2), -1.0,  t );  // symmetric
    _equal( sm.At(4,3), -2.0,  t );  // symmetric
 }


// =============================================================================
// symmetryTest — stiffness matrix must be symmetric for any operand
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::symmetryTest()
 {
    // scalar case — must be symmetric
    // NB: 'tensor diffusivity" was created earlier
    {
      NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                               "diffusivity",
                                               "fluid pressure",
                                               "fluid pressure" );
      setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );
      const size_t n{ sg_->Region("Model").Nodes() };
      SparseMatrix sm( n );
      calculateGlobalMatrix( sm, stiffness );
      for ( const auto& [i,j] : connected_pairs_ )
        _equal( sm.At(i,j), sm.At(j,i), tol_ );
    }

    // symmetric tensor case — must also be symmetric
    {
      NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                               "tensor diffusivity",
                                               "fluid pressure",
                                               "fluid pressure" );
      setElementTensor( { {1.0, 0.5, 0.5, 2.0},
                          {1.0, 0.5, 0.5, 2.0},
                          {1.0, 0.5, 0.5, 2.0},
                          {1.0, 0.5, 0.5, 2.0} },
                        "tensor diffusivity" );
      const size_t n{ sg_->Region("Model").Nodes() };
      SparseMatrix sm( n );
      calculateGlobalMatrix( sm, stiffness );
      for ( const auto& [i,j] : connected_pairs_ )
        _equal( sm.At(i,j), sm.At(j,i), tol_ );
    }
 }

// =============================================================================
// rowSumTest — row sums must be zero (conservation / no-flux property)
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::rowSumTest()
 {
    NumIntegral_dNT_lhsop_dN_dV<2U> stiffness( sg_->Database(),
                                             "diffusivity",
                                             "fluid pressure",
                                             "fluid pressure" );

    setElementScalar( { 1.0, 2.0, 3.0, 4.0 }, "diffusivity" );

    const size_t n{ sg_->Region("Model").Nodes() };
    SparseMatrix sm( n );
    calculateGlobalMatrix( sm, stiffness );

    if ( verbose_ )
      {
        cout << "\nrowSumTest: SparseMatrix after accumulation:\n";
        sm.Out();
        cout << "\nrowSumTest: Row sums:\n";
      }

    // known connected columns per row — avoids accessing non-existent entries
    const std::vector<std::vector<size_t>> connected{
        { 0U, 1U, 3U, 4U },        // node 0
        { 0U, 1U, 2U, 4U },        // node 1
        { 1U, 2U, 3U, 4U },        // node 2
        { 0U, 2U, 3U, 4U },        // node 3
        { 0U, 1U, 2U, 3U, 4U }     // node 4
    };

    // for a stiffness matrix with no Dirichlet BCs applied,
    // each row must sum to zero — this is the discrete conservation property
    for ( size_t i{0U}; i < n; ++i )
      {
        double row_sum{0.0};
        for ( const size_t j : connected[i] )
          row_sum += sm.At(i,j);

        if ( verbose_ )
          cout << "  row " << i << " sum = " << row_sum << "\n";

        _equal( row_sum, 0.0, tol_ * sm.At(i,i) );
      }
 }


// =============================================================================
// Helper methods
// =============================================================================

void NumIntegral_dNT_lhsop_dN_dV_Test::setElementScalar( const vector<double>& val,
                                                       const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_lhsop_dN_dV_Test::setElementScalar",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      (*it)->Store( key, ScalarVariable( PLAIN, val[i++] ) );
 }


void NumIntegral_dNT_lhsop_dN_dV_Test::setElementVector( const vector<array<double,2>>& val,
                                                       const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_lhsop_dN_dV_Test::setElementVector",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        VectorVariable<2U> vc;
        vc(0) = val[i][0];
        vc(1) = val[i][1];
        (*it)->Store( key, vc );
        ++i;
      }
 }


void NumIntegral_dNT_lhsop_dN_dV_Test::setElementTensor( const vector<array<double,4>>& val,
                                                       const char* var_name )
 {
    const csmp::Index key( sg_->Database().StorageKey( var_name ) );
    Region<2>& domain( sg_->Region("Model") );

    if ( val.size() != domain.Cells() )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_lhsop_dN_dV_Test::setElementTensor",
                             var_name, "Variable size does not match number of elements." );

    uint32_t i{0U};
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        // val[i] = {Dxx, Dxy, Dyx, Dyy}
        TensorVariable<2U> ts;
        ts(0,0) = val[i][0];  ts(0,1) = val[i][1];
        ts(1,0) = val[i][2];  ts(1,1) = val[i][3];
        (*it)->Store( key, ts );
        ++i;
      }
 }


void NumIntegral_dNT_lhsop_dN_dV_Test::calculateGlobalMatrix( SparseMatrix& sm,
                                                           MathOperatorLHS<2U>& oper )
 {
    Region<2>& domain( sg_->Region("Model") );
    for ( auto it = domain.CellsBegin(); it != domain.CellsEnd(); ++it )
      {
        oper.GetOperands(        *(*it) );
        oper.ComputeContribution(*(*it) );
        oper.AssignToGlobal(     *(*it), sm );
      }
 }

} // namespace csmp
