#include "IsoparametricQuadraticTetrahedron.h"
#include "CSMP_mathUtilities.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**

Default constructor initializes quadratic tetrahedron to use local coordinates
and a 4-point integration scheme where the Gauss points are located
near the corners.

*/
IsoparametricQuadraticTetrahedron::IsoparametricQuadraticTetrahedron()
  // CSMP_FEM_TYPE, isoparametric(y/n), uses_local_coordinates(y/n), order_of_shape_functions
  : FiniteElement( ISOPARAMETRIC_QUADRATIC_TETRAHEDRON, true, true, 2U ),
    DN(3,10),
    BEE(3,10),
    NXYZ(10,3),
    IP(4,3)
 {
   dim = 3;
   itp = 2;
   npf = 6;
   npe = 10;
   fpe = 4;
   spe = 6;
   epe = 4;
   nne = 6;
   cne = 6;
   XY.Resize(npe,dim);
   M.Resize(4,4);
   JAC.Resize(dim,dim);
   JINV.Resize(dim,dim);

   // base class vectors
   NRST.resize(npe);
   DNR.resize(npe);
   DNS.resize(npe);
   DNT.resize(npe);

   // initializing local node coordinates
   // (right-angle in origin tetrahedron of unit-length sides, volume = 1/6)
   // corner nodes
   NXYZ(0,0) = 0.0, NXYZ(0,1) = 0.0, NXYZ(0,2) = 0.0;
   NXYZ(1,0) = 1.0, NXYZ(1,1) = 0.0, NXYZ(1,2) = 0.0;
   NXYZ(2,0) = 0.0, NXYZ(2,1) = 1.0, NXYZ(2,2) = 0.0;
   NXYZ(3,0) = 0.0, NXYZ(3,1) = 0.0, NXYZ(3,2) = 1.0;
   // midside nodes at base
   NXYZ(4,0) = 0.5, NXYZ(4,1) = 0.0, NXYZ(4,2) = 0.0;
   NXYZ(5,0) = 0.5, NXYZ(5,1) = 0.5, NXYZ(5,2) = 0.0;
   NXYZ(6,0) = 0.0, NXYZ(6,1) = 0.5, NXYZ(6,2) = 0.0;
   // midside nodes on sides
   NXYZ(7,0) = 0.0, NXYZ(7,1) = 0.0, NXYZ(7,2) = 0.5;
   NXYZ(8,0) = 0.5, NXYZ(8,1) = 0.0, NXYZ(8,2) = 0.5;
   NXYZ(9,0) = 0.0, NXYZ(9,1) = 0.5, NXYZ(9,2) = 0.5;

   // intializing fixed size matrix with the shape function derivatives
   // at each node, see Mathematica Notebook IsoparametricQuadraticTetrahedron.nb
   DN(0,0) = -3.0, DN(1,0) = -3.0, DN(2,0) = -3.0; 
   DN(0,1) =  3.0, DN(1,1) =  0.0, DN(2,1) =  0.0; 
   DN(0,2) =  0.0, DN(1,2) =  3.0, DN(2,2) =  0.0;
   DN(0,3) =  0.0, DN(1,3) =  0.0, DN(2,3) =  3.0;
   DN(0,4) =  0.0, DN(1,4) = -2.0, DN(2,4) = -2.0;
   DN(0,5) =  2.0, DN(1,5) =  2.0, DN(2,5) =  0.0;
   DN(0,6) = -2.0, DN(1,6) =  0.0, DN(2,6) = -2.0;
   DN(0,7) = -2.0, DN(1,7) = -2.0, DN(2,7) =  0.0;
   DN(0,8) =  2.0, DN(1,8) =  0.0, DN(2,8) =  2.0;
   DN(0,9) =  0.0, DN(1,9) =  2.0, DN(2,9) =  2.0;

   // analytic location of integration points in r-s-t coordinates in unit-tetrahedron,
   // see Akin, FEA (course, Rice Uni, Houston), albeit with different numbering
   // see also Zienkiewicz and Taylor, Vol I, page 223
   gpe = 4;
   const double64 a = (5. + 3. * sqrt(5.)) / 20.;  // 0.585..
   const double64 b = (5. - sqrt(5.)) / 20.;       // 0.138..
   IP(0,0)= b; IP(0,1)= b; IP(0,2)= b; // next to node 0
   IP(1,0)= a; IP(1,1)= b; IP(1,2)= b; // node 1
   IP(2,0)= b; IP(2,1)= a; IP(2,2)= b; // node 2
   IP(3,0)= b; IP(3,1)= b; IP(3,2)= a; // node 3

   W.resize( gpe );
   const double64 evolume_in_param_space(1./6.);
   const double64 one_divided_by_ips(1./gpe);
   const double64 ip_weight(one_divided_by_ips * evolume_in_param_space);
   W[0] = W[1] = W[2] = W[3] = ip_weight;

   // QuadratureRules_4Points();
   UsesLocalCoordinates(true);
   Isoparametric(true);
   VolumeElement();
   ElementType(ISOPARAMETRIC_QUADRATIC_TETRAHEDRON);

 } // end (constructor)




IsoparametricQuadraticTetrahedron::~IsoparametricQuadraticTetrahedron()
 {
 }




/** Integration rules implemented by MN:

As the total volume of the tetrahedron element is 1/6, the weights have to add up to 1/6.
Location of the integration points in r-s-t coordinates see Zienkiewicz and Taylor page 177
note that L1 = 1 - r - s - t
Higher order rules for tetrahedrons check Keast_CMAME_1986 and Jinyun_CMAME_1984;

*/

// MN: Zienkiewicz and Taylor page 177; Rao Page 151; Suitable For Linear Elements (Polynomial of degree 1).
// Tested: NOT OK; gives zero for diagonal members of the stiffness matrix; 
void IsoparametricQuadraticTetrahedron::QuadratureRules_1Point()
{
	IP.Resize(1,3);
	gpe = 1;
	const double vol (1./6.);
	IP(0,0) = 0.25, IP(0,1) = 0.25, IP(0,2) = 0.25; W[0] = 1. * vol;
}

//MN: Zienkiewicz and Taylor page 177; Rao Page 151; Suitable for the polynomial of degree 2 (Jinyun_CMAEM_1984);
//MN: Tested OK;
void IsoparametricQuadraticTetrahedron::QuadratureRules_4Points() 
{
	IP.Resize(4,3);
	gpe = 4;  // Gauss points per element for numerical integration
	const double vol (1./6.);
	IP(0,0) = 0.13819660, IP(0,1) = 0.13819660, IP(0,2) = 0.13819660; W[0] = 1./4. * vol;
	IP(1,0) = 0.58541020, IP(1,1) = 0.13819660, IP(1,2) = 0.13819660; W[1] = 1./4. * vol;
	IP(2,0) = 0.13819660, IP(2,1) = 0.58541020, IP(2,2) = 0.13819660; W[2] = 1./4. * vol;
	IP(3,0) = 0.13819660, IP(3,1) = 0.13819660, IP(3,2) = 0.58541020; W[3] = 1./4. * vol;
}

//MN: Note that Negative weights may cause numerical instabilities 
//MN: Zienkiewicz and Taylor page 177 ; The weiths in Rao Page 151 are wrong;
//MN: Suitable for the polynomial of degree 3 (Jinyun_CMAME_1984);
//MN: Tested OK;
void IsoparametricQuadraticTetrahedron::QuadratureRules_5Points() 
 {
	 double vol (1./6.);
	 IP.Resize(5,3);
	 gpe = 5;
	 const double64 a(1./4.), b(1./6.), c(1./2.);

	 IP(0,0) = a, IP(0,1) = a, IP(0,2) = a; W[0] = -4./5. * vol;
	 IP(1,0) = b, IP(1,1) = b, IP(1,2) = b; W[1] = 9./20. * vol;
	 IP(2,0) = b, IP(2,1) = b, IP(2,2) = c; W[2] = 9./20. * vol;
	 IP(3,0) = b, IP(3,1) = c, IP(3,2) = b; W[3] = 9./20. * vol;
	 IP(4,0) = c, IP(4,1) = b, IP(4,2) = b; W[4] = 9./20. * vol;
}

//MN: Suitable for the polynomial of degree 4 (Keast_CMAME_1986);
//MN: The  volume of element has already considered in the weights;
//MN: Tested OK;
void IsoparametricQuadraticTetrahedron::QuadratureRules_11Points() 
 {
	 IP.Resize(11,3);
	 gpe = 11;

	 IP(0,0) = 0.25, IP(0,1) = 0.25, IP(0,2) = 0.25; W[0] = -0.0131555555555555550;

	 const double64 a(0.0714285714285714285), b(0.7857142857142857142);
	 IP(1,0) = b, IP(1,1) = a, IP(1,2) = a; W[1] = 0.0076222222222222222;
	 IP(2,0) = a, IP(2,1) = b, IP(2,2) = a; W[2] = 0.0076222222222222222;
	 IP(3,0) = a, IP(3,1) = a, IP(3,2) = b; W[3] = 0.0076222222222222222;
	 IP(4,0) = a, IP(4,1) = a, IP(4,2) = a; W[4] = 0.0076222222222222222;
	 
	 double64 c(0.399403576166799219), d(0.100596423833200785);
	 IP(5,0) = c, IP(5,1) = c, IP(5,2) = d; W[5] = 0.0248888888888888888;
	 IP(6,0) = c, IP(6,1) = d, IP(6,2) = c; W[6] = 0.0248888888888888888;
	 IP(7,0) = c, IP(7,1) = d, IP(7,2) = d; W[7] = 0.0248888888888888888;
	 IP(8,0) = d, IP(8,1) = c, IP(8,2) = c; W[8] = 0.0248888888888888888;
	 IP(9,0) = d, IP(9,1) = c, IP(9,2) = d; W[9] = 0.0248888888888888888;
	 IP(10,0)= d, IP(10,1)= d, IP(10,2)= c; W[10]= 0.0248888888888888888;
}

//MN: Suitable for the polynomial of degree 4 (Jinyun_CMAME_1984);
void IsoparametricQuadraticTetrahedron::QuadratureRules_14Points() 
 {
	 IP.Resize(14,3);
	 gpe = 14;
	 const double64 a (1./2.);
	 IP(0,0) = a, IP(0,1) = a, IP(0,2) = 0; W[0] = 0.00317460317460317450;
	 IP(1,0) = a, IP(1,1) = 0, IP(1,2) = a; W[1] = 0.00317460317460317450;
	 IP(2,0) = a, IP(2,1) = 0, IP(2,2) = 0; W[2] = 0.00317460317460317450;
	 IP(3,0) = 0, IP(3,1) = a, IP(3,2) = a; W[3] = 0.00317460317460317450;
	 IP(4,0) = 0, IP(4,1) = a, IP(4,2) = 0; W[4] = 0.00317460317460317450;
	 IP(5,0) = 0, IP(5,1) = 0, IP(5,2) = a; W[5] = 0.00317460317460317450;

	 const double64 b (0.100526765225204467), c (0.698419704324386603);
	 IP(6,0)  = c, IP(6,1) = b, IP(6,2) = b; W[6] = 0.0147649707904967828;
	 IP(7,0)  = b, IP(7,1) = c, IP(7,2) = b; W[7] = 0.0147649707904967828;
	 IP(8,0)  = b, IP(8,1) = b, IP(8,2) = c; W[8] = 0.0147649707904967828;
	 IP(9,0)  = b, IP(9,1) = b, IP(9,2) = b; W[9] = 0.0147649707904967828;
	 
	 const double64 d (0.314372873493192195), e (0.0568813795204234229);
	 IP(10,0) = d, IP(10,1) = d, IP(10,2) = d; W[10] = 0.0221397911142651221;
	 IP(11,0) = d, IP(11,1) = d, IP(11,2) = e; W[11] = 0.0221397911142651221;
	 IP(12,0) = d, IP(12,1) = e, IP(12,2) = d; W[12] = 0.0221397911142651221;
	 IP(13,0) = e, IP(13,1) = d, IP(13,2) = d; W[13] = 0.0221397911142651221;
}

//MN: Suitable for the polynomial of degree 6 (Jinyun_CMAME_1984);
void IsoparametricQuadraticTetrahedron::QuadratureRules_24Points() 
 {
	 IP.Resize(24,3);
	 gpe = 24;
	 const double64 a (0.214602871259151684), b (0.356191386222544953);
	 IP(0,0) = b, IP(0,1) = a, IP(0,2) = a; W[0] = 0.00665379170969464506;
	 IP(1,0) = a, IP(1,1) = b, IP(1,2) = a; W[1] = 0.00665379170969464506;
	 IP(2,0) = a, IP(2,1) = a, IP(2,2) = b; W[2] = 0.00665379170969464506;
	 IP(3,0) = a, IP(3,1) = a, IP(3,2) = a; W[3] = 0.00665379170969464506;
	 
	 const double64 c (0.322337890142275646), d (0.0329863295731730594);
	 IP(4,0) = d, IP(4,1) = c, IP(4,2) = c; W[4] = 0.00922619692394239843;
	 IP(5,0) = c, IP(5,1) = d, IP(5,2) = c; W[5] = 0.00922619692394239843;
	 IP(6,0) = c, IP(6,1) = c, IP(6,2) = d; W[6] = 0.00922619692394239843;
	 IP(7,0) = c, IP(7,1) = c, IP(7,2) = c; W[7] = 0.00922619692394239843;

	 const double64 e (0.0406739585346113397), f (0.877978124396165982);
	 IP(8,0)  = f, IP(8,1) = e,  IP(8,2) = e;  W[8] = 0.00167953517588677620;
	 IP(9,0)  = e, IP(9,1) = f,  IP(9,2) = e;  W[9] = 0.00167953517588677620;
	 IP(10,0) = e, IP(10,1) = e, IP(10,2) = f; W[10] = 0.00167953517588677620;
	 IP(11,0) = e, IP(11,1) = e, IP(11,2) = e; W[11] = 0.00167953517588677620;
	 
	 const double64 g (0.0636610018750175299), h (0.269672331458315867), i (0.603005664791649076);
	 IP(12,0) = i, IP(12,1) = h, IP(12,2) = g; W[12] = 0.00803571428571428248;
	 IP(13,0) = i, IP(13,1) = g, IP(13,2) = h; W[13] = 0.00803571428571428248;
	 IP(14,0) = i, IP(14,1) = g, IP(14,2) = g; W[14] = 0.00803571428571428248;
	 IP(15,0) = h, IP(15,1) = i, IP(15,2) = g; W[15] = 0.00803571428571428248;
	 IP(16,0) = h, IP(16,1) = g, IP(16,2) = i; W[16] = 0.00803571428571428248;
	 IP(17,0) = h, IP(17,1) = g, IP(17,2) = g; W[17] = 0.00803571428571428248;
	 IP(18,0) = g, IP(18,1) = g, IP(18,2) = h; W[18] = 0.00803571428571428248;
	 IP(19,0) = g, IP(19,1) = g, IP(19,2) = i; W[19] = 0.00803571428571428248;
	//* 
	 IP(20,0) = g, IP(20,1) = h, IP(20,2) = g; W[20] = 0.00803571428571428248;
	 IP(21,0) = g, IP(21,1) = i, IP(21,2) = g; W[21] = 0.00803571428571428248;
	 IP(22,0) = g, IP(22,1) = h, IP(22,2) = i; W[22] = 0.00803571428571428248;
	 IP(23,0) = g, IP(23,1) = i, IP(23,2) = h; W[23] = 0.00803571428571428248;
	 //*/
}

//MN: Suitable for the polynomial of degree 8 (Jinyun_CMAME_1984);
void IsoparametricQuadraticTetrahedron::QuadratureRules_45Points() 
 {
	 IP.Resize(45,3);
	 gpe = 45;

	 IP(0,0) = 0.25, IP(0,1) = 0.25, IP(0,2) = 0.25; W[0] = 0.0393270066412926145;

	 const double64 a (0.127470936566639015), b (0.617587190300082967);
	 IP(1,0) = b, IP(1,1) = a, IP(1,2) = a; W[1] = 0.00408131605934270525;
	 IP(2,0) = a, IP(2,1) = b, IP(2,2) = a; W[2] = 0.00408131605934270525;
	 IP(3,0) = a, IP(3,1) = a, IP(3,2) = b; W[3] = 0.00408131605934270525;
	 IP(4,0) = a, IP(4,1) = a, IP(4,2) = a; W[4] = 0.00408131605934270525;
	 
	 const double64 c (0.0320788303926322960), d (0.903763508822103123);
	 IP(5,0) = d, IP(5,1) = c, IP(5,2) = c; W[5] = 0.000658086773304341943;
	 IP(6,0) = c, IP(6,1) = d, IP(6,2) = c; W[6] = 0.000658086773304341943;
	 IP(7,0) = c, IP(7,1) = c, IP(7,2) = d; W[7] = 0.000658086773304341943;
	 IP(8,0) = c, IP(8,1) = c, IP(8,2) = c; W[8] = 0.000658086773304341943;

	 const double64 e (0.0497770956432810185), f (0.450222904356718978);
	 IP(9,0)  = e, IP(9,1) = e,  IP(9,2) = f;  W[9]  = 0.00438425882512284693;
	 IP(10,0) = e, IP(10,1) = f, IP(10,2) = e; W[10] = 0.00438425882512284693;
	 IP(11,0) = e, IP(11,1) = f, IP(11,2) = f; W[11] = 0.00438425882512284693;
	 IP(12,0) = f, IP(12,1) = f, IP(12,2) = e; W[12] = 0.00438425882512284693;
	 IP(13,0) = f, IP(13,1) = e, IP(13,2) = f; W[13] = 0.00438425882512284693;
	 IP(14,0) = f, IP(14,1) = e, IP(14,2) = e; W[14] = 0.00438425882512284693;

	 const double64 g (0.183730447398549945), h (0.316269552601450060);
	 IP(15,0) = g, IP(15,1) = g, IP(15,2) = h; W[15] = 0.0138300638425098166;
	 IP(16,0) = g, IP(16,1) = h, IP(16,2) = g; W[16] = 0.0138300638425098166;
	 IP(17,0) = g, IP(17,1) = h, IP(17,2) = h; W[17] = 0.0138300638425098166;
	 IP(18,0) = h, IP(18,1) = h, IP(18,2) = g; W[18] = 0.0138300638425098166;
	 IP(19,0) = h, IP(19,1) = g, IP(19,2) = h; W[19] = 0.0138300638425098166; 
	 IP(20,0) = h, IP(20,1) = g, IP(20,2) = g; W[20] = 0.0138300638425098166;

	 const double64 i (0.231901089397150906), j(0.0229177878448171174), k(0.513280033360881072);
	 IP(21,0) = i, IP(21,1) = i, IP(21,2) = j; W[21] = 0.00424043742468372453;
	 IP(22,0) = i, IP(22,1) = i, IP(22,2) = k; W[22] = 0.00424043742468372453;
	 IP(23,0) = i, IP(23,1) = j, IP(23,2) = i; W[23] = 0.00424043742468372453;
	 IP(24,0) = i, IP(24,1) = k, IP(24,2) = i; W[24] = 0.00424043742468372453;
	 IP(25,0) = i, IP(25,1) = j, IP(25,2) = k; W[25] = 0.00424043742468372453;
	 IP(26,0) = i, IP(26,1) = k, IP(26,2) = j; W[26] = 0.00424043742468372453;
	 IP(27,0) = j, IP(27,1) = i, IP(27,2) = i; W[27] = 0.00424043742468372453;
	 IP(28,0) = k, IP(28,1) = i, IP(28,2) = i; W[28] = 0.00424043742468372453;
	 IP(29,0) = j, IP(29,1) = i, IP(29,2) = k; W[29] = 0.00424043742468372453;
	 IP(30,0) = k, IP(30,1) = i, IP(30,2) = j; W[30] = 0.00424043742468372453;
	 IP(31,0) = j, IP(31,1) = k, IP(31,2) = i; W[31] = 0.00424043742468372453;
	 IP(32,0) = k, IP(32,1) = j, IP(32,2) = i; W[32] = 0.00424043742468372453;
	 
	 const double64 l (0.0379700484718286102), m(0.730313427807538396), n(0.193746475248804382);
	 IP(33,0) = l, IP(33,1) = l, IP(33,2) = m; W[33] = 0.00223873973961420164;
	 IP(34,0) = l, IP(34,1) = l, IP(34,2) = n; W[34] = 0.00223873973961420164;
	 IP(35,0) = l, IP(35,1) = m, IP(35,2) = l; W[35] = 0.00223873973961420164;
	 IP(36,0) = l, IP(36,1) = n, IP(36,2) = l; W[36] = 0.00223873973961420164;
	 IP(37,0) = l, IP(37,1) = m, IP(37,2) = n; W[37] = 0.00223873973961420164;
	 IP(38,0) = l, IP(38,1) = n, IP(38,2) = m; W[38] = 0.00223873973961420164;
	 IP(39,0) = m, IP(39,1) = l, IP(39,2) = l; W[39] = 0.00223873973961420164;
	 IP(40,0) = n, IP(40,1) = l, IP(40,2) = l; W[40] = 0.00223873973961420164;
	 IP(41,0) = m, IP(41,1) = l, IP(41,2) = n; W[41] = 0.00223873973961420164;
	 IP(42,0) = n, IP(42,1) = l, IP(42,2) = m; W[42] = 0.00223873973961420164;
	 IP(43,0) = m, IP(43,1) = n, IP(43,2) = l; W[43] = 0.00223873973961420164;
	 IP(44,0) = n, IP(44,1) = m, IP(44,2) = l; W[44] = 0.00223873973961420164;
}








/** Computes the interpolation function values at the specified integration
point 'ip'.

@section arguments Input Arguments

A reference to the parent Element, the number of the integration point.

The interpolation function values are returned into the third argument.

*/
void IsoparametricQuadraticTetrahedron::N_AtIntegrationPoint( size_t ip, std::vector<double64>& IPOL )
 {
    assert( ip < gpe );

    // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), IPOL );

 } // end N_AtIntegrationPoint






void
IsoparametricQuadraticTetrahedron::NodesOfSegment( size_t segm_id, std::vector<size_t>& snids ) const
 {
    snids.resize(3);
    if ( segm_id == 0 ) {
         snids[0] = 0;
         snids[1] = 1;
         snids[2] = 4;
      }
    else if ( segm_id == 1 ) {
         snids[0] = 1;
         snids[1] = 2;
         snids[2] = 5;
      }
    else if ( segm_id == 2 ) {
         snids[0] = 2;
         snids[1] = 0;
         snids[2] = 6;
      }
    else if ( segm_id == 3 ) {
         snids[0] = 0;
         snids[1] = 3;
         snids[2] = 7;
      }
    else if ( segm_id == 4 ) {
         snids[0] = 1;
         snids[1] = 3;
         snids[2] = 8;
      }
    else if ( segm_id == 5 ) {
         snids[0] = 2;
         snids[1] = 3;
         snids[2] = 9;
      }
    else
    std::cout <<"\nIsoparametricQuadraticTetrahedron::NodesOfSegment: Erratic segment id requested: "<< segm_id << std::endl;

 } // end NodesOfSegment





/**

The faces are numbered such that face 0 lies opposite of node 0, face 1
node 1 etc. These faces are of the type QuadraticTriangle (TRI_6),
as is also present in library.

The nodes of each face are numbered counter-clockwise, looking from the
outside of the tetrahedron at each face.

@section arguments Input Arguments

The number of the desired face (0...3, 0-n-1).

Local node numbers are returned into the unsigned integer vector 'fnids'.

@section application Application

When boundary conditions shall be applied it is necessary to determine
the properties associated with the nodes of that element face.
Furthermore, NodesOfFace() is used in the construction of faces if
these are made part of the mesh connectivity by the
MeshManager.
 */
void IsoparametricQuadraticTetrahedron::NodesOfFace( size_t face_id,
                                                     vector<size_t>& fnids ) const
 {
    fnids.resize(6);
    if      ( face_id == 0 ) { // opposite node 0
         fnids[0] = 1;
         fnids[1] = 2;
         fnids[2] = 3;
         fnids[3] = 5;
         fnids[4] = 9;
         fnids[5] = 8;
      }
    else if ( face_id == 1 ) { // opposite node 1
         fnids[0] = 0;
         fnids[1] = 3;
         fnids[2] = 2;
         fnids[3] = 6;
         fnids[4] = 7;
         fnids[5] = 9;
      }
    else if ( face_id == 2 ) { // opposite node 2
         fnids[0] = 0;
         fnids[1] = 1;
         fnids[2] = 3;
         fnids[3] = 4;
         fnids[4] = 8;
         fnids[5] = 7;
      }
    else if ( face_id == 3 ) { // opposite node 3
         fnids[0] = 0;
         fnids[1] = 2;
         fnids[2] = 1;
         fnids[3] = 4;
         fnids[4] = 6;
         fnids[5] = 5;
      }
    else std::cerr <<"\nQuadraticTetrahedron::NodesOfFace: Face ID not identified: "<< face_id << std::endl;
 }




/** Returns local node ids of the 4 nodes located at the corners of
the quadratic tetrahedral element.

This is an unsigned integer vector with the 4 local corner node ID numbers
for the element.

*/
void  IsoparametricQuadraticTetrahedron::CornerNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(4);
    ids[0] = 0;
    ids[1] = 1;
    ids[2] = 2;
    ids[3] = 3;
 }




/** Returns 6 local node ids of the nodes 4-9 located at the midsides of
the quadratic tetrahedral element.

Returns an unsigned integer vector with the 4 local midside-node ID numbers
for the element.

*/
void  IsoparametricQuadraticTetrahedron::MidSideNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(6);
    ids[0] = 4;
    ids[1] = 5;
    ids[2] = 6;
    ids[3] = 7;
    ids[4] = 8;
    ids[5] = 9;
 }


void  IsoparametricQuadraticTetrahedron::CounterClockwiseNodes( std::vector<size_t>& ids ) const
 {
    ids.resize(npe);
    ids[0] = 0;
    ids[1] = 4;
    ids[2] = 1;
    ids[3] = 5;
    ids[4] = 2;
    ids[5] = 6;
    ids[6] = 7;
    ids[7] = 8;
    ids[8] = 9;
    ids[9] = 3;
 }


double64 IsoparametricQuadraticTetrahedron::WeightAtIntegrationPoint( size_t i ) const { return W[i]; }


CSMP_FEM_TYPE IsoparametricQuadraticTetrahedron::ElementTypeOfFace( size_t ) const
 {
    return ISOPARAMETRIC_QUADRATIC_TRIANGLE;
 }


void  IsoparametricQuadraticTetrahedron::JacobianAtIntegrationPoint( size_t gauss_point )
 {
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    Jacobian( DNR, DNS, DNT ); // 3D
 }


/** 
    Returns the (const) values of the element interpolation functions at the element barycenter.
    
    @test O.K. SKM retested, Aug. 2014
*/
void  IsoparametricQuadraticTetrahedron::N_AtBaryCenter( std::vector<double64>& IPOL )
 {
    IPOL.resize(npe);
    IPOL[0] = -0.125;
    IPOL[1] = -0.125;
    IPOL[2] = -0.125;
    IPOL[3] = -0.125;
    IPOL[4] =  0.25;
    IPOL[5] =  0.25;
    IPOL[6] =  0.25;
    IPOL[7] =  0.25;
    IPOL[8] =  0.25;
    IPOL[9] =  0.25;
 }



/**

Segment length is approximated as the sum of the length of two linear
segments making up each edge of the quadratic tetrahedron. The first segment
connects corner node 0 with midpoint node 4 and midpoint node 4 with
corner node 1 and so forth.

The lengths of the 6 segments will be returned into the vector which is
the second method argument.

@section implementation Implementation

The segment lengths are computed as the sum of the straight line
segments which make up the element face.
*/
void  IsoparametricQuadraticTetrahedron::EdgeLengths( std::vector<double64>& len )
{
    double64 sum;
    len.resize(spe);

    // segment 1
    sum     = (XY(4,0)-XY(0,0)) * (XY(4,0)-XY(0,0));
    sum    += (XY(4,1)-XY(0,1)) * (XY(4,1)-XY(0,1));
    sum    += (XY(4,2)-XY(0,2)) * (XY(4,2)-XY(0,2));
    len[0]  = sqrt(sum);
    sum     = (XY(1,0)-XY(4,0)) * (XY(1,0)-XY(4,0));
    sum    += (XY(1,1)-XY(4,1)) * (XY(1,1)-XY(4,1));
    sum    += (XY(1,2)-XY(4,2)) * (XY(1,2)-XY(4,2));
    len[0] += sqrt(sum);

    // segment 2
    sum     = (XY(5,0)-XY(1,0)) * (XY(5,0)-XY(1,0));
    sum    += (XY(5,1)-XY(1,1)) * (XY(5,1)-XY(1,1));
    sum    += (XY(5,2)-XY(1,2)) * (XY(5,2)-XY(1,2));
    len[1]  = sqrt(sum);
    sum     = (XY(2,0)-XY(5,0)) * (XY(2,0)-XY(5,0));
    sum    += (XY(2,1)-XY(5,1)) * (XY(2,1)-XY(5,1));
    sum    += (XY(2,2)-XY(5,2)) * (XY(2,2)-XY(5,2));
    len[1] += sqrt(sum);

    // segment 3
    sum     = (XY(6,0)-XY(2,0)) * (XY(6,0)-XY(2,0));
    sum    += (XY(6,1)-XY(2,1)) * (XY(6,1)-XY(2,1));
    sum    += (XY(6,2)-XY(2,2)) * (XY(6,2)-XY(2,2));
    len[2]  = sqrt(sum);
    sum     = (XY(0,0)-XY(6,0)) * (XY(0,0)-XY(6,0));
    sum    += (XY(0,1)-XY(6,1)) * (XY(0,1)-XY(6,1));
    sum    += (XY(0,2)-XY(6,2)) * (XY(0,2)-XY(6,2));
    len[2] += sqrt(sum);

    // segment 4
    sum     = (XY(7,0)-XY(0,0)) * (XY(7,0)-XY(0,0));
    sum    += (XY(7,1)-XY(0,1)) * (XY(7,1)-XY(0,1));
    sum    += (XY(7,2)-XY(0,2)) * (XY(7,2)-XY(0,2));
    len[3]  = sqrt(sum);
    sum     = (XY(3,0)-XY(7,0)) * (XY(3,0)-XY(7,0));
    sum    += (XY(3,1)-XY(7,1)) * (XY(3,1)-XY(7,1));
    sum    += (XY(3,2)-XY(7,2)) * (XY(3,2)-XY(7,2));
    len[3] += sqrt(sum);

    // segment 5
    sum     = (XY(8,0)-XY(1,0)) * (XY(8,0)-XY(1,0));
    sum    += (XY(8,1)-XY(1,1)) * (XY(8,1)-XY(1,1));
    sum    += (XY(8,2)-XY(1,2)) * (XY(8,2)-XY(1,2));
    len[4]  = sqrt(sum);
    sum     = (XY(3,0)-XY(8,0)) * (XY(3,0)-XY(8,0));
    sum    += (XY(3,1)-XY(8,1)) * (XY(3,1)-XY(8,1));
    sum    += (XY(3,2)-XY(8,2)) * (XY(3,2)-XY(8,2));
    len[4] += sqrt(sum);

    // segment 6
    sum     = (XY(9,0)-XY(2,0)) * (XY(9,0)-XY(2,0));
    sum    += (XY(9,1)-XY(2,1)) * (XY(9,1)-XY(2,1));
    sum    += (XY(9,2)-XY(2,2)) * (XY(9,2)-XY(2,2));
    len[5]  = sqrt(sum);
    sum     = (XY(3,0)-XY(9,0)) * (XY(3,0)-XY(9,0));
    sum    += (XY(3,1)-XY(9,1)) * (XY(3,1)-XY(9,1));
    sum    += (XY(3,2)-XY(9,2)) * (XY(3,2)-XY(9,2));
    len[5] += sqrt(sum);

 } // end EdgeLengths




/**

Computes the element aspect ratio as the ratio between the longest
and the shortest boundary segment. This rule differs from that presented
for the straight-sided triangle which is is discussed in the
LinearTriangle subclass.

@section arguments Input Arguments

The Element is consulted for its global coordinates.
*/
double64  IsoparametricQuadraticTetrahedron::AspectRatio()
{
   static vector<double64> vec(spe);

   EdgeLengths( vec );

   double64 seg_max(vec[0]), seg_min(vec[0]);

   // find largest segment
   for ( size_t i=1; i<spe; i++ ) {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
     }

   return seg_max / seg_min;
}




/**

This method attempts to calculate the inner radius of the curved-sided
triangle by talking 1/2 of the perimeter and dividing this measure by
the area of the triangle.

The procedure is empirically based giving a good match if the triangles
are relatively even-sided and have straight edges.

@section arguments Input Arguments

The parent element is queried for its node coordinates.
*/
double64  IsoparametricQuadraticTetrahedron::InnerRadius()
{
   static  vector<double64> segms(spe);
   double64                 sum(0.0);

   EdgeLengths( segms );
   for ( size_t i=0; i<spe; i++ ) sum += segms[i];

   return Volume() / (sum/2.);
}





/** Projection function from rst->xyz
*/
void IsoparametricQuadraticTetrahedron::ParametricToPhysical( vector<double64>& rst,
                                                              vector<double64>& xyz )
{
    Nrst( rst[0],rst[1],rst[2], DNR );

    xyz.resize(dim);
    xyz[0]=xyz[1]=xyz[2]=0.;

    for( size_t i=0; i<npe; i++ ) {
         xyz[0]+=XY(i,0)*DNR[i];
         xyz[1]+=XY(i,1)*DNR[i];
         xyz[2]+=XY(i,2)*DNR[i];
      }
}



/** Projection function from xyz->rst

Starts from the 2/3 2/3 0 point in  the parametric space
*/
void IsoparametricQuadraticTetrahedron::PhysicalToParametric(
                                       vector<double64>& rSt,
                                       const vector<double64>& xyz )
{

    vector<double64> outxyz(dim);
    vector<double64> rstHatK(dim);

    std::vector<double64> distanceFromGivenPointLinf(dim,0.0);
    double64 distanceFromGivenPointL2;

    // Find largest and smallest segments in order to define precision
    vector<double64> vec(spe);
    EdgeLengths( vec );
    double64 seg_max(vec[0]), seg_min(vec[0]);
    for ( size_t i=1; i<spe; i++ )
    {
        if ( vec[i] > seg_max ) seg_max = vec[i];
        if ( vec[i] < seg_min ) seg_min = vec[i];
    }

    const double64 geometricTolerance = 0.005*seg_min;

    // First guess as BaryCenter
    rstHatK[0] = 0.25;
    rstHatK[1] = 0.25;
    rstHatK[2] = 0.25;

    ParametricToPhysical( rstHatK, outxyz );
    // out( outxyz );

    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
    distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                     distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

    if( distanceFromGivenPointL2 > geometricTolerance )
    //if( (distanceFromGivenPointLinf[0] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[1] > geometricTolerance) ||
    //    (distanceFromGivenPointLinf[2] > geometricTolerance)  )
    {

        vector<double64> rstHatK_PlusOne(dim);
        double64 minDistanceFromGivenPoint;

        const double64 constantMu               = 1.0;
        const size_t numberOfFirstIterrations   = 5;
        const size_t maxNumberOfIterrations     = 20;
        const size_t numberOfIterationsWhenJacobiIsNotConstant = maxNumberOfIterrations ;

        const double64 incrementR = 1.0/(numberOfFirstIterrations-1);
        const double64 incrementS = 1.0/(numberOfFirstIterrations-1);
        const double64 incrementT = 1.0/(numberOfFirstIterrations-1);

        minDistanceFromGivenPoint = distanceFromGivenPointL2;

        rstHatK_PlusOne[0] = 0.0;
        for(size_t i=0;i<numberOfFirstIterrations;i++)
        {
            rstHatK_PlusOne[1] = 0.0;
            for(size_t j=0;j<numberOfFirstIterrations;j++)
            {
                rstHatK_PlusOne[2] = 0.0;
                for(size_t k=0;k<numberOfFirstIterrations;k++)
                {
                    ParametricToPhysical( rstHatK_PlusOne, outxyz);

                    distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
                    distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
                    distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );


                    distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                                     distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                                     distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

                    if( minDistanceFromGivenPoint > distanceFromGivenPointL2 )
                    {
                        for(size_t i=0; i<dim; i++)
                            rstHatK[i] = rstHatK_PlusOne[i];

                        minDistanceFromGivenPoint = distanceFromGivenPointL2;
                    }
                    rstHatK_PlusOne[2] += incrementT;
                }
                rstHatK_PlusOne[1] += incrementS;
            }
            rstHatK_PlusOne[0] += incrementR;
        }

        ParametricToPhysical( rstHatK, outxyz);

        distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
        distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
        distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );

        distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                         distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                         distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

        size_t iteration = 1;

        /// Newton-Raphson iterations
        while ( ( distanceFromGivenPointL2 > geometricTolerance ) && ( iteration < maxNumberOfIterrations ) )
        //while ( ( ( distanceFromGivenPointLinf[0] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[1] > geometricTolerance ) ||
        //          ( distanceFromGivenPointLinf[2] > geometricTolerance )    )
        //      &&  ( iteration < maxNumberOfIterrations )                )
        {
            // Out of range check
            if( (rstHatK[0]< 0.) || (rstHatK[0]>1.) ||
                (rstHatK[1]< 0.) || (rstHatK[0]>1.) ||
                (rstHatK[2]< 0.) || (rstHatK[1]>1.)  )
            {
                rstHatK[2] = 0.0;
                rstHatK[2] = 0.0;
                rstHatK[2] = 2.0;

                break;
            }

            // Inverse Jacobian calculations
            if( iteration < numberOfIterationsWhenJacobiIsNotConstant )
            {
                dNr( rstHatK[0], rstHatK[1], rstHatK[2], DNR );
                dNs( rstHatK[0], rstHatK[1], rstHatK[2], DNS );
                dNt( rstHatK[0], rstHatK[1], rstHatK[2], DNT );
                Jacobian( DNR, DNS, DNT );
                // Check whether Jacobian is positive ( might be not true for the point outside the element )
                //const double64 detJ = JacobianDeterminant();
                //if( detJ > 0.0 )
                JacobianInverse();
            }

            // Newton iteration
            rstHatK_PlusOne[0] = rstHatK[0] - constantMu*(JINV(0,0)*(outxyz[0]-xyz[0]) + JINV(1,0)*(outxyz[1]-xyz[1]) +JINV(2,0)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[1] = rstHatK[1] - constantMu*(JINV(0,1)*(outxyz[0]-xyz[0]) + JINV(1,1)*(outxyz[1]-xyz[1]) +JINV(2,1)*(outxyz[2]-xyz[2]));
            rstHatK_PlusOne[2] = rstHatK[2] - constantMu*(JINV(0,2)*(outxyz[0]-xyz[0]) + JINV(1,2)*(outxyz[1]-xyz[1]) +JINV(2,2)*(outxyz[2]-xyz[2]));

            for(size_t i=0; i<dim; i++)
                rstHatK[i] = rstHatK_PlusOne[i];

            ParametricToPhysical( rstHatK, outxyz);
            distanceFromGivenPointLinf[0] = std::abs( outxyz[0]-xyz[0] );
            distanceFromGivenPointLinf[1] = std::abs( outxyz[1]-xyz[1] );
            distanceFromGivenPointLinf[2] = std::abs( outxyz[2]-xyz[2] );
            distanceFromGivenPointL2 = sqrt( distanceFromGivenPointLinf[0]*distanceFromGivenPointLinf[0] +
                                             distanceFromGivenPointLinf[1]*distanceFromGivenPointLinf[1] +
                                             distanceFromGivenPointLinf[2]*distanceFromGivenPointLinf[2] );

            iteration++;

        }

        JAC.Zero();

        if( iteration == maxNumberOfIterrations )
        {
            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric: Tolerance = "
                <<geometricTolerance<<endl;

            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric: Real Point:\t"
                <<"x = "<<xyz[0]<<" ;\t"
                <<"y = "<<xyz[1]<<" ;\t"
                <<"z = "<<xyz[2]<<"\n";
            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric: Found Point:\t"
               <<"x = "<<outxyz[0]<<" ;\t"
               <<"y = "<<outxyz[1]<<" ;\t"
               <<"z = "<<outxyz[2]<<"\n";
            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric: Distance:\t"
               <<"x = "<<distanceFromGivenPointLinf[0]<<" ;"
               <<"y = "<<distanceFromGivenPointLinf[1]<<" ;"
               <<"z = "<<distanceFromGivenPointLinf[2]<<endl;
            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric Parametric Point:\t"
               <<"r = "<<rstHatK[0]<<" ;\t"
               <<"s = "<<rstHatK[1]<<" ;\t"
               <<"t = "<<rstHatK[2]<<"\n";

            std::vector<double64> N(npe,0.0);
            Nrst(rstHatK[0], rstHatK[1], rstHatK[2], N );
            cout<<" IsoparametricQuadraticTetrahedron::PhysicalToParametric: Shape functions:\t"
                <<"N[0] = "<<N[0]<<" ;\t"
                <<"N[1] = "<<N[1]<<" ;\t"
                <<"N[2] = "<<N[2]<<" ;\t"
                <<"N[3] = "<<N[3]<<" ;\t"
                <<"N[4] = "<<N[4]<<" ;\t"
                <<"N[5] = "<<N[5]<<" ;\t"
                <<"N[6] = "<<N[6]<<" ;\t"
                <<"N[7] = "<<N[7]<<" ;\t"
                <<"N[8] = "<<N[8]<<" ;\t"
                <<"N[9] = "<<N[9]<<"\n";

            csmp::Exception( WARNING, "IsoparametricQuadraticTetrahedron::PhysicalToParametric",
                            "Newton-Raphson iteration not converged");
        }

    }

    for(size_t i=0; i<dim; i++)
        rSt[i] = rstHatK[i];

}









/** Returns the value of the element interpolation functions at the points
'rst' in local coordinates.

@section arguments Input Arguments

The floating point coordinates 'r', 's' and 't' of the point at which the
interpolation shall be carried out.

@param IPOL The third method argument is the vector into which the values of the
n=nodes interpolation functions at the point 'rst' will be returned.

@section implementation Implementation

See source code header file.

@section application Application

Method is used to compute property values at the integration points of
the element.

*/
void IsoparametricQuadraticTetrahedron::Nrst( double64 L2, double64 L3, double64 L4,
                                              vector<double64>& IPOL ) const
{
   // note that L1 = 1 - r - s - t, L2 = r, L3 = s, L4 = t
   // see Bathe, page 256 (midside nodes) and Huyakorn and Pinder page 91 (corner nodes)
   const double64 L1 = 1.0 - L2 - L3 - L4;
   IPOL.resize(npe);
   IPOL[0] = L1 * ( 2.0 * L1 - 1.0 );
   IPOL[1] = L2 * ( 2.0 * L2 - 1.0 );
   IPOL[2] = L3 * ( 2.0 * L3 - 1.0 );
   IPOL[3] = L4 * ( 2.0 * L4 - 1.0 );
   IPOL[4] = 4.0 * L1 * L2;
   IPOL[5] = 4.0 * L2 * L3;
   IPOL[6] = 4.0 * L3 * L1;
   IPOL[7] = 4.0 * L1 * L4;
   IPOL[8] = 4.0 * L2 * L4;
   IPOL[9] = 4.0 * L3 * L4;
}

void IsoparametricQuadraticTetrahedron::Nrst( double64 L2, double64 L3, double64 L4,
                                              double64* IPOL ) const
{
   // note that L1 = 1 - r - s - t, L2 = r, L3 = s, L4 = t
   // see Bathe, page 256 (midside nodes) and Huyakorn and Pinder page 91 (corner nodes)
   const double64 L1 = 1.0 - L2 - L3 - L4;
   IPOL[0] = L1 * ( 2.0 * L1 - 1.0 );
   IPOL[1] = L2 * ( 2.0 * L2 - 1.0 );
   IPOL[2] = L3 * ( 2.0 * L3 - 1.0 );
   IPOL[3] = L4 * ( 2.0 * L4 - 1.0 );
   IPOL[4] = 4.0 * L1 * L2;
   IPOL[5] = 4.0 * L2 * L3;
   IPOL[6] = 4.0 * L3 * L1;
   IPOL[7] = 4.0 * L1 * L4;
   IPOL[8] = 4.0 * L2 * L4;
   IPOL[9] = 4.0 * L3 * L4;
}



/** This method and the complementary method dNs() compute the shape function
derivates with respect to the local coordinate axis.

@section arguments Input Arguments

The method takes the local coordinates of the point at which the
shape function derivatives shall be evaluated as third argument.

@param DNR The shape function derivatives are returned into the second method
argument which is a Meschpp floating point vector of the size n=nodes
per element.

@section implementation Implementation

see inlined source code in the header file.

@section application Application

The shape function derivatives are needed in most integration
procedures for elements.

*/
void IsoparametricQuadraticTetrahedron::dNr( double64 r, double64 s, double64 t,
                                             vector<double64>& DNR ) const
{
   const double64 nrst = 1. - r - s - t;
   DNR.resize(npe);
   DNR[0] =  1. - 4. * nrst;
   DNR[1] = -1. + 4. * r;
   DNR[2] =  0.;
   DNR[3] =  0.;
   DNR[4] = -4. * r + 4. * nrst;
   DNR[5] =  4. * s;
   DNR[6] = -4. * s;
   DNR[7] = -4. * t;
   DNR[8] =  4. * t;
   DNR[9] =  0.;
}



void IsoparametricQuadraticTetrahedron::dNs( double64 r, double64 s, double64 t,
                                             vector<double64>& DNS ) const
{
   const double64 nrst = 1. - r - s - t;
   DNS.resize(npe);
   DNS[0] =  1. - 4. * nrst;
   DNS[1] =  0.;
   DNS[2] = -1. + 4. * s;
   DNS[3] =  0.;
   DNS[4] = -4. * r;
   DNS[5] =  4. * r;
   DNS[6] = -4. * s + 4. * nrst;
   DNS[7] = -4. * t;
   DNS[8] =  0.;
   DNS[9] =  4. * t;

}



void IsoparametricQuadraticTetrahedron::dNt( double64 r, double64 s, double64 t,
                                             vector<double64>& DNT ) const
{
   const double64 nrst = 1. - r - s - t;
   DNT.resize(npe);
   DNT[0] =  1. - 4. * nrst;
   DNT[1] =  0.;
   DNT[2] =  0.;
   DNT[3] = -1. + 4. * t;
   DNT[4] = -4. * r;
   DNT[5] =  0.;
   DNT[6] = -4. * s;
   DNT[7] =  4. * nrst - 4. * t;
   DNT[8] =  4. * r;
   DNT[9] =  4. * s;

}





/** Compute derivatives of shape functions at corresponding nodes with
respect to the global coordinate system.

@section arguments Input Arguments

The element is used to obtain the global shape of the tetrahedron and the
interpolation function derivatives at each node are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

The interpolation-function derivative matrix is returned into the
second method argument.

*/
void IsoparametricQuadraticTetrahedron::dN( DenseMatrix<DM_MIN>& DN10 )
  {
     DN10.Resize(dim,npe);
     DN10 = DN;

     M.Resize(dim,1);

     // Jacobian transformation to global coordinate system
     for ( size_t i=0; i<npe; i++ )
       {
          for ( size_t j=0; j<dim; j++ ) RST[j] = NXYZ(i,j);

          // here the global coordinates come in
          dNr( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNR );
          dNs( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNS );
          dNt( NXYZ(i,0), NXYZ(i,1), NXYZ(i,2), DNT );

          // compute Jacobian matrix, its determinant and inversex
          Jacobian( DNR, DNS, DNT );
          JacobianInverse();

          for ( size_t j=0; j<dim; j++ ) M(j,0) = DN10(j,i);

          // 3x3 * 3x1 = 3x1 gives the global DN entries
          JINV *= M;
          DN10(0,i) = JINV(0,0);
          DN10(1,i) = JINV(1,0);
          DN10(2,i) = JINV(2,0);
          JINV.Resize(dim,dim);
       }

  } // end dN






/** Computes derivatives of shape functions at a point XY given in global coordinates
within the element (also given in global coordinates).

@section arguments Input Arguments

The element is used to obtain the global shape of the element and the
interpolation function derivatives at the point XY are returned into the
matrix M of dimensions rows = spatial dimensions x columns = nodes.

@return The interpolation-function derivative matrix is returned into the
second method argument.

*/
double64 IsoparametricQuadraticTetrahedron::dN_At( DenseMatrix<DM_MIN>& DN2, const vector<double64>& xyz  )
  {

    DNR.resize(dim);

    PhysicalToParametric(DNR,xyz);

    RST[0]=DNR[0];
    RST[1]=DNR[1];
    RST[2]=DNR[2];

    // here the global coordinates come in
    dNr( RST[0], RST[1], RST[2], DNR );
    dNs( RST[0], RST[1], RST[2], DNS );
    dNt( RST[0], RST[1], RST[2], DNT );
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    DN2.Resize(dim,dim);
    DN2  = JINV;

    DenseMatrix<DM_MIN>DN(dim,npe);
    for(int i=0;i<npe;i++){
        DN(0,i)=DNR[i];
        DN(1,i)=DNS[i];
        DN(2,i)=DNT[i];
    }

    DN2*=DN;

    return detJ;

 } // end dN





/**

Computes Integral N dA = 1/2 sum_1...n Wi Ji Ni. This formulation also
gives a meaningful volume if the element boundaries are curved. For a
description of numerical integration  of isoparametric quadratic
triangular elements, refer to Cook et al. (1989) 3rd Ed., p. 183.

@section arguments Input Arguments

The finite element which supplies the nodal coordinates from which
the area is computed.

@return The volume (m3) of the finite element.

*/
double64  IsoparametricQuadraticTetrahedron::Volume()
{
    double64  volume(0.);
    // numerical integration:
    // looping over the 4 Gauss points calculating determinant
    // test-function products and applying uniform weights
    for ( size_t i=0U; i<gpe; i++ )
      {
         dNr( IP(i,0), IP(i,1), IP(i,2), DNR );
         dNs( IP(i,0), IP(i,1), IP(i,2), DNS );
         dNt( IP(i,0), IP(i,1), IP(i,2), DNT );

         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix
         Jacobian( DNR, DNS, DNT );
         volume += JacobianInverse() * W[i];
      }

    return volume;
}








/**

Computes the values of the shape functions at the global coordinate 'xy'.
Since a Jacobian matrix can only be obtained for points which lie inside
the local triangle, such a transformation cannot be performed here.
Because of this reason global area-coordinate functions are employed
to obtain the shape function values at the global point. These
functions are borrowed from the NaturalQuadratic triangle.

Note, that the use of global shape functions implies that the method
only gives highly accurate results if the element boundaries are
not curved.

@section arguments Input Arguments

The parent Element, the vector which will hold the shape function
values, and the coordinate vector 'xy' in the global coordinate system.
Note that if the global coordinates given by 'xy' do not lie within the
global extent of the quadratic triangular element, the shape function values
will no longer add to one and the interpolation will then be erroneous.

The second method argument, the vector FN, will hold the shape function
values as computed at the point 'xy'.

@section implementation Implementation

The inverse of the Jacobian matrix is found for the global point 'xy'.
Then 'xy' transposed is pre-multiplied with the Jacobian to find the
local coordinate pair 'rs' corresponding to 'xy'. Using these local
coordinates the shape function values are found.

@section application Application

To interpolate a property value withing the quadratic triangular
element.

@todo SKM: does not find the correct location of the point xyz in local coordinates.

 */
void IsoparametricQuadraticTetrahedron::N( vector<double64>& N, const vector<double64>& xyz )
 {
    DNR.resize(dim);

    PhysicalToParametric(DNR, xyz);

    Nrst(DNR[0],DNR[1],DNR[2], N );

 } // end








/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the gauss
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.
*/
double64 IsoparametricQuadraticTetrahedron::dN_AtIntegrationPoint( DenseMatrix<DM_MIN>& B, size_t gauss_point )
 {
    // 1. compute local test-function derivative matrix at gauss point
    // get local shape function derivatives at Gauss point
    assert( gauss_point < gpe );
    dNr( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNR );
    dNs( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNS );
    dNt( IP(gauss_point,0), IP(gauss_point,1), IP(gauss_point,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8], B(2,8) = DNT[8];
    B(0,9) = DNR[9], B(1,9) = DNS[9], B(2,9) = DNT[9];

    B = JINV * B;

    return detJ;
 }




/**

1. Pre-multiplies local interpolation function derivative matrix, B, with
node coordinates to obtain the Jacobian matrix, J at the desired node
point.

2. The Jacobian matrix J is inverted and multiplied with B to obtain
the global interpolation function derivative matrix at the
point i (i=0...i=gauss points-1).

@section arguments Input Arguments

The first argument is a reference to the parent Element object,
the second argument is the result matrix, the third argument indicates
the gauss point at which the derivative matrix is computed and the
fourth argument specifies the degrees of freedom per node, for which
the B matrix shall be transformed.

@return The global interpolation function derivative matrix is returned into
the second method argument. The method also returns the determinant
of the Jacobian matrix since it is often needed in integration
procedures.
*/
double64 IsoparametricQuadraticTetrahedron::dN_AtNode( DenseMatrix<DM_MIN>& B, size_t nd )
 {
    assert( nd < npe );
    dNr( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNR );
    dNs( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNS );
    dNt( NXYZ(nd,0), NXYZ(nd,1), NXYZ(nd,2), DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8], B(2,8) = DNT[8];
    B(0,9) = DNR[9], B(1,9) = DNS[9], B(2,9) = DNT[9];

    B = JINV * B;

    return detJ;
 }


double64 IsoparametricQuadraticTetrahedron::dN_AtBarycenter( DenseMatrix<DM_MIN>& B )
 {
    double64 one_fourth(1.0/4.0);
    dNr( one_fourth, one_fourth, one_fourth, DNR );
    dNs( one_fourth, one_fourth, one_fourth, DNS );
    dNt( one_fourth, one_fourth, one_fourth, DNT );

    // compute Jacobian matrix, its determinant and inversex
    Jacobian( DNR, DNS, DNT );
    double64 detJ = JacobianInverse();

    // compose matrix DN = 3 x 10 in global coordinates
    // by multiplication of JINV with local DN
    B.Resize(dim,npe);
    B(0,0) = DNR[0], B(1,0) = DNS[0], B(2,0) = DNT[0];
    B(0,1) = DNR[1], B(1,1) = DNS[1], B(2,1) = DNT[1];
    B(0,2) = DNR[2], B(1,2) = DNS[2], B(2,2) = DNT[2];
    B(0,3) = DNR[3], B(1,3) = DNS[3], B(2,3) = DNT[3];
    B(0,4) = DNR[4], B(1,4) = DNS[4], B(2,4) = DNT[4];
    B(0,5) = DNR[5], B(1,5) = DNS[5], B(2,5) = DNT[5];
    B(0,6) = DNR[6], B(1,6) = DNS[6], B(2,6) = DNT[6];
    B(0,7) = DNR[7], B(1,7) = DNS[7], B(2,7) = DNT[7];
    B(0,8) = DNR[8], B(1,8) = DNS[8], B(2,8) = DNT[8];
    B(0,9) = DNR[9], B(1,9) = DNS[9], B(2,9) = DNT[9];

    B = JINV * B;

    return detJ;
 }





/**

Method returns a vector with a size of 3, containing the consecutively
ordered local node numbers of the side of the element which lies at the
indicated model boundary. If the element is not on the model boundary an
error is reported and the vector is initialized to unspecified.

@section arguments Input Arguments

The parent Element, the target boundary of the the Model of elements,
and a Meschach vector which will hold the the local indices of the identified
nodes (0...5).

@param fnids The resulting local node id's are returned into the third method argument.

@section implementation Implementation

While only the element knows which nodes are located at the mode boundary,
the FiniteElement knows in which order these appear.

@section application Application

To assign Neumann boundary conditions with a PDE operator for surface
integrals.
 */
void  IsoparametricQuadraticTetrahedron::ConsecutiveNodesAtBoundary( const vector<size_t>& bnodes,
                                                                     vector<size_t>& fnids )
 {
    fnids.resize(bnodes.size());

     if ( bnodes.size() != 6 )
       throw csmp::Exception( ERROR, "IsoparametricQuadraticTetrahedron::ConsecutiveNodesAtBoundary",
               "Cannot resolve node sequence for element boundary",
               "Probably because element lies at two boundaries simultaneously" );


 } // end ConsecutiveNodesAtBoundary







/** Outputs supplied node values to user-specified output file.
 
Element gets displayed as  VTK_QUADRATIC_TETRAHEDRON=24
 
@section arguments Input Arguments

The parent element, the output file name to which the extension ".vtk"
will be appended automatically, the name of the variable,
a data matrix which will contain n-colums=nodes and n-rows = dimensions
of data (n-rows=1 = scalar data, n-rows=3 = vector data etc.).

@section application Application

Single elements are output to VTK as polygons in order to test the
quality of interpolation and the computed properties directly.

@section messages Messages

The method will indicate if there is a problem in opening the output
file.

@test O.K. SKM 15/7/14

*/
void IsoparametricQuadraticTetrahedron::OutputNodeDataToVTK( const char* file_name,
                                                             const char* var_name,
                                                             DenseMatrix<DM_MIN>& DATA ) const
  {
     char  outfile[NAME_STRING], elmt[30];
     strcpy( outfile, file_name );
     sprintf( elmt, "%lu", CurrentID() );
     strcat( outfile, elmt );
     strcat( outfile, ".vtk" );

     // 0. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       {
           cout <<"\nIsoparametricQuadraticTetrahedron::OutputNodeDataToVTK ";
           cout <<"Output file could not be opened."<< endl;
           return;
       }

     // 1. writing the file header
     // --------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;

     // 2. writing node coordinates
     // ---------------------------
     DenseMatrix<DM_MIN> COORD(XY);
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     ofs <<"POINTS " << npe <<" float"<< endl;
     for ( size_t i=0; i<npe; i++ )
       {
          for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
          ofs << endl;
       }
     ofs << endl;

     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< 1 <<" "<< 11 << endl;
     // NOTE: as in VTK spec sheet
     ofs << 10 <<" 0 1 2 3 4 5 6 7 8 9" << endl;
     ofs << endl;

     // 4. writing CELL_TYPES
     // ---------------------
     // VTK_QUADRATIC_TETRA
     ofs <<"CELL_TYPES "<< 1 << endl;
     ofs << 24 << endl; // VTK_TETRAEDER
     ofs << endl;

     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< npe << endl;
     ofs.setf( ios::scientific );

     if ( DATA.Rows() == 1 )
       {
           ofs <<"SCALARS "<< var_name <<" float"<< endl;
           ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
           // matrix DATA is 1x9
           for ( size_t i=0; i<DATA.Cols(); i++ ) ofs << DATA(0,i) <<" ";
           ofs << endl;
       }
     else
       {
          ofs <<"VECTORS "<< var_name <<" float"<< endl;
          // variables have always 3 components since view screen is 3D
          // matrix DATA is vec-dim x 10
          for ( size_t i=0; i<DATA.Cols(); i++ ) {
               for ( size_t j=0; j<DATA.Rows(); j++ ) ofs << DATA(j,i) <<"  ";
               ofs << endl;
            }
       }
     ofs << endl;
     ofs.close();
     cout <<"\nIsoparametricQuadraticTetrahedron::OutputNodeDataToVTK: file '"<< outfile <<"' written successfully."<< endl;

  } // end OutputNodeDataToVTK









/**

Method uses the four integration points of the element to define a set of linear
interpolation functions (just like in the linear tetrahedron). With these
functions variable values specified for the integration points are
linearily extrapolated to the elements nodes. Since the extrapolation
functions are linear and the distance between the integration points and
the nodes (in the local coordinate system which is used) is small, the
extrapolation is exact to numerical precision.

@section arguments Input Arguments

The first argument specifies how many variable values are supplied by
the IVAR vector for each integration point. This allows, for instance,
to extrapolate vector or tensor variables with this method which are
put into IVAR sequentially. Thus, IVAR has the dimensions 4 x vars-per-
integration point. The third method argument NVAR will hold the
extrapolated values for each node. Thus, its size should be
10 x vars-per-integration point.

@param IVAR the integration point values that are to be extrapolated to nodes

@param NVAR the extrapolated values for each node of the element are returned into
the third method argument NVAR.

@section application Application

Property values calculated at integration points are exact 
in contrast to those calculated at nodes. Since the element
interpolation function derivatives of this quadratic element are linear 
functions, linear interpolation of integration point values across the
element is better than interpolation using its quadratic interpolation functions
which may oscillate.
Such linear interpolation is facilitated by this method. 

*/
void  IsoparametricQuadraticTetrahedron::ExtrapolateIntegrationPointVariableToNodes(
                                                                     size_t nvars,
                                                                     const vector<double64>& IVAR,
                                                                     vector<double64>&       NVAR )
 const
{
   static double64  a[4], b[4], c[4], d[4], intpol[4], volume6;
   static bool      first_call(true);
   size_t i;
   int32  j;

   vector<double64>  sum(nvars);

   // 0. Decide which case is dealt with in terms of the integration points
   //    which are used (rr and ss contain the integr.p. locations)
   if ( IVAR.size() != (gpe*nvars) )
     throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
                           "Input vector must have 'nvars' x 4 entries");

   NVAR.resize( npe * nvars );

   // 1. Compute the local interpolation function coefficients for the tetrahedron which
   //    is defined by the four integration points.
   if ( first_call ) {
       // checking starting conditions
       if ( gpe != 4 )
         throw csmp::Exception( FATAL_ERROR, "IsoparametricQuadraticTetrahedron::ExtrapolateIntegrationPointVariableToNodes",
                               "This method expects four integration points on which extrapolation functions will be based on" );

        // test function coefficients based on corner nodes only
        for ( i=0, j=1; i<4; i++ )
          {
             // a(i)
             a[i]  = -IP(n(i,1),0) * (IP(n(i,3),1)*IP(n(i,2),2)-IP(n(i,3),2)*IP(n(i,2),1));
             a[i] -=  IP(n(i,2),0) * (IP(n(i,1),1)*IP(n(i,3),2)-IP(n(i,1),2)*IP(n(i,3),1));
             a[i] -=  IP(n(i,3),0) * (IP(n(i,2),1)*IP(n(i,1),2)-IP(n(i,2),2)*IP(n(i,1),1));
             // b(i)
             b[i]  = IP(n(i,3),1)*IP(n(i,2),2) - IP(n(i,3),2)*IP(n(i,2),1);
             b[i] += IP(n(i,1),1)*IP(n(i,3),2) - IP(n(i,1),2)*IP(n(i,3),1);
             b[i] += IP(n(i,2),1)*IP(n(i,1),2) - IP(n(i,2),2)*IP(n(i,1),1);
             // c(i)
             c[i]  = IP(n(i,3),2)*IP(n(i,2),0) - IP(n(i,3),0)*IP(n(i,2),2);
             c[i] += IP(n(i,1),2)*IP(n(i,3),0) - IP(n(i,1),0)*IP(n(i,3),2);
             c[i] += IP(n(i,2),2)*IP(n(i,1),0) - IP(n(i,2),0)*IP(n(i,1),2);
             // d(i)
             d[i]  = IP(n(i,3),0)*IP(n(i,2),1) - IP(n(i,3),1)*IP(n(i,2),0);
             d[i] += IP(n(i,1),0)*IP(n(i,3),1) - IP(n(i,1),1)*IP(n(i,3),0);
             d[i] += IP(n(i,2),0)*IP(n(i,1),1) - IP(n(i,2),1)*IP(n(i,1),0);
             //
             a[i] *=  static_cast<double64>(j);
             b[i] *=  static_cast<double64>(j);
             c[i] *=  static_cast<double64>(j);
             d[i] *=  static_cast<double64>(j);
             j    *= -1;
          }
        // computing local element volume x 6
        volume6 = a[0] + a[1] + a[2] + a[3];

        first_call = false;
     }

   // 2. For each node point compute the values of the linear extrapolation functions
   //    and use these to extrapolate the values of the variables at the nodes.
   for ( i=0; i<npe; i++ )
     {
        // compute interpolation function values at node i
        intpol[0] = (a[0] + b[0] * NXYZ(i,0) + c[0] * NXYZ(i,1) + d[0] * NXYZ(i,2)) / volume6;
        intpol[1] = (a[1] + b[1] * NXYZ(i,0) + c[1] * NXYZ(i,1) + d[1] * NXYZ(i,2)) / volume6;
        intpol[2] = (a[2] + b[2] * NXYZ(i,0) + c[2] * NXYZ(i,1) + d[2] * NXYZ(i,2)) / volume6;
        intpol[3] = (a[3] + b[3] * NXYZ(i,0) + c[3] * NXYZ(i,1) + d[3] * NXYZ(i,2)) / volume6;

        // carry out extrapolation
        fill( sum.begin(), sum.end(), static_cast<double64>(0.) );
        for ( size_t l=0; l<gpe; l++ )
          for ( size_t k=0; k<nvars; k++ ) sum[k] += intpol[l] * IVAR[l*nvars + k];

        // store result in output vector
        for ( size_t k=0; k<nvars; k++ ) NVAR[i*nvars + k] = sum[k];
     }

} // end ExtrapolateIntegrationPointVariableToNodes (STL vectors)




/// integration point location transformed into global coordinates
/// @warning the matrix XYZ must be uptodate
void  IsoparametricQuadraticTetrahedron::IntegrationPoint( size_t ip,
                                                           vector<double64>& xyz ) const
 {
    assert( ip < gpe );
    xyz.resize(3U); xyz[0]=xyz[1]=xyz[2]=0.;
     // local interpolation function values
    Nrst( IP(ip,0), IP(ip,1), IP(ip,2), NRST );

    for( size_t i=0U; i<npe; i++ ) {
          xyz[0] += XY(i,0) * NRST[i];
          xyz[1] += XY(i,1) * NRST[i];
          xyz[2] += XY(i,2) * NRST[i];
      }

 } // end IntegrationPoint




// AP 2006
void IsoparametricQuadraticTetrahedron::ReferenceCoordinates( DenseMatrix<DM_MIN>& matCoords ) const
{
    matCoords.Resize(npe, dim);
    matCoords = NXYZ;
}



void IsoparametricQuadraticTetrahedron::JacobianAt( const std::vector<double64>& rst )
{
   dNr(rst[0], rst[1], rst[2], DNR);
   dNs(rst[0], rst[1], rst[2], DNS);
   dNt(rst[0], rst[1], rst[2], DNT);

   Jacobian(DNR, DNS, DNT);
}


/**
    Computes the unit-normal to the target face at its barycentre.
    
    @note the qualification 'Barycenter' is important because the surface of the face may be warped
    in this isoparametric element.

    @author HA 2016
*/
void IsoparametricQuadraticTetrahedron::UnitNormalAtFaceBarycenter( size_t face, std::vector<double64>& nrmlAtbarycenter )
{
   assert(face < Faces());
   nrmlAtbarycenter.resize(3U);
   // if face lies opposite to node 0
   if (face == 0) { 
      vector<double64> rstAtBarycenter(3, 1. / 3.);  // defining coordinates of barycenter i at parametric space
      vector<double64> nrml_rst(3, 0.577350269189626); // normal vector in parametric space
      
      JacobianAt(rstAtBarycenter);  
      JacobianInverse();
      nrmlAtbarycenter = JINV*nrml_rst;
      double64 nrmlAtbarycenterLength = sqrt(nrmlAtbarycenter[0] * nrmlAtbarycenter[0] + 
         nrmlAtbarycenter[1] * nrmlAtbarycenter[1] + nrmlAtbarycenter[2] * nrmlAtbarycenter[2]);
      nrmlAtbarycenter[0] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[1] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[2] /= nrmlAtbarycenterLength;
      return;
   }
   

   // if face lies opposite to node 1
   if (face == 1) { 
      vector<double64> rstAtBarycenter(3, 1. / 3.);  // defining coordinates of barycenter i at parametric space
      rstAtBarycenter[0] = 0.;
      vector<double64> nrml_rst(3, 0.);  // normal vector in parametric space
      nrml_rst[0] = -1.;
      
      JacobianAt(rstAtBarycenter);
      JacobianInverse();
      nrmlAtbarycenter = JINV*nrml_rst;
      double64 nrmlAtbarycenterLength = sqrt(nrmlAtbarycenter[0] * nrmlAtbarycenter[0] +
         nrmlAtbarycenter[1] * nrmlAtbarycenter[1] + nrmlAtbarycenter[2] * nrmlAtbarycenter[2]);
      nrmlAtbarycenter[0] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[1] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[2] /= nrmlAtbarycenterLength;
      return;
   }


   // if face lies opposite to node 2
   if (face == 2) { 
      vector<double64> rstAtBarycenter(3, 1. / 3.);  // defining coordinates of barycenter i at parametric space
      rstAtBarycenter[1] = 0.;
      vector<double64> nrml_rst(3, 0.);  // normal vector in parametric space
      nrml_rst[1] = -1.;
      
      JacobianAt(rstAtBarycenter);
      JacobianInverse();
      nrmlAtbarycenter = JINV*nrml_rst;
      double64 nrmlAtbarycenterLength = sqrt(nrmlAtbarycenter[0] * nrmlAtbarycenter[0] +
         nrmlAtbarycenter[1] * nrmlAtbarycenter[1] + nrmlAtbarycenter[2] * nrmlAtbarycenter[2]);
      nrmlAtbarycenter[0] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[1] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[2] /= nrmlAtbarycenterLength;
      return;
   }


   // if face lies opposite to node 3
   if (face == 3) { 
      vector<double64> rstAtBarycenter(3, 1. / 3.);  // defining coordinates of barycenter i at parametric space
      rstAtBarycenter[2] = 0.;
      vector<double64> nrml_rst(3, 0.);  // normal vector in parametric space
      nrml_rst[2] = -1.;

      JacobianAt(rstAtBarycenter);
      JacobianInverse();
      nrmlAtbarycenter = JINV*nrml_rst;
      double64 nrmlAtbarycenterLength = sqrt(nrmlAtbarycenter[0] * nrmlAtbarycenter[0] +
         nrmlAtbarycenter[1] * nrmlAtbarycenter[1] + nrmlAtbarycenter[2] * nrmlAtbarycenter[2]);
      nrmlAtbarycenter[0] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[1] /= nrmlAtbarycenterLength;
      nrmlAtbarycenter[2] /= nrmlAtbarycenterLength;
      return;
   }
  
} // end UnitNormalAtFaceBarycenter



/**
    returns the barycenter coordinates of IsoparametricQuadraticTetrahedron in physical space
    
    @author HA 2016
*/
void IsoparametricQuadraticTetrahedron::FaceBarycenterCoordinates( size_t face, std::vector<double64>& barycenterCoord )
{
   assert(face < Faces());
   barycenterCoord.resize(3U);
   array<double64,10> N_barycenterFace = N_AtFaceBarycenter(face);
   // getting the coordinate at barycenter of face 0
   for (size_t i = 0; i < IntegrationPoints(); ++i ) {
        barycenterCoord[0] = N_barycenterFace[i] * XYZ(i, 0);
        barycenterCoord[1] = N_barycenterFace[i] * XYZ(i, 1);
        barycenterCoord[2] = N_barycenterFace[i] * XYZ(i, 2);
     }
   
}



/** 
     returns the value of interpolation functions for a quadratic tetrahedron at the barycenter of face i
    
    @author HA 2016
*/
array<double64,10> IsoparametricQuadraticTetrahedron::N_AtFaceBarycenter( size_t i ) const
{
   assert( i < Faces() );
   const size_t nodesPerElement(10U);
  
   if (i == 0)
     return array < double64, nodesPerElement >{ {-1.110223024625156e-16, -0.111111111111111, -0.111111111111111,
       -0.111111111111111, 1.480297366166875e-16, 1.480297366166875e-16,
       1.480297366166875e-16, 0.444444444444444, 0.444444444444444, 0.444444444444444} };

   if (i == 1)
     return array < double64, nodesPerElement >{ {-0.111111111111111, 0., -0.111111111111111,
       -0.111111111111111, 0., 0.4444444444444444,
       0.4444444444444444, 0., 0.444444444444444, 0.} };

   if (i == 2)
     return array < double64, nodesPerElement > { {-0.111111111111111, -0.111111111111111, 0.,
       -0.111111111111111, 0.444444444444444, 0.,
       0.444444444444444, 0., 0., 0.444444444444444}};

//   if (i == 3) the only option left
     return array < double64, nodesPerElement > { {-0.111111111111111, -0.111111111111111, -0.111111111111111,
       0, 0.444444444444444, 0.444444444444444,
       0., 0.444444444444444, 0., 0.}};

}  // end N_AtFaceBarycenter


} // end namespace csmp


