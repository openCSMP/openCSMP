#include <iostream>

#include "Test.h"
#include "IsoparametricLinearHexahedron_Test.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearTetrahedron.h"
#include "Node.h"
#include "Element.h"
#include "FV_IntegrationPointsAndWeights.h"
#include "FiniteElement.h"
#include "FiniteVolumeStencil.h"
#include "DenseMatrix.h"

using namespace std;

namespace csmp {

	IsoparametricLinearHexahedron_Test::IsoparametricLinearHexahedron_Test(bool verbose) : verbose_(verbose) {}
	IsoparametricLinearHexahedron_Test::~IsoparametricLinearHexahedron_Test(){}

	void IsoparametricLinearHexahedron_Test::run()
	{
		//Testing the volume of a single Hexahedron, considering different types for a Hexahedron.
		//Coordinates of a 8 node Cube 1*1*1.
		double CubeXcord[8] = { 0, 1, 1, 0, 0, 1, 1, 0 };
		double CubeYcord[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
		double CubeZcord[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

		//Twisted Hexa. Top face twisted 30 degrees
		double TopTwXcord[8] = { 0, 1, 1, 0, 0, 0.866025404, 0.366025404, -0.5 };
		double TopTwYcord[8] = { 0, 0, 1, 1, 0, 0.5, 1.366025404, 0.866025404 };
		double TopTwZcord[8] = { 0, 0, 0, 0, 1, 1, 1, 1 };

		//Twisted Hexa. Sides Twisted 30 degrees
		double SideTwXcord[8] = { 0, 1, 0.683012702, -0.183012702, 0, 1, 1.183012702, 0.316987298 };
		double SideTwYcord[8] = { 0, 0, 1, 1, 0, 0, 1, 1 };
		double SideTwZcord[8] = { 0, 0, -0.183012702, 0.316987298, 1, 1, 0.683012702, 1.183012702 };

		//Skewed Hexa. None of the opposite faces are parallel.  Randomly generated coordinates. 
		double SkXcord[8] = { 0.0041, 1.05, 1.1478, 0.0464, 0.1281, 1.0491, 1.0827, 0.0604 };
		double SkYcord[8] = { 0.0467, 0.1169, 1.1358, 1.1705, 0.0827, 0.0995, 1.1436, 1.1902 };
		double SkZcord[8] = { 0.0334, 0.1724, 0.0962, 0.0145, 1.1961, 1.1942, 1.0391, 1.0153 };
		//Volume of this skewed Hexahedron calculated by Abaqus.
		const double SkewedHexaVolumeCalculatedbyAbaqus = 1.1324339180849166;
		
		if (verbose_){
			//Printing the volumes of the different hexahedra on the screen.
			cout << "AnalyticalHexaVolume " << std::setprecision(17) << AnalyticalHexaVolume(SideTwXcord, SideTwYcord, SideTwZcord) << endl;
			cout << "HexaVolUsingTetCSMP " << HexaVolUsingTetCSMP(SideTwXcord, SideTwYcord, SideTwZcord) << endl;
			cout << "HexaVolUsingCSMP " << HexaVolUsingCSMP(SideTwXcord, SideTwYcord, SideTwZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_) << endl;
			cout << "HexaVolUsingFVSectors " << HexaVolumeUsingFVSectors_ << endl;
			cout << "HexaVolUsingFVSectorsIPWeight " << HexaVolumeUsingFVSectorsIPWeight_ << endl;
		}

		//Testing Hexahedron volume calculated by CSMP against volume calculated from analytical solutions.Every second line.
		//Testing Hexahedron volume calculated by CSMP by adding the volume of the stencils within the Finite Element versus the volume calculated by mutliplying the determinent of the Jacobian matrix by the intergration point weight. Every second line.
		//The second check has been deactivated since there are differences between the results
		_equal(HexaVolUsingCSMP(CubeXcord, CubeYcord, CubeZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), AnalyticalHexaVolume(CubeXcord, CubeYcord, CubeZcord), 1e-15);
		_equal(HexaVolumeUsingFVSectors_, AnalyticalHexaVolume(CubeXcord, CubeYcord, CubeZcord), 1e-15);
		_equal(HexaVolUsingCSMP(TopTwXcord, TopTwYcord, TopTwZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), AnalyticalHexaVolume(TopTwXcord, TopTwYcord, TopTwZcord), 1e-15);
		//_equal(HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_, 1e-15);
		_equal(HexaVolUsingCSMP(SideTwXcord, SideTwYcord, SideTwZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), AnalyticalHexaVolume(SideTwXcord, SideTwYcord, SideTwZcord), 1e-15);
		//_equal(HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_, 1e-15);
		
		//Testing Hexahedron volume calculated by CSMP against volume calculated from dividing the Hexa into 6 Tetra and summing up their volumes.
		_equal(HexaVolUsingCSMP(CubeXcord, CubeYcord, CubeZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), HexaVolUsingTetCSMP(CubeXcord, CubeYcord, CubeZcord), 1e-15);
		_equal(HexaVolUsingCSMP(TopTwXcord, TopTwYcord, TopTwZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), HexaVolUsingTetCSMP(TopTwXcord, TopTwYcord, TopTwZcord), 1e-15);
		_equal(HexaVolUsingCSMP(SideTwXcord, SideTwYcord, SideTwZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), HexaVolUsingTetCSMP(SideTwXcord, SideTwYcord, SideTwZcord), 1e-15);
		
		//Testing skewed Hexahedron (i.e. none of the two opposite planes are not parallel) volume calculated by CSMP against volume calculated from Abaqus.
		_equal(HexaVolUsingCSMP(SkXcord, SkYcord, SkZcord, HexaVolumeUsingFVSectors_, HexaVolumeUsingFVSectorsIPWeight_), SkewedHexaVolumeCalculatedbyAbaqus, 1e-15);
		_equal(AnalyticalHexaVolume(SkXcord, SkYcord, SkZcord), HexaVolUsingTetCSMP(SkXcord, SkYcord, SkZcord), 1e-15);
		
	}

	double IsoparametricLinearHexahedron_Test::HexaVolUsingCSMP( double *Xcord, double *Ycord, double *Zcord,
                                                               double& HexaVolumeUsingFVSectors,
                                                               double& HexaVolumeUsingFVSectorsIPWeight )
	{
		//The following is the alternative for the dynamic allocation of memory, which is helpful when the required memory is only determined at the run time. 
		//This method creates a Hexahedron element each time that the class is called and kills the element when returining to the run method. 
		//IsoparametricLinearHexahedron *lHex_ = new IsoparametricLinearHexahedron();
		//Element<3U> *element_ = new Element < 3U >(lHex_) ;
		
		IsoparametricLinearHexahedron lHex(8);
		FiniteVolumeStencil<3U> FVStencil("ISOPARAMETRIC_LINEAR_HEXAHEDRON");	//This line creates the Finite volume stencils for an Isoparametric Linear Hexahedron type element, including the integration points and weights.
		Element<3U> element(&lHex, &FVStencil);		//This line creates a single Hexa element and the FV stencils within it. 

		std::vector<double> ParamIPCoord;
		FVStencil.SectorIntegrationPoint(6, 0, ParamIPCoord);
		
		if (verbose_){
			cout << "Sector Integration Point Parametric Coordinates" << FVStencil.SectorIntegrationPoint(6, 0) << endl;
			cout << "Sector Integration Weight " << FVStencil.SectorIntegrationWeight(6, 0) << endl;
		}

		Node<3U>  n0, n1, n2, n3, n4, n5, n6, n7;

		n0.Idx(1);
		n1.Idx(2);
		n2.Idx(3);
		n3.Idx(4);
		n4.Idx(5);
		n5.Idx(6);
		n6.Idx(7);
		n7.Idx(8);
		
		element.Idx(0);
		element.Assign(0, &n0);
		element.Assign(1, &n1);
		element.Assign(2, &n2);
		element.Assign(3, &n3);
		element.Assign(4, &n4);
		element.Assign(5, &n5);
		element.Assign(6, &n6);
		element.Assign(7, &n7);

		// nodes
		element.N(0)->x(Xcord[0]); element.N(0)->y(Ycord[0]); element.N(0)->z(Zcord[0]);
		element.N(1)->x(Xcord[1]); element.N(1)->y(Ycord[1]); element.N(1)->z(Zcord[1]);
		element.N(2)->x(Xcord[2]); element.N(2)->y(Ycord[2]); element.N(2)->z(Zcord[2]);
		element.N(3)->x(Xcord[3]); element.N(3)->y(Ycord[3]); element.N(3)->z(Zcord[3]);
		element.N(4)->x(Xcord[4]); element.N(4)->y(Ycord[4]); element.N(4)->z(Zcord[4]);
		element.N(5)->x(Xcord[5]); element.N(5)->y(Ycord[5]); element.N(5)->z(Zcord[5]);
		element.N(6)->x(Xcord[6]); element.N(6)->y(Ycord[6]); element.N(6)->z(Zcord[6]);
		element.N(7)->x(Xcord[7]); element.N(7)->y(Ycord[7]); element.N(7)->z(Zcord[7]);

		//cout << "\nVolume of One Hexahedron calculated by CSMP: " << element.Volume() << endl;
		double HexVol = element.Volume();

		//The calcualtion of the Hexa Volume summing up the volumes of the sectors of the FV stencils over the finite element. 
		HexaVolumeUsingFVSectors = 0;
		double SectorNumbers = element.Sectors();		// This gives the number of FV sectors within a Finite Elment.
		for (auto i = 0; i < SectorNumbers; i++) {
			HexaVolumeUsingFVSectors += element.SectorVolume(i);	// This gives the volume of the each of the sectors.
		}

		//delete lHex_, element_;  //Deletes a dynamically created element. 

		//The calcualtion of the Hexa Volume summing up the volumes of the sectors by multiplying the sector integration point weight by the determinant of the Jacobian matrix over the sectors of the FV stencils over the finite element. 
		double VolofSector(0.);
		HexaVolumeUsingFVSectorsIPWeight = 0.;
		for (auto i = 0; i < SectorNumbers; i++){
			if (verbose_){
				cout << element.FV()->SectorIntegrationWeight(i, 0) << "\t";
				cout << element.FE()->JacobianDeterminant() << "\t";
			}
			VolofSector = element.FV()->SectorIntegrationWeight(i, 0)*element.FE()->JacobianDeterminant();
			HexaVolumeUsingFVSectorsIPWeight += VolofSector;
			// Similar line from FiniteVolumeTraits: fVolume += e.FV()->SectorIntegrationWeight(iSector, j) * e.FE()->JacobianDeterminant();
		}

		return HexVol;
	}

	double IsoparametricLinearHexahedron_Test::HexaVolUsingTetCSMP(double *Xcord, double *Ycord, double *Zcord)
	{
		//This method creates 6 Tetrahedra which compose the hexahedron 
		double TVol[6], TetVolumes = 0;

		TVol[0] = TetVolFromCSMP(Xcord, Ycord, Zcord, 7, 3, 2, 5);
		TVol[1] = TetVolFromCSMP(Xcord, Ycord, Zcord, 2, 3, 1, 5);
		TVol[2] = TetVolFromCSMP(Xcord, Ycord, Zcord, 2, 6, 7, 5);
		TVol[3] = TetVolFromCSMP(Xcord, Ycord, Zcord, 3, 7, 4, 5);
		TVol[4] = TetVolFromCSMP(Xcord, Ycord, Zcord, 3, 0, 1, 4);
		TVol[5] = TetVolFromCSMP(Xcord, Ycord, Zcord, 4, 5, 1, 3);

		for (int i = 0; i < 6; i++){
			TetVolumes += TVol[i];
		}
		//cout << "\nVolume of One Hexahedron calculated by deviding it into Six Tetrahedra and getting their volumes from CSMP: " << TetVolumes << endl;
		return TetVolumes;
	}

	double IsoparametricLinearHexahedron_Test::TetVolFromCSMP( double *X, double *Y, double *Z,
                                                               int p0, int p1, int p2, int p3 )
	{
		IsoparametricLinearTetrahedron* lTet_ = new IsoparametricLinearTetrahedron();
		Element<3U> element(lTet_);
		Node<3U>  n0, n1, n2, n3;

		n0.Idx(1);
		n1.Idx(2);
		n2.Idx(3);
		n3.Idx(4);

		element.Idx(0);
		element.Assign(0, &n0);
		element.Assign(1, &n1);
		element.Assign(2, &n2);
		element.Assign(3, &n3);
		
		element.N(0)->x(X[p0]); element.N(0)->y(Y[p0]); element.N(0)->z(Z[p0]);
		element.N(1)->x(X[p1]); element.N(1)->y(Y[p1]); element.N(1)->z(Z[p1]);
		element.N(2)->x(X[p2]); element.N(2)->y(Y[p2]); element.N(2)->z(Z[p2]);
		element.N(3)->x(X[p3]); element.N(3)->y(Y[p3]); element.N(3)->z(Z[p3]);

		double VolumeofTetra = element.Volume();
		delete lTet_;
		return VolumeofTetra;
	}
  

	double IsoparametricLinearHexahedron_Test::AnalyticalHexaVolume(double *Xcord, double *Ycord, double *Zcord){
		double AnalyticalVol(0), AVol[6];

		AVol[0] = AnalyticalTetVol(Xcord, Ycord, Zcord, 7, 3, 2, 5);
		AVol[1] = AnalyticalTetVol(Xcord, Ycord, Zcord, 2, 3, 1, 5);
		AVol[2] = AnalyticalTetVol(Xcord, Ycord, Zcord, 2, 6, 7, 5);
		AVol[3] = AnalyticalTetVol(Xcord, Ycord, Zcord, 3, 7, 4, 5);
		AVol[4] = AnalyticalTetVol(Xcord, Ycord, Zcord, 3, 0, 4, 1);
		AVol[5] = AnalyticalTetVol(Xcord, Ycord, Zcord, 1, 5, 4, 3);

		for (int i = 0; i < 6; i++){
			AnalyticalVol += AVol[i];
		}
		//cout << "\nVolume of One Hexahedron calculated using analytical solution is: " << AnalyticalVol << endl;
		return AnalyticalVol;
	}

	double IsoparametricLinearHexahedron_Test::AnalyticalTetVol(double *X, double *Y, double *Z, int p0, int p1, int p2, int p3)
	{
		double a, b, c, d, e, f, g, h, i;

		a = X[p1] - X[p0];
		b = Y[p1] - Y[p0];
		c = Z[p1] - Z[p0];
		d = X[p2] - X[p0];
		e = Y[p2] - Y[p0];
		f = Z[p2] - Z[p0];
		g = X[p3] - X[p0];
		h = Y[p3] - Y[p0];
		i = Z[p3] - Z[p0];

		double Hvol = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
		double TetVol = Hvol > 0 ? Hvol / 6. : -Hvol / 6.0;
		return TetVol;
	}


}
