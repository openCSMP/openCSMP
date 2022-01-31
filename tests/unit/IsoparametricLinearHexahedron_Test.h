//
//  IsoparametricLinearHexahedron_Test.h
//
//  Created by Hossein Agheshlui on 30/09/2015.
//  Copyright (c) 2015 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_ISOPARAMETRICLINEARHEXAHEDRON_TEST_H
#define CSMP_ISOPARAMETRICLINEARHEXAHEDRON_TEST_H

#include "Test.h"
#include "CSMP_definitions.h"
namespace csmp {

class IsoparametricLinearHexahedron_Test : public Test {
	public:
		explicit IsoparametricLinearHexahedron_Test( bool verbose );
		~IsoparametricLinearHexahedron_Test();
    
		virtual void run();
		double HexaVolUsingCSMP(double *Xcord, double *Ycord, double *Zcord, double &HexaVolumeUsingFVSectors, double &HexaVolumeUsingFVSectorsIPWeight_);
		double HexaVolUsingTetCSMP(double *Xcord, double *Ycord, double *Zcord);
		double HexaVolAnalyticalSolution(double *Xcord, double *Ycord, double *Zcord);
		double TetVolFromCSMP(double *Xcord, double *Ycord, double *Zcord, int, int, int, int);
		double AnalyticalHexaVolume(double *Xcord, double *Ycord, double *Zcord);
		double AnalyticalTetVol(double *X, double *Y, double *Z, int p0, int p1, int p2, int p3);
		//double HexaVolUsingFVSectors(double *X, double *Y, double *Z);

	private:
		double *Xcord, *Ycord, *Zcord;
		double HexaVolumeUsingFVSectors_;
		double HexaVolumeUsingFVSectorsIPWeight_;
		bool verbose_;

};

} // end csmp

#endif /* defined(CSMP_ISOPARAMETRICLINEARHEXAHEDRON_TEST_H) */
