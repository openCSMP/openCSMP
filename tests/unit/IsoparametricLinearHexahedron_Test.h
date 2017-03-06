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
		double64 HexaVolUsingCSMP(double64 *Xcord, double64 *Ycord, double64 *Zcord, double64 &HexaVolumeUsingFVSectors, double64 &HexaVolumeUsingFVSectorsIPWeight_);
		double64 HexaVolUsingTetCSMP(double64 *Xcord, double64 *Ycord, double64 *Zcord);
		double64 HexaVolAnalyticalSolution(double64 *Xcord, double64 *Ycord, double64 *Zcord);
		double64 TetVolFromCSMP(double64 *Xcord, double64 *Ycord, double64 *Zcord, int, int, int, int);
		double64 AnalyticalHexaVolume(double64 *Xcord, double64 *Ycord, double64 *Zcord);
		double64 AnalyticalTetVol(double64 *X, double64 *Y, double64 *Z, int p0, int p1, int p2, int p3);
		//double64 HexaVolUsingFVSectors(double64 *X, double64 *Y, double64 *Z);

	private:
		double64 *Xcord, *Ycord, *Zcord;
		double64 HexaVolumeUsingFVSectors_;
		double64 HexaVolumeUsingFVSectorsIPWeight_;
		bool verbose_;

};

} // end csmp

#endif /* defined(CSMP_ISOPARAMETRICLINEARHEXAHEDRON_TEST_H) */
