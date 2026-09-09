#include "ExponentialTransferFunction_Test.h"

#include "CSMP_definitions.h"
#include "ExponentialFractureMatrixTransfer.h"

using namespace std;

namespace csmp {

void ExponentialTransferFunction_Test::run()
{
    cout << "\n\nExponential Fracture Matrix Transfer Test...\n\n";
    const double TOLERANCE( 1.E-2), SWM( 0.5 ), SWI( 0.05 ), SWR( 0.05 ), SNR( 0.14 ), LAMBDA( 2.5 ),
                 PD( 3000 ), KM( 1.E-15 ), PHIM( 0.25 ), MUW( 1.E-3 ), MUN( 5.E-3 ), SWF( 0.5 );
    ExponentialFractureMatrixTransfer FMTfct( SWI, SWR, SNR, LAMBDA, PD, KM, PHIM, MUW, MUN );

    _equal( FMTfct.CurrentTransferRate( 0., SWM, SWF ), numeric_limits<double>::max(), TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 60., SWM, SWF ), 6.31873E-05, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 600., SWM, SWF ), 1.99816E-05, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 6000., SWM, SWF ), 6.31873E-05, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 12000., SWM, SWF ), 4.46802E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 18000., SWM, SWF ), 3.64812E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 24000., SWM, SWF ), 3.15937E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 30000., SWM, SWF ), 2.82582E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 36000., SWM, SWF ), 2.57961E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 42000., SWM, SWF ), 2.38826E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 48000., SWM, SWF ), 2.23401E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 54000., SWM, SWF ), 2.10624E-06, TOLERANCE );
    _equal( FMTfct.CurrentTransferRate( 60000., SWM, SWF ), 6.31873E-06, TOLERANCE );
}

} //csmp
