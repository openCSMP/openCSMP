#include "ErrorFunction_Test.h"

namespace csmp{

void ErrorFunction_Test::run()
{
  // TESTEE
  ErrorFunction errorFunction;

  // ERF
  _equal( errorFunction.Erf( 0.   ), 0.,        1.E-7 );
  _equal( errorFunction.Erf( 0.05 ), 0.0563720, 1.E-7 );
  _equal( errorFunction.Erf( 0.10 ), 0.1124629, 1.E-7 );
  _equal( errorFunction.Erf( 0.30 ), 0.3286268, 1.E-7 );
  _equal( errorFunction.Erf( 0.40 ), 0.4283924, 1.E-7 );
  _equal( errorFunction.Erf( 0.50 ), 0.5204999, 1.E-7 );
  _equal( errorFunction.Erf( 0.80 ), 0.7421010, 1.E-7 );
  _equal( errorFunction.Erf( 1.00 ), 0.8427008, 1.E-7 );
  _equal( errorFunction.Erf( 1.20 ), 0.9103140, 1.E-7 );
  _equal( errorFunction.Erf( 1.60 ), 0.9763484, 1.E-7 );
  _equal( errorFunction.Erf( 1.90 ), 0.9927904, 1.E-7 );
  _equal( errorFunction.Erf( 2.40 ), 0.9993115, 1.E-7 );
  _equal( errorFunction.Erf( 2.90 ), 0.9999589, 1.E-7 );
  _equal( errorFunction.Erf( 3.50 ), 0.9999993, 1.E-7 );

  // ERFC
  _equal( errorFunction.Erfc( 0.   ), 1.,        1.E-7 );
  _equal( errorFunction.Erfc( 0.05 ), 0.9436280, 1.E-7 );
  _equal( errorFunction.Erfc( 0.10 ), 0.8875371, 1.E-7 );
  _equal( errorFunction.Erfc( 0.30 ), 0.6713732, 1.E-7 );
  _equal( errorFunction.Erfc( 0.40 ), 0.5716076, 1.E-7 );
  _equal( errorFunction.Erfc( 0.50 ), 0.4795001, 1.E-7 );
  _equal( errorFunction.Erfc( 0.80 ), 0.2578990, 1.E-7 );
  _equal( errorFunction.Erfc( 1.00 ), 0.1572992, 1.E-7 );
  _equal( errorFunction.Erfc( 1.20 ), 0.0896860, 1.E-7 );
  _equal( errorFunction.Erfc( 1.60 ), 0.0236516, 1.E-7 );
  _equal( errorFunction.Erfc( 1.90 ), 0.0072096, 1.E-7 );
  _equal( errorFunction.Erfc( 2.40 ), 0.0006885, 1.E-7 );
  _equal( errorFunction.Erfc( 2.90 ), 0.0000411, 1.E-7 );
  _equal( errorFunction.Erfc( 3.50 ), 0.0000007, 1.E-7 );

}

} // csmp
