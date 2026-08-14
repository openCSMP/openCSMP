#include "EOS_CO2H2ONaCl_Spycher2004.h"

using namespace std;

namespace csmp {

/**
      //  Spycher et al (2003),
      //  Co2-H2O mixtures in the geological sequestration of co2.
      //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
      //  Geochimica et Cosmochimica acta, vol 67, No 16
*/
EOS_CO2H2ONaCl_Spycher04::EOS_CO2H2ONaCl_Spycher04()
    : R ( 8.314472 ), // gas constant in m3 Pa K-1 mol-1
      stoichio ( 2 ), // stoichiometric number of ions contained in the dissolved salt,
      // it is taken 2 for NaCl, (Spycher, 2004)
      a_h2oco2 ( 7.89e+0 ), //  Pa m6 K0.5 molñ2
      b_co2 ( 2.780E-5 ), //  m3/mol
      b_h2o ( 1.818e-5 ), //  m3/mol
      b_mix ( 2.78e-5 ), //  m3/mol

      p0(1.0e5),       // reference pressure P0 = 1 bar, Spychler 2004
      vH2o( 18.1e-6 ), // im m3/mol taken from literature Spychler 2003,
      vCo2( 32.6e-6 ), // im m3/mol taken from literature Spychler 2003,

      molarMassH2o ( 18.01528e-3),  // Kilograms per mole
      molarMassCo2 ( 44.010e-3 ),    // Kilograms per mole
      molarMassNacl ( 58.443e-3 ),    // Kilograms per mole

      // standard condition
      pSC( 100000. ), // Standard pressure: SPE
      tSC(15.) // standard temperature:SPE
{
}




EOS_CO2H2ONaCl_Spycher04::~EOS_CO2H2ONaCl_Spycher04()
{
}



/// acos function for complex numbers specifically implemented because it gives problems with GNU and Intel compilers
complex<double> EOS_CO2H2ONaCl_Spycher04::complex_acos(const complex<double>& x)
{
    complex<double> y = std::log(x + std::sqrt(x*x - 1.0));
    if (y.imag() < 0.0)
        return complex<double>(-y.imag(), y.real());
    return complex<double>(y.imag(), -y.real());
}



/**
    @brief returns the compressed volume of CO2.

    Spycher et al (2003),
    Co2-H2O mixtures in the geological sequestration of co2.
    I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    Geochimica et Cosmochimica acta, vol 67, No 16

    R.W.D. Nickalls (1993),
    A new approach to solving the cubic: Cardan's solution revealed,
    The mathematical gazette, Vol 77

    T in C
    P in Pa
    compressedVolumeCo2 in m3/mol

    Range T:  12 < T <100 C
    Range P:  P<600 bar

*/
double EOS_CO2H2ONaCl_Spycher04::CompressedVolumeCo2( double pressure, double temperature  )
{
    assert( temperature >= 10. ); // should be 12oC, below which hydrates may come into play
    assert( temperature <= 100. );
    assert( pressure <= 6.0e7 );

    complex<double> delta(0.,0.), lambda(0.,0.), theta(0.,0.), h(0.,0.), complex2RealDummy;

    const double temp(degreeCToKelvin(temperature)), coeff_a(1.);

    //  coeff_b = - R * T / P;
    double coeff_b = R;
    coeff_b *=  temp;
    coeff_b /=  pressure;
    coeff_b *=  -1.;

    //  coeff_c = -( R * T * b_co2 / P ) + a_mix / (P * sqrt( T ) ) - pow(b_co2,2);
    double coeff_c = R;
    coeff_c *=  temp;
    coeff_c *=  b_co2;
    coeff_c /=  pressure;
    coeff_c *=  -1.;
    double dummy1  = a_mix(temperature);
    dummy1  /=  pressure;
    double dummy2  = sqrt( temp );
    dummy1  /=  dummy2;
    coeff_c +=  dummy1;
    dummy1  =   b_co2 * b_co2;
    coeff_c -=  dummy1;

    //  coeff_d = - a_mix * b_co2 /  (P * sqrt( T ));
    double coeff_d = a_mix(temperature);
    coeff_d *=  b_co2;
    coeff_d /=  pressure;
    dummy1  =   sqrt( temp );
    coeff_d /=  dummy1;
    coeff_d *=  -1.;

    //  x_N     = - coeff_b / ( 3. * coeff_a );
    double x_N = coeff_b;
    x_N /=  3.;
    x_N /=  coeff_a;
    x_N *=  -1;

    //  y_N     = pow(x_N,3) * coeff_a +  pow(x_N,2) * coeff_b + x_N * coeff_c + coeff_d;
    double y_N = x_N * x_N * x_N;
    y_N     *=  coeff_a;
    dummy1  =   x_N * x_N;
    dummy1  *=  coeff_b;
    y_N     +=  dummy1;
    dummy1  =   x_N;
    dummy1  *=  coeff_c;
    y_N     +=  dummy1;
    y_N     +=  coeff_d;

    //  test = ( pow(coeff_b,2) - 3. * coeff_a * coeff_c) / (9. *  pow(coeff_a,2) );
    double test = coeff_b * coeff_b;
    dummy1  =   -3.;
    dummy1  *=  coeff_a;
    dummy1  *=  coeff_c;
    test    +=  dummy1;
    test    /=  9.;
    dummy1  =   coeff_a * coeff_a;
    test    /=  dummy1;

    // test =delta^2
    if ( test < 0. ) {
        const  complex<double>& worzel = test;
        delta = sqrt( worzel );
    }
    else
        // delta = sqrt( ( pow(coeff_b,2) - 3. * coeff_a * coeff_c) / (9. *  pow(coeff_a,2) ));
        delta = sqrt( test );

    //  lambda  = sqrt( 3. *  pow(delta,2) );
    lambda  = delta*delta;
    lambda *= 3.;
    lambda = sqrt( lambda );

    //  h  = 2. *  pow(delta,3);
    h  =  delta*delta*delta;
    h *=  2.;

    double yNSquare(y_N*y_N);
    complex2RealDummy = h*h;
    double hSquare = complex2RealDummy.real();

    valarray<double> root(0.0,1);

    if ( yNSquare > hSquare ) {
        //    aux1 = 0.5 * ( -y_N + sqrt( yNSquare - hSquare ) );
        double aux1 = yNSquare;
        aux1 -= hSquare;
        aux1 =  sqrt( aux1 );
        aux1 -= y_N;
        aux1 *= 0.5;

        //    aux2 = 0.5 * ( -y_N - sqrt( yNSquare - hSquare ) );
        double aux2 =  yNSquare;
        aux2 -= hSquare;
        aux2 =  sqrt( aux2 );
        aux2 *= -1.;
        aux2 -= y_N;
        aux2 *= 0.5;
      
        double alpha(numeric_limits<double>::quiet_NaN());
        if ( aux1 < 0. && aux2 < 0 ) {
            // alpha = x_N - pow( -aux1, 1.0/3. ) - pow( -aux2, 1.0/3. );
            alpha   =   x_N;
            dummy1  =   cbrt( -aux1 );
            alpha   -=  dummy1;
            dummy1  =   cbrt( -aux2 );
            alpha   -=  dummy1;
        }
        else if ( aux1 < 0. && aux2 > 0 ) {
            // alpha = x_N - pow( -aux1, 1.0/3. ) + pow( aux2, 1.0/3. );
            alpha   =   x_N;
            dummy1  =   cbrt( -aux1 );
            alpha   -=  dummy1;
            dummy1  =   cbrt( aux2 );
            alpha   +=  dummy1;
        }
        else if ( aux1 > 0. && aux2 < 0 ) {
            // alpha = x_N + pow( aux1, 1.0/3. ) - pow( -aux2, 1.0/3. );
            alpha   =   x_N;
            dummy1  =   cbrt( aux1 );
            alpha   +=  dummy1;
            dummy1  =   cbrt( -aux2 );
            alpha   -=  dummy1;
        }
        else {
            // alpha = x_N + pow( aux1, 1.0/3. ) + pow( aux2, 1.0/3. );
            alpha   =   x_N;
            dummy1  =   cbrt( aux1 );
            alpha   +=  dummy1;
            dummy1  =   cbrt( aux2 );
            alpha   +=  dummy1;
        }

    root[0] = alpha;
    
    }
    else if ( yNSquare < hSquare ) {
        // SKM FIX - pi as in Robert's original equation
        constexpr double PI(3.14159265359);
        //    theta  = acos ( -y_N / h ) / 3.;
        theta  =  -y_N;
        theta /=  h ;
        theta  =  complex_acos( theta );
        theta /=  3.;

        //    complex2RealDummy = delta * cos( theta);
        complex2RealDummy =   cos( theta);
        complex2RealDummy *=  delta;

        //    gamma1            = x_N + 2. * complex2RealDummy.real();
        double gamma1  =   complex2RealDummy.real();
        gamma1  *=  2.;
        gamma1  +=  x_N;

        //    complex2RealDummy = delta * cos( 2. * PI / 3. + theta);
        complex2RealDummy =   2.;
        complex2RealDummy *=  PI;
        complex2RealDummy /=  3.;
        complex2RealDummy +=  theta;
        complex2RealDummy =   cos( complex2RealDummy );
        complex2RealDummy *=  delta ;

        //    gamma2            = x_N + 2. * complex2RealDummy.real();
        double gamma2  = complex2RealDummy.real();
        gamma2  *=  2.;
        gamma2  +=  x_N;

        //    complex2RealDummy = delta * cos( 4. * PI / 3. + theta);
        complex2RealDummy =   4.;
        complex2RealDummy *=  PI;
        complex2RealDummy /=  3.;
        complex2RealDummy +=  theta;
        complex2RealDummy =   cos( complex2RealDummy );
        complex2RealDummy *=  delta ;

        //    gamma3            = x_N + 2. * complex2RealDummy.real();
        double gamma3  = complex2RealDummy.real();
        gamma3  *= 2.;
        gamma3  += x_N;

        root.resize(3);
        root[0] = gamma1;
        root[1] = gamma2;
        root[2] = gamma3;
    }
    else {
        root.resize(1);
        double beta1, beta2(2.);
        if (h != 0.) {
            // delta = pow( y_N / ( 2. * coeff_a), 1. / 3.);
            delta =   y_N;
            delta /=  2.;
            delta /=  coeff_a;
            delta  =  pow( delta, 1. / 3. );

            // beta1 = x_N + delta.real();
            beta1 = delta.real();
            beta1 +=  x_N;

            // beta2 = x_N - 2. * delta.real();
            beta2 *=  delta.real();
            beta2 *=  -1;
            beta2 +=  x_N;
        }
        else {
            beta1 = x_N;
            beta2 = x_N;
        }

        root.resize(2);
        root[0] = beta1;
        root[1] = beta2;
    }

    double V_gas = root.max();
    double V_liq = root.min();

    //  work1 = P * (V_gas - V_liq);
    double work1 = V_gas;
    work1 -=  V_liq;
    work1 *=  pressure ;

    //  work2 = R * T * log( (V_gas - b_co2) / ( V_liq - b_co2) ) + a_mix / ( sqrt(T) * b_co2 ) * log( ( V_gas + b_co2) / ( V_liq + b_co2) * V_liq / V_gas );
    double work2  =  V_gas;
    work2   -=  b_co2;
    dummy1  =   V_liq;
    dummy1  -=  b_co2;
    work2   /=  dummy1;
    work2   =   log( work2 );
    work2   *=  R;
    work2   *=  temp;
    dummy1  =   V_gas;
    dummy1  +=  b_co2;
    dummy1  *=  V_liq;
    dummy1  /=  V_gas;
    dummy2  =   V_liq ;
    dummy2  +=  b_co2;
    dummy1  /=  dummy2;
    dummy1  =   log( dummy1 );
    dummy1  *=  a_mix(temperature);
    dummy2  =   sqrt(temp);
    dummy1  /=  dummy2;
    dummy1  /=  b_co2;
    work2   +=  dummy1;

    if ( (work2 - work1) < 0. ) return V_liq;
    if ( (work2 - work1) > 0. ) return V_gas;

    // V_root = V_gas
    return V_gas;
}



/**
    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    // T in C
    // P in Pa
    // phaseVolumeCo2 in m3/mol
    // fugacityCo2 without unit

    // y_H2o=0

    // Range T:  12 < T <100 C
    // Range P:  P<600 bar
*/
double EOS_CO2H2ONaCl_Spycher04::FugacityCo2( double pressure,
                                                double temperature,
                                                double phaseVolumeCo2 )
{
    assert( temperature >= 10. );
    assert( temperature <= 100. );
    assert( pressure <= 6.0e7 );

    const double temp(degreeCToKelvin(temperature));

    //  A1 = log ( V_root / ( V_root - b_mix ) );
    double A1      =   phaseVolumeCo2;
    double dummy1  =   phaseVolumeCo2;
    dummy1  -=  b_mix;
    A1      /=  dummy1;
    A1      =   log( A1 );

    //  A2 = b_co2 / ( V_root - b_mix );
    double A2 =  b_co2;
    A2 /= dummy1;

    //  A3 = - 2 * ( y_co2 * a_co2  ) / ( R * pow( T, 1.5 ) * b_mix ) * log ( ( V_root + b_mix ) / V_root );
    double A3      =   -2.;
    // A3      *=  y_co2;
    A3      *=  a_Co2(temperature);
    dummy1  =   pow( temp, 1.5 );
    dummy1  *=  R;
    dummy1  *=  b_mix;
    A3      /=  dummy1;
    dummy1  =   phaseVolumeCo2;
    dummy1  +=  b_mix;
    double dummy2  =   dummy1;
    dummy1  /=  phaseVolumeCo2;
    dummy1  =   log( dummy1 );
    A3      *=  dummy1;

    //  A4 = ( a_mix * b_co2 / ( R * pow( T, 1.5 ) * pow( b_mix,2 )  ) * ( log ( (V_root + b_mix) / V_root ) - b_mix / (V_root + b_mix) ) );
    double A4      =   a_mix(temperature);
    A4      *=  b_co2;
    dummy1  =   R;
    double dummy3  =   pow( temp, 1.5 );
    dummy1  *=  dummy3;
    dummy3  =   pow( b_mix,2 );
    dummy1  *=  dummy3;
    A4      /=  dummy1;
    dummy1  =   dummy2;
    dummy1  /=  phaseVolumeCo2;
    dummy1  =   log ( dummy1 );
    dummy3  =   b_mix;
    dummy3  /=  dummy2;
    dummy1  -=  dummy3;
    A4      *=  dummy1;

    //  A5 = - log( P * V_root / ( R * T )  );
    double A5 =    pressure;
    A5 *=   phaseVolumeCo2;
    A5 /=   R;
    A5 /=   temp;
    A5 =    log( A5 );
    A5 *=   -1;

    double phi_co2 =   A1;
    phi_co2 +=  A2;
    phi_co2 +=  A3;
    phi_co2 +=  A4;
    phi_co2 +=  A5;
    phi_co2 =   exp( phi_co2  );

    return phi_co2;
}




/**

    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    // T in C
    // P in Pa
    // phaseVolumeCo2 in m3/mol
    // fugacityH2o without unit

    // y_H2o=0

    // Range T:  12 < T <100 C
    // Range P:  P<600 bar
*/
double EOS_CO2H2ONaCl_Spycher04::FugacityH2o( double pressure,
                                                double temperature,
                                                double phaseVolumeCo2 )
{
    assert( temperature >= 10. );
    assert( temperature <= 100. );
    assert( pressure <= 6.0e7 );

    const double temp(degreeCToKelvin(temperature));

    //  B1 = log ( V_root / ( V_root - b_mix ) );
    double B1      =   phaseVolumeCo2;
    double dummy1  =   phaseVolumeCo2;
    dummy1  -=  b_mix;
    B1      /=  dummy1;
    B1      =   log(B1);

    //  B2 = b_h2o / ( V_root - b_mix );
    double B2 =  b_h2o;
    B2 /= dummy1;

    //  B3 = - 2 * ( y_co2 * a_h2oco2 ) / ( R * pow( T, 1.5 ) * b_mix ) * log ( ( V_root + b_mix ) / V_root) ;
    double B3 =  -2.;
    // B3      *=  1.0 (y_co2);
    B3      *=  a_h2oco2;
    dummy1  =   pow( temp, 1.5 );
    dummy1  *=  R;
    dummy1  *=  b_mix;
    B3      /=  dummy1;
    dummy1  =   phaseVolumeCo2;
    dummy1  +=  b_mix;
    double dummy2  =   dummy1;
    dummy1  /=  phaseVolumeCo2;
    dummy1  =   log( dummy1 );
    B3      *=  dummy1;

    //  B4 = ( a_mix * b_h2o / ( R * pow( T, 1.5 ) *  pow( b_mix,2 ) ) ) * ( log ( (V_root + b_mix) / V_root ) - b_mix / (V_root + b_mix) );
    double B4      =   a_mix(temperature);
    B4      *=  b_h2o;
    dummy1  =   R;
    double dummy3  =   pow( temp, 1.5 );
    dummy1  *=  dummy3;
    dummy3  =   pow( b_mix,2 );
    dummy1  *=  dummy3;
    B4      /=  dummy1;
    dummy1  =   dummy2;
    dummy1  /=  phaseVolumeCo2;
    dummy1  =   log ( dummy1 );
    dummy3  =   b_mix;
    dummy3  /=  dummy2;
    dummy1  -=  dummy3;
    B4      *=  dummy1;

    //  B5 = - log( P * V_root / ( R * T )  );
    double B5 =  pressure;
    B5 *= phaseVolumeCo2;
    B5 /= R;
    B5 /= temp;
    B5 =  log( B5 );
    B5 *= -1;

    //  phi_h2o = B1 + B2 + B3 + B4 + B5;
    double phi_h2o =   B1;
    phi_h2o +=  B2;
    phi_h2o +=  B3;
    phi_h2o +=  B4;
    phi_h2o +=  B5;
    phi_h2o =   exp( phi_h2o );

    return phi_h2o;
}



/**
    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    // T in C
    // thermEquilConstH2o without unit

    // range T 10ñ110∞C

    //thermEquilConstH2o for P0=1 bar

*/
double EOS_CO2H2ONaCl_Spycher04::thermEquilConstH2o( double temperature )
{
    //  kH2o    = -2.209 + 3.097e-2 * (T) - 1.098e-4 * pow( T, 2. ) + 2.048e-7 * pow( T, 3. );
    //  kH2o    = pow(10, kH2o);
    double kH2o    =   -2.209;
    double dummy1  =   temperature;
    dummy1  *=  3.097e-2;
    kH2o    +=  dummy1;
    dummy1  =   temperature;
    dummy1  =   pow( dummy1, 2. );
    dummy1  *=  1.098e-4;
    kH2o    -=  dummy1;
    dummy1  =   temperature;
    dummy1  =   pow( dummy1, 3. );
    dummy1  *=  2.048e-7;
    kH2o    +=  dummy1;
    kH2o    =   pow(10, kH2o);

    return kH2o;
}




/**
    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    // T in C
    // thermEquilConstH2o without unit

    // range T 12ñ110∞C

    //thermEquilConstCo2G for P0=1 bar
*/
double EOS_CO2H2ONaCl_Spycher04::thermEquilConstCo2G( double temperature )
{
    //  kCo2G   = 1.198 + 1.304e-2 * (T) - 5.446e-5 * pow( T, 2. );
    //  kCo2G   = pow(10, kCo2G);
    double kCo2G   =   1.198;
    double dummy1  =   temperature;
    dummy1  *=  1.304e-2;
    kCo2G   +=  dummy1;
    dummy1  =   temperature;
    dummy1  =   pow( dummy1, 2. );
    dummy1  *=  5.446e-5;
    kCo2G   -=  dummy1;
    kCo2G   =   pow( 10, kCo2G );

    return kCo2G;
}




/**
    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    // T in C
    // thermEquilConstCo2L without unit

    // range T 12ñ31∞C

*/
double EOS_CO2H2ONaCl_Spycher04::thermEquilConstCo2L( double temperature )
{
    //thermEquilConstCo2L for P0=1 bar
    //  kCo2L   = 1.169 + 1.368e-2 * (T) - 5.380e-5 * pow( T, 2.);
    //  kCo2L   = pow( 10, kCo2L);
    double kCo2L   =   1.169;
    double dummy1  =   temperature;
    dummy1  *=  1.368e-2;
    kCo2L   +=  dummy1;
    dummy1  =   temperature;
    dummy1  =   pow( dummy1, 2. );
    dummy1  *=  5.380e-5;
    kCo2L   -=  dummy1;
    kCo2L   =   pow( 10, kCo2L);

    return kCo2L;
}




/**
    //  Drummond 1981
    //  Boiling and mixing of hydrothermal fluids: chemical effects on mineral precipitation
    //  PhD thesis, Pennsylvania State University

    // we can find it also in
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // T in C
    // activityCoefficient without unit, we should transform it to fraction molar to be used
    // mSalt molality in mol/Kg

    // range T : 20 ñ 400∞C
    // range msalt: 0ñ6.5 m

*/
double EOS_CO2H2ONaCl_Spycher04::activityCoefficientDrummond1981(double pressure, double temperature, double mSalt )
 {
    double temp(degreeCToKelvin(temperature)); // conversion from C to K

    //  activityCoeff = (-1.0312 + 1.2806e-3 * T + 255.9 / T ) * mSalt - ( 4.445e-1 - 1.606e-3 * T ) * mSalt / (mSalt + 1) ;
    //  activityCoeff = exp ( activityCoeff );
    double dummy1 =   -1.0312;
    double dummy2 =   1.2806e-3;
    dummy2        *=  temp;
    dummy1        +=  dummy2;
    dummy2        =   255.9;
    dummy2        /=  temp;
    dummy1        +=  dummy2;
    dummy1        *=  mSalt;
    double activityCoeff =   dummy1;
    dummy1        =   -1.606e-3;
    dummy1        *=  temp;
    dummy1        +=  4.445e-1;
    dummy1        *=  mSalt;
    dummy2        =   mSalt;
    dummy2        +=  1;
    dummy1        /=  dummy2;
    activityCoeff -=  dummy1;
    activityCoeff =   exp ( activityCoeff );// activity in terms of molality gamm_m

    // conversion to activity in terms of mole fraction gamm_x
    double xCo2 = x_Co2(pressure,temperature,mSalt);

    double mCO2 = molalCo2FromBrineMoleFractionCo2(xCo2,mSalt);

    double aux = mCO2+mSalt;
    aux /= 55.508;
    aux += 1.;

    double dummy = mCO2;
    dummy /= 55.508;
    dummy += 1.;

    return activityCoeff * aux/dummy;
}



/** Computes the solubility of CO2 in brine
    
    @ref Duan and Sun 2003
    
    An improved model calculating CO2 solubility in pure water and aqueous NaCl solutions from 257 to 533 K and from 0 to  2000 bar
    Chem. Geol. 193, 257-272

    we can find it also in
      Spycher et al (2004),
      CO2-H2O Mixtures in the Geological Sequestration of CO2.
      II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    T      in C
    P      in Pa (total pressure: PCo2 + PH2o)
    mSalt  molality in mol/Kg
    activityCoefficient without unit

    range T:        from   273 to 533 K
    range P:        from 0 to 2000 bar
    ionic strength: from 0 to 4.3 m
    up to 6m NaCl molality
    up to 4m CaCl molality

    This version is just valid for NaCl [mNaCl](mK=mCa=mMg=mSo4=0)
*/
double EOS_CO2H2ONaCl_Spycher04::activityCoefficientDuanSun2003( double temperature,
                                                                   double pressure,
                                                                   double mSalt )
{
    const double
            aa1( -0.411370585 ),
            aa2( 6.07632013e-4 ),
            aa3( 97.5347708 ),
            aa4( -0.0237622469 ),
            aa5( 0.0170656236 ),
            aa6( 1.41335834e-5 ),
            b1( 3.36389723e-4 ),
            b2( -1.98298980e-5 ),
            b3( 2.12220830e-3 ),
            b4( -5.24873303e-3 ),
            p( pressure * 1.e-5 ), //conversion pressure from Pa to bar
            temp(degreeCToKelvin(temperature)); // conversion of temperature from C to K

    double lambda  = aa1 + aa2 * temp + aa3 / temp + aa4 * p / temp + aa5 * p / (630 - temp) + aa6 * temp * log( p );
    double zeta    = b1 + b2 * temp  + b3 * p / temp + b4 * p / ( 630 - temp );

    double activityCoeff = 2 * lambda * mSalt + zeta * pow( mSalt, 2);
    activityCoeff          = exp( activityCoeff );

    return activityCoeff;
}





/**

    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // T in C
    // mSalt molality in mol/Kg
    // activityCoefficient without unit, en mole fraction ready to use


    // range T : 0 ñ 350 C
    // range NaCl molality: 0ñ1.95 m
*/
double EOS_CO2H2ONaCl_Spycher04::activityCoefficientBattistelliEtal1997( double temperature,
                                                                 double mSalt )
{
    const double
            aa1( 1.19784e-1 ),
            aa2( -7.17823e-4 ),
            aa3( 4.93854e-6 ),
            aa4( -1.03826e-8 ),
            aa5( 1.08233e-11 );

    double ks( aa1 );
    ks += aa2 * temperature;
    ks += aa3 * pow( temperature, 2 );
    ks += aa4 * pow( temperature, 3 );
    ks += aa5 * pow( temperature, 4 );

    double activityCoeff = ks * mSalt;
    activityCoeff = pow(10, activityCoeff );

    return activityCoeff;

}




/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // T in C
    // mSalt molality in mol/Kg

    // range T : 40 ñ 160 C
    // range NaCl molality: 0ñ6.0 m
    // just the NaCl was taken into account
*/
double EOS_CO2H2ONaCl_Spycher04::activityCoefficientRumpf1994( double pressure ,double temperature,double mSalt )
{
    // activityCoefficient need conversion to mole fraction to be used (we added the conversion)
    const double
            aa1( 0.254 ),
            aa2( -76.82 ),
            aa3( -10656 ),
            aa4( 6312e3 ),
            gamma( 0.0028 ),
            temp(degreeCToKelvin(temperature));

    double B0 = aa1;
    B0 += aa2 / temp;
    B0 += aa3 / pow( temp, 2 );
    B0 += aa4 / pow( temp, 3 );

    double activityCoeff = 2 * mSalt *B0 + 3 * pow(mSalt,2) * gamma;
    activityCoeff = exp (activityCoeff );

    // conversion to molar fraction
    double xCo2 = x_Co2(pressure,temperature,mSalt);

    double mCO2 = molalCo2FromBrineMoleFractionCo2(xCo2,mSalt);

    double aux = mCO2+mSalt;
    aux /= 55.508;
    aux += 1.;

    double dummy = mCO2;
    dummy /= 55.508;
    dummy += 1.;

    return activityCoeff * aux/dummy;
}




/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // T in C
    // P in Pa
    // SpycherA in bar-1

*/
double EOS_CO2H2ONaCl_Spycher04::calculateSpycherA( double pressure,
                                            double temperature,
                                            double kH2o,
                                            double phiH2o )
{
    const double pTotal(pressure),
                   temp(degreeCToKelvin(temperature)); // conversion from C to K

    //  A = kH2o / ( phiH2o * pTotal * 1.0e-5) * exp( ( P - p0 ) * vH2o / (R * T) );
    double A       =   kH2o;
    double dummy1  =   phiH2o;
    dummy1  *=  pTotal;
    dummy1  *=  1.0e-5;
    A       /=  dummy1;
    dummy1  =   pressure;
    dummy1  -=  p0;
    dummy1  *=  vH2o;
    dummy1  /=  R;
    dummy1  /=  temp;
    dummy1  =   exp(dummy1);
    A       *=  dummy1;

    return A;
}



/** Spycher et al (2004),
    CO2-H2O Mixtures in the Geological Sequestration of CO2.
    II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    T in C
    P in Pa
    phaseVolumeCo2 in m3/mol
    SpycherB in bar
    activityCoeff (molar fraction gamm_s)
*/
double EOS_CO2H2ONaCl_Spycher04::calculateSpycherB( double pressure,
                                                      double temperature,
                                                      double phaseVolumeCo2,
                                                      double  phiCo2,
                                                      double kCo2L,
                                                      double kCo2G,
                                                      double activityCoeffCO2 )
{
    const double pTotal(pressure), temp(degreeCToKelvin(temperature));

    if ( ( temp < 304 ) && (phaseVolumeCo2 < 94 * 1e-6) )
      {
        //    B = phiCo2 * pTotal* 1.0e-5 / ( 55.508 * activityCoeff * kCo2L ) * exp( -( P - p0 ) * vCo2 / (R * T) );
        double B       =   phiCo2;
        B       *=  pTotal;
        B       *=  1.0e-5;
        double dummy1  =   55.508;
        dummy1  *=  activityCoeffCO2;
        dummy1  *=  kCo2L;
        B       /=  dummy1;
        dummy1  =   -pressure;
        dummy1  +=  p0;
        dummy1  *=  vCo2;
        dummy1  /=  R;
        dummy1  /=  temp;
        dummy1  =   exp(dummy1);
        B       *=  dummy1;
        return B;
      }
  
    // else
    //    B = phiCo2 * pTotal * 1.0e-5/ ( 55.508 * activityCoeff * kCo2G ) * exp( -( P - p0 ) * vCo2 / (R * T) );
    double B       =   phiCo2;
    B       *=  pTotal;
    B       *=  1.0e-5;
    double dummy1  =   55.508;
    dummy1  *=  activityCoeffCO2;
    dummy1  *=  kCo2G;
    B       /=  dummy1;
    dummy1  =   -pressure;
    dummy1  +=  p0;
    dummy1  *=  vCo2;
    dummy1  /=  R;
    dummy1  /=  temp;
    dummy1  =   exp(dummy1);
    B       *=  dummy1;

    return B;
}




/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar
*/
double EOS_CO2H2ONaCl_Spycher04::molarFracH2oCarbon( double spycherA,
                                             double spycherB,
                                             double mSalt )
{
    // spycherA in bar-1
    // SpycherB in bar
    // msalt in mol/Kg

    // approximation XH20= aH20
    //molarFracH2oCarbon = yH20
    const double A(spycherA),
                   B(spycherB);

    //  yH2o = ( 1 - B ) * 55.508 / ( ( 1/ A - B ) * ( stoichio * mSalt + 55.508 ) + stoichio * mSalt * B ) ;
    double yH2o    =   1;
    yH2o    /=  A;
    yH2o    -=  B;
    double dummy1   =  stoichio;
    dummy1  *=  mSalt;
    dummy1  +=  55.508;
    yH2o    *=  dummy1;
    dummy1  =   stoichio;
    dummy1  *=  mSalt;
    dummy1  *=  B;
    yH2o    +=  dummy1;
    yH2o    =   1./ yH2o;
    yH2o    *=  1- B;
    yH2o    *=  55.508;

    return yH2o;
}



/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar
*/
double EOS_CO2H2ONaCl_Spycher04::molarFracCO2Brine( double yH2o,
                                                      double spycherB )
{
    // spycherA in bar-1
    // SpycherB in bar

    // molarFracCO2Brine XC02

    const double B(spycherB);

    //  xCo2 =  B * (1 - yH2o);
    double xCo2(1.);
    xCo2  -=  yH2o;
    xCo2  *=  B;

    return xCo2;
}



/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // in the case of pure water => molality salt : msalt=0
*/
double EOS_CO2H2ONaCl_Spycher04::molalCo2FromPureWaterMolarFracCo2( double MolarFracCo2 )
{
    double molalCo2 = 55.508;
    molalCo2 *= MolarFracCo2;
    molalCo2 /= (1 - MolarFracCo2);

    return molalCo2;
}





/**  returns molalCo2

    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // determining the molality of Co2 in saline water from a molality of Co2 in pure water
    // activity coefficient in this case is the same as of Duan

*/
double EOS_CO2H2ONaCl_Spycher04::molalCo2FromPureWaterMolalCo2( double molalCo2Pure,
                                                        double activityCoeff )
{
    double molalCo2 = molalCo2Pure;
    molalCo2  /= activityCoeff;

    return molalCo2;
}




/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // in this case we have saline water => molality salt : msalt
*/
double EOS_CO2H2ONaCl_Spycher04::molalCo2FromBrineMoleFractionCo2( double xCo2,
                                                           double mSalt )
{
    //  molalityCo2 = xCo2 * ( stoichio * mSalt + 55.508 )  / ( 1 - xCo2 );
    double molalityCo2 =   stoichio;
    molalityCo2 *=  mSalt;
    molalityCo2 +=  55.508;
    molalityCo2 *=  xCo2;
    double dummy1      =   1.;
    dummy1      -=  xCo2;
    molalityCo2 /=  dummy1;

    return molalityCo2;
}



/**
    returns molarFracCo2
 
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar
*/
double EOS_CO2H2ONaCl_Spycher04::molarFracCo2FromBrineMolalCo2( double molalCo2,
                                                                  double mSalt)
{
   return molalCo2 / ( molalCo2 + 55.508 + stoichio * mSalt );
}



/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // Xsalt=f(mco2,ms)
    // water contains NaCl and Co2

*/
double EOS_CO2H2ONaCl_Spycher04::molalNaClToMolarFracNaClInCo2SatAqueousPhase(double molalCo2,
                                                                                double mSalt)
{
    //  massSalt = stoichio * mSalt / ( 55.508 + stoichio * mSalt + molalityCo2 );
    double xSalt = stoichio;
    xSalt    *=  mSalt;
    double dummy1 = stoichio;
    dummy1   *=  mSalt;
    dummy1   +=  55.508;
    dummy1   +=  molalCo2;
    xSalt    /=  dummy1;

    return xSalt;
}



/**
    //  Spycher et al (2004),
    //  CO2-H2O Mixtures in the Geological Sequestration of CO2.
    //  II. Partitioning in Chloride Brines at 12-100oC and up to 600 bar

    // XH20=f(XC02,xSalt)
    // water contains NaCl and Co2

*/
double EOS_CO2H2ONaCl_Spycher04::molarFracH2oBrine ( double xCo2,
                                                       double xSalt )
{
    return 1. - xCo2 - xSalt;
}


/// returns yCo2Brine calculated from 1-yH2oBrine
double EOS_CO2H2ONaCl_Spycher04::molarFracCo2Carbon( double yH2oBrine )
{
   return 1. - yH2oBrine;
}



double EOS_CO2H2ONaCl_Spycher04::massFracCo2InCarbonicPhase( double molarFracCo2InCarbonicPhase,
                                                               double molarFracH2oInCarbonicPhase )
{
    double bulkMassFrac  = molarFracCo2InCarbonicPhase * molarMassCo2;
    bulkMassFrac += molarFracH2oInCarbonicPhase * molarMassH2o;

    return    (molarFracCo2InCarbonicPhase * molarMassCo2 / bulkMassFrac)*100.;// mass fraction in % weight;
}



double EOS_CO2H2ONaCl_Spycher04::massFracH2oInCarbonicPhase( double molarFracCo2InCarbonicPhase,
                                                               double molarFracH2oInCarbonicPhase )
{
    double bulkMassFrac  = molarFracCo2InCarbonicPhase * molarMassCo2;
    bulkMassFrac += molarFracH2oInCarbonicPhase * molarMassH2o;

    return    (molarFracH2oInCarbonicPhase * molarMassH2o / bulkMassFrac)*100.;// mass fraction in % weight;
}



double EOS_CO2H2ONaCl_Spycher04::massFracCo2inAqueousPhase( double molarFracCo2InAqueousPhase,
                                                    double molarFracH2oInAqueousPhase,
                                                    double molarFracNaclInAqueousPhase )
{
    double bulkMassFrac  = molarFracCo2InAqueousPhase * molarMassCo2;

    bulkMassFrac += molarFracH2oInAqueousPhase * molarMassH2o;
    bulkMassFrac += molarFracNaclInAqueousPhase * molarMassNacl;

    return (molarFracCo2InAqueousPhase * molarMassCo2 / bulkMassFrac)*100;// mass fraction in % weight
}

double EOS_CO2H2ONaCl_Spycher04::massFracH2oInAqueousPhase( double molarFracCo2InAqueousPhase,
                                                    double molarFracH2oInAqueousPhase,
                                                    double molarFracNaclInAqueousPhase )
{
    double bulkMassFrac  = molarFracCo2InAqueousPhase * molarMassCo2;

    bulkMassFrac += molarFracH2oInAqueousPhase * molarMassH2o;
    bulkMassFrac += molarFracNaclInAqueousPhase * molarMassNacl;

    return (molarFracH2oInAqueousPhase * molarMassH2o / bulkMassFrac)*100;// mass fraction in % weight
}



double EOS_CO2H2ONaCl_Spycher04::massFracNaClInAqueousPhase( double molarFracCo2InAqueousPhase,
                                                     double molarFracH2oInAqueousPhase,
                                                     double molarFracNaclInAqueousPhase )
{
    double bulkMassFrac  = molarFracCo2InAqueousPhase * molarMassCo2;
    bulkMassFrac += molarFracH2oInAqueousPhase * molarMassH2o;
    bulkMassFrac += molarFracNaclInAqueousPhase * molarMassNacl;

    return (molarFracNaclInAqueousPhase * molarMassNacl / bulkMassFrac)*100;// mass fraction in % weight
}




/**
    // THERMODYNAMIC MODEL (after Duan and Sun)

    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16


    //  calculation of a_Co2=a_mix

    //  a_Co2   = (7.54 * 1.0e7 - 4.13 * 1.0e4 * T ) ;

    // Range T:  283 < T < 380

    // a_Co2 unit is Pa m6 K0.5 molñ2


    //assumption infinite H2O diluation
    // yH20=0 and yCo2=1

*/
double EOS_CO2H2ONaCl_Spycher04::a_Co2(double temperature)
{
    double dummy1,aCo2(0.),
             temp(degreeCToKelvin(temperature)); //conversion from C to K

    aCo2    =   7.54;
    aCo2   *=  1.0e7;
    dummy1  =   4.13;
    dummy1 *=  1.0e4;
    dummy1 *=  temp;
    aCo2   -=  dummy1;
    aCo2   *=  1.0e5;
    aCo2   *=  1.0e-12; // conversion from bar cm6 K0.5 molñ2 to Pa m6 K0.5 molñ2

    return   aCo2;
}


/// activity of mixture
double EOS_CO2H2ONaCl_Spycher04::a_mix(double temperature)
{
    return a_Co2(temperature);
}


// thermodynamic properties as a function of P, T and msalt


/// compressed volume molar CO2
double EOS_CO2H2ONaCl_Spycher04::V_Co2(double pressure,double temperature )
{
    return CompressedVolumeCo2 (pressure,temperature);
}


// molar fractions

/// molar fraction CO2 in carbonic phase, yCO2
double EOS_CO2H2ONaCl_Spycher04::y_Co2(double pressure,double temperature, double msalt  )
{
    double yH2o=y_H2o(pressure,temperature,msalt);

    return molarFracCo2Carbon(yH2o);
}


/**
    computes mole fraction CO2 in aqueous phase, xCo2 at equilibrium
*/
double EOS_CO2H2ONaCl_Spycher04::x_Co2(double pressure,double temperature, double msalt )
{
    double yH2o=y_H2o(pressure,temperature,msalt);

    double phaseVolumeCo2=V_Co2(pressure,temperature);

    double kCo2L  = thermEquilConstCo2L( temperature );

    double kCo2G = thermEquilConstCo2G( temperature );

    double phiCo2 = FugacityCo2( pressure,temperature,phaseVolumeCo2 );

// to use the activity Coefficient of Battistelli Etal 1997 as it is ready to use
// activityCoeff = activityCoefficientBattistelliEtal1997( temperature,msalt );
// or activityCoeff = activityCoefficientRumpf1994( pressure,temperature,msalt );

    double activityCoeff = activityCoefficientDuanSun2003( temperature,pressure,msalt );

    double spycherB= calculateSpycherB( pressure,temperature,phaseVolumeCo2,phiCo2,kCo2L,kCo2G,activityCoeff );

    return molarFracCO2Brine(yH2o,spycherB);
}



/// molality CO2 in aqueuse phase, mco2
double EOS_CO2H2ONaCl_Spycher04::m_Co2(double pressure,double temperature, double msalt )
{
    double xco2 = x_Co2(pressure,temperature,msalt);

    return molalCo2FromBrineMoleFractionCo2(xco2,msalt);
}



/// molar fraction H20 in carbonic phase, yH2o
double EOS_CO2H2ONaCl_Spycher04::y_H2o(double pressure,double temperature, double msalt  )

{
     double phaseVolumeCo2=V_Co2(pressure,temperature);

     double kH2o = thermEquilConstH2o( temperature );

     double kCo2L  = thermEquilConstCo2L( temperature );

     double kCo2G = thermEquilConstCo2G( temperature );

     double phiCo2 = FugacityCo2( pressure,temperature,phaseVolumeCo2 );

     double phiH2o = FugacityH2o( pressure,temperature,phaseVolumeCo2 );

     // we will use the  activity Coefficient of Battistelli Etal 1997 as it is ready to use

   //  activityCoeff = activityCoefficientBattistelliEtal1997( temperature,msalt );

     double activityCoeff = activityCoefficientDuanSun2003( temperature,pressure,msalt );

 //    activityCoeff = activityCoefficientRumpf1994( pressure,temperature,msalt );

     double spycherA = calculateSpycherA( pressure,temperature,kH2o,phiH2o );

     double spycherB= calculateSpycherB( pressure,temperature,phaseVolumeCo2,phiCo2,kCo2L,kCo2G,activityCoeff );

    return molarFracH2oCarbon(spycherA, spycherB, msalt );
}




/// molar fraction H20 in aqueuse phase, xH2o
double EOS_CO2H2ONaCl_Spycher04::x_H2o(double pressure,double temperature, double msalt  )
{
    double xSalt=X_s(pressure,temperature,msalt);

    double xCo2=x_Co2(pressure,temperature,msalt );

    return molarFracH2oBrine( xCo2,xSalt );
}



/// molar fraction salt in aqueous phase, Xs
double EOS_CO2H2ONaCl_Spycher04::X_s(double pressure,double temperature, double msalt )
{
    double xCo2=x_Co2(pressure,temperature,msalt );

    double mCo2=molalCo2FromBrineMoleFractionCo2(xCo2,msalt); // CO2 molality

    return molalNaClToMolarFracNaClInCo2SatAqueousPhase(mCo2,msalt);
}


// mas fractions in % weight


/// mass fraction CO2 in carbonic phase Y_CO2 @todo DOES THIS MAKE SENSE?
double EOS_CO2H2ONaCl_Spycher04::Y_Co2(double pressure,double temperature, double msalt  )
{
    double yCo2=y_Co2(pressure,temperature,msalt );

    double yH20=y_H2o(pressure,temperature,msalt );

    return massFracCo2InCarbonicPhase(yCo2,yH20);
}



/// mass fraction H20 in carbonic phase, Y_H2O
double EOS_CO2H2ONaCl_Spycher04::Y_H2o(double pressure,double temperature, double msalt  )
{
    double yCo2=y_Co2(pressure,temperature,msalt );

    double yH20=y_H2o(pressure,temperature,msalt );

    return massFracH2oInCarbonicPhase(yCo2,yH20);
}



/// mass fraction CO2 in aqueous phase, X_CO2
double EOS_CO2H2ONaCl_Spycher04::X_Co2(double pressure,double temperature, double msalt  )
{
    double xCo2=x_Co2(pressure,temperature,msalt );

    double xH20=x_H2o(pressure,temperature,msalt );

    double xs=X_s(pressure,temperature,msalt );

    return massFracCo2inAqueousPhase(xCo2,xH20,xs);
}




/// mass fraction H20 in aqueous phase, X_H2O
double EOS_CO2H2ONaCl_Spycher04::X_H2o(double pressure,double temperature, double msalt  )
{
    double xCo2=x_Co2(pressure,temperature,msalt );

    double xH20=x_H2o(pressure,temperature,msalt );

    double xs=X_s(pressure,temperature,msalt );

    return massFracH2oInAqueousPhase(xCo2,xH20,xs);
}



/// mass fraction of salt in aqueous phase after equilibration
double EOS_CO2H2ONaCl_Spycher04::X_salt(double pressure,double temperature, double msalt )
{
    double xCo2=x_Co2(pressure,temperature,msalt );

    double xH20=x_H2o(pressure,temperature,msalt );

    double xs=X_s(pressure,temperature,msalt );

    return massFracNaClInAqueousPhase(xCo2,xH20,xs);
}


/// mole fraction of salt in aqueous phase after equilibration
double EOS_CO2H2ONaCl_Spycher04::x_salt(double molalityCO2, double msalt )
{
   return  molalNaClToMolarFracNaClInCo2SatAqueousPhase(molalityCO2,msalt);
}



//  --------------------
//  Transport Properties
//  --------------------





/** @brief Partial molar volume of dissolved co2 in brine

    //  Garcia 2001,
    //  Density of aqueous solutions of co2,
    //  LBNL Report 49023,
    //  Lawrence Berkeley National Laboratory,
    //  Berkely, CA

    // valid for   0 < T < 300 ∞C
    //  Temperature in          ∞C
    //  volume Partial Molar CO2 is  the apparent molar volume of dissolved CO2 (m3/mol)

*/
double EOS_CO2H2ONaCl_Spycher04::volumePartialMolarCo2( double temperature )
{
    const double temperature2(temperature * temperature);
    double vPartialmolar = ( 37.51 - 9.585e-2 * temperature + 8.740e-4 * temperature2 - 5.044e-7 * temperature2*temperature );

    return vPartialmolar* 1e-6;// conversion from cm^3/mol to m^3/mol
}






/** @brief  Density of brine phase

    Rowe and Chou 1970
    //  Pressure-volume-temperature-concentration relation of aqueous NaCl solutions
    //  J. Chem. Eng. Data 15, 61-66
    //
    //  from cgs to SI by
    //
    //  Kestin et al 1981
    //  Tables of the Dynamic and Kinematic Viscosity of Aqueous NaCl Solutions in
    //  the Temperature Range 20-150 C and the Pressure Range 0.1-35 MPa
    //  J. Phys. Chem. Ref. Data, Vol. 10, No. 1,
    //
    //  Density in              kg/m^3
    //  Pressure in             Pa
    //  Temperature in          ∞C
    //  Salt Concentration in   mass fraction

    // valid for a range of pressure: 1 < P < 350 bar
    // valid for a range of temperature: 20 < T < 150 ∞C
    // valid for a range of salt molality : 0 < m_salt < 6 molality
    //   equivalent to salt mass fraction : 0 < massFracSalt < 0.95347 %weight
*/
double EOS_CO2H2ONaCl_Spycher04::densityBrine( double pressure,
                                                 double temperature,
                                                 double mSalt )
{
    double
            A(0.),
            B(0.),
            C(0.),
            D(0.),
            E(0.),
            F(0.),
            G(0.),
            H(0.),
            densBrine(0.),
            dummy1(0.),
            press ( pressure * 1e-6 ), //converting from Pa to MPa
            temp(degreeCToKelvin(temperature)), //converting from ∞C to K
            massFracSalt(0.0);

    massFracSalt = molalNaClToMassFracNaClInAqueousPhase(mSalt)*1.e-2; //conversion from salt molality to salt mass fraction

    //  density brine = 1/(A - B * press - C * pow(press, 2.) + massFracSalt * D  + pow(massFracSalt, 2.) * E - massFracSalt * F * press  - pow( massFracSalt, 2.) * G * press - 0.5 * H  * pow( press, 2.))
    const double temp2(temp * temp);
    const double press2(press * press);
    A         = 1.006741e2 / temp2 - 1.127522 / temp + 5.916365e-3 - 1.035794e-5 * temp + 9.270048e-9 * temp2;
    B         = 1.042948 / temp2 - 1.1933677e-2 / temp + 5.307535e-5 - 1.0688768e-7 * temp + 8.492739e-11 * temp2;
    C         = 1.23268e-9 - 6.861928e-12 * temp;
    D         = -2.5166e-3 + 1.11766e-5 * temp - 1.70552e-8 * temp2;
    E         = 2.84851e-3 - 1.54305e-5 * temp + 2.23982e-8 * temp2;
    F         = -1.5106e-5 + 8.4605e-8 * temp - 1.2715e-10 * temp2;
    G         = 2.7676e-5 - 1.5694e-7 * temp + 2.3102e-10 * temp2;
    H         = 6.4633e-8 - 4.1671e-10 * temp + 6.8599e-13 * temp2;
    const double  massFracSalt2(massFracSalt * massFracSalt);
    dummy1    = A - B * press - C * press2 + massFracSalt * D  + massFracSalt2 * E - massFracSalt * F * press  - massFracSalt2 * G * press - 0.5 * H  * press2;
    densBrine = 1 / dummy1;

    return    densBrine;
}



/** @brief Density of co2 saturated brine phase
    //  Garcia 2001,
    //  Density of aqueous solutions of co2,
    //  LBNL Report 49023,
    //  Lawrence Berkeley National Laboratory,
    //  Berkely, CA

    //  Density in              kg/m^3
    //  densBrine in            kg/m^3
    //  vPartialmolar in        m^3/mol
    //  xCo2 in                 mol fraction

*/
double EOS_CO2H2ONaCl_Spycher04::densityAqueousPhase( double vPartialmolar,
                                                        double densBrine,
                                                        double xCo2 )
{
    double xH2o       = 1 - xCo2;
    double itsDensAq  = 1 + xCo2 * molarMassCo2 / ( xH2o * molarMassH2o);

    itsDensAq  /= ( xCo2 * vPartialmolar / ( molarMassH2o * xH2o) + 1 / densBrine );

    return itsDensAq;
}




/// Density of carbonic phase
double EOS_CO2H2ONaCl_Spycher04::densityCarbonicPhase( double phaseVolumeCo2 )
{
    //  phaseVolumeCo2 in              m^3 /mol
    //  densityCarbonicPhase in        kg/m^3
    // double itsDensCo2    =  molarMassCo2 / phaseVolumeCo2;
    assert( molarMassCo2 > 0. );

    return  molarMassCo2 / phaseVolumeCo2;
}






/** @brief  Brine compressibility

    //  Rowe and Chou 1970
    //  Pressure-volume-temperature-concentration relation of aqueous NaCl solutions
    //  J. Chem. Eng. Data 15, 61-66
    //
    //  [Brine compressibility] = 1/ kPa
    //  [pressure]              = kPa
    //  [density]               = kg/m^3

*/
double EOS_CO2H2ONaCl_Spycher04::compressibilityBrine( double densBrine,
                                                         double pressure,
                                                         double densBrineRef )
{
    const double pressureRef( 101.325 ), press( pressure * 1e-3);

    double comprBrine = ( densBrine - densBrineRef ) / ( densBrine * ( press  - pressureRef ) );

    comprBrine *= 1.0e-03; //  1/kPa -> 1/Pa

    return comprBrine;

}


/** @brief Compressibility of carbonic phase

    //  Mathias, Hardisty, Trudell and Zimmerman 2009
    //  Screening and selection of sites for CO2 sequestration based on pressure buildup
    //  Int. J. Greenhouse Gas Control 3
    //
    // isothermal compressibility Pa-1

    //  [temperature] = C
    //  [pressure]    = Pa
    //  [volume]      = molar volume at pressure p  here , unit is m3/mol:
*/
double EOS_CO2H2ONaCl_Spycher04::compressibilityCarbonicPhase( double temperature,
                                                                 double phaseVolumeCo2 )
{
    double pressureVolumeGradient( 0. ),
            a ( a_mix(temperature) ),
            b ( b_mix ),
            phaseVolCo2(0.),
            betaCarb( 0. ),
            temp(degreeCToKelvin(temperature)),
            Aux1( 0. ),
            Aux2( 0. ),
            Aux3( 0. );

    phaseVolCo2 = phaseVolumeCo2;

    //  pressureVolumeGradient *=  - GasConstant * temperature  / pow( ( phaseVolCo2 -  b), 2) ;
    pressureVolumeGradient  =   -1.;
    pressureVolumeGradient  *=  R;
    pressureVolumeGradient  *=  temp;
    Aux1                    =   phaseVolCo2;
    Aux1                    -=  b;
    Aux1                    *=  Aux1;
    pressureVolumeGradient  /=  Aux1;

    cout << "pressureVolumeGradient1 " << pressureVolumeGradient <<endl;

    //  pressureVolumeGradient +=  a * ( 2*phaseVolCo2 + b ) / ( pow( temperature, 0.5) * pow( phaseVolCo2, 2 ) * pow ( ( phaseVolCo2 + b ), 2) );
    Aux2                    =   2.*a;
    Aux2                    *=  phaseVolCo2;
    Aux2                    =+  b*a;
    Aux3                    =   sqrt(temp);
    Aux3                   *=  phaseVolCo2 * phaseVolCo2;
    Aux3                   *=  pow ( ( phaseVolCo2 + b ), 2);
    Aux2                   /=  Aux3;
    pressureVolumeGradient +=  Aux2;

    betaCarb  =   -1;
    betaCarb  /=  phaseVolCo2;
    betaCarb  /=  pressureVolumeGradient;

    return betaCarb;

}



/// Compressibility of carbonic phase
double EOS_CO2H2ONaCl_Spycher04::compressibilityCarbonicPhaseZ( double pressure,
                                                        double temperature,
                                                        double phaseVolumeCo2 )
{
    //  Pressure in             Pa
    //  Temperature in          ∞C
    //  phaseVolumeCo2 in       m^3
    double temp(degreeCToKelvin(temperature)); //conversion to K

    // zFactor
    return pressure * phaseVolumeCo2 / ( R * temp );
}




/** @brief Viscosity of aqueos phase

    // Kestin, Khalifa and Correia 1981
    // Tables of the dynamic and kinemaic viscosity of aqueous NaCl solutions in the temperature range 20 - 150 C
    // and the pressure range 0.1-35 MPa
    // J. Phys. Chem. Ref. Data 10 (1), 71-87
    //
    // and
    //
    // Kestin, Khalifa, Abe, Grimes, Sookiazan, and Wakeham 1978
    // Effect of pressure on the viscosity of aqueous NaCl solutions in the temperature wange 20 - 150 C
    // J. Chem. Eng. Data, Vol 23, No4, 328-336
    //

    //  Viscosity in            Pa s
    //  Pressure in             Pa
    //  Temperature in          ∞C
    //  Salt molality in   mol/Kg

    // valid for a range of pressure: 1 < P < 350 bar
    // valid for a range of temperature: 20 < T < 150 ∞C
    // valid for a range of salt molality : 0 < msalt < 6 molality


    // Attention: pressure coeeficients \betas ' s in GPa^(-1 )

*/
double EOS_CO2H2ONaCl_Spycher04::viscosityBrine( double pressure,
                                         double temperature,
                                         double mSalt )
{
    const double
            aa1( 3.324e-2 ),
            aa2( 3.624e-3 ),
            aa3( -1.879e-4 ),
            bb1(-3.96e-2 ),
            bb2(1.02e-2 ),
            bb3(-7.02e-4 ),
            cc1( 1.2378 ),
            cc2( -1.303e-3 ),
            cc3( 3.06e-6 ),
            cc4( 2.55e-8 ),
            dd1( 6.044 ),
            dd2( 2.8e-3 ),
            dd3( 3.6e-5 ),
            beta1( -1.297 ),
            beta2( 5.74e-2 ),
            beta3( -6.97e-4 ),
            beta4( 4.47e-6),
            beta5( -1.05e-8);

    double
            cs(0.),
            betaW(0.),
            betaStar(0.),
            viscosBeta(0.),
            muW(0.),
            muR(0.),
            muZero(0.),
            viscosBrine(0.),
            aux4(0.),
            AA(0.),
            BB(0.),
            press( pressure * 1e-6);

    const double temperature2(temperature * temperature);
    const double mSalt2(mSalt * mSalt);

    cs          = dd1 + dd2 * temperature + dd3 * temperature2;
    betaW       = beta1 + beta2 * temperature + beta3 * temperature2 + beta4 * temperature2 * temperature + beta5 * temperature2 * temperature2;
    betaStar    = ( 2.5 * ( mSalt / cs ) - 2.0 * pow( ( mSalt / cs ), 2 ) + 0.5 * pow( ( mSalt / cs ), 3 ) );
    viscosBeta  = ( 0.545 + 2.8e-3 * temperature - betaW ) * betaStar + betaW;
    aux4        = ( cc1 * ( 20 - temperature )  + cc2 * pow( ( 20. - temperature ), 2 ) + cc3 * pow( ( 20. - temperature ), 3 ) + cc4 * pow( ( 20. - temperature ), 4 ) )  / ( 96 + temperature );
    muW         = pow( 10, aux4 ) * 1002;
    BB          = bb1 * mSalt + bb2 * mSalt2 + bb3 * mSalt2 * mSalt;
    AA          = aa1 * mSalt + aa2 * mSalt2 + aa3 * mSalt2 * mSalt;
    muR         = AA + BB *  aux4;
    muR         = pow( 10, muR );
    muZero      = muR * muW ;

    viscosBrine = muZero * ( 1 +  1e-3 * viscosBeta * press ); //  to account for GPa^(-1 )

    return viscosBrine*1e-6; //conversion from micro Pa s to Pa s
}





/** @brief Viscosity of carbon dioxide

    //  Fenghour and Wakeham 1998
    //  The viscosity of carbon dioxide
    //  J. Phys. Chem. Ref. Data, Vol. 27, No. 1
    //
    //  and
    //
    //  Vesovic et al 1990
    //  The transport properties of carbon dioxide
    //  J. Phys. Chem. Ref. Data, Vol. 19, No. 3
    //
    //  [temperature] = ∞C
    //  [density]     = kg/m^3
    //  [viscosity]   = Pa s

*/
double EOS_CO2H2ONaCl_Spycher04::viscosityCarbonicPhase( double temperature,
                                                           double phaseVolumeCo2 )
{
    const double
            aa0(0.235156),
            aa1(-0.491266),
            aa2(5.211155e-2),
            aa3(5.347906e-2),
            aa4(-1.537102e-2),
            dd1(0.4071119e-2),
            dd2(0.7198037e-4),
            dd3(0.2411697e-16),
            dd4(0.2971072e-22),
            dd5(-0.1627888e-22);
            //ee1(5.5934e-3),
            //ee2(6.1757e-5),
            //ee3(0.0),
            //ee4(2.6430e-11);

    const double molarvolume(phaseVolumeCo2), temp(degreeCToKelvin(temperature));

    double rho      = densityCarbonicPhase(molarvolume);
    double T_red    = temp / 251.196;
    double logT_red = log(T_red);
    double dummy1 = aa0 * 1. + aa1 * logT_red + aa2 * logT_red*logT_red + aa3 * pow( logT_red, 3. ) + aa4 * pow( logT_red, 4. );
    dummy1          = exp( dummy1 );
    double eta_0       = 1.00697 * pow( temp, 0.5 ) / dummy1;  //Eq. (3)
    double eta_excess  = dd1 * rho + dd2 * rho*rho + dd3 * pow(rho, 6. ) / pow(T_red, 3.) + dd4 * pow(rho, 8.) + dd5 * pow(rho, 8) / T_red;  //Eq.(8)

    //    eta_crit    = ee1 * rho + ee2 * pow(rho, 2.) + ee3 * pow(rho, 3.) + ee4 * pow(rho, 4.);
    //    viscosCo2   = eta_0 +  eta_excess + eta_crit;
    double viscosCo2 = eta_0 +  eta_excess ;
    viscosCo2 *= 1.0e-6; // convert muPa s -> Pa s

    return viscosCo2;
}



/**
    @brief Molecular diffusion coeffisient of co2 into brine
 
    Ratcliff and Holdcroft 1963
    Diffusivities of gases in aqueous electrolyte solutions
    Trans. Inst. Chem. Eng. 41, 315-319
 
    and
 
    Al-Rawajfeh 2004
    Modelling and simulation of co2 release in multiple-effect distillers for seawater desalination
    PhD Dissertation, Martin-Luther university of Halle-Wittenberg
    Department of engineering sciences, isntitute of theral process engineering
    Halle (Saale), Germany, pp. 8-10
 
      and
 
    McLachlan and Danckwerts 1972
    Desorption of carbon dioxide from aqueous potash solutions with and without the addition of arsenite as a catalyst
    Trans. Inst. Chem. Eng. 50, 300-309

    [temperature] = ∞C
    [molecularDiffCoeffCo2intoBrine] = m^2/s
    [viscoBrine] = Pa s
*/
double EOS_CO2H2ONaCl_Spycher04::molecularDiffCoeffCo2intoBrine( double temperature,
                                                         double viscoBrine )
{
    double temp(degreeCToKelvin(temperature));

    double D0 =  pow( 10.,  -4.1764 + 712.52 /  temp - 2.5907e5 / ( temp*temp ) );

    double Db = D0 / pow( 10, 0.87 * log10( viscoBrine *1e6 /1002 ) ); //conversion from Pa s to micro Pa s

    return Db * 1e-4; // conversion from cm^2/s to m^2/s
}




/**
    //  Spycher et al (2003),
    //  Co2-H2O mixtures in the geological sequestration of co2.
    //  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
    //  Geochimica et Cosmochimica acta, vol 67, No 16

    //  Pressure in             Pa
    //  Temperature in          ∞C
*/
double EOS_CO2H2ONaCl_Spycher04::equilKH2o( double temperature,
                                    double pressure,
                                    double kH2o )
{
    double temp(degreeCToKelvin(temperature)); //conversion from ∞C to K

    return kH2o * exp( (pressure - p0) * vH2o / (R * temp) );
}



/**
//  Spycher et al (2003),
//  Co2-H2O mixtures in the geological sequestration of co2.
//  I. assessment and calculation of mutual solubilities from 12 to 100 C  and up to 600 bar,
//  Geochimica et Cosmochimica acta, vol 67, No 16

//  Pressure in             Pa
//  Temperature in          ∞C
*/
double EOS_CO2H2ONaCl_Spycher04::equilKCo2( double temperature,
                                    double pressure,
                                    double kCo2G )
{
    double temp(degreeCToKelvin(temperature)); //conversion from ∞C to K

    return kCo2G * exp( (pressure - p0) * vCo2 / (R * temp) );
}



/// volume molar dissolved CO2 in brine
double EOS_CO2H2ONaCl_Spycher04::Vdiss_Co2(double temperature)
{
    return volumePartialMolarCo2(temperature);
}



/// density brine
double EOS_CO2H2ONaCl_Spycher04::Rho_brine(double pressure,double temperature, double msalt )
{
    return densityBrine(pressure,temperature,msalt);
}



/// viscosity brine
double EOS_CO2H2ONaCl_Spycher04::mu_brine(double pressure,double temperature, double msalt )
{
    return viscosityBrine(pressure,temperature,msalt);
}



/// compressibility brine
double EOS_CO2H2ONaCl_Spycher04::C_brine(double pressure,double temperature, double msalt )
{
    double densbrine=Rho_brine(pressure,temperature,msalt);

    // pressure =p0 = 1 bar refernce pressure
    double densBrineRef=Rho_brine(p0,temperature,msalt);

    // compressibility
    return compressibilityBrine(densbrine,pressure,densBrineRef);
}


/// density Aqueous Phase (contains dissolved CO2)
double EOS_CO2H2ONaCl_Spycher04::Rho_AqueousPhase(double pressure,double temperature, double msalt )
{
    double vPartialmolar=Vdiss_Co2(temperature);

    double densbrine= Rho_brine(pressure,temperature,msalt);

    double xCo2=x_Co2(pressure,temperature,msalt);

    // density
    return densityAqueousPhase(vPartialmolar,densbrine,xCo2);
}



/**
    viscosity Aqueous Phase (contain dissolved CO2)
    @note we ignore The effect of dissolved CO2 in the viscosity of aqeuse phase
*/
double EOS_CO2H2ONaCl_Spycher04::mu_AqueousPhase(double pressure,double temperature, double msalt )
{
   return mu_brine(pressure,temperature,msalt);
}


/// compressibility Aqueous Phase (contain dissolved CO2)
double EOS_CO2H2ONaCl_Spycher04::C_AqueousPhase(double pressure,double temperature, double msalt )
{
    // we ignore The effect of dissolved CO2 in the viscosity of aqeuse phase
    return C_brine(pressure,temperature,msalt);
}



/// density Carbonic Phase
double EOS_CO2H2ONaCl_Spycher04::Rho_CarbonicPhase(double pressure,double temperature )
{
    const double phaseVolumeCo2=V_Co2(pressure,temperature);

    return densityCarbonicPhase(phaseVolumeCo2);
}



/// viscosity Carbonic Phase
double EOS_CO2H2ONaCl_Spycher04::mu_CarbonicPhase(double pressure,double temperature )
{
    double phaseVolumeCo2=V_Co2(pressure,temperature);

    return viscosityCarbonicPhase(temperature,phaseVolumeCo2);
}



/// compressibility of the carbonic phase
double EOS_CO2H2ONaCl_Spycher04::C_CarbonicPhase(double pressure,double temperature )
{
//    double comp(0.0),phaseVolumeCo2(0.);
//    phaseVolumeCo2=V_Co2(pressure,temperature);
    // problem with compressibilityCarbonicPhase
//    comp= compressibilityCarbonicPhase(temperature,phaseVolumeCo2);
    const double intP(1.);

    return -(Bg(pressure+intP,temperature) - Bg(pressure,temperature) ) / (Bg(pressure,temperature) * intP);
}




/// Z compressbility factor Carbonic Phase
double EOS_CO2H2ONaCl_Spycher04::Z_CarbonicPhase(double pressure, double temperature )
{
    double phaseVolumeCo2=V_Co2(pressure,temperature);
    // Zcomp=
    return compressibilityCarbonicPhaseZ(pressure,temperature,phaseVolumeCo2);
}



/// molecular Diffusivity Coefficient Co2 into Brine
double EOS_CO2H2ONaCl_Spycher04::D_Co2(double pressure, double temperature, double msalt )
{
    double visBrine=viscosityBrine(pressure,temperature,msalt);
    return molecularDiffCoeffCo2intoBrine(temperature,visBrine);
}



/// equilibrium KCO2
double EOS_CO2H2ONaCl_Spycher04::KCo2(double pressure, double temperature)
{
    double kco2const=thermEquilConstCo2G(temperature);
    return equilKCo2(temperature,pressure,kco2const);
}



/// equilibrium KH2O
double EOS_CO2H2ONaCl_Spycher04::KH2o(double pressure, double temperature )
{
    double kh2oconst=thermEquilConstH2o(temperature);
    return equilKH2o(temperature,pressure,kh2oconst);
}




//  Reservoir Properties


/**

    // Hassanzadeh et al 2008,
    // Predicting PVT data for CO2ñbrine mixtures for
    // black-oil simulation of CO2 geological storage

*/
double EOS_CO2H2ONaCl_Spycher04::solutionAqueousCarbonicRatio( double pressure, double temperature ,double msalt )
{
    // Rs volume of dissolved Co2 @ Sc/
    //  rs    =   Rho_AqueousPhase(pSC,tSC,msalt);
    //  rs    *=  X_Co2(pressure,temperature,msalt);
    //  rs    /=  Rho_CarbonicPhase(pSC,tSC);
    //  rs    /=  aux;

    double rs = Rho_brine(pSC,tSC,msalt);
    rs    *=  x_Co2(pressure,temperature,msalt);
    rs    /=  Rho_CarbonicPhase(pSC,tSC);
    rs    /=  1-x_Co2(pressure,temperature,msalt);

    return rs;
}



/// vCo2 @RC / vCo2 @SC unit => Rm^3/Sm^3
double EOS_CO2H2ONaCl_Spycher04::GasFormationVolumeFactor( double pressure, double temperature )
{
    double bg = V_Co2(pressure,temperature );
    bg         /= V_Co2( pSC, tSC );

    return bg;
}




/**
    // Hassanzadeh et al 2008,
    // Predicting PVT data for CO2ñbrine mixtures for
    // black-oil simulation of CO2 geological storage

    // Vwater @RC / Vwater @SC  unit => Rm^3/Sm^3

*/
double EOS_CO2H2ONaCl_Spycher04::WaterFormationVolumeFactor( double pressure, double temperature,double msalt )
{
    //  bw    =   Rho_AqueousPhase(pSC,tSC,msalt);
    //  bw    /=  Rho_AqueousPhase(pressure,temperature,msalt);
    double bw =   Rho_brine(pSC,tSC,msalt);
    bw    /=  Rho_brine(pressure,temperature,msalt);
    bw    /=  1-X_Co2(pressure,temperature,msalt)*1e-2;

    return bw;
}




/// Reservoir solution Aqueous Carbonic Ratio (gas-water ratio)
double EOS_CO2H2ONaCl_Spycher04::Rs( double pressure, double temperature,double msalt )
{
   return solutionAqueousCarbonicRatio(pressure,temperature,msalt);
}



/// gas (CO2) formation volume factor, Bg (CO2)
double EOS_CO2H2ONaCl_Spycher04::Bg( double pressure, double temperature )
{
    return GasFormationVolumeFactor(pressure,temperature);
}



/// formation factor, Bw (brine
double EOS_CO2H2ONaCl_Spycher04::Bw( double pressure, double temperature,double msalt )
{
    return WaterFormationVolumeFactor(pressure,temperature,msalt);
}



// conversions


/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double EOS_CO2H2ONaCl_Spycher04::massFracNaClToMolalNaClInAqueousPhase( double massFracSalt)
{
    // mass fraction in % weight
    return 1. * massFracSalt / ( molarMassNacl * ( 100. - massFracSalt ) );
}



/// Mass Fraction salt (massFracNaCl) in % weight substance in weight solvent NOT in ppm
double EOS_CO2H2ONaCl_Spycher04::molalNaClToMassFracNaClInAqueousPhase( double mSalt)
{
    double dummy = mSalt * molarMassNacl;
    return dummy / (1 + dummy) * 100.; // in % weight
}



double EOS_CO2H2ONaCl_Spycher04::massFracNaClToMolarFracNaClInAqueousPhase( double massFracSalt)
{
    double molarFrac =  massFracSalt * molarMassH2o;

    double dummy = massFracSalt * (molarMassH2o - molarMassNacl) + 100. * molarMassNacl;

    return molarFrac / dummy;
}




double EOS_CO2H2ONaCl_Spycher04::molalNaClToMolarFracNaClInAqueousPhase( double mSalt)
{
    double molarFrac = mSalt / (mSalt + 55.508);

    return molarFrac;
}



double EOS_CO2H2ONaCl_Spycher04::molalNaClToPpmInAqueousPhase( double mSalt )
{
    double
            conversion( molalNaClToMassFracNaClInAqueousPhase( mSalt) );

    conversion  *=  1.e-2*1.e6; //  ppm \in [0, 1.e6] ,1.e-2 because mass fraction is in %

    return      conversion;
}



double EOS_CO2H2ONaCl_Spycher04::ppmNaClToMolalNaClInAqueousPhase( double ppmSalt )
{
    return massFracNaClToMolalNaClInAqueousPhase(ppmSalt*1.e-6*1.e2); //conversion from massfraction to molality
}

// NON MEMBER FUNCTIONS

double psiToPa( double pressureInPsi )
{
    return pressureInPsi * 6894.75729;
}


double paToPsi( double pressureInPa )
{
    return pressureInPa * 0.000145037738;
}



double paTobar( double pressureInPa )
{
    return  pressureInPa * 1.0e-5;
}



double barTopa( double pressureInbar )
{
    return  pressureInbar * 1.0e+5;
}



double degreeCToKelvin( double temperatureInC )
{
   return  temperatureInC + 273.15;
}



double KelvinTodegreeC( double temperatureInK )
{
   return  temperatureInK - 273.15;
}



static double temp( double depth )
 {
    const double temp_surface(20); // 20 C
    const double therm_grad(0.03); //3 C/100m
    return temp_surface + depth * therm_grad;
 }


static double pres( double depth )
 {
    const double pres_grad(10000.); // 0.1 bar/m
    const double pres_surface(100000); // 1 bar
    return pres_surface + depth * pres_grad;
 }


/// plotting routine (data are in text file for spreadsheet)
void EOS_CO2H2ONaCl_Spycher04::plot_brine()
{
    // needed for plotting
    const double Tmax(100.);
    const double Tmin(20.);
    const double Pmax(35000000.); //350 bar
    const double Pmin(100000.); // 1 bar
    const double msaltmin(0.);
    const double msaltmax(1.5); //87664.5 ppm
    const double nTP(100);
    const double nmsalt(4);
    const double max_depth(2600.); //2600 m

    ofstream fout,fout_density_Brine_depth,fout_viscosity_Brine_depth,fout_compressibility_Brine_depth;
    vector<double> T(nTP),P(nTP),msalt(nmsalt),depth(nTP);

    fout.open("density_brine.txt");
    fout_density_Brine_depth.open("density_brine_depth.txt");
    fout_viscosity_Brine_depth.open("viscosity_brine_depth.txt");
    fout_compressibility_Brine_depth.open("compressibility_brine_depth.txt");

    fout_density_Brine_depth<< "Depth ";
    fout_viscosity_Brine_depth<< "Depth ";
    fout_compressibility_Brine_depth<< "Depth ";

    for (int n=0; n<nTP;n++)
    {
        T[n]=Tmin+n*(Tmax-Tmin)/(nTP-1);
        P[n]=Pmin+n*(Pmax-Pmin)/(nTP-1);
        depth[n]=n*(max_depth)/(nTP-1);
    }

    for (int nm=0; nm<nmsalt;nm++)
    {

        msalt[nm]=msaltmin+nm*(msaltmax-msaltmin)/(nmsalt-1);
        fout_density_Brine_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_viscosity_Brine_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_compressibility_Brine_depth<< "\t" << "msalt= " <<msalt[nm];
    }

    fout_density_Brine_depth<<"\n";
    fout_viscosity_Brine_depth<< "\n";
    fout_compressibility_Brine_depth<<"\n";

    for (int n=0; n<nTP;n++)
    {
        fout<<T[n]<<"\t"<<paTobar(P[n]);
        fout_density_Brine_depth<<depth[n];
        fout_viscosity_Brine_depth<<depth[n];
        fout_compressibility_Brine_depth<<depth[n];

        for (int nm=0; nm<nmsalt;nm++)
        {
            fout<<"\t"<<Rho_brine(P[n],T[n],msalt[nm]);
            fout_density_Brine_depth<<"\t"<<Rho_brine(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_viscosity_Brine_depth<<"\t"<<mu_brine(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_compressibility_Brine_depth<<"\t"<<C_brine(pres(depth[n]),temp(depth[n]),msalt[nm]);
        }
        fout<<"\n";
        fout_density_Brine_depth<<"\n";
        fout_viscosity_Brine_depth<<"\n";
        fout_compressibility_Brine_depth<<"\n";
    }
}







void EOS_CO2H2ONaCl_Spycher04::plot_AqueousPhase()
{
    // needed for plotting
    const double Tmax(100.);
    const double Tmin(20.);
    const double Pmax(35000000.); //350 bar
    const double Pmin(100000.); // 1 bar
    const double msaltmin(0.);
    const double msaltmax(1.5); //87664.5 ppm
    const double nTP(100);
    const double nmsalt(4);
    const double max_depth(2600.); //2600 m

    ofstream fout,fout_density_AqueousPhase_depth,fout_density_difference_depth,
            fout_Rs_depth,fout_D_Co2_depth,fout_Bw_depth,
            fout_density_AqueousPhase_T50_msalt0,fout_density_AqueousPhase_T100_msalt0,
            fout_viscosity_AqueousPhase_T50_msalt0,fout_viscosity_AqueousPhase_T100_msalt0,
            fout_viscosity_AqueousPhase_T50_SeaWater,fout_density_AqueousPhase_T50_SeaWater;
    vector<double> T(nTP), P(nTP), msalt(nmsalt), depth(nTP), PP(1000);

    fout.open("density_AqueousPhase.txt");
    fout_density_AqueousPhase_depth.open("density_AqueousPhase_depth.txt");
    fout_density_difference_depth.open("density_difference_depth.txt");
    fout_D_Co2_depth.open("D_Co2_depth.txt");
    fout_Rs_depth.open("Rs_depth.txt");
    fout_Bw_depth.open("Bw_depth.txt");

    fout_density_AqueousPhase_T50_msalt0.open("H2O_density_T50_msalt0.txt");
    fout_density_AqueousPhase_T100_msalt0.open("H2O_density_T100_msalt0.txt");
    fout_viscosity_AqueousPhase_T50_msalt0.open("H2O_viscosity_T50_msalt0.txt");
    fout_viscosity_AqueousPhase_T100_msalt0.open("H2O_viscosity_T100_msalt0.txt");
    fout_viscosity_AqueousPhase_T50_SeaWater.open("H2O_viscosity_T50_SeaWater.txt");
    fout_density_AqueousPhase_T50_SeaWater.open("H2O_density_T50_SeaWater.txt");

    fout_density_AqueousPhase_depth<< "Depth ";
    fout_density_difference_depth<< "Depth ";
    fout_D_Co2_depth<< "Depth ";
    fout_Rs_depth<< "Depth ";
    fout_Bw_depth<< "Depth ";


    for (int n=0; n<nTP;n++)
    {
        T[n]=Tmin+n*(Tmax-Tmin)/(nTP-1);
        P[n]=Pmin+n*(Pmax-Pmin)/(nTP-1);
        depth[n]=n*(max_depth)/(nTP-1);
    }

    for (int nm=0; nm<nmsalt;nm++)
    {

        msalt[nm]=msaltmin+nm*(msaltmax-msaltmin)/(nmsalt-1);
        fout_density_AqueousPhase_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_density_difference_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_D_Co2_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_Rs_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_Bw_depth<< "\t" << "msalt= " <<msalt[nm];

    }

    fout_density_AqueousPhase_depth<<"\n";
    fout_density_difference_depth<<"\n";
    fout_D_Co2_depth<<"\n";
    fout_Rs_depth<<"\n";
    fout_Bw_depth<<"\n";

    for (int n=0; n<nTP;n++)
    {
        fout_density_AqueousPhase_depth<<depth[n];
        fout_density_difference_depth<<depth[n];
        fout_D_Co2_depth<<depth[n];
        fout_Rs_depth<<depth[n];
        fout_Bw_depth<<depth[n];

        for (int nm=0; nm<nmsalt;nm++)
        {
            fout<<"\t"<<Rho_AqueousPhase(P[n],T[n],msalt[nm]);
            fout_density_AqueousPhase_depth<<"\t"<<Rho_AqueousPhase(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_density_difference_depth<<"\t"<<100*(Rho_AqueousPhase(pres(depth[n]),temp(depth[n]),msalt[nm])-Rho_brine(pres(depth[n]),temp(depth[n]),msalt[nm]))/Rho_brine(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_D_Co2_depth<<"\t"<<D_Co2(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_Rs_depth<<"\t"<<Rs(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_Bw_depth<<"\t"<<Bw(pres(depth[n]),temp(depth[n]),msalt[nm]);
        }
        fout_density_AqueousPhase_depth<<"\n";
        fout_density_difference_depth<<"\n";
        fout_D_Co2_depth<<"\n";
        fout_Rs_depth<<"\n";
        fout_Bw_depth<<"\n";
    }

    for (int np=0; np<1000;np++)
    {
        PP[np]=Pmin+np*(barTopa(710)-Pmin)/(1000-1);
        fout_density_AqueousPhase_T50_msalt0<< paTobar(PP[np])<<"\t"<<Rho_AqueousPhase(PP[np],50,0.0)<<endl;
        fout_density_AqueousPhase_T100_msalt0<<paTobar(PP[np])<<"\t"<<Rho_AqueousPhase(PP[np],100,0.0)<<endl;
        fout_viscosity_AqueousPhase_T50_msalt0<<paTobar(PP[np])<<"\t"<<mu_AqueousPhase(PP[np],50,0.0)*1000<<endl; // in cp
        fout_viscosity_AqueousPhase_T100_msalt0<<paTobar(PP[np])<<"\t"<<mu_AqueousPhase(PP[np],100,0.0)*1000<<endl; // in cp
        fout_viscosity_AqueousPhase_T50_SeaWater<<paTobar(PP[np])<<"\t"<<mu_AqueousPhase(PP[np],50,0.62)*1000<<endl; // sea water mass fraction 3.5% molality 0.6
        fout_density_AqueousPhase_T50_SeaWater<< paTobar(PP[np])<<"\t"<<Rho_AqueousPhase(PP[np],50,0.62)<<endl;

    }
}




void EOS_CO2H2ONaCl_Spycher04::plot_CarbonicPhase()
{
    // needed for plotting
    const double Tmax(100.);
    const double Tmin(20.);
    const double Pmax(35000000.); //350 bar
    const double Pmin(100000.); // 1 bar
    const double nTP(100);
    const double nmsalt(4);
    const double max_depth(2600.); //2600 m

    ofstream fout,fout_density_CarbonicPhase_depth,fout_viscosity_CarbonicPhase_depth,
            fout_compressibility_CarbonicPhase_depth,fout_Zfacor_CarbonicPhase_depth,
            fout_KCo2_depth,fout_KH2o_depth,fout_Bg_depth,
            fout_density_CarbonicPhase_T50,fout_density_CarbonicPhase_T100,
            fout_viscosity_CarbonicPhase_T50,fout_viscosity_CarbonicPhase_T100,
            fout_viscosity_CarbonicPhase_T320K,fout_compressibility_CarbonicPhase_T50,fout_Bg_T50;
    vector<double> T(nTP),P(nTP),msalt(nmsalt),depth(nTP),PP(1000);

    fout.open("density_CarbonicPhase.txt");
    fout_density_CarbonicPhase_depth.open("density_CarbonicPhase_depth.txt");
    fout_viscosity_CarbonicPhase_depth.open("viscosity_CarbonicPhase_depth.txt");
    fout_compressibility_CarbonicPhase_depth.open("compressibility_CarbonicPhase_depth.txt");
    fout_Zfacor_CarbonicPhase_depth.open("ZFactor_CarbonicPhase_depth.txt");
    fout_KCo2_depth.open("KCo2_depth.txt");
    fout_KH2o_depth.open("KH2o_depth.txt");
    fout_Bg_depth.open("Bg_depth.txt");
    fout_Bg_T50.open("Bg_T50.txt");

    fout_density_CarbonicPhase_T50.open("CO2_density_T50.txt");
    fout_density_CarbonicPhase_T100.open("CO2_density_T100.txt");
    fout_viscosity_CarbonicPhase_T50.open("CO2_viscosity_T50.txt");
    fout_viscosity_CarbonicPhase_T100.open("CO2_viscosity_T100.txt");
    fout_viscosity_CarbonicPhase_T320K.open("CO2_viscosity_T320K.txt");

    fout_compressibility_CarbonicPhase_T50.open("CO2_compressibility_T50.txt");

    fout_density_CarbonicPhase_depth<< "Depth "<<"\t" << "Density ";
    fout_viscosity_CarbonicPhase_depth<< "Depth "<<"\t" << "viscosity ";
    fout_compressibility_CarbonicPhase_depth<< "Depth "<<"\t" << "compressibility ";
    fout_Zfacor_CarbonicPhase_depth<< "Depth "<<"\t" << "ZFactor ";
    fout_KCo2_depth<< "Depth "<<"\t" << "KCo2 ";
    fout_KH2o_depth<< "Depth "<<"\t" << "KH2o ";
    fout_Bg_depth<< "Depth "<<"\t" << "Bg ";

    fout_density_CarbonicPhase_T50<< "Pressure "<<"\t" << " CO2_density"<<endl;
    fout_density_CarbonicPhase_T100<< "Pressure "<<"\t" << " CO2_density"<<endl;
    fout_viscosity_CarbonicPhase_T50<< "Pressure "<<"\t" << " CO2_viscosity"<<endl;
    fout_viscosity_CarbonicPhase_T100<< "Pressure "<<"\t" << " CO2_viscosity"<<endl;
    fout_viscosity_CarbonicPhase_T320K<< "Pressure "<<"\t" << " CO2_viscosity"<<endl;
    fout_compressibility_CarbonicPhase_T50<< "Pressure "<<"\t" << " CO2_compressibility"<<endl;
    fout_Bg_T50<< "pressure "<<"\t" << "Bg ";

    fout_density_CarbonicPhase_depth<<"\n";
    fout_viscosity_CarbonicPhase_depth<< "\n";
    fout_compressibility_CarbonicPhase_depth<<"\n";
    fout_Zfacor_CarbonicPhase_depth<<"\n";
    fout_KCo2_depth<<"\n";
    fout_KH2o_depth<<"\n";
    fout_Bg_depth<<"\n";
    fout_compressibility_CarbonicPhase_T50<<"\n";
    fout_Bg_T50<<"\n";

    for (int n=0; n<nTP;n++)
    {
        T[n]=Tmin+n*(Tmax-Tmin)/(nTP-1);
        P[n]=Pmin+n*(Pmax-Pmin)/(nTP-1);
        depth[n]=n*(max_depth)/(nTP-1);
    }

    for (int n=0; n<nTP;n++)
    {
        fout<<T[n]<<"\t"<<paTobar(P[n]);
        fout_density_CarbonicPhase_depth<<depth[n];
        fout_viscosity_CarbonicPhase_depth<<depth[n];
        fout_compressibility_CarbonicPhase_depth<<depth[n];
        fout_Zfacor_CarbonicPhase_depth<<depth[n];
        fout_KCo2_depth<<depth[n];
        fout_KH2o_depth<<depth[n];
        fout_Bg_depth<<depth[n];

        fout<<"\t"<<Rho_CarbonicPhase(P[n],T[n]);
        fout_density_CarbonicPhase_depth<<"\t"<<Rho_CarbonicPhase(pres(depth[n]),temp(depth[n]));
        fout_viscosity_CarbonicPhase_depth<<"\t"<<mu_CarbonicPhase(pres(depth[n]),temp(depth[n]));
        fout_compressibility_CarbonicPhase_depth<<"\t"<<C_CarbonicPhase(pres(depth[n]),temp(depth[n]));
        fout_Zfacor_CarbonicPhase_depth<<"\t"<<Z_CarbonicPhase(pres(depth[n]),temp(depth[n]));
        fout_KCo2_depth<<"\t"<<KCo2(pres(depth[n]),temp(depth[n]));
        fout_KH2o_depth<<"\t"<<KH2o(pres(depth[n]),temp(depth[n]));
        fout_Bg_depth<<"\t"<<Bg(pres(depth[n]),temp(depth[n]));

        fout<<"\n";
        fout_density_CarbonicPhase_depth<<"\n";
        fout_viscosity_CarbonicPhase_depth<<"\n";
        fout_compressibility_CarbonicPhase_depth<<"\n";
        fout_Zfacor_CarbonicPhase_depth<<"\n";
        fout_KCo2_depth<<"\n";
        fout_KH2o_depth<<"\n";
        fout_Bg_depth<<"\n";
    }

    for (int np=0; np<1000;np++)
    {
        PP[np]=Pmin+np*(barTopa(710)-Pmin)/(1000-1);
        fout_density_CarbonicPhase_T50<< paTobar(PP[np])<<"\t"<<Rho_CarbonicPhase(PP[np],50)<<endl;
        fout_density_CarbonicPhase_T100<<paTobar(PP[np])<<"\t"<<Rho_CarbonicPhase(PP[np],100)<<endl;
        fout_viscosity_CarbonicPhase_T50<<paTobar(PP[np])<<"\t"<<mu_CarbonicPhase(PP[np],50)*1000<<endl; // in cp
        fout_viscosity_CarbonicPhase_T100<<paTobar(PP[np])<<"\t"<<mu_CarbonicPhase(PP[np],100)*1000<<endl; // in cp
        fout_viscosity_CarbonicPhase_T320K<<paTobar(PP[np])<<"\t"<<mu_CarbonicPhase(PP[np],KelvinTodegreeC(320))*1000<<endl; // in cp
        fout_compressibility_CarbonicPhase_T50<< paTobar(PP[np])<<"\t"<<C_CarbonicPhase(PP[np],50)*100000<<endl; //in bar-1
        fout_Bg_T50<< paTobar(PP[np])<<"\t"<<Bg(PP[np],50)<<endl;
    }
}





void EOS_CO2H2ONaCl_Spycher04::plot_thermodynamics()
{
    // needed for plotting
    const double Tmax(100.);
    const double Tmin(20.);
    const double Pmax(35000000.); //350 bar
    const double Pmin(100000.); // 1 bar
    const double msaltmin(0.);
    const double msaltmax(1.5); //87664.5 ppm
    const double nTP(100);
    const double nmsalt(4);
    const double max_depth(2600.); //2600 m

    ofstream fout_y_Co2_depth,fout_y_H2o_depth,fout_x_Co2_depth,
            fout_x_H2o_depth,fout_x_s_depth,
            fout_y_H2o_T50,fout_y_H2o_T100,fout_x_Co2_T50,fout_x_Co2_T100,fout_x_Co2_T50_SeaWater,
            fout_y_H2o_T50_SeaWater,fout_x_H2o_T50_SeaWater,fout_x_s_T50_SeaWater,
            fout_m_Co2_T40_SeaWater,fout_m_Co2_T60_SeaWater;
    vector<double> T(nTP),P(nTP),msalt(nmsalt),depth(nTP),PP(1000);

    //    fout.open("density_brine.txt");
    fout_y_Co2_depth.open("y_Co2_depth.txt");
    fout_y_H2o_depth.open("y_H2o_depth.txt");
    fout_x_Co2_depth.open("x_Co2_depth.txt");
    fout_x_H2o_depth.open("x_H2o_depth.txt");
    fout_x_s_depth.open("x_s_depth.txt");
    fout_y_H2o_T50.open("y_H2o_T50.txt");
    fout_y_H2o_T100.open("y_H2o_T100.txt");
    fout_x_Co2_T50.open("x_Co2_T50.txt");
    fout_x_Co2_T100.open("x_Co2_T100.txt");
    fout_y_H2o_T50_SeaWater.open("y_H2o_T50_SeaWater.txt");
    fout_x_Co2_T50_SeaWater.open("x_Co2_T50_SeaWater.txt");
    fout_x_H2o_T50_SeaWater.open("x_H2o_T50_SeaWater.txt");
    fout_x_s_T50_SeaWater.open("x_s_T50_SeaWater.txt");
    fout_m_Co2_T40_SeaWater.open("m_Co2_T40_SeaWater.txt");
    fout_m_Co2_T60_SeaWater.open("m_Co2_T60_SeaWater.txt");

    //    fout_density_Brine_depth<< "Depth ";
    fout_y_Co2_depth<< "Depth ";
    fout_y_H2o_depth<< "Depth ";
    fout_x_Co2_depth<< "Depth ";
    fout_x_H2o_depth<< "Depth ";
    fout_y_H2o_T50<< "pressure " << endl;
    fout_y_H2o_T100<< "pressure "<< endl;
    fout_x_Co2_T50<< "pressure "<< endl;
    fout_x_Co2_T100<< "pressure "<< endl;

    for (int n=0; n<nTP;n++)
    {
        T[n]=Tmin+n*(Tmax-Tmin)/(nTP-1);
        P[n]=Pmin+n*(Pmax-Pmin)/(nTP-1);
        depth[n]=n*(max_depth)/(nTP-1);
    }

    for (int nm=0; nm<nmsalt;nm++)
    {

        msalt[nm]=msaltmin+nm*(msaltmax-msaltmin)/(nmsalt-1);
        fout_y_Co2_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_y_H2o_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_x_Co2_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_x_H2o_depth<< "\t" << "msalt= " <<msalt[nm];
        fout_x_s_depth<< "\t" << "msalt= " <<msalt[nm];
    }

    fout_y_Co2_depth<<"\n";
    fout_y_H2o_depth<< "\n";
    fout_x_Co2_depth<<"\n";
    fout_x_H2o_depth<< "\n";
    fout_x_s_depth<<"\n";

    for (int n=0; n<nTP;n++)
    {
        //        fout<<T[n]<<"\t"<<paTobar(P[n]);
        fout_y_Co2_depth<<depth[n];
        fout_y_H2o_depth<<depth[n];
        fout_x_Co2_depth<<depth[n];
        fout_x_H2o_depth<<depth[n];
        fout_x_s_depth<<depth[n];

        for (int nm=0; nm<nmsalt;nm++)
        {
            //            fout<<"\t"<<Rho_brine(P[n],T[n],msalt[nm]);
            fout_y_Co2_depth<<"\t"<<y_Co2(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_y_H2o_depth<<"\t"<<y_H2o(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_x_Co2_depth<<"\t"<<x_Co2(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_x_H2o_depth<<"\t"<<x_H2o(pres(depth[n]),temp(depth[n]),msalt[nm]);
            fout_x_s_depth<<"\t"<<X_s(pres(depth[n]),temp(depth[n]),msalt[nm]);

        }
        //        fout<<"\n";
        fout_y_Co2_depth<<"\n";
        fout_y_H2o_depth<<"\n";
        fout_x_Co2_depth<<"\n";
        fout_x_H2o_depth<<"\n";
        fout_x_s_depth<<"\n";
    }


    for (int np=0; np<1000;np++)
    {
        PP[np]=Pmin+np*(barTopa(710)-Pmin)/(1000-1);
        fout_y_H2o_T50<< paTobar(PP[np])<<"\t"<<y_H2o(PP[np],50,0.0)*1000<<endl;
        fout_y_H2o_T100<<paTobar(PP[np])<<"\t"<<y_H2o(PP[np],100,0.0)*1000<<endl;
        fout_x_Co2_T50<<paTobar(PP[np])<<"\t"<<x_Co2(PP[np],50,0.0)*100<<endl;
        fout_x_Co2_T100<<paTobar(PP[np])<<"\t"<<x_Co2(PP[np],100,0.0)*100<<endl;
        fout_x_Co2_T50_SeaWater<<paTobar(PP[np])<<"\t"<<x_Co2(PP[np],50,0.62)*100<<endl;
        fout_y_H2o_T50_SeaWater<<paTobar(PP[np])<<"\t"<<y_H2o(PP[np],50,0.62)*1000<<endl;
        fout_x_H2o_T50_SeaWater<<paTobar(PP[np])<<"\t"<<x_H2o(PP[np],50,0.62)*100<<endl;
        fout_x_s_T50_SeaWater<<paTobar(PP[np])<<"\t"<<X_s(PP[np],50,0.62)*100<<endl;
        fout_m_Co2_T40_SeaWater<<paTobar(PP[np])<<"\t"<<m_Co2(PP[np],40,0.62)<<endl;
        fout_m_Co2_T60_SeaWater<<paTobar(PP[np])<<"\t"<<m_Co2(PP[np],60,0.62)<<endl;
    }
}

} //end csmp



