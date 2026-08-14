#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>
#include <sstream>

#include "StatisticalDistributionGenerator.h"

using namespace std;

namespace csmp {

/**
 Normally (Gaussian) distributed random numbers are generated, using the Box-Muller
 transformation.  This transformation takes two uniformly distributed deviates
 within the unit circle, and transforms them into two independently
 distributed normal deviates.  Utilizes the internal rand() function; this can
 easily be changed to use a better and faster RNG.

 The parameters passed to the function are the mean and standard deviation of
 the desired distribution.  The default values used, when no arguments are
 passed, are 0 and 1 - the standard normal distribution.

 Two functions are provided:

 The first uses the so-called polar version of the B-M transformation, using
 multiple calls to a uniform RNG to ensure the initial deviates are within the
 unit circle.  This avoids making any costly trigonometric function calls.

 The second makes only a single set of calls to the RNG, and calculates a
 position within the unit circle with two trigonometric function calls.

 The polar version is generally superior in terms of speed; however, on some
 systems, the optimization of the math libraries may result in better
 performance of the second.  Try it out on the target system to see which
 works best for you.  On my test machine (Athlon 3800+), the non-trig version
 runs at about 3x10^6 calls/s; while the trig version runs at about
 1.8x10^6 calls/s (-O2 optimization).

 Example calls:
 RandomNormalNonTrigonometric();          //returns normal deviate with mean = 0.0, std. deviation = 1.0
 RandomNormalNonTrigonometric(5.2,3.0);   //returns deviate with mean = 5.2, std. deviation = 3.0

 Dependencies - requires <cmath> for the sqrt(), sin(), and cos() calls, and a defined value for PI.

  @note "Polar" version without trigonometric calls.

*/
double StatisticalDistributionGenerator::RandomNormalNonTrigonometric( double mu, double sigma )
{
    static bool   deviateAvailable = false;     //        flag
    static double storedDeviate;               //        deviate from previous calculation
    double polar, rsquared, var1, var2;
    // If no deviate has been stored, the polar Box-Muller transformation is
    // performed, producing two independent normally-distributed random
    // deviates.  One is stored for the next round, and one is returned.
    if (!deviateAvailable) {
        // choose pairs of uniformly distributed deviates, discarding those
        // that don't fall within the unit circle
        do {
            var1 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            var2 = 2.0 * (double(rand()) / double(RAND_MAX)) - 1.0;
            rsquared = var1 * var1 + var2 * var2;
        } while ( rsquared >= 1.0 || rsquared == 0.0);
        // calculate polar tranformation for each deviate
        polar = sqrt(-2.0 * log(rsquared) / rsquared);
        // store first deviate and set flag
        storedDeviate = var1 * polar;
        deviateAvailable = true;
        // return second deviate
        return var2 * polar * sigma + mu;
    }
    // If a deviate is available from a previous call to this function, it is
    // returned, and the flag is set to false.
    else {
        deviateAvailable = false;
        return storedDeviate * sigma + mu;
    }
}


/**
 Standard version with trigonometric calls.
*/
double StatisticalDistributionGenerator::RandomNormalTrigonometric( double mu, double sigma )
{
    static bool   deviateAvailable = false;      //        flag
    static double storedDeviate;                //        deviate from previous calculation
    double      dist, angle;
    // If no deviate has been stored, the standard Box-Muller transformation is
    // performed, producing two independent normally-distributed random
    // deviates.  One is stored for the next round, and one is returned.
    if (!deviateAvailable) {
        // choose a pair of uniformly distributed deviates, one for the
        // distance and one for the angle, and perform transformations
        dist = sqrt(-2.0 * log(double(rand()) / double(RAND_MAX)) );
        angle = 2.0 * atan(1) * 4 * (double(rand()) / double(RAND_MAX));
        // calculate and store first deviate and set flag
        storedDeviate = dist*cos(angle);
        deviateAvailable = true;
        // calcaulate return second deviate
        return dist * sin(angle) * sigma + mu;
    }
    // If a deviate is available from a previous call to this function, it is
    // returned, and the flag is set to false.
    else {
        deviateAvailable = false;
        return storedDeviate * sigma + mu;
    }
}

double StatisticalDistributionGenerator::RandomNormal( double mu, double sigma )
{
    return RandomNormalNonTrigonometric( mu, sigma );
}

double StatisticalDistributionGenerator::RandomNormalLimit( double mu, double sigma, double lower, double upper )
{
    double rnd;
    do {
        rnd = RandomNormal( mu, sigma );
    } while ( rnd < lower || rnd > upper );
    return rnd;
}

double StatisticalDistributionGenerator::RandomLogNormal( double mu, double sigma )
{
    //generate normal distributed random variable with mean 0 and standard deviation 1
    double normalRnd = RandomNormalNonTrigonometric();

    double meanlog = log(mu) - 1 / 2 * log(1 + (sigma / mu) * (sigma / mu));
    double stdvlog = sqrt(log(1 + (sigma / mu) * (sigma / mu)));

    //generate lognormal distributed random variable
    double logNormalRnd = exp(meanlog + stdvlog * normalRnd);
    return logNormalRnd;
}

double StatisticalDistributionGenerator::RandomLogNormalLimit( double mu, double sigma, double lower, double upper )
{
    double rnd;
    do {
        rnd = RandomLogNormal( mu, sigma );
    } while ( rnd < lower || rnd > upper );
    return rnd;
}




/** ----------------------------------------------------------------------------
 Calculates the eigenvalues and normalized eigenvectors of a symmetric
 2x2 and 3x3 matrix A using the Jacobi algorithm.
 The upper triangular part of A is destroyed during the calculation,
 the diagonal elements are read but not destroyed, and the lower
 triangular elements are not referenced at all.
 ----------------------------------------------------------------------------
 Parameters:
    A: The symmetric input matrix
    V: Storage buffer for eigenvectors
    d: Storage buffer for eigenvalues
 ----------------------------------------------------------------------------
*/
void StatisticalDistributionGenerator::EigenDecomposition( Matrix A, Matrix& V, vector<double>& d )
{
    const size_t n = A.Rows();
    double sd, so;                  // Sums of diagonal resp. off-diagonal elements
    double s, c, t;                 // sin(phi), cos(phi), tan(phi) and temporary storage
    double g, h, z, theta;          // More temporary storage
    double thresh;

    // initialize V to the identitity matrix
    for (size_t i=0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            V(i,j) = V(j,i) = 0.;
        }
        V(i,i) = 1.0;
    }

    // initialize d to diag(A)
    for (size_t i = 0; i < n; ++i)
        d[i] = A(i,i);

    // calculate (tr(A)^2)
    sd = 0.;
    for (size_t i = 0; i < n; ++i)
        sd += fabs(d[i]);
    sd = sd * sd;

    // main iteration loop
    for (size_t nIter = 0; nIter < 50; nIter++) {
        // test for convergence
        so = 0.;
        for (size_t p = 0; p < n; ++p)
            for (size_t q = p + 1; q < n; ++q)
                so += fabs(A(p,q));
            if (so == 0.) { // method converged
                // sort eigenvalues and corresponding vectors.
                for (size_t i = 0; i < n - 1; i++) {
                    size_t k = i;
                    double p = d[i];
                    for (size_t j = i + 1; j < n; j++) {
                        if (d[j] < p) {
                            k = j;
                            p = d[j];
                        }
                    }
                    if (k != i) {
                        d[k] = d[i];
                        d[i] = p;
                        for (size_t j = 0; j < n; j++) {
                            p = V(j,i);
                            V(j,i) = V(j,k);
                            V(j,k) = p;
                        }
                    }
                }
                return;
            }
        if (nIter < 4)
            thresh = 0.2 * so / (n * n);
        else
            thresh = 0.0;
        // do sweep
        for (size_t p = 0; p < n; ++p) {
            for (size_t q = p+1; q < n; ++q) {
                g = 100.0 * fabs(A(p,q));
                if ( nIter > 4  &&  fabs(d[p]) + g == fabs(d[p]) &&  fabs(d[q]) + g == fabs(d[q]) ) {
                  A(p,q) = 0.0;
                }
                else if ( fabs(A(p,q)) > thresh ) {
                    // calculate Jacobi transformation
                    h = d[q] - d[p];
                    if ( fabs(h) + g == fabs(h) ) {
                        t = A(p,q) / h;
                    }
                    else {
                        theta = 0.5 * h / A(p,q);
                        if (theta < 0.)
                            t = -1. / (sqrt(1. + (theta * theta)) - theta);
                        else
                            t = 1. / (sqrt(1. + (theta * theta)) + theta);
                    }
                    c = 1. / sqrt(1. + (t * t));
                    s = t * c;
                    z = t * A(p,q);

                    // apply Jacobi transformation
                    A(p,q) = 0.;
                    d[p] -= z;
                    d[q] += z;
                    for (size_t r = 0; r < p; r++) {
                        t = A(r,p);
                        A(r,p) = c*t - s*A(r,q);
                        A(r,q) = s*t + c*A(r,q);
                    }
                    for (size_t r = p + 1; r < q; r++) {
                        t = A(p,r);
                        A(p,r) = c*t - s*A(r,q);
                        A(r,q) = s*t + c*A(r,q);
                    }
                    for (size_t r = q + 1; r < n; r++) {
                        t = A(p,r);
                        A(p,r) = c*t - s*A(q,r);
                        A(q,r) = s*t + c*A(q,r);
                    }
                    // update eigenvectors
                    for (size_t r = 0; r < n; r++) {
                        t = V(r,p);
                        V(r,p) = c*t - s*V(r,q);
                        V(r,q) = s*t + c*V(r,q);
                    }
                }
            }
        }
    }
    cout << "StatisticalDistributionGenerator::EigenDecomposition" << endl;
    cout << "Jacobi's eigen decomposition method doesn't converge" << endl;
}



void StatisticalDistributionGenerator::Histogram( vector<double>& prop , const size_t no_bins )
{
    double min = *min_element(prop.begin(),prop.end());
    double max = *max_element(prop.begin(),prop.end());
    vector<double> temp_bins;
    double del_bin = (max - min) / no_bins;
    for (size_t i = 0; i < no_bins + 1; ++i) {
        temp_bins.push_back(min + del_bin * i);
    }
    for(size_t i = 0; i < no_bins; ++i) {
        auto count = temp_bins.at(i) + del_bin / 2;
        hist_bins_.push_back( count );
        // Siroos: count_left  = count_if(prop.begin(), prop.end(), bind2nd(less<double>(), hist_bins_.at(i) - del_bin / 2));
        long count_left  = count_if( prop.begin(), prop.end(), [count]( auto const& x ){ return x < count; } );
        long count_right = ( i == no_bins - 1 ) ?
                count_if( prop.begin(), prop.end(), [count]( auto const& x ){ return x <= count; } ) :
                count_if( prop.begin(), prop.end(), [count]( auto const& x ){ return x  < count; } );

        hist_freqs_.push_back( static_cast<size_t>(abs(count_right - count_left)) );
    }
    temp_bins.clear();
}



void StatisticalDistributionGenerator::OutputHistogramToFile( string fileName )
{
    string extension = ".txt";
    ofstream file(fileName.c_str() + extension);
    for (size_t it = 0; it < hist_freqs_.size(); ++it)
    {
        file << setw(20) << hist_bins_.at(it) << setw(20) << hist_freqs_.at(it);
        file << endl;
    }
    file.close();

}

vector<double> StatisticalDistributionGenerator::Bins() const
{
    return hist_bins_;
}

vector<size_t> StatisticalDistributionGenerator::Frequencies() const 
{
    return hist_freqs_;
}

} // csmp
