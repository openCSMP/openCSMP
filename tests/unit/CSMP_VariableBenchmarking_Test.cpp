//
//  CSMP_VariableBenchmarking_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 1/02/2017.
//
//  Revised: tests fit-for-purposeness of all CSMP variable types
//  for production use on large meshes (1M nodes).
//

#include "CSMP_VariableBenchmarking_Test.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

using namespace std;

namespace csmp {

void VariableBenchmarking_Test::run()
{
    // ========================================================================
    // SETUP
    // ========================================================================
    
    constexpr size_t    N              = 1'000'000;  ///< representative large mesh node count
    constexpr size_t    N_RUNS         = 5;          ///< median over this many runs
    constexpr long long MAX_NS_PER_OP  = 10LL;       ///< maximum acceptable nanoseconds per operation

    cerr << "\nVariableBenchmarking_Test::run: testing fit-for-purposeness of CSMP variable types"
         << "\n  N = " << N << " nodes, median of " << N_RUNS << " runs"
         << "\n  acceptance threshold: " << MAX_NS_PER_OP << " ns/op\n" << endl;

    // Random number generators
    default_random_engine              intgen;
    uniform_int_distribution<>         uniform_int{0, 8};  // 9 VARIABLE_FLAGs
    auto random_flag = [&]() { return static_cast<VARIABLE_FLAG>(uniform_int(intgen)); };

    default_random_engine              doublegen;
    uniform_real_distribution<double>  uniform_double{0., 1.};
    auto random_double = [&]() { return uniform_double(doublegen); };

    // Pre-generate random flags and values to avoid measuring RNG overhead
    vector<VARIABLE_FLAG> flags(N);
    vector<double>        values(N);
    for ( size_t i{0U}; i < N; ++i ) {
        flags[i]  = random_flag();
        values[i] = random_double();
    }

    // ========================================================================
    // BENCHMARK HELPER
    // Returns median nanoseconds per operation over N_RUNS
    // ========================================================================
    
    auto benchmark = [&]( const string& name, auto func ) -> long long
    {
        vector<long long> timings(N_RUNS);
        for ( size_t run{0U}; run < N_RUNS; ++run ) {
            auto t0 = chrono::high_resolution_clock::now();
            func();
            auto t1 = chrono::high_resolution_clock::now();
            timings[run] = chrono::duration_cast<chrono::nanoseconds>(t1-t0).count();
        }
        sort( timings.begin(), timings.end() );
        const long long median_total  = timings[N_RUNS/2];
        const long long ns_per_op     = median_total / static_cast<long long>(N);
        
        cerr << "  " << left << setw(40) << name
             << setw(8) << right << ns_per_op << " ns/op  "
             << ( ns_per_op <= MAX_NS_PER_OP ? "PASS" : "FAIL" ) << endl;
        
        return ns_per_op;
    };

    // Checksum to prevent dead-code elimination by the compiler
    double checksum = 0.0;

    // ========================================================================
    // 1. SCALAR VARIABLE
    // ========================================================================
    
    cerr << "ScalarVariable:" << endl;
    {
        vector<ScalarVariable> scalars(N);
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i )
            scalars[i] = makeScalar( flags[i], values[i] );
        
        const long long ns = benchmark( "assignment (makeScalar)", [&]() {
            for ( size_t i{0U}; i < N; ++i )
                scalars[i] = makeScalar( flags[i], values[i] );
        });
        
        for ( const auto& s : scalars ) checksum += s();
        _test( ns <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 2. VECTOR VARIABLE (2D)
    // ========================================================================
    
    cerr << "\nVectorVariable<2>:" << endl;
    {
        vector<VectorVariable<2>> vecs(N);
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i ) {
            VectorVariable<2> v;
            v(0) = values[i];
            v(1) = values[(i+1) % N];
            vecs[i] = v;
        }
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                vecs[i](0) = values[i];
                vecs[i](1) = values[(i+1) % N];
            }
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                VectorVariable<2> v( vecs[i] );
                checksum += v(0);
            }
        });
        
        for ( const auto& v : vecs ) checksum += v(0) + v(1);
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 3. VECTOR VARIABLE (3D)
    // ========================================================================
    
    cerr << "\nVectorVariable<3>:" << endl;
    {
        vector<VectorVariable<3>> vecs(N);
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i ) {
            vecs[i](0) = values[i];
            vecs[i](1) = values[(i+1) % N];
            vecs[i](2) = values[(i+2) % N];
        }
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                vecs[i](0) = values[i];
                vecs[i](1) = values[(i+1) % N];
                vecs[i](2) = values[(i+2) % N];
            }
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                VectorVariable<3> v( vecs[i] );
                checksum += v(0);
            }
        });
        
        for ( const auto& v : vecs ) checksum += v(0) + v(1) + v(2);
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 4. TENSOR VARIABLE (2D)
    // ========================================================================
    
    cerr << "\nTensorVariable<2>:" << endl;
    {
        vector<TensorVariable<2>> tensors(N);
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i ) {
            tensors[i](0,0) = values[i];
            tensors[i](0,1) = values[(i+1) % N];
            tensors[i](1,0) = values[(i+2) % N];
            tensors[i](1,1) = values[(i+3) % N];
        }
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                tensors[i](0,0) = values[i];
                tensors[i](0,1) = values[(i+1) % N];
                tensors[i](1,0) = values[(i+2) % N];
                tensors[i](1,1) = values[(i+3) % N];
            }
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                TensorVariable<2> t( tensors[i] );
                checksum += t(0,0);
            }
        });
        
        for ( const auto& t : tensors ) checksum += t(0,0) + t(1,1);
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 5. TENSOR VARIABLE (3D)
    // ========================================================================
    
    cerr << "\nTensorVariable<3>:" << endl;
    {
        vector<TensorVariable<3>> tensors(N);
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i )
            for ( uint32_t r{0U}; r < 3U; ++r )
                for ( uint32_t c{0U}; c < 3U; ++c )
                    tensors[i](r,c) = values[(i + r*3 + c) % N];
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i )
                for ( uint32_t r{0U}; r < 3U; ++r )
                    for ( uint32_t c{0U}; c < 3U; ++c )
                        tensors[i](r,c) = values[(i + r*3 + c) % N];
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                TensorVariable<3> t( tensors[i] );
                checksum += t(0,0);
            }
        });
        
        for ( const auto& t : tensors ) checksum += t(0,0) + t(2,2);
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 6. ARRAY VARIABLE
    // ========================================================================
    
    cerr << "\nArrayVariable (depth=4):" << endl;
    {
        constexpr uint32_t depth = 4U;
        vector<ArrayVariable> arrays(N, ArrayVariable(depth));
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i )
            for ( uint32_t d{0U}; d < depth; ++d )
                arrays[i](d) = values[(i+d) % N];
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i )
                for ( uint32_t d{0U}; d < depth; ++d )
                    arrays[i](d) = values[(i+d) % N];
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                ArrayVariable a( arrays[i] );
                checksum += a(0);
            }
        });
        
        for ( const auto& a : arrays ) checksum += a[0];
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // ========================================================================
    // 7. FLAGGED ARRAY VARIABLE
    // ========================================================================
    
    cerr << "\nFlaggedArrayVariable (depth=4):" << endl;
    {
        constexpr uint32_t depth = 4U;
        vector<FlaggedArrayVariable> farrays(N, FlaggedArrayVariable(depth));
        
        // Warm-up
        for ( size_t i{0U}; i < N; ++i )
            for ( uint32_t d{0U}; d < depth; ++d )
                farrays[i](d) = values[(i+d) % N];
        
        const long long ns_assign = benchmark( "component assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i )
                for ( uint32_t d{0U}; d < depth; ++d )
                    farrays[i](d) = values[(i+d) % N];
        });
        
        const long long ns_flag = benchmark( "flag assignment", [&]() {
            for ( size_t i{0U}; i < N; ++i )
                for ( uint32_t d{0U}; d < depth; ++d )
                    farrays[i].Flag( d, flags[i] );
        });
        
        const long long ns_copy = benchmark( "copy construction", [&]() {
            for ( size_t i{0U}; i < N; ++i ) {
                FlaggedArrayVariable fa( farrays[i] );
                checksum += fa(0);
            }
        });
        
        for ( const auto& fa : farrays ) checksum += fa[0];
        _test( ns_assign <= MAX_NS_PER_OP );
        _test( ns_flag   <= MAX_NS_PER_OP );
        _test( ns_copy   <= MAX_NS_PER_OP );
    }

    // Prevent dead-code elimination
    if ( checksum < 0.0 ) cerr << "checksum: " << checksum << endl;

    cerr << "\nVariableBenchmarking_Test::run: finished." << endl;
    cout << endl;

} // end run

} // end csmp
