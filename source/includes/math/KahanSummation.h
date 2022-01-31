//
//  KahanSummation.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 29/8/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef KahanSummation_h
#define KahanSummation_h

/**

Kahan summation algorithm, also known as compensated summation,[1] significantly reduces the numerical error in the total obtained by adding a sequence of finite-precision floating-point numbers, compared to the obvious approach. This is done by keeping a separate running compensation (a variable to accumulate small errors).

In particular, simply summing n numbers in sequence has a worst-case error that grows proportional to n, and a root mean square error that grows as
n
{\sqrt {n}} for random inputs (the roundoff errors form a random walk).[2] With compensated summation, the worst-case error bound is effectively independent of n, so a large number of values can be summed with an error that only depends on the floating-point precision.[2]

See Wikipedia

- alternatively, order the values by number and then start the summation with the smallest values !

*/

/* pseudo code

function KahanSum(input)
    var sum = 0.0                    // Prepare the accumulator.
    var c = 0.0                      // A running compensation for lost low-order bits.

    for i = 1 to input.length do     // The array input has elements indexed input[1] to input[input.length].
        var y = input[i] - c         // c is zero the first time around.
        var t = sum + y              // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y            // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
        sum = t                      // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
    next i                           // Next time around, the lost low part will be added to y in a fresh attempt.

    return sum
*/

#endif /* KahanSummation_h */
