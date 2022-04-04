
/*
 * This is an executable test containing a number of loops to measure
 * the performance of a compiler. Arrays' length is LEN_1D by default
 * and if you want a different array length, you should replace every
 * LEN_1D by your desired number which must be a multiple of 40. If you
 * want to increase the number of loop calls to have a longer run time
 * you have to manipulate the constant value iterations. There is a dummy
 * function called in each loop to make all computations appear required.
 * The time to execute this function is included in the time measurement
 * for the output but it is neglectable.
 *
 *  The output includes three columns:
 *    Loop:        The name of the loop
 *    Time(Sec):     The time in seconds to run the loop
 *    Checksum:    The checksum calculated when the test has run
 *
 * In this version of the codelets arrays are static type.
 *
 * All functions/loops are taken from "TEST SUITE FOR VECTORIZING COMPILERS"
 * by David Callahan, Jack Dongarra and David Levine except those whose
 * functions' name have 4 digits.
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "common.h"
#include "array_defs.h"
#include <xmmintrin.h>
#include <immintrin.h>
#include <string.h>
#include "tsvc_baseline.h"
#include "tsvc_sse.h"
#include "tsvc_avx.h"

// array definitions
__attribute__((aligned(ARRAY_ALIGNMENT))) real_t flat_2d_array[LEN_2D*LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t x[LEN_1D];

__attribute__((aligned(ARRAY_ALIGNMENT))) real_t a[LEN_1D],b[LEN_1D],c[LEN_1D],d[LEN_1D],e[LEN_1D],
        aa[LEN_2D][LEN_2D],bb[LEN_2D][LEN_2D],cc[LEN_2D][LEN_2D],tt[LEN_2D][LEN_2D];

__attribute__((aligned(ARRAY_ALIGNMENT))) int indx[LEN_1D];

real_t* __restrict__ xx;
real_t* yy;

real_t s123(struct args_t * func_args)
{

//    induction variable recognition
//    induction variable under an if
//    not vectorizable, the condition cannot be speculated

    initialise_arrays(__func__);
    s123_sse(func_args);
    return calc_checksum(__func__);
}

real_t s124(struct args_t * func_args)
{

//    induction variable recognition
//    induction variable under both sides of if (same value)

    initialise_arrays(__func__);
    s124_sse(func_args);
    return calc_checksum(__func__);
}
real_t s161(struct args_t * func_args)
{

//    control flow
//    tests for recognition of loop independent dependences
//    between statements in mutually exclusive regions.

    initialise_arrays(__func__);
    s161_sse(func_args);
    return calc_checksum(__func__);
}
real_t s1161(struct args_t * func_args)
{

//    control flow
//    tests for recognition of loop independent dependences
//    between statements in mutually exclusive regions.

    initialise_arrays(__func__);
    s1161_sse(func_args);
    return calc_checksum(__func__);
}
real_t s253(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansio assigned under if

    initialise_arrays(__func__ );
    s253_sse(func_args);
    return calc_checksum(__func__ );
}
real_t s258(struct args_t * func_args)
{

//    scalar and array expansion
//    wrap-around scalar under an if

    initialise_arrays(__func__);
    s258_sse(func_args);
    return calc_checksum(__func__);
}
real_t s271(struct args_t * func_args)
{

//    control flow
//    loop with singularity handling

    initialise_arrays(__func__);
    s271_sse(func_args);
    return calc_checksum(__func__);
}
real_t s273(struct args_t * func_args)
{

//    control flow
//    simple loop with dependent conditional

    initialise_arrays(__func__);
    s273_sse(func_args);
    return calc_checksum(__func__);
}
real_t s274(struct args_t * func_args)
{

//    control flow
//    complex loop with dependent conditional
    initialise_arrays(__func__ );
    s274_sse(func_args);
    return calc_checksum(__func__ );
}
real_t s277(struct args_t * func_args)
{

//    control flow
//    test for dependences arising from guard variable computation.

    initialise_arrays(__func__);
    s277_sse(func_args);
    return calc_checksum(__func__);
}
real_t s278(struct args_t * func_args)
{

//    control flow
//    if/goto to block if-then-else

    initialise_arrays(__func__);
    s278_sse(func_args);
    return calc_checksum(__func__);
}
real_t s2711(struct args_t * func_args)
{

//    control flow
//    semantic if removal

    initialise_arrays(__func__);
    s2711_sse(func_args);
    return calc_checksum(__func__);
}
real_t s2712(struct args_t * func_args)
{

//    control flow
//    if to elemental min

    initialise_arrays(__func__);
    s2712_sse(func_args);
    return calc_checksum(__func__);
}
real_t s314(struct args_t * func_args)
{

//    reductions
//    if to max reduction

    initialise_arrays(__func__);
    return s314_sse(func_args);
}
real_t s315(struct args_t * func_args)
{

//    reductions
//    if to max with index reductio 1 dimension

    initialise_arrays(__func__);
    return s315_sse(func_args);
}
real_t s316(struct args_t * func_args)
{

//    reductions
//    if to min reduction

    initialise_arrays(__func__);
    return s316_sse(func_args);
}
real_t s3111(struct args_t * func_args)
{

//    reductions
//    conditional sum reduction

    initialise_arrays(__func__);
    return s3111_sse(func_args);
}
real_t s3113(struct args_t * func_args)
{

//    reductions
//    maximum of absolute value

    initialise_arrays(__func__);
    return s3113_sse(func_args);
}
real_t s341(struct args_t * func_args)
{

//    packing
//    pack positive values
//    not vectorizable, value of j in unknown at each iteration

    initialise_arrays(__func__);
    s341_sse(func_args);
    return calc_checksum(__func__);
}
real_t s342(struct args_t * func_args)
{

//    packing
//    unpacking
//    not vectorizable, value of j in unknown at each iteration

    initialise_arrays(__func__);
    s342_sse(func_args);
    return calc_checksum(__func__);
}
real_t s343(struct args_t * func_args)
{

//    packing
//    pack 2-d array into one dimension
//    not vectorizable, value of k in unknown at each iteration

    initialise_arrays(__func__);
    s343_sse(func_args);
    return calc_checksum(__func__);
}
real_t s443(struct args_t * func_args)
{

//    non-logical if's
//    arithmetic if

    initialise_arrays(__func__);
    s443_sse(func_args);
    return calc_checksum(__func__);
}
real_t vif(struct args_t * func_args)
{

//    control loops
//    vector if

    initialise_arrays(__func__);
    vif_sse(func_args);
    return calc_checksum(__func__);
}
typedef real_t(*test_function_t)(struct args_t *);
void time_function(test_function_t vector_func, void * arg_info)
{
    struct args_t func_args = {.arg_info=arg_info};

    double result = vector_func(&func_args);

    double tic=func_args.t1.tv_sec+(func_args.t1.tv_usec/1000000.0);
    double toc=func_args.t2.tv_sec+(func_args.t2.tv_usec/1000000.0);

    double taken = toc-tic;

    printf("%10.3f\t%f\n", taken, result);
}
int main(int argc, char ** argv){
    int n1 = 1;
    int n3 = 1;
    int* ip;
    real_t s1,s2;
    init(&ip, &s1, &s2);
    printf("Loop \tTime(sec) \tChecksum\n");
    time_function(&s123, NULL);
    time_function(&s124, NULL);
    time_function(&s161, NULL);
    time_function(&s1161, NULL);
    time_function(&s253, NULL);
    time_function(&s258, NULL);
    time_function(&s271, NULL);
    time_function(&s273, NULL);
    time_function(&s274, NULL);
    time_function(&s277, NULL);
    time_function(&s278, NULL);
    time_function(&s2711, NULL);
    time_function(&s2712, NULL);
    time_function(&s314, NULL);
    time_function(&s315, NULL);
    time_function(&s316, NULL);
    time_function(&s3111, NULL);
    time_function(&s3113, NULL);
    time_function(&s341, NULL);
    time_function(&s342, NULL);
    time_function(&s343, NULL);
    time_function(&s443, NULL);
    time_function(&vif, NULL);
    return EXIT_SUCCESS;
}
