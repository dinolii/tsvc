
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
#include "tsvc_sse.h"

void s123_sse(struct args_t * func_args)
{
//    induction variable recognition
//    induction variable under an if
//    not vectorizable, the condition cannot be speculated

    gettimeofday(&func_args->t1, NULL);

    int j;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        j = -1;
        int i = 0;
        int upper_bound = (LEN_1D/2)/vf*vf;
        for (; i < upper_bound; i+=vf) {
            j++;
            a[j] = b[i] + d[i] * e[i];
            if(c[i] > (real_t)0. && c[i+1] > (real_t) 0. && c[i+2] > (real_t)0. && c[i+3] > (real_t)0.){
                j++;
                a[j] = c[i] + d[i] * e[i];
                j++;
                a[j] = b[i+1] + d[i+1] * e[i+1];
                j++;
                a[j] = c[i+1] + d[i+1] * e[i+1];
                j++;
                a[j] = b[i+2] + d[i+2] * e[i+2];
                j++;
                a[j] = c[i+2] + d[i+2] * e[i+2];
                j++;
                a[j] = b[i+3] + d[i+3] * e[i+3];
                j++;
                a[j] = c[i+3] + d[i+3] * e[i+3];
            }
            else if(!(c[i] > (real_t)0.) && !(c[i+1] > (real_t) 0.) && !(c[i+2] > (real_t)0.) && !(c[i+3] > (real_t)0.)){
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                _mm_store_ps(&a[j], _mm_add_ps(_mm_load_ps(&b[i]), _mm_mul_ps(d_, e_)));
                j += 3;
            }
            else{
                if (c[i] > (real_t)0.) {
                    j++;
                    a[j] = c[i] + d[i] * e[i];
                }
                j++;
                a[j] = b[i+1] + d[i+1] * e[i+1];
                if (c[i+1] > (real_t)0.) {
                    j++;
                    a[j] = c[i+1] + d[i+1] * e[i+1];
                }
                j++;
                a[j] = b[i+2] + d[i+2] * e[i+2];
                if (c[i+2] > (real_t)0.) {
                    j++;
                    a[j] = c[i+2] + d[i+2] * e[i+2];
                }
                j++;
                a[j] = b[i+3] + d[i+3] * e[i+3];
                if (c[i+3] > (real_t)0.) {
                    j++;
                    a[j] = c[i+3] + d[i+3] * e[i+3];
                }
            }
        }
        for(; i < (LEN_1D/2); i++){
            j++;
            a[j] = b[i] + d[i] * e[i];
            if (c[i] > (real_t)0.) {
                j++;
                a[j] = c[i] + d[i] * e[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s124_sse(struct args_t * func_args)
{

//    induction variable recognition
//    induction variable under both sides of if (same value)

    gettimeofday(&func_args->t1, NULL);

    int j;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        j = -1;
        int i = 0;
        int upper_bound = (LEN_1D / vf * vf);
        for (; i < upper_bound; i+=vf) {
            if(b[i] > (real_t)0. && b[i+1] > (real_t)0. && b[i+2] > (real_t)0. && b[i+3] > (real_t)0.){
                j++;
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                _mm_store_ps(&a[j], _mm_add_ps(b_, _mm_mul_ps(d_, e_)));
                j += 3;
            }
            else if(!(b[i] > (real_t)0.) && !(b[i+1] > (real_t)0.) && !(b[i+2] > (real_t)0.) && !(b[i+3] > (real_t)0.)){
                j++;
                __m128 c_ = _mm_load_ps(&b[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                _mm_store_ps(&a[j], _mm_add_ps(c_, _mm_mul_ps(d_, e_)));
                j += 3;
            }
            else{
                if (b[i] > (real_t)0.) {
                    j++;
                    a[j] = b[i] + d[i] * e[i];
                } else {
                    j++;
                    a[j] = c[i] + d[i] * e[i];
                }
                if (b[i+1] > (real_t)0.) {
                    j++;
                    a[j] = b[i+1] + d[i+1] * e[i+1];
                } else {
                    j++;
                    a[j] = c[i+1] + d[i+1] * e[i+1];
                }
                if (b[i+2] > (real_t)0.) {
                    j++;
                    a[j] = b[i+2] + d[i+2] * e[i+2];
                } else {
                    j++;
                    a[j] = c[i+2] + d[i+2] * e[i+2];
                }
                if (b[i+3] > (real_t)0.) {
                    j++;
                    a[j] = b[i+3] + d[i+3] * e[i+3];
                } else {
                    j++;
                    a[j] = c[i+3] + d[i+3] * e[i+3];
                }
            }
        }
        for(; i < LEN_1D; i++){
            if (b[i] > (real_t)0.) {
                j++;
                a[j] = b[i] + d[i] * e[i];
            } else {
                j++;
                a[j] = c[i] + d[i] * e[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s161_sse(struct args_t * func_args)
{

//    control flow
//    tests for recognition of loop independent dependences
//    between statements in mutually exclusive regions.

    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for(int nl = 0; nl < iterations/2; nl++){
        int i =0;
        int upper_bound = (LEN_1D - 1)/vf*vf;
        for(; i < upper_bound; i+=vf){
            if(b[i] < (real_t)0. && b[i+1] < (real_t)0. && b[i+2] < (real_t)0. && b[i+3] < (real_t)0.){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                _mm_store_ps(&c[i+1], _mm_add_ps(a_, _mm_mul_ps(d_, d_)));
            }
            else if(!(b[i] < (real_t)0.) && !(b[i+1] < (real_t)0.) && !(b[i+2] < (real_t)0.) && !(b[i+3] < (real_t)0.)){
                __m128 c_ = _mm_load_ps(&c[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                _mm_store_ps(&a[i], _mm_add_ps(c_, _mm_mul_ps(d_, e_)));
            }
            else{
                if(b[i] < (real_t)0.){
                    c[i+1] = a[i] + d[i] * d[i];
                }
                else{
                    a[i] = c[i] + d[i] * e[i];
                }
                if(b[i+1] < (real_t)0.){
                    c[i+2] = a[i+1] + d[i+1] * d[i+1];
                }
                else{
                    a[i+1] = c[i+1] + d[i+1] * e[i+1];
                }
                if(b[i+2] < (real_t)0.){
                    c[i+3] = a[i+2] + d[i+2] * d[i+2];
                }
                else{
                    a[i+2] = c[i+2] + d[i+2] * e[i+2];
                }
                if(b[i+3] < (real_t)0.){
                    c[i+4] = a[i+3] + d[i+3] * d[i+3];
                }
                else{
                    a[i+3] = c[i+3] + d[i+3] * e[i+3];
                }
            }
        }
        for(; i < LEN_1D - 1; ++i){
            if(b[i] < (real_t)0.){
                c[i+1] = a[i] + d[i] * d[i];
            }
            else{
                a[i] = c[i] + d[i] * e[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s1161_sse(struct args_t * func_args){
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for(int nl = 0; nl < iterations; nl++){
        int i =0;
        int upper_bound = (LEN_1D - 1)/vf*vf;
        for(; i < upper_bound; i+=vf){
            if(c[i] < (real_t)0. && c[i+1] < (real_t)0. && c[i+2] < (real_t)0. && c[i+3] < (real_t)0.){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                _mm_store_ps(&b[i], _mm_add_ps(a_, _mm_mul_ps(d_, d_)));
            }
            else if(!(c[i] < (real_t)0.) && !(c[i+1] < (real_t)0.) && !(c[i+2] < (real_t)0.) && !(c[i+3] < (real_t)0.)){
                __m128 c_ = _mm_load_ps(&c[i]);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                _mm_store_ps(&a[i], _mm_add_ps(c_, _mm_mul_ps(d_, e_)));
            }
            else{
                if(c[i] < (real_t)0.){
                    b[i] = a[i] + d[i] * d[i];
                }
                else{
                    a[i] = c[i] + d[i] * e[i];
                }
                if(c[i+1] < (real_t)0.){
                    b[i+1] = a[i+1] + d[i+1] * d[i+1];
                }
                else{
                    a[i+1] = c[i+1] + d[i+1] * e[i+1];
                }
                if(c[i+2] < (real_t)0.){
                    b[i+2] = a[i+2] + d[i+2] * d[i+2];
                }
                else{
                    b[i+2] = c[i+2] + d[i+2] * e[i+2];
                }
                if(c[i+3] < (real_t)0.){
                    c[i+3] = a[i+3] + d[i+3] * d[i+3];
                }
                else{
                    a[i+3] = c[i+3] + d[i+3] * e[i+3];
                }
            }
        }
        for(; i < LEN_1D - 1; ++i){
            if(c[i] < (real_t)0.){
                b[i] = a[i] + d[i] * d[i];
            }
            else{
                a[i] = c[i] + d[i] * e[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);

}
void s253_sse(struct args_t * func_args)
{

//    scalar and array expansion
//    scalar expansio assigned under if

    gettimeofday(&func_args->t1, NULL);

    real_t s;
    int vf = 4;
    int upper_bound = LEN_1D / vf * vf;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        for (; i < upper_bound; i+=vf) {
            __m128 a_ = _mm_load_ps(&a[i]);
            __m128 b_ = _mm_load_ps(&b[i]);
            if(a[i] > b[i] && a[i+1] > b[i+1] && a[i+2] > b[i+2] && a[i+3] > b[i+3]){
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 s_ = _mm_sub_ps(a_, _mm_mul_ps(b_, d_));
                __m128 c_ = _mm_load_ps(&c[i]);
                __m128 temp = _mm_add_ps(c_, s_);
                _mm_store_ps(&c[i], temp);
                _mm_store_ps(&a[i], s_);
            }
            else if(!(a[i+0] > b[i+0]) && !(a[i+1] > b[i+1]) && !(a[i+2] > b[i+2]) && !(a[i+3] > b[i+3])){

            }
            else{
                if (a[i] > b[i]) {
                    s = a[i] - b[i] * d[i];
                    c[i] += s;
                    a[i] = s;
                }
                if (a[i+1] > b[i+1]) {
                    s = a[i+1] - b[i+1] * d[i+1];
                    c[i+1] += s;
                    a[i+1] = s;
                }
                if (a[i+2] > b[i+2]) {
                    s = a[i+2] - b[i+2] * d[i+2];
                    c[i+2] += s;
                    a[i+2] = s;
                }
                if (a[i+3] > b[i+3]) {
                    s = a[i+3] - b[i+3] * d[i+3];
                    c[i+3] += s;
                    a[i+3] = s;
                }
            }
        }
        for(; i < LEN_1D; i++){
            if (a[i] > b[i]) {
                s = a[i] - b[i] * d[i];
                c[i] += s;
                a[i] = s;
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s258_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);

    real_t s;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        s = 0.;
        int i = 0;
        int upper_bound = LEN_2D / vf * vf;
        for (; i < upper_bound; i += vf) {
            __m128 aa_ = _mm_load_ps(&(aa[0][i]));
            __m128 c_ = _mm_load_ps(&c[i]);
            __m128 d_ = _mm_load_ps(&d[i]);
            __m128 s_ = _mm_load1_ps(&s);
            __m128 one_ = _mm_set1_ps(1.);
            if(a[i] > 0 && a[i+1] > 0 && a[i+2] > 0 && a[i+3] > 0){
                s_ = _mm_mul_ps(d_, d_);
                s = d[i+3] * d[i+3];
                _mm_store_ps(&b[i], _mm_add_ps(_mm_mul_ps(s_, c_), d_));
                _mm_store_ps(&e[i], _mm_mul_ps(_mm_add_ps(s_, one_), aa_));
            }
            else if(!(a[i] > 0) && !(a[i+1] > 0) && !(a[i+2] > 0) && !(a[i+3] > 0)){
                _mm_store_ps(&b[i], _mm_add_ps(_mm_mul_ps(s_, c_), d_));
                _mm_store_ps(&e[i], _mm_mul_ps(_mm_add_ps(s_, one_), aa_));
            }
            else{
                if (a[i] > 0.) {
                    s = d[i] * d[i];
                }
                b[i] = s * c[i] + d[i];
                e[i] = (s + (real_t)1.) * aa[0][i];
                if (a[i+1] > 0.) {
                    s = d[i+1] * d[i+1];
                }
                b[i+1] = s * c[i+1] + d[i+1];
                e[i+1] = (s + (real_t)1.) * aa[0][i+1];
                if (a[i+2] > 0.) {
                    s = d[i+2] * d[i+2];
                }
                b[i+2] = s * c[i+2] + d[i+2];
                e[i+2] = (s + (real_t)1.) * aa[0][i+2];
                if (a[i+3] > 0.) {
                    s = d[i+3] * d[i+3];
                }
                b[i+3] = s * c[i+3] + d[i+3];
                e[i+3] = (s + (real_t)1.) * aa[0][i+3];
            }
        }
        for(; i < LEN_2D; ++i){
            if (a[i] > 0.) {
                s = d[i] * d[i];
            }
            b[i] = s * c[i] + d[i];
            e[i] = (s + (real_t)1.) * aa[0][i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s271_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < 4*iterations; nl++) {
        int i = 0;
        int upper_bound = (LEN_1D / vf * vf);
        for (; i < upper_bound; i+=vf) {
            if (b[i] > (real_t)0. && b[i+1] > (real_t)0. && b[i+2] > (real_t)0. && b[i+3] > (real_t)0.) {
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 c_ = _mm_load_ps(&c[i]);
                _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(b_, c_)));
            }
            else if(!(b[i] > (real_t)0.) && !(b[i+1] > (real_t)0.) && !(b[i+2] > (real_t)0.) && !(b[i+3] > (real_t)0.)){

            }
            else{
                if (b[i] > (real_t)0.) {
                    a[i] += b[i] * c[i];
                }
                if (b[i+1] > (real_t)0.) {
                    a[i+1] += b[i+1] * c[i+1];
                }
                if (b[i+2] > (real_t)0.) {
                    a[i+2] += b[i+2] * c[i+2];
                }
                if (b[i+3] > (real_t)0.) {
                    a[i+3] += b[i+3] * c[i+3];
                }
            }
        }
        for(;i<LEN_1D;i++){
            if (b[i] > (real_t)0.) {
                a[i] += b[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s273_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D / vf * vf;
        for (; i < upper_bound; i+=vf) {
            __m128 a_ = _mm_load_ps(&a[i]);
            __m128 d_ = _mm_load_ps(&d[i]);
            __m128 e_ = _mm_load_ps(&e[i]);
            _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(d_, e_)));
            if(a[i] < (real_t)0. && a[i+1] < (real_t)0. && a[i+2] < (real_t)0. && a[i+3] < (real_t)0.){
                __m128 b_ = _mm_load_ps(&b[i]);
                _mm_store_ps(&b[i], _mm_add_ps(b_, _mm_mul_ps(d_, e_)));
            }
            else if(!(a[i] < (real_t)0.) && !(a[i+1] < (real_t)0.) && !(a[i+2] < (real_t)0.) && !(a[i+3] < (real_t)0.)){

            }
            else{
                if (a[i] < (real_t)0.)
                    b[i] += d[i] * e[i];
                if (a[i+1] < (real_t)0.)
                    b[i+1] += d[i+1] * e[i+1];
                if (a[i+2] < (real_t)0.)
                    b[i+2] += d[i+2] * e[i+2];
                if (a[i+3] < (real_t)0.)
                    b[i+3] += d[i+3] * e[i+3];
            }
            __m128 c_ = _mm_load_ps(&c[i]);
            a_ = _mm_load_ps(&a[i]);
            _mm_store_ps(&c[i], _mm_add_ps(c_, _mm_mul_ps(a_, d_)));
        }
        for(; i < LEN_1D; i++){
            a[i] += d[i] * e[i];
            if (a[i] < (real_t)0.)
                b[i] += d[i] * e[i];
            c[i] += a[i] * d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s274_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for (;i < upper_bound; i+=vf) {
            __m128 c_ = _mm_load_ps(&c[i]);
            __m128 d_ = _mm_load_ps(&d[i]);
            __m128 e_ = _mm_load_ps(&e[i]);
            __m128 mul = _mm_mul_ps(e_, d_);
            __m128 res = _mm_add_ps(c_, mul);
            _mm_store_ps(&a[i], res);
            if(a[i] > (real_t)0. && a[i+1] > (real_t)0. && a[i+2] > (real_t)0. && a[i+3] > (real_t)0.){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                _mm_store_ps(&b[i], _mm_add_ps(a_, b_));
            }
            else if(!(a[i] > (real_t) 0.) && !(a[i+1] > (real_t)0.) && !(a[i+2] > (real_t)0.) && !(a[i+3] > (real_t)0.)){
                _mm_store_ps(&a[i], mul);
            }
            else{
                if(a[i] > (real_t)0.){
                    b[i] = a[i] + b[i];
                }
                else{
                    a[i] = d[i] * e[i];
                }
                if(a[i+1] > (real_t)0.){
                    b[i+1] = a[i+1] + b[i+1];
                }
                else{
                    a[i+1] = d[i+1] * e[i+1];
                }
                if(a[i+2] > (real_t)0.){
                    b[i+2] = a[i+2] + b[i+2];
                }
                else{
                    a[i+2] = d[i+2] * e[i+2];
                }
                if(a[i+3] > (real_t)0.){
                    b[i+3] = a[i+3] + b[i+3];
                }
                else{
                    a[i+3] = d[i+3] * e[i+3];
                }
            }
        }
        for(;i<LEN_1D;i++){
            a[i] = c[i] + e[i] * d[i];
            if (a[i] > (real_t)0.) {
                b[i] = a[i] + b[i];
            } else {
                a[i] = d[i] * e[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s277_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        int upper_bound = (LEN_1D - 1) / vf * vf;
        for(; i < upper_bound; i+=vf){
            if(a[i] >= (real_t)0. && a[i+1] >= (real_t)0. && a[i+2] >= (real_t)0. && a[i+3] >= (real_t)0.){

            }
            else if(!(a[i] >= (real_t)0.) && !(a[i+1] >= (real_t)0.) && !(a[i+2] >= (real_t)0.) && !(a[i+3] >= (real_t)0.)){
                if(b[i] >= (real_t)0. && b[i+1] >= (real_t)0. && b[i+2] >= (real_t)0. && b[i+3] >= (real_t)0.){
                    __m128 c_ = _mm_load_ps(&c[i]);
                    __m128 d_ = _mm_load_ps(&d[i]);
                    __m128 e_ = _mm_load_ps(&e[i]);
                    _mm_store_ps(&b[i+1], _mm_add_ps(c_, _mm_mul_ps(d_, e_)));
                }
                else if(!(b[i] >= (real_t)0.) && !(b[i+1] >= (real_t)0.) && !(b[i+2] >= (real_t)0.) && !(b[i+3] >= (real_t)0.)){
                    __m128 a_ = _mm_load_ps(&a[i]);
                    __m128 c_ = _mm_load_ps(&c[i]);
                    __m128 d_ = _mm_load_ps(&d[i]);
                    _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(c_, d_)));
                }
                else{
                    if(b[i] >= (real_t)0.){
                        b[i+1] = c[i] + d[i] * e[i];
                    }
                    else{
                        a[i] += c[i] * d[i];
                    }
                    if(b[i+1] >= (real_t)0.){
                        b[i+2] = c[i+1] + d[i+1] * e[i+1];
                    }
                    else{
                        a[i+1] += c[i+1] * d[i+1];
                    }
                    if(b[i+2] >= (real_t)0.){
                        b[i+3] = c[i+2] + d[i+2] * e[i+2];
                    }
                    else{
                        a[i+2] += c[i+2] * d[i+2];
                    }
                    if(b[i+3] >= (real_t)0.){
                        b[i+4] = c[i+3] + d[i+3] * e[i+3];
                    }
                    else{
                        a[i+3] += c[i+3] * d[i+3];
                    }
                }

            }
            else{
                if(a[i] >= (real_t)0.){

                }
                else{
                    if(b[i] >= (real_t)0.){
                        b[i+1] = c[i] + d[i] * e[i];
                    }
                    else{
                        a[i] += c[i] * d[i];
                    }
                }
                if(a[i+1] >= (real_t)0.){

                }
                else{
                    if(b[i+1] >= (real_t)0.){
                        b[i+2] = c[i+1] + d[i+1] * e[i+1];
                    }
                    else{
                        a[i+1] += c[i+1] * d[i+1];
                    }
                }
                if(a[i+2] >= (real_t)0.){

                }
                else{
                    if(b[i+2] >= (real_t)0.){
                        b[i+3] = c[i+2] + d[i+2] * e[i+2];
                    }
                    else{
                        a[i+2] += c[i+2] * d[i+2];
                    }
                }
                if(a[i+3] >= (real_t)0.){

                }
                else{
                    if(b[i+3] >= (real_t)0.){
                        b[i+4] = c[i+3] + d[i+3] * e[i+3];
                    }
                    else{
                        a[i+3] += c[i+3] * d[i+3];
                    }
                }
            }
        }
        for (; i < LEN_1D-1; i++) {
            if (a[i] >= (real_t)0.) {
                goto L20;
            }
            if (b[i] >= (real_t)0.) {
                goto L30;
            }
            a[i] += c[i] * d[i];
            L30:
            b[i+1] = c[i] + d[i] * e[i];
            L20:
            ;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s278_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for(;i<upper_bound;i+=vf){
            if(a[i] > (real_t)0. && a[i+1] > (real_t)0. && a[i+2] > (real_t)0. && a[i+3] > (real_t)0.){
                __m128 c_ = _mm_load_ps(&c[i]);
                __m128 neg_c = _mm_sub_ps(_mm_setzero_ps(), c_);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                __m128 res_ = _mm_add_ps(neg_c, _mm_mul_ps(d_, e_));
                _mm_store_ps(&c[i], res_);
                __m128 b_ = _mm_load_ps(&b[i]);
                _mm_store_ps(&a[i], _mm_add_ps(b_, _mm_mul_ps(res_, d_)));

            }
            else if(!(a[i] > (real_t)0.) && !(a[i+1] > (real_t)0.) && !(a[i+2] > (real_t)0.) && !(a[i+3] > (real_t)0.)){
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 neg_b = _mm_sub_ps(_mm_setzero_ps(), b_);
                __m128 d_ = _mm_load_ps(&d[i]);
                __m128 e_ = _mm_load_ps(&e[i]);
                __m128 res = _mm_add_ps(neg_b, _mm_mul_ps(d_, e_));
                _mm_store_ps(&b[i], res);
                __m128 c_ = _mm_load_ps(&c[i]);
                b_ = _mm_load_ps(&b[i]);
                _mm_store_ps(&a[i], _mm_add_ps(b_, _mm_mul_ps(c_, d_)));
            }
            else{
                if(a[i] > (real_t)0.){
                    c[i] = -c[i] + d[i] * e[i];
                }
                else{
                    b[i] = -b[i] + d[i] * e[i];
                    a[i] = b[i] + c[i] * d[i];
                }
                if(a[i+1] > (real_t)0.){
                    c[i+1] = -c[i+1] + d[i+1] * e[i+1];
                }
                else{
                    b[i+1] = -b[i+1] + d[i+1] * e[i+1];
                    a[i+1] = b[i+1] + c[i+1] * d[i+1];
                }
                if(a[i+2] > (real_t)0.){
                    c[i+2] = -c[i+2] + d[i+2] * e[i+2];
                }
                else{
                    b[i+2] = -b[i+2] + d[i+2] * e[i+2];
                    a[i+2] = b[i+2] + c[i+2] * d[i+2];
                }
                if(a[i+3] > (real_t)0.){
                    c[i+3] = -c[i+3] + d[i+3] * e[i+3];
                }
                else{
                    b[i+3] = -b[i+3] + d[i+3] * e[i+3];
                    a[i+3] = b[i+3] + c[i+3] * d[i+3];
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if (a[i] > (real_t)0.) {
                goto L20;
            }
            b[i] = -b[i] + d[i] * e[i];
            goto L30;
            L20:
            c[i] = -c[i] + d[i] * e[i];
            L30:
            a[i] = b[i] + c[i] * d[i];
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s2711_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < 4*iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for (; i < upper_bound; i+=vf) {
            if(b[i] != (real_t)0.0 && b[i+1] != (real_t)0.0 && b[i+2]!=(real_t)0.0 && b[i+3]!=(real_t)0.0){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 c_ = _mm_load_ps(&c[i]);
                _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(b_, c_)));
            }
            else if(!(b[i] != (real_t)0.0) && !(b[i+1] != (real_t)0.0) && !(b[i+2]!=(real_t)0.0) && !(b[i+3]!=(real_t)0.0)){

            }
            else{
                if (b[i] != (real_t)0.0) {
                    a[i] += b[i] * c[i];
                }
                if (b[i+1] != (real_t)0.0) {
                    a[i+1] += b[i+1] * c[i+1];
                }
                if (b[i+2] != (real_t)0.0) {
                    a[i+2] += b[i+2] * c[i+2];
                }
                if (b[i+3] != (real_t)0.0) {
                    a[i+3] += b[i+3] * c[i+3];
                }
            }
        }
        for(;i<LEN_1D;i++){
            if (b[i] != (real_t)0.0) {
                a[i] += b[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s2712_sse(struct args_t * func_args)
{

//    control flow
//    if to elemental min
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < 4*iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D / vf * vf;
        for (; i < upper_bound; i+=vf) {
            if(a[i] > b[i] && a[i+1] > b[i+1] && a[i+2] > b[i+2] && a[i+3] > b[i+3]){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 c_ = _mm_load_ps(&c[i]);
                _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(b_, c_)));
            }
            else if(!(a[i] > b[i]) && !(a[i+1] > b[i+1]) && !(a[i+2] > b[i+2]) & !(a[i+3] > b[i+3])){

            }
            else{
                if (a[i] > b[i]) {
                    a[i] += b[i] * c[i];
                }
                if (a[i+1] > b[i+1]) {
                    a[i+1] += b[i+1] * c[i+1];
                }
                if (a[i+2] > b[i+2]) {
                    a[i+2] += b[i+2] * c[i+2];
                }
                if (a[i+3] > b[i+3]) {
                    a[i+3] += b[i+3] * c[i+3];
                }
            }
        }
        for(; i<LEN_1D;i++){
            if (a[i] > b[i]) {
                a[i] += b[i] * c[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
real_t s314_sse(struct args_t * func_args)
{

//    reductions
//    if to max reduction

    gettimeofday(&func_args->t1, NULL);

    real_t x;
    int vf = 4;
    for (int nl = 0; nl < iterations*5; nl++) {
        x = a[0];
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for (; i < upper_bound; i+=vf) {
            if (a[i] > x && a[i+1] > x && a[i+2] > x && a[i+3] > x) {
                if(a[i] > a[i+1] && a[i] > a[i+2] && a[i] > a[i+3]){
                    x = a[i];
                }
                else if(a[i+1] > a[i] && a[i+1] > a[i+2] && a[i+1] > a[i+3]){
                    x = a[i+1];
                }
                else if(a[i+2] > a[i] && a[i+2] > a[i+1] && a[i+2] > a[i+3]){
                    x = a[i+2];
                }
                else if(a[i+3] > a[i] && a[i+3] > a[i+1] && a[i+3] > a[i+2]){
                    x = a[i+3];
                }
            }
            else if(!(a[i] > x) && !(a[i+1] > x) && !(a[i+2] > x) && !(a[i+3] > x)){

            }
            else{
                if (a[i] > x) {
                    x = a[i];
                }
                if (a[i+1] > x) {
                    x = a[i+1];
                }
                if (a[i+2] > x) {
                    x = a[i+2];
                }
                if (a[i+3] > x) {
                    x = a[i+3];
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if (a[i] > x) {
                x = a[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, x);
    }

    gettimeofday(&func_args->t2, NULL);
    return x;
}
real_t s315_sse(struct args_t * func_args)
{

//    reductions
//    if to max with index reductio 1 dimension

    gettimeofday(&func_args->t1, NULL);

    for (int i = 0; i < LEN_1D; i++)
        a[i] = (i * 7) % LEN_1D;

    real_t x, chksum;
    int index;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        x = a[0];
        index = 0;
        int i = 0;
        int upper_bound = LEN_1D / vf * vf;
        for(; i < upper_bound; i+= vf){
            if (a[i] > x && a[i+1] > x && a[i+2] > x && a[i+3] > x) {
                if(a[i] > a[i+1] && a[i] > a[i+2] && a[i] > a[i+3]){
                    x = a[i];
                    index = i;
                }
                else if(a[i+1] > a[i] && a[i+1] > a[i+2] && a[i+1] > a[i+3]){
                    x = a[i+1];
                    index = i+1;
                }
                else if(a[i+2] > a[i] && a[i+2] > a[i+1] && a[i+2] > a[i+3]){
                    x = a[i+2];
                    index = i + 2;
                }
                else if(a[i+3] > a[i] && a[i+3] > a[i+1] && a[i+3] > a[i+2]){
                    x = a[i+3];
                    index = i + 3;
                }
            }
            else if(!(a[i] > x) && !(a[i+1] > x) && !(a[i+2] > x) && !(a[i+3] > x)){

            }
            else{
                if (a[i] > x) {
                    x = a[i];
                    index = i;
                }
                if (a[i+1] > x) {
                    x = a[i+1];
                    index = i + 1;
                }
                if (a[i+2] > x) {
                    x = a[i+2];
                    index = i + 2;
                }
                if (a[i+3] > x) {
                    x = a[i+3];
                    index = i + 3;
                }
            }
        }
        for (; i < LEN_1D; ++i) {
            if (a[i] > x) {
                x = a[i];
                index = i;
            }
        }
        chksum = x + (real_t) index;
        dummy(a, b, c, d, e, aa, bb, cc, chksum);
    }

    gettimeofday(&func_args->t2, NULL);
    return index + x + 1;
}
real_t s316_sse(struct args_t * func_args)
{

//    reductions
//    if to min reduction

    gettimeofday(&func_args->t1, NULL);

    real_t x;
    int vf = 4;
    for (int nl = 0; nl < iterations*5; nl++) {
        x = a[0];
        int i = 1;
        int upper_bound = LEN_1D/vf*vf;
        for(; i < upper_bound; i+=vf){
            if (a[i] < x && a[i+1] < x && a[i+2] < x && a[i+3] < x) {
                if(a[i] < a[i+1] && a[i] < a[i+2] && a[i] < a[i+3]){
                    x = a[i];
                }
                else if(a[i+1] < a[i] && a[i+1] < a[i+2] && a[i+1] < a[i+3]){
                    x = a[i+1];
                }
                else if(a[i+2] < a[i] && a[i+2] < a[i+1] && a[i+2] < a[i+3]){
                    x = a[i+2];
                }
                else if(a[i+3] < a[i] && a[i+3] < a[i+1] && a[i+3] < a[i+2]){
                    x = a[i+3];
                }
            }
            else if(!(a[i] < x) && !(a[i+1] < x) && !(a[i+2] < x) && !(a[i+3] < x)){

            }
            else{
                if (a[i] < x) {
                    x = a[i];
                }
                if (a[i+1] < x) {
                    x = a[i+1];
                }
                if (a[i+2] < x) {
                    x = a[i+2];
                }
                if (a[i+3] < x) {
                    x = a[i+3];
                }
            }
        }
        for (; i < LEN_1D; ++i) {
            if (a[i] < x) {
                x = a[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, x);
    }

    gettimeofday(&func_args->t2, NULL);
    return x;
}
real_t s3111_sse(struct args_t * func_args)
{

//    reductions
//    conditional sum reduction
    gettimeofday(&func_args->t1, NULL);

    real_t sum;
    int vf = 4;
    for (int nl = 0; nl < iterations/2; nl++) {
        sum = 0.;
        int i = 0;
        int upper_bound = LEN_1D/vf * vf;
        for (; i < upper_bound; i+=vf) {
            if (a[i] > (real_t)0. && a[i + 1] > (real_t)0. && a[i+2] > (real_t)0. && a[i+3] > (real_t)0.) {
                sum += a[i];
                sum += a[i+1];
                sum += a[i+2];
                sum += a[i+3];
            }
            else if(!(a[i] > (real_t)0.) && !(a[i + 1] > (real_t)0.) && !(a[i+2] > (real_t)0.) && !(a[i+3] > (real_t)0.)){

            }
            else{
                if (a[i] > (real_t)0.) {
                    sum += a[i];
                }
                if (a[i+1] > (real_t)0.) {
                    sum += a[i+1];
                }
                if (a[i+2] > (real_t)0.) {
                    sum += a[i+2];
                }
                if (a[i+3] > (real_t)0.) {
                    sum += a[i+3];
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if (a[i] > (real_t)0.) {
                sum += a[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, sum);
    }

    gettimeofday(&func_args->t2, NULL);
    return sum;
}
real_t s3113_sse(struct args_t * func_args)
{

//    reductions
//    maximum of absolute value

    gettimeofday(&func_args->t1, NULL);
    int vf =4;
    real_t max;
    for (int nl = 0; nl < iterations*4; nl++) {
        max = ABS(a[0]);
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for(; i < upper_bound; i+=vf){
            if((ABS(a[i])) > max && (ABS(a[i+1])) > max && (ABS(a[i+2])) > max && (ABS(a[i+3])) > max){
                if((ABS(a[i])) > (ABS(a[i+1])) && (ABS(a[i])) > (ABS(a[i+2])) && (ABS(a[i])) > (ABS(a[i+3]))){
                    max = (ABS(a[i]));
                }
                else if((ABS(a[i+1])) > (ABS(a[i])) && (ABS(a[i+1])) > (ABS(a[i+2])) && (ABS(a[i+1])) > (ABS(a[i+3]))){
                    max = (ABS(a[i+1]));
                }
                else if((ABS(a[i+2])) > (ABS(a[i])) && (ABS(a[i+2])) > (ABS(a[i+1])) && (ABS(a[i+2])) > (ABS(a[i+3]))){
                    max = (ABS(a[i+2]));
                }
                else if((ABS(a[i+3])) > (ABS(a[i])) && (ABS(a[i+3])) > (ABS(a[i+1])) && (ABS(a[i+3])) > (ABS(a[i+2]))){
                    max = (ABS(a[i+3]));
                }
            }
            else if(!((ABS(a[i])) > max) && !((ABS(a[i+1])) > max) && !((ABS(a[i+2])) > max) && !((ABS(a[i+3])) > max)){

            }
            else{
                if ((ABS(a[i])) > max) {
                    max = ABS(a[i]);
                }
                if ((ABS(a[i+1])) > max) {
                    max = ABS(a[i+1]);
                }
                if ((ABS(a[i+2])) > max) {
                    max = ABS(a[i+2]);
                }
                if ((ABS(a[i+3])) > max) {
                    max = ABS(a[i+3]);
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if ((ABS(a[i])) > max) {
                max = ABS(a[i]);
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, max);
    }

    gettimeofday(&func_args->t2, NULL);
    return max;
}
void s341_sse(struct args_t * func_args)
{

//    packing
//    pack positive values
//    not vectorizable, value of j in unknown at each iteration

    gettimeofday(&func_args->t1, NULL);

    int j;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        j = -1;
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for (; i < upper_bound; i+=vf) {
            if(b[i] > (real_t)0. && b[i+1] > (real_t)0. && b[i+2] > (real_t)0. && b[i+3] > (real_t)0.){
                j++;
                _mm_store_ps(&a[j], _mm_load_ps(&b[i]));
                j+=3;
            }
            else if(!(b[i] > (real_t)0.) && !(b[i+1] > (real_t)0.) && !(b[i+2] > (real_t)0.) && !(b[i+3] > (real_t)0.)){

            }
            else{
                if (b[i] > (real_t)0.) {
                    j++;
                    a[j] = b[i];
                }
                if (b[i+1] > (real_t)0.) {
                    j++;
                    a[j] = b[i+1];
                }
                if (b[i+2] > (real_t)0.) {
                    j++;
                    a[j] = b[i+2];
                }
                if (b[i+3] > (real_t)0.) {
                    j++;
                    a[j] = b[i+3];
                }
            }
        }
        for(;i<LEN_1D;i++){
            if (b[i] > (real_t)0.) {
                j++;
                a[j] = b[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s342_sse(struct args_t * func_args)
{

//    packing
//    unpacking
//    not vectorizable, value of j in unknown at each iteration

    gettimeofday(&func_args->t1, NULL);

    int j = 0;
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        j = -1;
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for(;i<upper_bound;i+=vf){
            if (a[i] > (real_t)0. && a[i+1] > (real_t)0. && a[i+2] > (real_t)0. && a[i+3] > (real_t)0.) {
                j++;
                _mm_store_ps(&a[i], _mm_load_ps(&b[j]));
                j+=3;
            }
            else if(!(a[i] > (real_t)0.) && !(a[i+1] > (real_t)0.) && !(a[i+2] > (real_t)0.) && !(a[i+3] > (real_t)0.)){

            }
            else{
                if (a[i] > (real_t)0.) {
                    j++;
                    a[i] = b[j];
                }
                if (a[i+1] > (real_t)0.) {
                    j++;
                    a[i+1] = b[j];
                }
                if (a[i+2] > (real_t)0.) {
                    j++;
                    a[i+2] = b[j];
                }
                if (a[i+3] > (real_t)0.) {
                    j++;
                    a[i+3] = b[j];
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if (a[i] > (real_t)0.) {
                j++;
                a[i] = b[j];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s343_sse(struct args_t * func_args)
{

//    packing
//    pack 2-d array into one dimension
//    not vectorizable, value of k in unknown at each iteration

    gettimeofday(&func_args->t1, NULL);

    int k;
    int vf = 4;
    for (int nl = 0; nl < 10*(iterations/LEN_2D); nl++) {
        k = -1;
        int i = 0;
        int upper_bound = LEN_2D / vf * vf;
        for (; i < upper_bound; i+=vf) {
            for (int j = 0; j < LEN_2D; j++) {
                if (bb[j][i] > (real_t)0. && bb[j][i+1] > (real_t)0. && bb[j][i+2] > (real_t)0. && bb[j][i+3] > (real_t)0.) {
                    k++;
                    _mm_store_ps(&flat_2d_array[k], _mm_load_ps(&aa[j][i]));
                    k+=3;
                }
                else if(!(bb[j][i] > (real_t)0.) && !(bb[j][i+1] > (real_t)0.) && !(bb[j][i+2] > (real_t)0.) && !(bb[j][i+3] > (real_t)0.)){

                }
                else{
                    if (bb[j][i] > (real_t)0.) {
                        k++;
                        flat_2d_array[k] = aa[j][i];
                    }
                    if (bb[j][i+1] > (real_t)0.) {
                        k++;
                        flat_2d_array[k] = aa[j][i+1];
                    }
                    if (bb[j][i+2] > (real_t)0.) {
                        k++;
                        flat_2d_array[k] = aa[j][i+2];
                    }
                    if (bb[j][i+3] > (real_t)0.) {
                        k++;
                        flat_2d_array[k] = aa[j][i+3];
                    }
                }
            }
        }
        for(; i < LEN_2D;i++){
            for (int j = 0; j < LEN_2D; j++) {
                if (bb[j][i] > (real_t)0.) {
                    k++;
                    flat_2d_array[k] = aa[j][i];
                }
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void s443_sse(struct args_t * func_args)
{
    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < 2*iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D/vf*vf;
        for (; i < upper_bound; i+=vf) {
            if(d[i] <= (real_t)0. && d[i+1] <= (real_t)0. && d[i+2] <= (real_t)0. && d[i+3] <= (real_t)0.){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                __m128 c_ = _mm_load_ps(&c[i]);
                _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(b_, c_)));
            }
            else if(!(d[i] <= (real_t)0.) && !(d[i+1] <= (real_t)0.) && !(d[i+2] <= (real_t)0.) && !(d[i+3] <= (real_t)0.)){
                __m128 a_ = _mm_load_ps(&a[i]);
                __m128 b_ = _mm_load_ps(&b[i]);
                _mm_store_ps(&a[i], _mm_add_ps(a_, _mm_mul_ps(b_, b_)));
            }
            else{
                if (d[i] <= (real_t)0.){
                    a[i] += b[i] * c[i];
                }
                else{
                    a[i] += b[i] * b[i];
                }
                if (d[i+1] <= (real_t)0.){
                    a[i+1] += b[i+1] * c[i+1];
                }
                else{
                    a[i+1] += b[i+1] * b[i+1];
                }
                if (d[i+2] <= (real_t)0.){
                    a[i+2] += b[i+2] * c[i+2];
                }
                else{
                    a[i+2] += b[i+2] * b[i+2];
                }
                if (d[i+3] <= (real_t)0.){
                    a[i+3] += b[i+3] * c[i+3];
                }
                else{
                    a[i+3] += b[i+3] * b[i+3];
                }

            }
        }
        for (; i < LEN_1D; i++) {
            if (d[i] <= (real_t)0.) {
                goto L20;
            } else {
                goto L30;
            }
            L20:
            a[i] += b[i] * c[i];
            goto L50;
            L30:
            a[i] += b[i] * b[i];
            L50:
            ;
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
void vif_sse(struct args_t * func_args)
{

//    control loops
//    vector if

    gettimeofday(&func_args->t1, NULL);
    int vf = 4;
    for (int nl = 0; nl < iterations; nl++) {
        int i = 0;
        int upper_bound = LEN_1D / vf * vf;
        for (; i < LEN_1D; i+=vf) {
            if (b[i] > (real_t)0. && b[i+1] > (real_t)0. && b[i+2] > (real_t)0. && b[i+3] > (real_t)0.) {
                _mm_store_ps(&a[i], _mm_load_ps(&b[i]));
            }
            else if(!(b[i] > (real_t)0.) && !(b[i+1] > (real_t)0.) && !(b[i+2] > (real_t)0.) && !(b[i+3] > (real_t)0.)){

            }
            else{
                if (b[i] > (real_t)0.) {
                    a[i] = b[i];
                }
                if (b[i+1] > (real_t)0.) {
                    a[i+1] = b[i+1];
                }
                if (b[i+2] > (real_t)0.) {
                    a[i+2] = b[i+2];
                }
                if (b[i+3] > (real_t)0.) {
                    a[i+3] = b[i+3];
                }
            }
        }
        for (; i < LEN_1D; i++) {
            if (b[i] > (real_t)0.) {
                a[i] = b[i];
            }
        }
        dummy(a, b, c, d, e, aa, bb, cc, 0.);
    }

    gettimeofday(&func_args->t2, NULL);
}
