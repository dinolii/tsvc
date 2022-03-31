

#ifndef TSVC_2_MASTER_TSVC_AVX_H
#define TSVC_2_MASTER_TSVC_AVX_H
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
void s123_avx(struct args_t * func_args);
void s124_avx(struct args_t * func_args);
void s161_avx(struct args_t * func_args);
void s1161_avx(struct args_t * func_args);
void s253_avx(struct args_t * func_args);
void s258_avx(struct args_t * func_args);
void s271_avx(struct args_t * func_args);
void s273_avx(struct args_t * func_args);
void s274_avx(struct args_t * func_args);
void s277_avx(struct args_t * func_args);
void s278_avx(struct args_t * func_args);
void s2711_avx(struct args_t * func_args);
void s2712_avx(struct args_t * func_args);
real_t s314_avx(struct args_t * func_args);
real_t s315_avx(struct args_t * func_args);
real_t s316_avx(struct args_t * func_args);
real_t s3111_avx(struct args_t * func_args);
real_t s3113_avx(struct args_t * func_args);
void s341_avx(struct args_t * func_args);
void s342_avx(struct args_t * func_args);
void s343_avx(struct args_t * func_args);
void s443_avx(struct args_t * func_args);
void vif_avx(struct args_t * func_args);
#endif //TSVC_2_MASTER_TSVC_AVX_H
