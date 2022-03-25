

#ifndef TSVC_2_MASTER_TSVC_SSE_H
#define TSVC_2_MASTER_TSVC_SSE_H
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
void s123_sse(struct args_t * func_args);
void s124_sse(struct args_t * func_args);
void s161_sse(struct args_t * func_args);
void s1161_sse(struct args_t * func_args);
void s253_sse(struct args_t * func_args);
void s258_sse(struct args_t * func_args);
void s271_sse(struct args_t * func_args);
void s273_sse(struct args_t * func_args);
void s274_sse(struct args_t * func_args);
void s277_sse(struct args_t * func_args);
void s278_sse(struct args_t * func_args);
void s2711_sse(struct args_t * func_args);
void s2712_sse(struct args_t * func_args);
real_t s314_sse(struct args_t * func_args);
real_t s315_sse(struct args_t * func_args);
real_t s316_sse(struct args_t * func_args);
real_t s3111_sse(struct args_t * func_args);
real_t s3113_sse(struct args_t * func_args);
void s341_sse(struct args_t * func_args);
void s342_sse(struct args_t * func_args);
void s343_sse(struct args_t * func_args);
void s443_sse(struct args_t * func_args);
void vif_sse(struct args_t * func_args);
#endif //TSVC_2_MASTER_TSVC_SSE_H
