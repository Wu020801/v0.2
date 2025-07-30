/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once
#include "commonDef.h"

#define B_BITS 10
#define QUARTER (1 << (B_BITS - 2))

#define MAKE_CONTEXT(lg_pmps, mps, cycno) \
  (((uint16_t)(cycno) << 12) | ((uint16_t)(mps) << 0) | (uint16_t)(lg_pmps << 1))

typedef union context_t {
	struct {
		unsigned MPS : 1;       // 1  bit
		unsigned LG_PMPS : 11;  // 11 bits
		unsigned cycno : 2;     // 2  bits
	};
	uint16_t v;
} context_t;

typedef struct geom_ctx_set_t {
	context_t ctx_length_eq0;
	context_t ctx_length_prefix;
	context_t ctx_length_suffix;

	context_t ctx_geom_sign[3];
	context_t ctx_geom_eq0[3];
	context_t ctx_geom_prefix[3];
	context_t ctx_geom_suffix[3];
}geom_ctx_set_t;

