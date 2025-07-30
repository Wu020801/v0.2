/* The copyright in this software is being made available under the BSD
* License, included below. This software may be subject to other third party
* and contributor rights, including patent rights, and no such rights are
* granted under this license.
*
* Copyright (c) 2019-2033, Audio Video coding Standard Workgroup of China
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
*  * Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*  * Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*  * Neither the name of Audio Video coding Standard Workgroup of China
*    nor the names of its contributors maybe used to endorse or promote products
*    derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
* THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "decBacTop.h"
#include "decBacCore.h"
#include "comRom.h"

///< \in TLibDecoder \{

/**
 * Class TDecBacTop
 * entropy decoder
 */

//////////////////////////////////////////////////////////////////////////
// Public class functions
//////////////////////////////////////////////////////////////////////////

int TDecBacTop::parseRunlength() {
  int val = 0;
  if (m_bac->biari_decode_symbol(p_aec, &p_aec->geometry_syn_ctx.ctx_length_eq0))
    return 0;
  else {
    val = parseExpGolomb(2, &p_aec->geometry_syn_ctx.ctx_length_prefix,
						 &p_aec->geometry_syn_ctx.ctx_length_suffix);
    return (1 + val);
  }
  return val + 1;
}

int64_t TDecBacTop::parseResidual(const int ctx_id, const unsigned int& golombNum) {

  int64_t val = 0;
  if (m_bac->biari_decode_symbol(p_aec,
	  &p_aec->geometry_syn_ctx.ctx_geom_eq0[ctx_id])) {
	  return 0;
  }
  int sign_bit = m_bac->biari_decode_symbol(p_aec,
 				 &p_aec->geometry_syn_ctx.ctx_geom_sign[ctx_id]);
  val = parseExpGolomb(golombNum, &p_aec->geometry_syn_ctx.ctx_geom_prefix[ctx_id],
					   &p_aec->geometry_syn_ctx.ctx_geom_suffix[ctx_id]) + 1;

  val = (sign_bit == 1) ? val : -val;
  return val;
}

int TDecBacTop::parseExpGolomb(int k, context_t* p_ctxPrefix, context_t* p_ctxSufffix) {
	unsigned int l;
	int symbol = 0;
	int binary_symbol = 0;
	do {
		l = m_bac->biari_decode_symbol(p_aec, p_ctxPrefix);
		if (l == 1) {
			symbol += (1 << k);
			k++;
		}
	} while (l != 0);
	while (k--)  //next binary part
		if (m_bac->biari_decode_symbol(p_aec, p_ctxSufffix) == 1) {
			binary_symbol |= (1 << k);
		}
	return static_cast<unsigned int>(symbol + binary_symbol);
}

bool TDecBacTop::decodeTerminationFlag() {
  return !!m_bac->biari_decode_final(p_aec);
}

TDecBacTop::TDecBacTop() {
  m_bac = nullptr;
  reset();
}

TDecBacTop::~TDecBacTop() {
  if (m_bac) {
    delete m_bac;
    m_bac = nullptr;
  }
}

void TDecBacTop::initBac() {
  m_bac = new TDecBacCore;
  m_bac->dec_sbac_init(&m_bitStream);
}

void TDecBacTop::reset() {
  if (m_bac) {
    delete m_bac;
    m_bac = nullptr;
  }
  memset(&m_bitStream, 0, sizeof(m_bitStream));
  aec.init();
  p_aec = &aec;
  std::fill(begin(memoryChannel), end(memoryChannel), 15);
}

void TDecBacTop::setBitstreamBuffer(TComBufferChunk& buffer) {
  m_bac->com_bsr_init(&m_bitStream, (uint8_t*)buffer.addr, buffer.ssize, NULL);
  m_bac->init_geometry_contexts(p_aec);
  m_bac->aec_start_decoding(p_aec, m_bitStream.beg, 0,
	  buffer.ssize);
}

///< \{
