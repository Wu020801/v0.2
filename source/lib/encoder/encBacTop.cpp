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

#include "encBacTop.h"
#include "comRom.h"
#include <algorithm>

///< \in TLibEncoder \{

/**
 * Class TEncBacTop
 * entropy encoder
 */

//////////////////////////////////////////////////////////////////////////
// Public class functions
//////////////////////////////////////////////////////////////////////////

void TEncBacTop::encodeRunlength(int32_t& length) {
  const bool isZero = length == 0;
  m_bac->biari_encode_symbol_aec(p_aec, isZero, &p_aec->p_geometry_ctx_set->ctx_length_eq0);
  if (isZero)
    return;
  encodeExpGolomb(length - 1, 2, &p_aec->p_geometry_ctx_set->ctx_length_prefix,
	  &p_aec->p_geometry_ctx_set->ctx_length_suffix);

  m_bitStream.cur = p_aec->p;
}

void TEncBacTop::encodeResidual(int32_t& delta, const uint32_t ctx_id) {
	const bool isZero = delta == 0;
	m_bac->biari_encode_symbol_aec(p_aec, isZero, &p_aec->p_geometry_ctx_set->ctx_geom_eq0[ctx_id]);
	if (isZero)
		return;
	int sign_bit = (delta > 0) ? 1 : 0;
	m_bac->biari_encode_symbol_aec(p_aec, sign_bit, &p_aec->p_geometry_ctx_set->ctx_geom_sign[ctx_id]);
	int32_t abs_delta = (delta > 0) ? delta : (-delta);
	encodeExpGolomb(abs_delta - 1, 1, &p_aec->p_geometry_ctx_set->ctx_geom_prefix[ctx_id],
					&p_aec->p_geometry_ctx_set->ctx_geom_suffix[ctx_id]);
	m_bitStream.cur = p_aec->p;
}

void TEncBacTop::encodeExpGolomb(unsigned int symbol, int k, context_t* p_ctxPrefix,
								 context_t* p_ctxSufffix) {
	while (1) {
		if (symbol >= (1u << k)) {
			m_bac->biari_encode_symbol_aec(p_aec, 1, p_ctxPrefix);
			symbol -= (1u << k);
			k++;
		}
		else {
			m_bac->biari_encode_symbol_aec(p_aec, 0, p_ctxPrefix);
			while (k--)
				m_bac->biari_encode_symbol_aec(p_aec, (symbol >> k) & 1, p_ctxSufffix);
			break;
		}
	}
}

void TEncBacTop::encodeTerminationFlag(bool terminateFlag) {
	m_bac->biari_encode_symbol_final_aec(p_aec, terminateFlag ? 1 : 0);
	m_bac->aec_done(p_aec);

	m_bitStream.cur = p_aec->p;
}

void TEncBacTop::encodeFinish() {
  if (m_bac)
    m_bac->enc_sbac_finish(&m_bitStream);

  m_bac->Demulate(&m_bitStream);

  m_bac->com_bsw_deinit(&m_bitStream);  ///< flush the buffer
}

TEncBacTop::TEncBacTop() {
  m_bac = nullptr;
  reset();
}

TEncBacTop::~TEncBacTop() {
  if (m_bac) {
    delete m_bac;
    m_bac = nullptr;
  }
}

void TEncBacTop::setBitstreamBuffer(TComBufferChunk& buffer) {
	buffer.allocateBuffSize(buffersize);
	m_bac->com_bsw_init(&m_bitStream, (uint8_t*)buffer.addr, (uint8_t*)buffer.addr2, buffer.bsize, NULL);
	m_bac->init_geometry_contexts(p_aec);
	m_bac->aec_start(p_aec, m_bitStream.beg, m_bitStream.end, 1);
}

void TEncBacTop::initBac() {
  m_bac = new TEncBacCore;
  m_bac->enc_sbac_init();
}

void TEncBacTop::reset() {
  if (m_bac) {
    delete m_bac;
    m_bac = nullptr;
  }
  m_bitStream.init();
  aec.init();
  p_aec = &aec;
}

uint64_t TEncBacTop::getBitStreamLength() {
  return m_bitStream.cur - m_bitStream.beg;
}

///< \{
