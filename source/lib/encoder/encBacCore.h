#pragma once
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
 *    nor the names of its contributors maybe used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
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

#include "contextModel.h"
#include "comBacCore.h"
#include "comBitStream.h"

///< \in TLibEncoder \{

/**
 * Class TEncBacCore
 * entropy encoder engine
 */
static const uint32_t NUM_FLUSH_BITS = 24;

class TEncBacCore {
public:
  uint32_t range;
  uint32_t code;
  int left_bits;
  uint32_t stacked_ff;
  uint32_t pending_byte;
  uint32_t is_pending_byte;

public:
  TEncBacCore();
  ///< entropy encoder (with context)
  void enc_sbac_init();
  void aec_start(aec_t* p_aec, uint8_t* p_bs_start, uint8_t* p_bs_end, int b_writing);
  void bitstr_put_one_bit(aec_t* p_aec, uint32_t b);
  void bitstr_end_stream(aec_t* p_aec);
  void aec_done(aec_t* p_aec);
  void init_geometry_contexts(aec_t* p_aec);

  int aec_get_written_bits(aec_t* p_aec);
  int aec_get_shift(uint32_t v);
  void bitstr_flush_bits(aec_t* p_aec);
  void bitstt_put_one_bit_and_remainder(aec_t* p_aec, const int b);
  void biari_encode_symbol_aec(aec_t* p_aec, uint8_t symbol, context_t* p_ctx);
  void biari_encode_symbol_final_aec(aec_t* p_aec, uint8_t symbol);
  void biari_encode_symbol_eq_prob_aec(aec_t* p_aec, uint8_t symbol);
  void biari_encode_symbols_eq_prob_aec(aec_t* p_aec, uint32_t val, int len);
  void enc_sbac_finish(COM_BS* bsw);
  void sbac_write_unary_sym_ep(uint32_t sym, COM_BS* bs, aec_t* p_aec);

  ///< bypass encoder
  void com_bsw_init(COM_BS* bs, uint8_t* buf, uint8_t* buftmp, int size, COM_BS_FN_FLUSH fn_flush);
  void com_bsw_deinit(COM_BS* bs);  ///< write out all the remaining bits in the bitstream buffer
  int com_bsw_write1(COM_BS* bs, int val);
  int com_bsw_write(COM_BS* bs, uint32_t val, int len);
  void com_bsw_write_ue(COM_BS* bs, uint64_t val);
  void com_bsw_write_se(COM_BS* bs, int64_t val);
  void com_bsw_write_byte_align(COM_BS* bs);

  void Demulate(COM_BS* bs);  ///< prevent start code inside buffer chunk data

};  ///< END CLASS TEncBacCore

///< \{
