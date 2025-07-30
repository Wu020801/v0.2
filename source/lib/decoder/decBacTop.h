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

#include "decBacCore.h"
#include "commonDef.h"
#include "contextModel.h"
// #include "common/TComBufferChunk.h"

#define CHECK_ALL_CTX 0

///< \in TLibDecoder \{

/**
 * Class TDecBacTop
 * entropy decoder
 */

class TDecBacTop {
private:
  COM_BS m_bitStream;
  TDecBacCore* m_bac;
  aec_t aec;
  aec_t* p_aec;
  geom_ctx_set_t m_ctxBackup;  // use to store the context for later restoration
  uint8_t memoryChannel[1024];
  uint8_t temp[1024];

public:
  ///< decoding syntax

  ///< attribute related syntax
  int parseRunlength();
  int parseExpGolomb(int k, context_t* p_ctxPrefix, context_t* p_ctxSufffix);
  int64_t parseResidual(const int ctx_id, const unsigned int& golombNum);

  TDecBacTop();
  ~TDecBacTop();
  void reset();
  void initBac();
  bool decodeTerminationFlag();
  void setBitstreamBuffer(TComBufferChunk& buffer);

  void saveContext() {
    m_ctxBackup = p_aec->geometry_syn_ctx;
    for (int i = 0; i < 1024; i++)
      temp[i] = memoryChannel[i];
  }

  void restoreContext() {
    p_aec->geometry_syn_ctx = m_ctxBackup;
    for (int i = 0; i < 1024; i++)
      memoryChannel[i] = temp[i];
  }

};  ///< END CLASS TDecBacTop

///< \{
