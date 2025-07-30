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

#include "comBitStream.h"

///< \in TLibCommon

/**
 * Implementation of TComBitstream
 * bitstream buffer operations
 */

//////////////////////////////////////////////////////////////////////////
// Public class functions
//////////////////////////////////////////////////////////////////////////

TComBitstream::TComBitstream() {
  init();
}

TComBitstream::~TComBitstream() {
  reset();
}

void TComBitstream::reset() {
  if (addr)
    free((unsigned char*)addr);
  if (addr2)
    free((unsigned char*)addr2);
  init();
}

void TComBitstream::allocateBuffSize(const size_t buffSize) {
  ///< allocate bitstream buffer
  if (addr == nullptr) {
    addr = (unsigned char*)calloc(1, buffSize);
    assert(checkCond(addr != NULL, "Error: cannot allocate bitstream buffer!"));
  }
  if (addr2 == nullptr) {
    addr2 = (unsigned char*)calloc(1, buffSize);
    assert(checkCond(addr2 != NULL, "Error: cannot allocate bitstream buffer!"));
  }
  ssize = (int)buffSize;
  bsize = (int)buffSize;
  err = 0;
}

void* TComBitstream::getBitStreamBuffer() {
  return addr;
}

void TComBufferChunk::writeToBitstream(ofstream* outBitstream, uint64_t length) {
	outBitstream->write((char*)getBitStreamBuffer(), length);
}

int TComBufferChunk::readFromBitstream(ifstream& inBitstream, int buffersize) {
	return 1;
}

bool TComBufferChunk::readFromBitstream(std::vector<char>& inBuffer) {
	size_t buffersize = (1 << 28);
	allocateBuffSize(buffersize);
	char* pucBuffer = (char*)getBitStreamBuffer();
	if (pucBuffer != nullptr) {
		std::memcpy(pucBuffer, inBuffer.data(), inBuffer.size());
		buffersize = initParsingConvertPayloadToRBSP(inBuffer.size(), (unsigned char*)addr, (unsigned char*)addr2);;
		allocateBuffSize(buffersize);
		return true;
	}
	else {
		return false;
	}
}

size_t TComBufferChunk::initParsingConvertPayloadToRBSP(const size_t uiBytesRead, unsigned char* pBuffer,
	unsigned char* pBuffer2) {
	unsigned int uiZeroCount = 0;
	unsigned int uiBytesReadOffset = 0;
	unsigned int uiBitsReadOffset = 0;
	const unsigned char* pucRead = pBuffer;
	unsigned char* pucWrite = pBuffer2;
	unsigned int uiWriteOffset = uiBytesReadOffset;
	unsigned char ucCurByte = pucRead[uiBytesReadOffset];

	for (uiBytesReadOffset = 0; uiBytesReadOffset < uiBytesRead; uiBytesReadOffset++) {
		ucCurByte = pucRead[uiBytesReadOffset];
		if (2 <= uiZeroCount && 0x02 == pucRead[uiBytesReadOffset]) {
			pucWrite[uiWriteOffset] = ((pucRead[uiBytesReadOffset] >> 2) << (uiBitsReadOffset + 2));
			uiBitsReadOffset += 2;
			uiZeroCount = 0;
			if (uiBitsReadOffset >= 8) {
				uiBitsReadOffset = 0;
				continue;
			}
			if (uiBytesReadOffset >= uiBytesRead) {
				break;
			}
		}
		else if (2 <= uiZeroCount && 0x01 == pucRead[uiBytesReadOffset]) {
			uiBitsReadOffset = 0;
			pucWrite[uiWriteOffset] = pucRead[uiBytesReadOffset];
		}
		else {
			pucWrite[uiWriteOffset] = (pucRead[uiBytesReadOffset] << uiBitsReadOffset);
		}

		if (uiBytesReadOffset + 1 < uiBytesRead) {
			pucWrite[uiWriteOffset] |= (pucRead[uiBytesReadOffset + 1] >> (8 - uiBitsReadOffset));
		}
		uiWriteOffset++;

		if (0x00 == ucCurByte) {
			uiZeroCount++;
		}
		else {
			uiZeroCount = 0;
		}
	}

	// th just clear the remaining bits in the buffer
	for (unsigned int ui = uiWriteOffset; ui < uiBytesRead; ui++) {
		pucWrite[ui] = 0;
	}
	memcpy(pBuffer, pBuffer2, uiWriteOffset);
	return uiBytesRead;
}
//////////////////////////////////////////////////////////////////////////
// Private class functions
//////////////////////////////////////////////////////////////////////////

void TComBitstream::init() {
  addr = addr2 = pddr = nullptr;
  bsize = ssize = err = 0;
  memset(ndata, 0, 4 * sizeof(int));
  memset(pdata, 0, 4 * sizeof(void*));
  memset(ts, 0, sizeof(long long));
}

///< \}
