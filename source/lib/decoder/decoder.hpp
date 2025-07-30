#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cinttypes>
#include <algorithm>
#include "mesh.hpp"
#include "decBacTop.h"

struct decoderParams {
	bool flag_texture = false;
	bool flag_recon_dec = false;
	std::string input = {};
	std::string output = {};

	// subdivision
	int subdivIterCount = 0;
};

struct decoderDracoParams {
	std::string input = {};
	std::string output = {};

	int pos_quantization_bits = 11;
	int tex_coords_quantization_bits = 10;
};

struct decoderHPMParams {
	char* input = nullptr;
	char* output_yuv = nullptr;
	char* output_map = nullptr;
};

class mcemDecoder {
private:
	TDecBacTop m_decBac;            ///< pointer to bac
	TComBufferChunk m_bufferChunk;  ///< bitstream buffer chunk

public:
	void decode(decoderParams& params,
				decoderDracoParams&	   paramsDraco, 
				decoderHPMParams&      paramsHPM,
				std::vector<char>& bitstreamBase,
				std::vector<char>& bitstreamGeom);
	void decodeMetadata(std::vector<char>&	bitstreamBase, 
						decoderParams&		params);
	void decodeBaseMesh(MeshBundle&         MB,
						TriMesh<MeshType>&  reconMesh,
						std::vector<char>&	bitstreamBase,	
						decoderParams&		params,
						decoderDracoParams& paramsDraco);
	int decodeDracoFile(decoderDracoParams&	paramsDraco);
	bool decodeDraco(TriMesh<MeshType>&      reconMesh,
					 std::vector<char>&      bitstream);

	void decodeGeometryResidual(MeshBundle&        MB, 
								std::vector<char>& bitstreamGeom,
								decoderParams&	   params);
	bool decodeTexture(decoderHPMParams& paramsHPM);
};

extern "C" {
	int decoder_hpm(int argc, char** argv, int *w, int *h);
}
