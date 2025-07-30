#pragma once

#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "mesh.hpp"
#include "image.hpp"
#include <UVAtlas.h>

#include "encBacTop.h"
#include "comBitStream.h"

#define DEBUG_ENC_LOG 1

//============================================================================

struct encoderDracoParams {
	int input_pos_bitdepth = 16;
	bool is_point_cloud = false;
	int pos_quantization_bits = 16;
	int tex_coords_quantization_bits = 12;
	bool tex_coords_deleted = false;
	int normals_quantization_bits = 8;
	bool normals_deleted = false;
	int generic_quantization_bits = 8;
	bool generic_deleted = false;
	int compression_level = 7;
	bool preserve_polygons = false;
	bool use_metadata = false;
	std::string input;
	std::string output;
};

struct encoderParams {
	// general
	bool flag_texture = false;
	bool flag_recon_enc = false;
	std::string input_mesh;
	std::string output_bin;
	std::string output_mesh;
	std::string output_geometry;
	std::string recon_mesh;
	int input_pos_bitdepth = 16;
	uint32_t inputWidth = 2048;
	uint32_t inputHeight = 2048;

	MeshType scaleUV;
	Vector2D<MeshType> bboxMinUV = { 0.0, 0.0 };
	Vector2D<MeshType> bboxMaxUV = { 1.0, 1.0 };

	// decimation
	double		 targetTriangleRatio = 0.01;

	// reparametrization
	DirectX::UVATLAS uvOptions = DirectX::UVATLAS_DEFAULT;
	size_t           maxCharts = size_t();
	float            maxStretch = 0.16667F;
	float            gutter = 2.F;

	// subdivfit
	bool         removeDuplicateVert = true;
	uint32_t	 subdivFitSubdivIterCount = 0;

	// texture transfer
	int32_t       textureTransferPaddingBoundaryIterationCount = 2;
	int32_t       textureTransferPaddingDilateIterationCount = 2;
	double        textureTransferPaddingSparseLinearThreshold = 0.05;
};

struct encoderHPMParams {
	char* q = nullptr;
	char* w = nullptr;
	char* h = nullptr;
	char* input = nullptr;
	char* output = nullptr;
	char* config = nullptr;
	char* yuv = nullptr;
};

class mcemEncoder {
private:
	TEncBacTop m_encBac;            ///< pointer to bac
	TComBufferChunk m_bufferChunk;  ///< bitstream buffer chunk

public:
	void encode(encoderParams& params, 
				encoderDracoParams&	   paramsDraco, 
				encoderHPMParams&      paramsHPM);
	void cleanMesh(TriMesh<MeshType>& inputMesh, encoderParams const& params);

	bool simplifyInput(TriMesh<MeshType>&     input,
					   TriMesh<MeshType>&     decimate,
					   TriMesh<MeshType>&     mapped,
					   TriMesh<MeshType>&     reference,
					   encoderParams& params);
	void computeXYZmapping(TriMesh<MeshType>           inputMesh,
						   TriMesh<MeshType>           reconMesh,
						   std::vector<uint32_t>&      reconXYZ2inputXYZ,
						   uint32_t					   subdivFitSubdivIterCount,
						   bool                        isSubdiv = true) const;

	void encodeMetadata(encoderParams&			  params,
						ofstream&                 bitstreamFile);
	void encodeBaseMesh(MeshBundle&               MB,
						TriMesh<MeshType>&        reconMesh,
					    encoderParams&			  params,
						encoderDracoParams &	  paramsDraco,
						bool                      flag_lossless_geometry);
	bool encodeDracoFile(encoderDracoParams&      params, 
						 TriMesh<MeshType>&       reconMesh);
	bool encodeDraco(TriMesh<MeshType>&           baseMesh, 
					 TriMesh<MeshType>&           reconMesh,
					 std::vector<char>&           bitstream,
					 encoderDracoParams&		  params);

	void computeGeometryResidual(MeshBundle&             MB,
								TriMesh<MeshType> const& reconMesh);
	void encodeGeometryResidual(MeshBundle&    MB,
								encoderParams& params);

	bool encodeTexture(Image3<uint8_t>& rgb_data, encoderParams& params, 
						encoderHPMParams& paramsHPM);
};

extern "C" {
	int encoder_hpm(int argc, char** argv);
}
