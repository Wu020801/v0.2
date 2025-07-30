
#include "decoder.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "draco/io/stdio_file_reader.h"
#include "draco/io/file_reader_factory.h"
#include "draco/io/stdio_file_writer.h"
#include "draco/io/file_writer_factory.h"
bool draco::StdioFileReader::registered_in_factory_ = draco::FileReaderFactory::RegisterReader(draco::StdioFileReader::Open);
bool draco::StdioFileWriter::registered_in_factory_ = draco::FileWriterFactory::RegisterWriter(draco::StdioFileWriter::Open);

char** load_hpm_options(decoderHPMParams& params_hpm, int* argc_hpm) {

    int i = *argc_hpm;
    char** argv_hpm = new char* [20];
    if (params_hpm.input != nullptr) {
        argv_hpm[i] = "-i";
        argv_hpm[++i] = params_hpm.input;
    }else{
		printf("not enough input for hpm.");
    }
    
    if (params_hpm.output_yuv != nullptr) {
        argv_hpm[++i] = "-o";
        argv_hpm[++i] = params_hpm.output_yuv;
    }
    *argc_hpm = i + 1;
   
    return argv_hpm;
}

void YUV420ToRGB(unsigned char* yuvBuffer, unsigned char* rgbBuffer, int width, int height) {
    int frameSize = width * height;
    int chromaSize = frameSize / 4;

    unsigned char* yPlane = yuvBuffer;
    unsigned char* uPlane = yuvBuffer + frameSize;
    unsigned char* vPlane = uPlane + chromaSize;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int yIndex = y * width + x;
            int uvIndex = (y / 2) * (width / 2) + (x / 2);

            int Y = yPlane[yIndex];
            int U = uPlane[uvIndex] - 128;
            int V = vPlane[uvIndex] - 128;

            int R = Y + 1.5748 * V;
            int G = Y - 0.1873 * U - 0.4681 * V;
            int B = Y + 1.8556 * U;

            R = std::min(std::max(R, 0), 255);
            G = std::min(std::max(G, 0), 255);
            B = std::min(std::max(B, 0), 255);

            rgbBuffer[3 * yIndex + 0] = R;
            rgbBuffer[3 * yIndex + 1] = G;
            rgbBuffer[3 * yIndex + 2] = B;
        }
    }
}

bool mcemDecoder::decodeTexture(decoderHPMParams& paramsHPM) {
	char** argv_hpm = new char*[20];
	int argc_hpm = 1;
	argv_hpm = load_hpm_options(paramsHPM, &argc_hpm);

	int w = 0;
	int h = 0;
	decoder_hpm(argc_hpm, argv_hpm, &w, &h);
	std::cout << "texture map yuv saved to:" << paramsHPM.output_yuv << std::endl;

	char* yuvFilename = paramsHPM.output_yuv;
	char* pngFilename = paramsHPM.output_map;

	std::ifstream yuvFile(yuvFilename, std::ios::binary);
	if (!yuvFile.is_open()) {
		std::cerr << "Error opening yuv file." << std::endl;
		return 1;
	}

	int frameSize = w * h * 3 / 2;
	std::vector<unsigned char> yuvBuffer(frameSize);
	std::vector<unsigned char> rgbBuffer(w * h * 3);

	yuvFile.read(reinterpret_cast<char*>(yuvBuffer.data()), frameSize);
	if (!yuvFile) {
		std::cerr << "Error reading YUV file." << std::endl;
		return 1;
	}

	YUV420ToRGB(yuvBuffer.data(), rgbBuffer.data(), w, h);
	if (stbi_write_png(pngFilename, w, h, 3, rgbBuffer.data(), w * 3)) {
		std::cout << "texture map png file saved to: " << pngFilename << std::endl;
	}
	else {
		std::cerr << "Error saving PNG file." << std::endl;
	}
	return true;
}

void mcemDecoder::decodeMetadata(std::vector<char>&	bitstreamBase,
								 decoderParams&		params) {
	assert(bitstreamBase.size() >= 2);
	int count = 0;
	int32_t tempVal = 0;
	std::memcpy(&tempVal, bitstreamBase.data() + count, 2);
	params.subdivIterCount = tempVal;
	count += 2;
	bitstreamBase.erase(bitstreamBase.begin(), bitstreamBase.begin() + count);
}

void mcemDecoder::decodeBaseMesh(MeshBundle&        MB,
								TriMesh<MeshType>&  reconMesh,
								std::vector<char>&	bitstreamBase,
								decoderParams&		params,
								decoderDracoParams&	paramsDraco) {
	// decode metadata
	decodeMetadata(bitstreamBase, params);

	decodeDraco(reconMesh, bitstreamBase);
	if (params.subdivIterCount > 0) {
		subdivBaseMesh(MB, reconMesh, params.subdivIterCount);
	}		
}

void mcemDecoder::decodeGeometryResidual(MeshBundle&         MB,
										 std::vector<char>&	 bitstreamGeom,
										 decoderParams&      params) {
	if (!m_bufferChunk.readFromBitstream(bitstreamGeom)) {
		std::cout << "Fail in readFromBitstream" << std::endl;
	}
	m_decBac.setBitstreamBuffer(m_bufferChunk);
	m_decBac.initBac();
	init_aec_context_tab();

	for (uint32_t idx = 0; idx < MB.disp.size(); idx++) {
		for (uint32_t dim = 0; dim < 3; dim++) {
			int32_t res = m_decBac.parseResidual(dim, 1);
			MB.disp[idx][dim] = res;
		}
	}
	assert(m_decBac.decodeTerminationFlag());
	m_bufferChunk.reset();
	m_decBac.reset();
}

void mcemDecoder::decode(decoderParams& params, 
						 decoderDracoParams& paramsDraco, 
						 decoderHPMParams& paramsHPM, 
						 std::vector<char>& bitstreamBase,
						 std::vector<char>& bitstreamGeom) {
	TriMesh<MeshType> reconMesh;
	MeshBundle        MB;
	std::vector<uint32_t> mtlPos;
	double userTimeTotal = 0.0;
	double userTimeGeometry = 0.0;
	double userTimeTexture = 0.0;
	clock_t userTimeTotalBegin = clock();
	clock_t userTimeGeometryBegin = clock();
	decodeBaseMesh(MB, reconMesh, bitstreamBase, params, paramsDraco);
	MB.disp.clear();
	MB.disp.resize(reconMesh.XYZ.size());
	if (bitstreamGeom.size() > 0) {
		decodeGeometryResidual(MB, bitstreamGeom, params);
		reconMesh.reconGeometry(MB.disp);
	}
	userTimeGeometry = (double)(clock() - userTimeGeometryBegin) / CLOCKS_PER_SEC;

	clock_t userTimeTextureBegin = clock();
	if (params.flag_texture) {
		decodeTexture(paramsHPM);
	}
	userTimeTexture = (double)(clock() - userTimeTextureBegin) / CLOCKS_PER_SEC;
	userTimeTotal = (double)(clock() - userTimeTotalBegin) / CLOCKS_PER_SEC;

	if (params.flag_recon_dec) {
		reconMesh.writeOBJ(params.output, mtlPos, true);
	}
	std::cout << "Total decode time (user):" << userTimeTotal << " sec." << std::endl;
	std::cout << "Total geometry decode time (user):" << userTimeGeometry << " sec." << std::endl;
	std::cout << "Total texture decode time (user):" << userTimeTexture << " sec." << std::endl;
}
   
