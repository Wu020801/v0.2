
#include "decoder.hpp"

void Usage() {
    printf("Usage: avs_mesh_decoder \n");
    printf("Main options:\n");
    printf("  -h | -?                   show help.\n");
    printf("  -i <input>               *input bitstream path.\n");
    printf("  -oMesh <output>           output mesh file path.\n");
    printf("  -oYuv <output>            output texture map yuv file path.\n");
    printf("  -oMap <output>            output texture map png file path.\n");
}

char* removeExtension(char* filename) {
    char* file = filename;
    std::string strFilename(file);
    size_t dotPos = strFilename.rfind('.');
    std::string result = strFilename.substr(0, dotPos);

    char* newFilename = new char[result.size() + 1];
    std::strcpy(newFilename, result.c_str());
    return newFilename;
}

bool decodeBitstreams(const std::string& inputPath, std::vector<char>& bufferBase, 
					  std::vector<char>& bufferGeom, char* & inputMapBin) {
    std::ifstream inputFile(inputPath, std::ios::binary);
    uint32_t size1 = 0;
    inputFile.read(reinterpret_cast<char*>(&size1), sizeof(size1));
	bufferBase.resize(size1);
    inputFile.read(bufferBase.data(), size1);

	uint32_t size2 = 0;
	inputFile.read(reinterpret_cast<char*>(&size2), sizeof(size2));
	bufferGeom.resize(size2);
	inputFile.read(bufferGeom.data(), size2);

    //if (inputMeshBin.empty()) {
    //    inputMeshBin = removeExtension(inputPath);
    //    inputMeshBin += "_mesh.bin";
    //}
    //std::ofstream outputBitstream1(inputMeshBin, std::ios::binary);
    //outputBitstream1.write(buffer1.data(), buffer1.size());
    //outputBitstream1.close();
	std::cout << "Bitstream size basemesh:" << size1 << " bytes." << endl;
	std::cout << "Bitstream size geometry:" << size2 << " bytes." << endl;

	std::vector<char> bufferTex((std::istreambuf_iterator<char>(inputFile)), std::istreambuf_iterator<char>());
	inputFile.close();
    if (bufferTex.empty()) {
        return false;
    }
    else {
        if (inputMapBin == nullptr) {
            std::string input = removeExtension(inputPath);
            size_t totalSize = input.length() + std::strlen("_map.bin") + 1;
            inputMapBin = new char[totalSize];
            std::strcpy(inputMapBin, input.c_str());
            std::strcat(inputMapBin, "_map.bin");
        }
        std::ofstream outputBitstream2(inputMapBin, std::ios::binary);
        outputBitstream2.write(bufferTex.data(), bufferTex.size());
        outputBitstream2.close();
		std::cout << "Bitstream size texture:" << bufferTex.size() << " bytes." << endl;
    }
	return true;
}

int main(int argc, char** argv) {
	std::cout << "avs mesh decoder version 0.2" << std::endl;
	decoderParams params;
	decoderDracoParams params_draco;
	decoderHPMParams params_hpm;
    const int argc_check = argc - 1;
    std::string inputBin;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp("-h", argv[i]) || !strcmp("-?", argv[i])) {
            Usage();
            return 0;
        }
        else if (!strcmp("-i", argv[i]) && i < argc_check) {
            inputBin = argv[++i];
        }
        else if (!strcmp("-oMap", argv[i]) && i < argc_check) {
            params_hpm.output_map = argv[++i];
        }
        else if (!strcmp("-oYuv", argv[i]) && i < argc_check) {
            params_hpm.output_yuv = argv[++i];
        }
        else if (!strcmp("-oMesh", argv[i]) && i < argc_check) {
            params.output = argv[++i];
			params_draco.output = params.output;
			params.flag_recon_dec = true;
        }
    }
    if (inputBin.empty()) {
        Usage();
        return -1;
    }

	std::vector<char> bitstreamBase;
	std::vector<char> bitstreamGeom;
	params.flag_texture = decodeBitstreams(inputBin, bitstreamBase, bitstreamGeom, params_hpm.input);

	if (params.flag_texture) {
		if (params_hpm.input == nullptr) {
			std::cout << "Error: has no input map bitstream!" << std::endl;
			return -1;
		}
		else {
			if (params_hpm.output_map == nullptr) {
				char* mappath = {};
				mappath = removeExtension(params_hpm.input);
				strcat(mappath, ".png");
				params_hpm.output_map = mappath;
			}
			if (params_hpm.output_yuv == nullptr) {
				char* yuvpath = {};
				yuvpath = removeExtension(params_hpm.input);
				strcat(yuvpath, ".yuv");
				params_hpm.output_yuv = yuvpath;
			}
		}
	}

	mcemDecoder mcemDecoder;
	mcemDecoder.decode(params, params_draco, params_hpm, bitstreamBase, bitstreamGeom);
	std::cout << "Decode finished! " << std::endl;
}
   
