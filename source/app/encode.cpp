
#include "encoder.hpp"

int StringToInt(const std::string& s) {
    char* end;
    return strtol(s.c_str(), &end, 10);  
}

float StringToFloat(const std::string& s) {
	char* end;
	return strtof(s.c_str(), &end);
}

void Usage() {
    printf("Usage: avs_mesh_encoder\n");
    printf("  -h | -?                           show help.\n");
    printf("  -iMesh <input mesh>              *input mesh file name.\n");
    printf("  -iMap <input map>                 input texture map file name.\n");
    printf("  -oMesh <output mesh>              output mesh bitstream.\n");
    printf("  -oMap <output map>                output texture map bitstream.\n");
    printf("  -outputPath <output bitstream>    output total bitstream.\n");
	printf("  -recMesh <reconstructed mesh>     output reconstructed mesh file.\n");
    printf("  -encodeMap <flag>                 encode texture map(1) or not(0). default = 0 \n");
	printf("  -inputPosBitdepth                 input position bitdepth \n");
	printf("  -removeDuplicate                  disable in lossless mode \n");
	printf("  -targetTriangleRatio              target triangle ratio in decimation \n"); 
	printf("  -subdivIterCount                  iteration count in subdivision fitting \n");

    printf("---------------------draco parameters---------------------------\n");
    printf("  -point_cloud                      forces the input to be encoded as a point "
												"cloud.\n");
    printf("  -draco_qp <value>                 quantization bits for the position "
												 "attribute, default = 11.\n");
    printf("  -draco_qt <value>                 quantization bits for the texture coordinate "
												"attribute, default = 10.\n");
    printf("  -draco_qn <value>                 quantization bits for the normal vector "
												"attribute, default = 8.\n");
    printf("  -draco_qg <value>                 quantization bits for any generic attribute, "
												"default = 8.\n");
    printf("  -draco_cl <value>                 compression level [0-10], most = 10, least = 0, "
												"default = 7.\n");
    printf("  --skip ATTRIBUTE_NAME             skip a given attribute (NORMAL, TEX_COORD, "
												"GENERIC)\n");
    printf("  --metadata                        use metadata to encode extra information in "
												"mesh files.\n");

    printf("  -preserve_polygons                encode polygon info as an attribute.\n");
    printf("Use negative quantization values to skip the specified attribute\n");

    printf("---------------------hpm parameters---------------------------\n");
    printf("  -hpm_cfg <hpm config>             hpm encode config path.\n");
    printf("  -hpm_q <qp>                       hpm qp.\n");
    printf("  -hpm_w <width>                    texture map width.\n");
    printf("  -hpm_h <height>                   texture map height.\n");
    printf("  -hpm_yuv <yuv file>               texture map yuv path.\n");

}

bool load_options(int argc, char** argv, encoderParams& params, 
				  encoderDracoParams& params_draco, encoderHPMParams& params_hpm) {
    const int argc_check = argc - 1;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp("-h", argv[i]) || !strcmp("-?", argv[i])) {
            Usage();
            return 0;
        }
        else if (!strcmp("-iMesh", argv[i]) && i < argc_check) {
			params.input_mesh = argv[++i];
			params_draco.input = params.input_mesh;
        }
        else if (!strcmp("-iMap", argv[i]) && i < argc_check) {
            params_hpm.input = argv[++i];
        }
        else if (!strcmp("-oMesh", argv[i]) && i < argc_check) {
			params.output_mesh = argv[++i];
			params_draco.output = params.output_mesh;
        }
        else if (!strcmp("-oMap", argv[i]) && i < argc_check) {
            params_hpm.output = argv[++i];
        }
		else if (!strcmp("-recMesh", argv[i]) && i < argc_check) {
			params.recon_mesh = argv[++i];
			params.flag_recon_enc = true;
		}
		else if (!strcmp("-outputPath", argv[i]) && i < argc_check) {
			params.output_bin = argv[++i];
		}
		else if (!strcmp("-encodeMap", argv[i]) && i < argc_check) {
			int encodeMap = StringToInt(argv[++i]);
			if (encodeMap == 1) {
				params.flag_texture = true;
			}
		}
		else if (!strcmp("-inputPosBitdepth", argv[i])) {
			params.input_pos_bitdepth = StringToInt(argv[++i]);
			params_draco.input_pos_bitdepth = params.input_pos_bitdepth;
		}
		else if (!strcmp("-removeDuplicate", argv[i])) {
			params.removeDuplicateVert = StringToInt(argv[++i]);
		}
		else if (!strcmp("-targetTriangleRatio", argv[i])) {
			params.targetTriangleRatio = StringToFloat(argv[++i]);
		}
		else if (!strcmp("-subdivIterCount", argv[i])) {
			params.subdivFitSubdivIterCount = StringToInt(argv[++i]);
		}
		else if (!strcmp("-hpm_cfg", argv[i]) && i < argc_check) {
			params_hpm.config = argv[++i];
		}
		else if (!strcmp("-hpm_q", argv[i]) && i < argc_check) {
			params_hpm.q = argv[++i];
			int qp = std::stoi(params_hpm.q);
			if (qp < 0 || qp > 63) {
				printf(
					"Error: The maximum number of quantization bits for the video encoder "
					" is 63.\n");
				return false;
			}
		}
		else if (!strcmp("-hpm_w", argv[i]) && i < argc_check) {
			params_hpm.w = argv[++i];
			params.inputWidth = StringToInt(params_hpm.w);
		}
		else if (!strcmp("-hpm_h", argv[i]) && i < argc_check) {
			params_hpm.h = argv[++i];
			params.inputHeight = StringToInt(params_hpm.h);
		}
        else if (!strcmp("-hpm_yuv", argv[i]) && i < argc_check) {
            params_hpm.yuv = argv[++i];
        }
        else if (!strcmp("-point_cloud", argv[i])) {
            params_draco.is_point_cloud = true;
        }
        else if (!strcmp("-draco_qp", argv[i]) && i < argc_check) {
			params_draco.pos_quantization_bits = StringToInt(argv[++i]);
        }
        else if (!strcmp("-draco_qt", argv[i]) && i < argc_check) {
			params_draco.tex_coords_quantization_bits = StringToInt(argv[++i]);
        }
        else if (!strcmp("-draco_qn", argv[i]) && i < argc_check) {
			params_draco.normals_quantization_bits = StringToInt(argv[++i]);
        }
        else if (!strcmp("-draco_qg", argv[i]) && i < argc_check) {
			params_draco.generic_quantization_bits = StringToInt(argv[++i]);
        }
        else if (!strcmp("-draco_cl", argv[i]) && i < argc_check) {
			params_draco.compression_level = StringToInt(argv[++i]);
        }
        else if (!strcmp("--skip", argv[i]) && i < argc_check) {
            if (!strcmp("NORMAL", argv[i + 1])) {
				params_draco.normals_quantization_bits = -1;
            }
            else if (!strcmp("TEX_COORD", argv[i + 1])) {
				params_draco.tex_coords_quantization_bits = -1;
            }
            else if (!strcmp("GENERIC", argv[i + 1])) {
				params_draco.generic_quantization_bits = -1;
            }
            else {
                printf("Error: Invalid attribute name after --skip\n");
                return -1;
            }
            ++i;
        }
        else if (!strcmp("--metadata", argv[i])) {
			params_draco.use_metadata = true;
        }
        else if (!strcmp("-preserve_polygons", argv[i])) {
            params_draco.preserve_polygons = true;
        }
    }
    if (argc < 3 || params.input_mesh.empty()) {
        std::cout << "argc < 3 or options.input is empty" << std::endl;
        Usage();
        return false;
    }
	if (params_draco.output.empty()) {
		params_draco.output = removeExtension(params.input_mesh) + "_draco.bin";
	}
	if (!params.output_bin.empty()) {
		params.output_geometry = removeExtension(params.output_bin) + "_geom.bin";
	}
	if (params.flag_texture) {
		if (params_hpm.input == nullptr || params_hpm.config == nullptr) {
			std::cout << "Error: -encodeMap = 1,but has no input png or hpm config!" << std::endl;
			return -1;
		}
		if (params_hpm.yuv == nullptr) {
			size_t out_length = removeExtension(params_draco.input).size() + strlen(".yuv") + 1;
			char* out = (char*)malloc(out_length);
			if (out != nullptr) {
				strcpy(out, removeExtension(params_draco.input).c_str());
				strcat(out, ".yuv");
				params_hpm.yuv = out;
			}
		}
		if (params_hpm.output == nullptr) {
			size_t out_length = removeExtension(params_draco.input).size() + strlen("_hpm.bin") + 1;
			char* out = (char*)malloc(out_length);
			if (out != nullptr) {
				strcpy(out, removeExtension(params_draco.input).c_str());
				params_hpm.output = strcat(out, "_hpm.bin");
			}
		}
	}
	if (params.output_bin.empty() && !params.input_mesh.empty()) {
		size_t out_length = removeExtension(params.input_mesh).size() + strlen(".bin") + 1;
		char* out = (char*)malloc(out_length);
		if (out != nullptr) {
			strcpy(out, removeExtension(params.input_mesh).c_str());
			params.output_bin = strcat(out, ".bin");
		}
	}
	// print
	std::cout << "input mesh: " << params.input_mesh << std::endl;
	std::cout << "input map : " << params_hpm.input << std::endl;
	// std::cout << "resQuantBit: " << params.resQuantBit << std::endl;
	std::cout << "input pos bitdepth: " << params.input_pos_bitdepth << std::endl;
	std::cout << "targetTriangleRatio: " << params.targetTriangleRatio << std::endl;
	std::cout << "subdivFitSubdivIterCount: " << params.subdivFitSubdivIterCount << std::endl;
	std::cout << "draco_qp: " << params_draco.pos_quantization_bits << std::endl;
	std::cout << "draco_qt: " << params_draco.tex_coords_quantization_bits << std::endl;
	std::cout << "hpm_q: " << StringToInt(params_hpm.q) << std::endl;
	std::cout << "hpm_w: " << params.inputWidth << std::endl;
	std::cout << "hpm_h: " << params.inputHeight << std::endl;
	return true;
}

void encodeBitstreams(bool flag_texture, const std::string& bitstreamDracoPath, const std::string& bitstreamGeomPath, 
					  const std::string& bitstreamTexPath, const std::string& outputBinPath) {
	ofstream outputFile;
	outputFile.open(outputBinPath, fstream::binary | fstream::out);

	std::ifstream bitstream1(bitstreamDracoPath, std::ios::binary);
    std::vector<char> buffer1((std::istreambuf_iterator<char>(bitstream1)), std::istreambuf_iterator<char>());
    bitstream1.close();
    uint32_t size1 = buffer1.size();
    outputFile.write(reinterpret_cast<const char*>(&size1), sizeof(size1));
    outputFile.write(buffer1.data(), buffer1.size());

	std::ifstream bitstream2(bitstreamGeomPath, std::ios::binary);
	std::vector<char> buffer2((std::istreambuf_iterator<char>(bitstream2)), std::istreambuf_iterator<char>());
	bitstream2.close();
	uint32_t size2 = buffer2.size();
	outputFile.write(reinterpret_cast<const char*>(&size2), sizeof(size2));
	outputFile.write(buffer2.data(), buffer2.size());
	std::cout << "Bitstream size basemesh:" << size1 << " bytes." << endl;
	std::cout << "Bitstream size geometry:" << size2 << " bytes." << endl;

	if (flag_texture) {
		std::ifstream bitstream3(bitstreamTexPath, std::ios::binary);
		std::vector<char> buffer3((std::istreambuf_iterator<char>(bitstream3)), std::istreambuf_iterator<char>());
		bitstream3.close();
		uint32_t size3 = buffer3.size();
		outputFile.write(buffer3.data(), buffer3.size());
		std::cout << "Bitstream size texture:" << size3 << " bytes." << endl;
	}
	outputFile.close();
}

int main(int argc, char** argv) {
	std::cout << "avs mesh encoder version 0.2" << std::endl;
	
	encoderParams params;
	encoderDracoParams paramsDraco;
	encoderHPMParams paramsHPM;
    char* outputPath = {};

    if (!load_options(argc, argv, params, paramsDraco, paramsHPM)) {
        return -1;
    }

	mcemEncoder mcemEncoder;
	mcemEncoder.encode(params, paramsDraco, paramsHPM);

	if (!params.output_bin.empty()) {
		encodeBitstreams(params.flag_texture, paramsDraco.output, params.output_geometry,
						 paramsHPM.output, params.output_bin);
		std::cout << "Bitstream saved to:" << params.output_bin << std::endl;
	}
	std::cout << "Encode finished! " << std::endl;
	return 0;
}