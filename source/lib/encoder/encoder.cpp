
#include "encoder.hpp"
#include "misc.hpp"
#include "geometryDecimate.hpp"
#include "geometryParametrization.hpp"
#include "transferColor.hpp"
#include "textureParametrization.hpp"
#include "contextModel.h"

#include "draco/io/stdio_file_reader.h"
#include "draco/io/file_reader_factory.h"
#include "draco/io/stdio_file_writer.h"
#include "draco/io/file_writer_factory.h"

bool draco::StdioFileReader::registered_in_factory_ = draco::FileReaderFactory::RegisterReader(draco::StdioFileReader::Open);
bool draco::StdioFileWriter::registered_in_factory_ = draco::FileWriterFactory::RegisterWriter(draco::StdioFileWriter::Open);

char** load_hpm_options(encoderHPMParams& params_hpm, int* argc_hpm) {
	std::vector<std::string> args;
	if (params_hpm.config) {
		args.push_back("--config");
		args.push_back(params_hpm.config);
	}
	else {
		std::cerr << "have no cfg!" << std::endl;
		exit(-1);
	}

	if (params_hpm.yuv) {
		args.push_back("-i");
		args.push_back(params_hpm.yuv);
	}
	else {
		std::cerr << "have no yuv path!" << std::endl;
		exit(-1);
	}

	if (params_hpm.output) {
		args.push_back("-o");
		args.push_back(params_hpm.output);
	}
	else {
		std::cerr << "have no output path!" << std::endl;
		exit(-1);
	}

	if (params_hpm.q) {
		args.push_back("-q");
		args.push_back(params_hpm.q);
	}
	else {
		std::cerr << "hpm qp error!" << std::endl;
		exit(-1);
	}

	if (params_hpm.w && params_hpm.h) {
		args.push_back("-h");
		args.push_back(params_hpm.h);
		args.push_back("-w");
		args.push_back(params_hpm.w);
	}

	*argc_hpm = args.size() + 1;
	char** argv_hpm = new char*[*argc_hpm];
	for (int i = 0; i < *argc_hpm - 1; ++i) {
		argv_hpm[i + 1] = new char[args[i].size() + 1];
		std::strcpy(argv_hpm[i + 1], args[i].c_str());
		//std::cerr << "argv_hpm[" << i << "] = " << argv_hpm[i+1] << std::endl;
	}
	return argv_hpm;
}

bool mcemEncoder::encodeTexture(Image3<uint8_t>& rgb_data, encoderParams& params, 
								encoderHPMParams& paramsHPM) {

	std::cout << "transfer rgb to yuv420" << std::endl;
	uint32_t width = params.inputWidth;
	uint32_t height = params.inputHeight;
	std::vector<unsigned char> y_plane(width * height);
	std::vector<unsigned char> u_plane(width * height / 4);
	std::vector<unsigned char> v_plane(width * height / 4);
	rgb_data.rgb2yuv420(y_plane, u_plane, v_plane);

	std::ofstream file(paramsHPM.yuv, std::ios::binary);
	file.write(reinterpret_cast<const char*>(y_plane.data()), y_plane.size());
	file.write(reinterpret_cast<const char*>(u_plane.data()), u_plane.size());
	file.write(reinterpret_cast<const char*>(v_plane.data()), v_plane.size());
	file.close();

	char** argv_hpm = new char*[20];
	int argc_hpm = 1;
	argv_hpm = load_hpm_options(paramsHPM, &argc_hpm);
	encoder_hpm(argc_hpm, argv_hpm);

	return true;
}

bool mcemEncoder::simplifyInput(TriMesh<MeshType>&     input,
								TriMesh<MeshType>&     decimate,
								TriMesh<MeshType>&     mapped,
								TriMesh<MeshType>&     reference,
								encoderParams& params) {

	reference = input;
	GeometrySimplify<MeshType> geometrySimplify;

	int faceCount = std::ceil(input.FaceXYZ.size() * params.targetTriangleRatio);
	if (!geometrySimplify.simplify(reference, decimate, mapped, faceCount)) {
		std::cout << "Error: simplify failed \n";
		exit(1);
	}

	for (int32_t i = 0; i < decimate.XYZ.size(); i++) {
		for (int j = 0; j < 3; j++) {
			decimate.XYZ[i][j] = std::round(decimate.XYZ[i][j]);
		}
	}
	if (!geometrySimplify.removeDuplVert(decimate)) {
		std::cout << "Error: removeDuplVert failed \n";
		exit(1);
	}
	return true;
}

void mcemEncoder::cleanMesh(TriMesh<MeshType>& inputMesh, encoderParams const& params) {

	if (params.removeDuplicateVert) {
		std::vector<Vector3D<MeshType>> XYZ1;
		std::vector<Vector3D<uint32_t>> FaceXYZ1;
		std::vector<uint32_t>           uniqueOriginalIdx, original2unique;
		removeDuplicateVert(inputMesh.XYZ, inputMesh.FaceXYZ, XYZ1, FaceXYZ1, uniqueOriginalIdx, original2unique);
#if DEBUG_ENC_LOG
		std::cout << "cleanMeshes removed " << inputMesh.XYZ.size() - XYZ1.size() << " vertices" << std::endl;
		std::cout << "cleanMeshes removed " << inputMesh.FaceXYZ.size() - FaceXYZ1.size() << " faces" << std::endl;
#endif
		std::swap(inputMesh.XYZ, XYZ1);
		std::swap(inputMesh.FaceXYZ, FaceXYZ1);

		removeDegenerateTriangles(inputMesh);
#if DEBUG_ENC_LOG
		std::cout << "cleanMeshes removed " << FaceXYZ1.size() - inputMesh.FaceXYZ.size() << " faces" << std::endl;
#endif
	}
}

void
mcemEncoder::computeXYZmapping(TriMesh<MeshType>           inputMesh,
							   TriMesh<MeshType>           reconMesh,
							   std::vector<uint32_t>&      reconXYZ2inputXYZ,
							   uint32_t					   subdivFitSubdivIterCount,
							   bool                        isSubdiv) const {
	if (isSubdiv) {
		inputMesh.subdivMidpoint(subdivFitSubdivIterCount);
	}

	std::map<Vector3D<MeshType>, uint32_t> xyz2idx;
	for (uint32_t i = 0; i < inputMesh.XYZ.size(); ++i) { xyz2idx[inputMesh.XYZ[i]] = i; }

	reconXYZ2inputXYZ.resize(reconMesh.XYZ.size());
	for (uint32_t i = 0; i < reconMesh.XYZ.size(); ++i) {
		auto it = xyz2idx.find(reconMesh.XYZ[i]);
		if (it != xyz2idx.end()) {
			reconXYZ2inputXYZ[i] = it->second;
		}
		else {
			assert(false);
		}
	}
}

void mcemEncoder::encodeMetadata(encoderParams&	params, 
								 ofstream&		bitstreamFile) {
	bitstreamFile.write(reinterpret_cast<const char*>(&params.subdivFitSubdivIterCount), 2);
}

void mcemEncoder::encodeBaseMesh(MeshBundle&        MB,
								TriMesh<MeshType>&  reconMesh,
								encoderParams&		params,
								encoderDracoParams&	paramsDraco,
								bool                flag_lossless_geometry) {
	auto& base = MB.base;
	if (!flag_lossless_geometry) {
		// quantize uv
		uint8_t input_tex_bitdepth = ceilLog2(std::max(params.inputWidth, params.inputHeight) + 1);
		params.scaleUV = ((1 << input_tex_bitdepth) - 1) / (params.bboxMaxUV - params.bboxMinUV).Linfnorm();
		if (params.flag_texture) {
			quantize(base.UV, input_tex_bitdepth, params.scaleUV, params.bboxMinUV, params.bboxMaxUV, false);
		}
	}

	std::vector<char>        bitstream;
	ofstream bitstreamFile;       ///< file bitstream
	bitstreamFile.open(params.output_mesh, fstream::binary | fstream::out);
	if (!bitstreamFile) {
		cerr << "Error: failed to open bitstream file " << params.output_mesh << " for writing!"
			<< endl;
	}
	encodeMetadata(params, bitstreamFile);
	// encode base mesh using draco
	encodeDraco(base, reconMesh, bitstream, paramsDraco);
	bitstreamFile.write(bitstream.data(), bitstream.size());
	bitstreamFile.close();

	if (!flag_lossless_geometry) {
		// mapping subdivided-reconMesh with subdivided-baseMesh
		std::vector<uint32_t> reconXYZ2inputXYZ;
		subdivBaseMesh(MB, reconMesh, params.subdivFitSubdivIterCount);
		computeXYZmapping(base, reconMesh, reconXYZ2inputXYZ, params.subdivFitSubdivIterCount);

		// rearrange MB.subdivFit using the mapping relation for geometry residual
		auto subdivFit0 = reconMesh;
		std::swap(subdivFit0, MB.subdivFit);
		for (uint32_t i = 0; i < MB.subdivFit.XYZ.size(); ++i) {
			MB.subdivFit.XYZ[i] = subdivFit0.XYZ[reconXYZ2inputXYZ[i]];
		}
	}
}

void mcemEncoder::computeGeometryResidual(MeshBundle& MB, TriMesh<MeshType> const& reconMesh) {
	MB.disp.resize(reconMesh.XYZ.size());

	for (uint32_t i = 0; i < reconMesh.XYZ.size(); ++i) {
		auto const& xyz0 = reconMesh.XYZ[i];
		auto const& xyz1 = MB.subdivFit.XYZ[i];
		auto disp = xyz1 - xyz0;
		for (int j = 0; j < 3; j++) {
			MB.disp[i][j] = std::round(disp[j]);
		}
	}
}

void mcemEncoder::encodeGeometryResidual(MeshBundle& MB, encoderParams& params){
	ofstream bitstreamFile;       ///< file bitstream
	bitstreamFile.open(params.output_geometry, fstream::binary | fstream::out);
	if (!bitstreamFile) {
		cerr << "Error: failed to open bitstream file " << params.output_geometry << " for writing!"
			<< endl;
	}
	init_aec_context_tab();
	m_encBac.setBitstreamBuffer(m_bufferChunk);
	m_encBac.initBac();

	for (uint32_t idx = 0; idx < MB.disp.size(); idx++) {	
		for (uint32_t dim = 0; dim < 3; dim++) {
			int32_t res = MB.disp[idx][dim];
			m_encBac.encodeResidual(res, dim);
		}
	}
	m_encBac.encodeTerminationFlag();
	m_encBac.encodeFinish();
	m_bufferChunk.writeToBitstream(&bitstreamFile, m_encBac.getBitStreamLength());
	m_bufferChunk.reset();
	m_encBac.reset();
	bitstreamFile.close();
}

void mcemEncoder::encode(encoderParams&      params,
						 encoderDracoParams& paramsDraco, 
						 encoderHPMParams&	 paramsHPM) {
	// read files
	TriMesh<MeshType> inputMesh;
	inputMesh.readOBJ(params.input_mesh);
	Image3<uint8_t> inputTexture;
	if (params.flag_texture) {
		inputTexture.readImage(paramsHPM.input, params.inputWidth, params.inputHeight);	
	}
	else {
		inputMesh.UV.clear();
		inputMesh.FaceUV.clear();
	}
	double userTimeTotal = 0.0;
	double userTimeGeometry = 0.0;
	double userTimeTexture = 0.0;
	clock_t userTimeTotalBegin = clock();
	cleanMesh(inputMesh, params);

	TriMesh<MeshType> decimatedMesh, decimatedReparametrizedMesh;
	TriMesh<MeshType> reconMesh;
	Image3<uint8_t> reconTexture;
	MeshBundle        MB;
	std::vector<uint32_t> mtlPos;
	bool flag_lossless_geometry = false;
	if (params.targetTriangleRatio == 1.0 && params.subdivFitSubdivIterCount == 0) {
		clock_t userTimeGeometryBegin1 = clock();
		flag_lossless_geometry = true;
		MB.base = inputMesh;
		MB.subdivFit = inputMesh;
		encodeBaseMesh(MB, reconMesh, params, paramsDraco, flag_lossless_geometry);
		userTimeGeometry = (double)(clock() - userTimeGeometryBegin1) / CLOCKS_PER_SEC;
		reconTexture = inputTexture;
	}
	else { // decimation
		inputMesh.computeNormals();
		MB.reference = inputMesh;
		simplifyInput(inputMesh, decimatedMesh, MB.mapped, MB.reference, params);

		// reparameterization
		if (params.flag_texture) {
			TextureParametrization parametrizer;
			parametrizer.generate(decimatedMesh, decimatedReparametrizedMesh, params);
		}
		else {
			decimatedReparametrizedMesh = decimatedMesh;
		}
		// subdivision
		SubdivFit subdivFit;
		SubdivFitParameters subdivfitParams;
		subdivfitParams.subdivFitSubdivIterCount = params.subdivFitSubdivIterCount;
		subdivfitParams.removeDuplicateVert = params.removeDuplicateVert;
		subdivFit.generate(MB.reference, decimatedReparametrizedMesh, MB.mapped, MB.base, MB.subdivFit, subdivfitParams);
		
		clock_t userTimeGeometryBegin2 = clock();
		encodeBaseMesh(MB, reconMesh, params, paramsDraco, flag_lossless_geometry);
		computeGeometryResidual(MB, reconMesh);
		encodeGeometryResidual(MB, params);
		reconMesh.reconGeometry(MB.disp);
		userTimeGeometry = (double)(clock() - userTimeGeometryBegin2) / CLOCKS_PER_SEC;

		if (params.flag_texture) {
			// color transfer
			TransferColor transferColor;
			transferColor.textureTransfer(inputMesh, inputTexture, reconMesh, reconTexture, params);
		}
	}
	
	clock_t userTimeTextureBegin = clock();
	if (params.flag_texture) {
		encodeTexture(reconTexture, params, paramsHPM);
	}
	userTimeTexture = (double)(clock() - userTimeTextureBegin) / CLOCKS_PER_SEC;
	userTimeTotal = (double)(clock() - userTimeTotalBegin) / CLOCKS_PER_SEC;

	if (params.flag_recon_enc) {
		reconMesh.writeOBJ(params.recon_mesh, mtlPos, true);
	}
	std::cout << "Output number of faces:" << reconMesh.FaceXYZ.size() << std::endl;
	std::cout << "Total encode time (user):" << userTimeTotal << " sec." << std::endl;
	std::cout << "Total geometry encode time (user):" << userTimeGeometry << " sec." << std::endl;
	std::cout << "Total texture encode time (user):" << userTimeTexture << " sec." << std::endl;
}
