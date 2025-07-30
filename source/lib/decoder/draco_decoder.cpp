#pragma once
#include "draco/compression/decode.h"
#include "draco/core/cycle_timer.h"
#include "draco/io/file_utils.h"
#include "draco/io/obj_encoder.h"
#include "draco/io/parser_utils.h"
#include "draco/io/ply_encoder.h"
#include "draco/io/stl_encoder.h"
#include "decoder.hpp"

template<typename T>
void
convert(draco::Mesh const& meshDraco, TriMesh<T>& meshSrc) {
	const auto* xyzAtt = meshDraco.GetNamedAttribute(draco::GeometryAttribute::POSITION);
	const auto* uvAtt = meshDraco.GetNamedAttribute(draco::GeometryAttribute::TEX_COORD);
	const auto* normalAtt = meshDraco.GetNamedAttribute(draco::GeometryAttribute::NORMAL);

	if (xyzAtt) {
		for (draco::AttributeValueIndex i(0); i < static_cast<uint32_t>(xyzAtt->size()); ++i) {
			Vector3D<int64_t> xyz;
			if (!xyzAtt->ConvertValue<int64_t, 3>(i, xyz.data())) { return; }
			meshSrc.XYZ.push_back(xyz);
		}
	}
	if (uvAtt) {
		for (draco::AttributeValueIndex i(0); i < static_cast<uint32_t>(uvAtt->size()); ++i) {
			Vector2D<int64_t> uv;
			if (!uvAtt->ConvertValue<int64_t, 2>(i, uv.data())) { return; }
			meshSrc.UV.push_back(uv);
		}
	}
	if (normalAtt) {
		for (draco::AttributeValueIndex i(0); i < static_cast<uint32_t>(normalAtt->size()); ++i) {
			Vector3D<int64_t> normal;
			if (!normalAtt->ConvertValue<int64_t, 3>(i, normal.data())) { return; }
			meshSrc.Normal.push_back(normal);
		}
	}

	for (draco::FaceIndex f(0); f < meshDraco.num_faces(); ++f) {
		const auto& face = meshDraco.face(f);
		uint32_t    i = xyzAtt->mapped_index(face[0]).value();
		uint32_t    j = xyzAtt->mapped_index(face[1]).value();
		uint32_t    k = xyzAtt->mapped_index(face[2]).value();
		meshSrc.FaceXYZ.emplace_back(i, j, k);
		if (uvAtt && uvAtt->size() > 0) {
			i = uvAtt->mapped_index(face[0]).value();
			j = uvAtt->mapped_index(face[1]).value();
			k = uvAtt->mapped_index(face[2]).value();
			meshSrc.FaceUV.emplace_back(i, j, k);
		}
		if (normalAtt && normalAtt->size() > 0) {
			i = normalAtt->mapped_index(face[0]).value();
			j = normalAtt->mapped_index(face[1]).value();
			k = normalAtt->mapped_index(face[2]).value();
			meshSrc.FaceNormal.emplace_back(i, j, k);
		}
	}
}

//============================================================================

int ReturnError(const draco::Status& status) {
    printf("Failed to decode the input file %s\n", status.error_msg());
    return -1;
}

int mcemDecoder::decodeDracoFile(decoderDracoParams &options) {
    std::vector<char> data;
    if (!draco::ReadFileToBuffer(options.input, &data)) {
        printf("Failed opening the input file.\n");
        return -1;
    }

    if (data.empty()) {
        printf("Empty input file.\n");
        return -1;
    }

    // Create a draco decoding buffer. Note that no data is copied in this step.
    draco::DecoderBuffer buffer;
    buffer.Init(data.data(), data.size());

    draco::CycleTimer timer;
    // Decode the input data into a geometry.
    std::unique_ptr<draco::PointCloud> pc;
    draco::Mesh* mesh = nullptr;
    auto type_statusor = draco::Decoder::GetEncodedGeometryType(&buffer);
    if (!type_statusor.ok()) {
        return ReturnError(type_statusor.status());
    }
    const draco::EncodedGeometryType geom_type = type_statusor.value();
    if (geom_type == draco::TRIANGULAR_MESH) {
        timer.Start();
        draco::Decoder decoder;
        auto statusor = decoder.DecodeMeshFromBuffer(&buffer);
        if (!statusor.ok()) {
            return ReturnError(statusor.status());
        }
        std::unique_ptr<draco::Mesh> in_mesh = std::move(statusor).value();
        timer.Stop();
        if (in_mesh) {
            mesh = in_mesh.get();
            pc = std::move(in_mesh);
        }
    }
    else if (geom_type == draco::POINT_CLOUD) {
        // Failed to decode it as mesh, so let's try to decode it as a point cloud.
        timer.Start();
        draco::Decoder decoder;
        auto statusor = decoder.DecodePointCloudFromBuffer(&buffer);
        if (!statusor.ok()) {
            return ReturnError(statusor.status());
        }
        pc = std::move(statusor).value();
        timer.Stop();
    }

    if (pc == nullptr) {
        printf("Failed to decode the input file.\n");
        return -1;
    }

    if (options.output.empty()) {
        // Save the output model into a ply file.
        options.output = options.input + ".ply";
    }

    // Save the decoded geometry into a file.
    const std::string extension = draco::parser::ToLower(
        options.output.size() >= 4
        ? options.output.substr(options.output.size() - 4)
        : options.output);

    if (extension == ".obj") {
        draco::ObjEncoder obj_encoder;
        if (mesh) {
            if (!obj_encoder.EncodeToFile(*mesh, options.output)) {
                printf("Failed to store the decoded mesh as OBJ.\n");
                return -1;
            }
        }
        else {
            if (!obj_encoder.EncodeToFile(*pc, options.output)) {
                printf("Failed to store the decoded point cloud as OBJ.\n");
                return -1;
            }
        }
    }
    else if (extension == ".ply") {
        draco::PlyEncoder ply_encoder;
        if (mesh) {
            if (!ply_encoder.EncodeToFile(*mesh, options.output)) {
                printf("Failed to store the decoded mesh as PLY.\n");
                return -1;
            }
        }
        else {
            if (!ply_encoder.EncodeToFile(*pc, options.output)) {
                printf("Failed to store the decoded point cloud as PLY.\n");
                return -1;
            }
        }
    }
    else if (extension == ".stl") {
        draco::StlEncoder stl_encoder;
        if (mesh) {
            draco::Status s = stl_encoder.EncodeToFile(*mesh, options.output);
            if (s.code() != draco::Status::OK) {
                printf("Failed to store the decoded mesh as STL.\n");
                return -1;
            }
        }
        else {
            printf("Can't store a point cloud as STL.\n");
            return -1;
        }
    }
    else {
        printf("Invalid output file extension. Use .obj .ply or .stl.\n");
        return -1;
    }
    printf("Decoded geometry saved to %s (%" PRId64 " ms to decode)\n",
        options.output.c_str(), timer.GetInMs());
    return 0;

}

bool mcemDecoder::decodeDraco(TriMesh<MeshType>&           reconMesh,
							 std::vector<char>&           bitstream) {
	draco::Decoder decoder;
	draco::DecoderBuffer decBuffer;
	decBuffer.Init((const char*)bitstream.data(), bitstream.size());
	auto type = draco::Decoder::GetEncodedGeometryType(&decBuffer);
	if (!type.ok()) {
		printf("Draco Failed GetEncodedGeometryType: %s.\n", type.status().error_msg());
		exit(-1);
	}
	auto statusType = draco::Decoder::GetEncodedGeometryType(&decBuffer);
	if (type.value() == draco::TRIANGULAR_MESH) {
		draco::Decoder decoder;
		auto           status = decoder.DecodeMeshFromBuffer(&decBuffer);
		if (!status.ok()) {
			printf("Draco Failed DecodeMeshFromBuffer: %s.\n", status.status().error_msg());
			exit(-1);
		}
		std::unique_ptr<draco::Mesh> decMesh = std::move(status).value();
		if (decMesh) {
			convert(*(decMesh.get()), reconMesh);
		}
		else {
			printf("Draco Failed no in mesh \n");
			exit(-1);
		}
	}
	else {
		printf("Draco Failed no mesh type not supported.\n");
		exit(-1);
	}
	return true;
}