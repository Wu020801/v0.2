#pragma once
#include "draco/compression/config/compression_shared.h"
#include "draco/compression/encode.h"
#include "draco/compression/expert_encode.h"
#include "draco/core/cycle_timer.h"
#include "draco/io/file_utils.h"
#include "draco/io/mesh_io.h"
#include "draco/io/point_cloud_io.h"

#include "draco/io/obj_encoder.h"
#include "draco/io/parser_utils.h"

#include "encoder.hpp"

//============================================================================

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

template<typename T>
void
convert(TriMesh<T> const& meshSrc, draco::Mesh& meshDraco, int input_pos_bitdepth) {
	auto const& xyzCount = meshSrc.XYZ.size();
	auto const& uvCount = meshSrc.UV.size();
	auto const& normalCount = meshSrc.Normal.size();
	auto const& faceCount = meshSrc.FaceXYZ.size();
	meshDraco.SetNumFaces(faceCount);
	meshDraco.set_num_points(3 * faceCount);
	
	int                      xyzAttIdx, uvAttIdx, normalAttIdx;
	const bool               use_identity_mapping = false;
	
	if (input_pos_bitdepth < 32) {
		draco::GeometryAttribute ga;
		ga.Init(draco::GeometryAttribute::POSITION, nullptr, 3, draco::DT_INT32, false, sizeof(int32_t) * 3, 0);
		xyzAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, xyzCount);

		if (uvCount > 0) {
			draco::GeometryAttribute ga;
			ga.Init(draco::GeometryAttribute::TEX_COORD, nullptr, 2, draco::DT_INT32, false, sizeof(int32_t) * 2, 0);
			uvAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, uvCount);
		}

		if (normalCount > 0) {
			draco::GeometryAttribute ga;
			ga.Init(draco::GeometryAttribute::NORMAL, nullptr, 3, draco::DT_INT32, false, sizeof(int32_t) * 3, 0);
			normalAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, normalCount);
		}

		draco::AttributeValueIndex xyzAttValIdx(0), uvAttValIdx(0), normalAttValIdx(0);
		for (uint32_t i = 0; i < xyzCount; ++i) {
			Vector3D<int32_t> xyz(meshSrc.XYZ[i]);
			meshDraco.attribute(xyzAttIdx)->SetAttributeValue(xyzAttValIdx++, xyz.data());
		}
		if (uvCount > 0) {
			for (uint32_t i = 0; i < uvCount; ++i) {
				Vector2D<int32_t> uv(meshSrc.UV[i]);
				meshDraco.attribute(uvAttIdx)->SetAttributeValue(uvAttValIdx++, uv.data());
			}
		}
		if (normalCount > 0) {
			for (uint32_t i = 0; i < normalCount; ++i) {
				Vector3D<int32_t> normal(meshSrc.Normal[i]);
				meshDraco.attribute(normalAttIdx)->SetAttributeValue(normalAttValIdx++, normal.data());
			}
		}
	}
	else {
		draco::GeometryAttribute ga;
		ga.Init(draco::GeometryAttribute::POSITION, nullptr, 3, draco::DT_INT64, false, sizeof(int64_t) * 3, 0);
		xyzAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, xyzCount);

		if (uvCount > 0) {
			draco::GeometryAttribute ga;
			ga.Init(draco::GeometryAttribute::TEX_COORD, nullptr, 2, draco::DT_INT64, false, sizeof(int64_t) * 2, 0);
			uvAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, uvCount);
		}

		if (normalCount > 0) {
			draco::GeometryAttribute ga;
			ga.Init(draco::GeometryAttribute::NORMAL, nullptr, 3, draco::DT_INT64, false, sizeof(int64_t) * 3, 0);
			normalAttIdx = meshDraco.AddAttribute(ga, use_identity_mapping, normalCount);
		}

		draco::AttributeValueIndex xyzAttValIdx(0), uvAttValIdx(0), normalAttValIdx(0);
		for (uint32_t i = 0; i < xyzCount; ++i) {
			Vector3D<int64_t> xyz(meshSrc.XYZ[i]);
			meshDraco.attribute(xyzAttIdx)->SetAttributeValue(xyzAttValIdx++, xyz.data());
		}
		if (uvCount > 0) {
			for (uint32_t i = 0; i < uvCount; ++i) {
				Vector2D<int64_t> uv(meshSrc.UV[i]);
				meshDraco.attribute(uvAttIdx)->SetAttributeValue(uvAttValIdx++, uv.data());
			}
		}
		if (normalCount > 0) {
			for (uint32_t i = 0; i < normalCount; ++i) {
				Vector3D<int64_t> normal(meshSrc.Normal[i]);
				meshDraco.attribute(normalAttIdx)->SetAttributeValue(normalAttValIdx++, normal.data());
			}
		}
	}


	for (uint32_t i = 0; i < faceCount; ++i) {
		draco::Mesh::Face face;
		Vector3D<int32_t> faceXYZ(meshSrc.FaceXYZ[i]);
		for (uint32_t j = 0; j < 3; ++j) {
			face[j] = 3 * i + j;
			meshDraco.attribute(xyzAttIdx)->SetPointMapEntry(face[j], draco::AttributeValueIndex(faceXYZ[j]));
		}
		if (uvCount > 0) {
			Vector3D<int32_t> faceUV(meshSrc.FaceUV[i]);
			for (uint32_t j = 0; j < 3; ++j) {
				meshDraco.attribute(uvAttIdx)->SetPointMapEntry(face[j], draco::AttributeValueIndex(faceUV[j]));
			}
		}
		if (normalCount > 0) {
			Vector3D<int32_t> normalUV(meshSrc.FaceNormal[i]);
			for (uint32_t j = 0; j < 3; ++j) {
				meshDraco.attribute(normalAttIdx)->SetPointMapEntry(face[j], draco::AttributeValueIndex(normalUV[j]));
			}
		}
		meshDraco.SetFace(draco::FaceIndex(i), face);
	}
	meshDraco.DeduplicatePointIds();
}

//============================================================================

void PrintOptions(const draco::PointCloud& pc, const encoderDracoParams& options) {
    printf("Encoder options:\n");
    printf("  Compression level = %d\n", options.compression_level);
    if (options.pos_quantization_bits == 0) {
        printf("  Positions: No quantization\n");
    }
    else {
        printf("  Positions: Quantization = %d bits\n",
            options.pos_quantization_bits);
    }

    if (pc.GetNamedAttributeId(draco::GeometryAttribute::TEX_COORD) >= 0) {
        if (options.tex_coords_quantization_bits == 0) {
            printf("  Texture coordinates: No quantization\n");
        }
        else {
            printf("  Texture coordinates: Quantization = %d bits\n",
                options.tex_coords_quantization_bits);
        }
    }
    else if (options.tex_coords_deleted) {
        printf("  Texture coordinates: Skipped\n");
    }

    if (pc.GetNamedAttributeId(draco::GeometryAttribute::NORMAL) >= 0) {
        if (options.normals_quantization_bits == 0) {
            printf("  Normals: No quantization\n");
        }
        else {
            printf("  Normals: Quantization = %d bits\n",
                options.normals_quantization_bits);
        }
    }
    else if (options.normals_deleted) {
        printf("  Normals: Skipped\n");
    }

    if (pc.GetNamedAttributeId(draco::GeometryAttribute::GENERIC) >= 0) {
        if (options.generic_quantization_bits == 0) {
            printf("  Generic: No quantization\n");
        }
        else {
            printf("  Generic: Quantization = %d bits\n",
                options.generic_quantization_bits);
        }
    }
    else if (options.generic_deleted) {
        printf("  Generic: Skipped\n");
    }
    printf("\n");
}

bool EncodePointCloudToFile(const draco::PointCloud& pc, const std::string& file,
							draco::ExpertEncoder* encoder) {
    draco::CycleTimer timer;
    // Encode the geometry.
    draco::EncoderBuffer buffer;
    timer.Start();
    const draco::Status status = encoder->EncodeToBuffer(&buffer);
    if (!status.ok()) {
        printf("Failed to encode the point cloud.\n");
        printf("%s\n", status.error_msg());
        return false;
    }
    timer.Stop();
    // Save the encoded geometry into a file.
    if (!draco::WriteBufferToFile(buffer.data(), buffer.size(), file)) {
        printf("Failed to write the output file.\n");
        return false;
    }
    printf("Encoded point cloud saved to %s (%" PRId64 " ms to encode).\n",
        file.c_str(), timer.GetInMs());
    printf("\nEncoded size = %zu bytes\n\n", buffer.size());
    return true;
}

bool EncodeMeshToFile(const draco::Mesh& mesh, const std::string& file,
					  draco::ExpertEncoder* encoder) {
    draco::CycleTimer timer;
    // Encode the geometry.
    draco::EncoderBuffer buffer;
    timer.Start();
    const draco::Status status = encoder->EncodeToBuffer(&buffer);
    if (!status.ok()) {
        printf("Failed to encode the mesh.\n");
        printf("%s\n", status.error_msg());
        return false;
    }
    timer.Stop();
    // Save the encoded geometry into a file.
    if (!draco::WriteBufferToFile(buffer.data(), buffer.size(), file)) {
        printf("Failed to create the output file.\n");
        return false;
    }
    printf("Encoded mesh saved to %s (%" PRId64 " ms to encode).\n", file.c_str(),
        timer.GetInMs());
    printf("\nEncoded size = %zu bytes\n\n", buffer.size());
    return true;
}

int ReturnError(const draco::Status& status) {
	printf("Failed to decode the input file %s\n", status.error_msg());
	return -1;
}

bool mcemEncoder::encodeDracoFile(encoderDracoParams& options, TriMesh<MeshType>& reconMesh) {
    std::unique_ptr<draco::PointCloud> pc;
    draco::Mesh* mesh = nullptr;
    if (!options.is_point_cloud) {
        draco::Options load_options;
        load_options.SetBool("use_metadata", options.use_metadata);
        load_options.SetBool("preserve_polygons", options.preserve_polygons);
        auto maybe_mesh = draco::ReadMeshFromFile(options.input, load_options);
        if (!maybe_mesh.ok()) {
            printf("Failed loading the input mesh: %s.\n",
                maybe_mesh.status().error_msg());
            return false;
        }
        mesh = maybe_mesh.value().get();
        pc = std::move(maybe_mesh).value();
    }
    else {
        auto maybe_pc = draco::ReadPointCloudFromFile(options.input);
        if (!maybe_pc.ok()) {
            printf("Failed loading the input point cloud: %s.\n",
                maybe_pc.status().error_msg());
            return false;
        }
        pc = std::move(maybe_pc).value();
    }

    if (options.pos_quantization_bits < 0) {
        printf("Error: Position attribute cannot be skipped.\n");
        return false;
    }

    // Delete attributes if needed. This needs to happen before we set any
    // quantization settings.
    if (options.tex_coords_quantization_bits < 0) {
        if (pc->NumNamedAttributes(draco::GeometryAttribute::TEX_COORD) > 0) {
            options.tex_coords_deleted = true;
        }
        while (pc->NumNamedAttributes(draco::GeometryAttribute::TEX_COORD) > 0) {
            pc->DeleteAttribute(
                pc->GetNamedAttributeId(draco::GeometryAttribute::TEX_COORD, 0));
        }
    }
    if (options.normals_quantization_bits < 0) {
        if (pc->NumNamedAttributes(draco::GeometryAttribute::NORMAL) > 0) {
            options.normals_deleted = true;
        }
        while (pc->NumNamedAttributes(draco::GeometryAttribute::NORMAL) > 0) {
            pc->DeleteAttribute(
                pc->GetNamedAttributeId(draco::GeometryAttribute::NORMAL, 0));
        }
    }
    if (options.generic_quantization_bits < 0) {
        if (pc->NumNamedAttributes(draco::GeometryAttribute::GENERIC) > 0) {
            options.generic_deleted = true;
        }
        while (pc->NumNamedAttributes(draco::GeometryAttribute::GENERIC) > 0) {
            pc->DeleteAttribute(
                pc->GetNamedAttributeId(draco::GeometryAttribute::GENERIC, 0));
        }
    }
#ifdef DRACO_ATTRIBUTE_INDICES_DEDUPLICATION_SUPPORTED
    // If any attribute has been deleted, run deduplication of point indices again
    // as some points can be possibly combined.
    if (options.tex_coords_deleted || options.normals_deleted ||
        options.generic_deleted) {
        pc->DeduplicatePointIds();
    }
#endif

    // Convert compression level to speed (that 0 = slowest, 10 = fastest).
    const int speed = 10 - options.compression_level;

    draco::Encoder encoder;
    // Setup encoder options.
    if (options.pos_quantization_bits > 0) {
        encoder.SetAttributeQuantization(draco::GeometryAttribute::POSITION,
            options.pos_quantization_bits);
    }
    if (options.tex_coords_quantization_bits > 0) {
        encoder.SetAttributeQuantization(draco::GeometryAttribute::TEX_COORD,
            options.tex_coords_quantization_bits);
    }
    if (options.normals_quantization_bits > 0) {
        encoder.SetAttributeQuantization(draco::GeometryAttribute::NORMAL,
            options.normals_quantization_bits);
    }
    if (options.generic_quantization_bits > 0) {
        encoder.SetAttributeQuantization(draco::GeometryAttribute::GENERIC,
            options.generic_quantization_bits);
    }
    encoder.SetSpeedOptions(speed, speed);

    if (options.output.empty()) {
        // Create a default output file by attaching .drc to the input file name.
        options.output = options.input + ".drc";
    }

    PrintOptions(*pc, options);

    const bool input_is_mesh = mesh && mesh->num_faces() > 0;

    // Convert to ExpertEncoder that allows us to set per-attribute options.
    std::unique_ptr<draco::ExpertEncoder> expert_encoder;
    if (input_is_mesh) {
        expert_encoder.reset(new draco::ExpertEncoder(*mesh));
    }
    else {
        expert_encoder.reset(new draco::ExpertEncoder(*pc));
    }
    expert_encoder->Reset(encoder.CreateExpertEncoderOptions(*pc));

    // Check if there is an attribute that stores polygon edges. If so, we disable
    // the default prediction scheme for the attribute as it actually makes the
    // compression worse.
    const int poly_att_id =
        pc->GetAttributeIdByMetadataEntry("name", "added_edges");
    if (poly_att_id != -1) {
        expert_encoder->SetAttributePredictionScheme(
            poly_att_id, draco::PredictionSchemeMethod::PREDICTION_NONE);
    }

	bool ret = false;

    if (input_is_mesh) {
        ret = EncodeMeshToFile(*mesh, options.output, expert_encoder.get());
    }
    else {
        ret = EncodePointCloudToFile(*pc, options.output, expert_encoder.get());
    }

    if (ret && options.compression_level < 10) {
        printf(
            "For better compression, increase the compression level up to '-cl 10' "
            ".\n\n");
    }

	// draco decoder
	std::vector<char> data;
	if (!draco::ReadFileToBuffer(options.output, &data)) {
		printf("Failed opening the input file.\n");
		return -1;
	}
	draco::DecoderBuffer buffer;
	buffer.Init(data.data(), data.size());
	// Decode the input data into a geometry.
	auto type_statusor = draco::Decoder::GetEncodedGeometryType(&buffer);
	if (!type_statusor.ok()) {
		return ReturnError(type_statusor.status());
	}
	const draco::EncodedGeometryType geom_type = type_statusor.value();
	if (geom_type == draco::TRIANGULAR_MESH) {
		draco::Decoder decoder;
		auto statusor = decoder.DecodeMeshFromBuffer(&buffer);
		if (!statusor.ok()) {
			return ReturnError(statusor.status());
		}
		std::unique_ptr<draco::Mesh> in_mesh = std::move(statusor).value();
		std::unique_ptr<draco::Mesh> decmesh = std::move(statusor).value();
		if (in_mesh) {
			mesh = in_mesh.get();
			pc = std::move(in_mesh);
		}
		// reconMesh
		convert(*mesh, reconMesh);
	}
	else if (geom_type == draco::POINT_CLOUD) {
		// Failed to decode it as mesh, so let's try to decode it as a point cloud.
		draco::Decoder decoder;
		auto statusor = decoder.DecodePointCloudFromBuffer(&buffer);
		if (!statusor.ok()) {
			return ReturnError(statusor.status());
		}
		pc = std::move(statusor).value();
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
	else {
		printf("Invalid output file extension. Use .obj .ply or .stl.\n");
		return -1;
	}

	return 0;
}

bool mcemEncoder::encodeDraco(TriMesh<MeshType>&           baseMesh,
								TriMesh<MeshType>&         reconMesh,
								std::vector<char>&         bitstream,
								encoderDracoParams&		   params) {
	draco::Mesh meshDraco;
	convert(baseMesh, meshDraco, params.input_pos_bitdepth);

	draco::Encoder encoder;
	encoder.SetSpeedOptions(10 - params.compression_level, 10 - params.compression_level);
	encoder.SetEncodingMethod(draco::MeshEncoderMethod::MESH_EDGEBREAKER_ENCODING);
	encoder.SetAttributeQuantization(draco::GeometryAttribute::POSITION, params.pos_quantization_bits);
	encoder.SetAttributeQuantization(draco::GeometryAttribute::TEX_COORD, params.tex_coords_quantization_bits);
	if (params.normals_quantization_bits > 0) {
		encoder.SetAttributeQuantization(draco::GeometryAttribute::NORMAL, params.normals_quantization_bits);
	}

	draco::EncoderBuffer encBuffer;
	const draco::Status status = encoder.EncodeMeshToBuffer(meshDraco, &encBuffer);
	if (!status.ok()) {
		printf("Draco Failed to encode the mesh: %s\n", status.error_msg());
		exit(-1);
	}
	// printf("EncodeMeshToBuffer => buffer size = %zu \n", encBuffer.size());
	fflush(stdout);
	bitstream.resize(encBuffer.size());
	std::copy(encBuffer.data(), encBuffer.data() + encBuffer.size(), bitstream.data());

	// reconMesh
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