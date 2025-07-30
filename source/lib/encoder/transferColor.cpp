/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#include "transferColor.hpp"
#include "kdtree.hpp"
#include "encoder.hpp"

#include <nanoflann/nanoflann.hpp>
#include <unordered_map>

//============================================================================

template<typename T>
inline static T
clip(const T& n, const T& lower, const T& upper) {
  return (std::max)(lower, (std::min)(n, upper));
}

Vector3D<MeshType>
closestPointColor(Vector3D<MeshType> const& xyz,
                  TriMesh<MeshType> const&  mesh,
                  Image3<uint8_t> const&    texture,
                  std::vector<uint32_t>&    xyz2uv,
                  AdjInfo const&            vertex2triangle,
                  KdTree<MeshType> const&   kdtree,
                  MeshType*                 minSqDist = nullptr) {
  MeshType dist2;
  if (minSqDist == nullptr) minSqDist = &dist2;
  uint32_t xyzIdx;
  kdtree.query(xyz.data(), 1, &xyzIdx, minSqDist);
  auto        uvClosest = mesh.UV[xyz2uv[xyzIdx]];
  auto const& faceIdx   = vertex2triangle.getNeighbors(xyzIdx);
  for (auto const& fidx : faceIdx) {
    auto const&        faceXYZ = mesh.FaceXYZ[fidx];
    auto const&        xyz0    = mesh.XYZ[faceXYZ[0]];
    auto const&        xyz1    = mesh.XYZ[faceXYZ[1]];
    auto const&        xyz2    = mesh.XYZ[faceXYZ[2]];
    Vector3D<MeshType> bc;
    auto const&        xyzClosest = closestPointInTriangle(xyz, xyz0, xyz1, xyz2, &bc);
    auto const&        sqDist     = (xyzClosest - xyz).L2normSq();
    if (sqDist < *minSqDist) {
      *minSqDist         = sqDist;
      auto const& faceUV = mesh.FaceUV[fidx];
      auto const& uv0    = mesh.UV[faceUV[0]];
      auto const& uv1    = mesh.UV[faceUV[1]];
      auto const& uv2    = mesh.UV[faceUV[2]];
      uvClosest          = bc[0] * uv0 + bc[1] * uv1 + bc[2] * uv2;
    }
  }
  return texture.interp(uvClosest[0], uvClosest[1]);
}

//============================================================================

void
TransferColor::textureTransfer(TriMesh<MeshType> const&    inputMesh,
                               Image3<uint8_t> const&      inputTexture,
                               TriMesh<MeshType>&          reconMesh,
                               Image3<uint8_t>&            outputTexture,
							   encoderParams& params) {
  auto       subdivInputMesh = inputMesh;
  auto	     outputMesh = reconMesh;
  // dequantize uv
  dequantize(subdivInputMesh.UV, params.scaleUV, params.bboxMinUV);
  dequantize(outputMesh.UV, params.scaleUV, params.bboxMinUV);

  subdivInputMesh.subdivMidpoint(params.subdivFitSubdivIterCount); 
  std::vector<std::vector<uint8_t>> occupancy;
  textureTransferByMeshes(subdivInputMesh, inputTexture, outputMesh, outputTexture, occupancy, params);
  outputTexture.dilatePadding(occupancy, 2 * params.textureTransferPaddingDilateIterationCount);
}

//============================================================================

void
TransferColor::textureTransferByMeshes(TriMesh<MeshType> const&           inputMesh,
                                       Image3<uint8_t> const&             inputTexture,
                                       TriMesh<MeshType>&                 outputMesh,
                                       Image3<uint8_t>&                   outputTexture,
                                       std::vector<std::vector<uint8_t>>& occupancy,
									   encoderParams const&        params) {
  printf("Texture transfer by meshes\n");
  fflush(stdout);
  std::vector<uint32_t> xyz2uv;
  inputMesh.computeXYZ2UV(xyz2uv);
  AdjInfo vertex2triangle;
  computeVertex2triangle(inputMesh.XYZ.size(), inputMesh.FaceXYZ, vertex2triangle);
  KdTree<MeshType> kdtree(3, inputMesh.XYZ, 10);

  if (outputTexture.height() != params.inputHeight || outputTexture.width() != params.inputWidth) {
    outputTexture.resize(params.inputHeight, params.inputWidth);
  }

  auto const otHeight   = outputTexture.height();
  auto const otWidth    = outputTexture.width();
  auto const otHeight_1 = otHeight - 1;
  auto const otWidth_1  = otWidth - 1;

  occupancy.resize(otHeight);
  for (int v = 0; v < otHeight; v++) occupancy[v].resize(otWidth, 0);
  std::vector<std::vector<uint32_t>> faceMap(otHeight, std::vector<uint32_t>(otWidth));

  for (uint32_t fidx = 0; fidx < outputMesh.FaceXYZ.size(); fidx++) {
    auto const& faceUV  = outputMesh.FaceUV[fidx];
    auto const& faceXYZ = outputMesh.FaceXYZ[fidx];

    std::array<Vector2D<MeshType>, 3> uv3 = {
      outputMesh.UV[faceUV[0]], outputMesh.UV[faceUV[1]], outputMesh.UV[faceUV[2]]};
    std::array<Vector3D<MeshType>, 3> xyz3 = {
      outputMesh.XYZ[faceXYZ[0]], outputMesh.XYZ[faceXYZ[1]], outputMesh.XYZ[faceXYZ[2]]};


    MeshType const area = std::abs((uv3[1] - uv3[0]) ^ (uv3[2] - uv3[0]));
    if (area == 0.0) continue;
    MeshType const iarea = 1.0 / area;
    MeshType const umin  = std::min({uv3[0][0], uv3[1][0], uv3[2][0]});
    MeshType const umax  = std::max({uv3[0][0], uv3[1][0], uv3[2][0]});
    MeshType const vmin  = std::min({uv3[0][1], uv3[1][1], uv3[2][1]});
    MeshType const vmax  = std::max({uv3[0][1], uv3[1][1], uv3[2][1]});
    auto const     u0    = uint32_t(std::ceil(umin * otWidth_1));
    auto const     u1    = uint32_t(std::floor(umax * otWidth_1));
    auto const     v0    = uint32_t(std::ceil(vmin * otHeight_1));
    auto const     v1    = uint32_t(std::floor(vmax * otHeight_1));

    for (uint32_t iv = v0; iv <= v1; ++iv) {
      auto const v   = MeshType(iv) / otHeight_1;
      auto const row = otHeight_1 - iv;
	  assert(row <= otHeight_1);
      for (uint32_t iu = u0; iu <= u1; ++iu) {
        auto const& col = iu;
		assert(col <= otWidth_1);
        //if (occupancy[row][col]) continue;
        auto const         u = MeshType(iu) / otWidth_1;
        Vector2D<MeshType> uv(u, v);
        MeshType const     w0 = (uv3[2] - uv3[1]) ^ (uv - uv3[1]) * iarea;
        MeshType const     w1 = (uv3[0] - uv3[2]) ^ (uv - uv3[2]) * iarea;
        MeshType const     w2 = (uv3[1] - uv3[0]) ^ (uv - uv3[0]) * iarea;
        if (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0) {
          auto const& xyz   = w0 * xyz3[0] + w1 * xyz3[1] + w2 * xyz3[2];
          auto const& color = closestPointColor(xyz, inputMesh, inputTexture, xyz2uv, vertex2triangle, kdtree);
          outputTexture.set(row, col, color.round());
          occupancy[row][col] = 255;
          faceMap[row][col]   = fidx;
        }
      }
    }
  }
  std::array<std::array<int32_t, 2>, 4> shift = {{{-1, -0}, {0, -1}, {1, 0}, {0, 1}}};
  for (uint32_t it = 0; it < params.textureTransferPaddingBoundaryIterationCount; it++) {
    auto const checkVal = uint8_t(255 - it);
    for (uint32_t row = 0; row < otHeight; ++row) {
      auto const v = MeshType(otHeight_1 - row) / otHeight_1;
      for (uint32_t col = 0; col < otWidth; ++col) {
        if (occupancy[row][col]) continue;
        auto const         u = MeshType(col) / otWidth_1;
        Vector2D<MeshType> uv(u, v);
        Vector3D<MeshType> color(0.0);
        auto               minSqDist = std::numeric_limits<MeshType>::max();
        uint32_t           count     = 0;
        for (auto const& [drow, dcol] : shift) {
          auto const& row1 = row + drow;
          auto const& col1 = col + dcol;
          if (row1 < 0 || row1 > otHeight_1 || col1 < 0 || col1 > otWidth_1 || occupancy[row1][col1] != checkVal) {
            continue;
          }

          auto const                        fidx    = faceMap[row1][col1];
          auto const&                       faceUV  = outputMesh.FaceUV[fidx];
          auto const&                       faceXYZ = outputMesh.FaceXYZ[fidx];
          std::array<Vector2D<MeshType>, 3> uv3     = {
            outputMesh.UV[faceUV[0]], outputMesh.UV[faceUV[1]], outputMesh.UV[faceUV[2]]};
          std::array<Vector3D<MeshType>, 3> xyz3 = {
            outputMesh.XYZ[faceXYZ[0]], outputMesh.XYZ[faceXYZ[1]], outputMesh.XYZ[faceXYZ[2]]};

          MeshType const area = std::abs((uv3[1] - uv3[0]) ^ (uv3[2] - uv3[0]));
          if (area == 0.0) continue;
          MeshType const iarea = 1.0 / area;
          MeshType const w0    = (uv3[2] - uv3[1]) ^ (uv - uv3[1]) * iarea;
          MeshType const w1    = (uv3[0] - uv3[2]) ^ (uv - uv3[2]) * iarea;
          MeshType const w2    = (uv3[1] - uv3[0]) ^ (uv - uv3[0]) * iarea;
          auto const&    xyz   = w0 * xyz3[0] + w1 * xyz3[1] + w2 * xyz3[2];
          MeshType       sqDist;
          color += closestPointColor(xyz, inputMesh, inputTexture, xyz2uv, vertex2triangle, kdtree, &sqDist);
          if (sqDist < minSqDist) {
            minSqDist         = sqDist;
            faceMap[row][col] = fidx;
          }
          count++;
        }
        if (count) {
          color /= count;
          outputTexture.set(row, col, color.round());
          occupancy[row][col] = checkVal - 1;
        }
      }
    }
  }
}

