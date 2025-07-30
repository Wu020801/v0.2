/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include "vector.hpp"
#include "triangle.hpp"
#include "misc.hpp"

#define MeshType float

//============================================================================

const double kSubdivWeights[2] = { 0.5, 0.5 };

struct SubdivLoD {
	SubdivLoD() = default;
	~SubdivLoD() = default;

	uint32_t xyzCount = 0;
	uint32_t uvCount = 0;
	uint32_t faceXYZCount = 0;
	uint32_t faceUVCount = 0;
};
//============================================================================

template<typename T>
class TriMesh {
public:
	TriMesh() = default;
	~TriMesh() = default;

	template<typename S>
	TriMesh(TriMesh<S> const& mesh) {
		*this = mesh;
	}

	template<typename S>
	TriMesh& operator=(TriMesh<S> const& mesh) {
		Mtllib = mesh.Mtllib;
		FaceXYZ = mesh.FaceXYZ;
		FaceUV = mesh.FaceUV;
		FaceNormal = mesh.FaceNormal;
		XYZ.resize(mesh.XYZ.size());
		for (uint32_t i = 0; i < XYZ.size(); i++) XYZ[i] = mesh.XYZ[i];
		UV.resize(mesh.UV.size());
		for (uint32_t i = 0; i < UV.size(); i++) UV[i] = mesh.UV[i];
		Normal.resize(mesh.Normal.size());
		for (uint32_t i = 0; i < Normal.size(); i++) Normal[i] = mesh.Normal[i];
		return *this;
	}

	TriMesh& setFromFaces(const TriMesh& mesh, const size_t firstFace, const size_t lastFace);

	void          removeDegeneratedVertices();
	void          computeNormals(bool const normalize = true, uint32_t const mode = 0);
	void          computeTriangleNormals(std::vector<Vector3D<T>>& normals, bool const normalize = true) const;
	void          computeUV2XYZ(std::vector<uint32_t>& uv2xyz) const;
	void          computeXYZ2UV(std::vector<uint32_t>& xyz2uv) const;
	void          invertOrientation();
	void          append(TriMesh const& mesh);
	inline double computeFaceArea(const size_t faceIdx) {
		auto& v1 = XYZ[FaceXYZ[faceIdx][0]];
		auto& v2 = XYZ[FaceXYZ[faceIdx][1]];
		auto& v3 = XYZ[FaceXYZ[faceIdx][2]];
		return 0.5 * ((v2 - v1) ^ (v3 - v1)).norm();
	}

	void subdivMidpoint(uint32_t const          iterCount,
						std::vector<SubdivLoD>* lod = nullptr,
						std::vector<uint64_t>*  midpoint2endpoints3D = nullptr,
						std::vector<uint64_t>*  midpoint2endpoints2D = nullptr,
						std::vector<uint32_t>*  faceIndex = nullptr);
	void reconGeometry(std::vector<Vector3D<T>>& disp);
	bool readOBJ(const std::string& filePath);
	bool writeOBJ(const std::string& filePath, const std::vector<uint32_t>& mtlPos = {}, bool isInt = false) const;
	bool read(const std::string& filePath);
	bool write(const std::string& filePath) const;
	void clear();
	void print() const;

	std::vector<Vector3D<T>>        XYZ;
	std::vector<Vector2D<T>>        UV;
	std::vector<Vector3D<T>>        Normal;
	std::vector<Vector3D<uint32_t>> FaceXYZ;
	std::vector<Vector3D<uint32_t>> FaceUV;
	std::vector<Vector3D<uint32_t>> FaceNormal;
	std::vector<Vector3D<T>>        RGB;
	std::string                     Mtllib;
};

//============================================================================

template<typename T, typename S>
bool compareData(T const& data1, S const& data2);

//============================================================================

struct MeshBundle {
	MeshBundle() = default;
	~MeshBundle() = default;
	TriMesh<MeshType>               base;
	TriMesh<MeshType>               reference;
	TriMesh<MeshType>               mapped;
	TriMesh<MeshType>               subdivFit;
	std::vector<Vector3D<MeshType>> disp;
	std::vector<SubdivLoD>          LoD;
	std::vector<uint64_t>           midpoint2endpoints;
	std::vector<uint32_t>           faceIndex;
};

//============================================================================

class AdjInfo {
public:
  AdjInfo()  = default;
  ~AdjInfo() = default;

  void clear() {
    BeginEnd.clear();
    NeighborCount.clear();
    Neighbors.clear();
  }

  uint32_t size() const { return NeighborCount.size(); }

  void resize(uint32_t const sz) {
    BeginEnd.resize(sz + 1, 0);
    NeighborCount.resize(sz, 0);
  }

  void reserve(uint32_t const sz) {
    BeginEnd.reserve(sz + 1);
    NeighborCount.reserve(sz);
    Neighbors.reserve(6 * sz);
  }

  void computeBeginEnd() {
    assert((size() + 1) == BeginEnd.size() && BeginEnd.front() == 0);

    for (uint32_t i = 0; i < size(); i++) { BeginEnd[i + 1] = BeginEnd[i] + NeighborCount[i]; }
    Neighbors.resize(BeginEnd.back());
  }

  void addNeighbor(uint32_t const i, uint32_t const neighbor) {
    Neighbors[BeginEnd[i] + NeighborCount[i]++] = neighbor;
  }

  std::vector<uint32_t> getNeighbors(uint32_t const i) const {
    auto const& begin = BeginEnd[i];
    auto const& end   = BeginEnd[i + 1];
    return std::vector<uint32_t>(Neighbors.begin() + begin, Neighbors.begin() + end);
  }

  std::vector<uint32_t> BeginEnd;
  std::vector<uint32_t> NeighborCount;
  std::vector<uint32_t> Neighbors;
};

//============================================================================
//
//struct Material {
//  std::string      name             = "material000";
//  std::string      texture          = "texture";
//  Vector3D<double> ambiant          = Vector3D<double>(1.0);
//  Vector3D<double> diffuse          = Vector3D<double>(1.0);
//  Vector3D<double> specular         = Vector3D<double>(1.0);
//  double           specularExponent = 0.1;
//  double           transparency     = 1.0;
//  uint32_t         illumination     = 2;
//
//  bool write(const std::string&        mtlPath,
//             std::vector<std::string>& texPath,
//             uint32_t                  mtlCount          = 1,
//             const bool                mtlDefaultSetting = true) {
//    if (mtlPath.empty()) return false;
//    std::ofstream ofs(mtlPath);
//    if (ofs.is_open()) {
//      for (uint32_t i = 0; i < mtlCount; ++i) {
//        ofs << "newmtl " << name << i << "\n";
//        if (mtlDefaultSetting) {
//          ofs << "Ka " << ambiant << "\n";
//          ofs << "Kd " << diffuse << "\n";
//          ofs << "Ks " << specular << "\n";
//          ofs << "Tr " << transparency << "\n";
//          ofs << "illum " << illumination << "\n";
//          ofs << "Ns " << specularExponent << "\n";
//        }
//        /*if (!texPath.empty()) {
//          auto tp = fillInNum(texPath, i);
//          texture = std::filesystem::path(tp).filename().string();
//        } */
//        ofs << "map_Kd " << texPath[i] << "\n\n";
//      }
//      ofs.close();
//      return true;
//    }
//    return false;
//  }
//};

//============================================================================

template<typename T, typename S>
bool
compareData(T const& data1, S const& data2) {
  if (data1.size() != data2.size()) return false;
  if (data1.size() == 0) return true;
  if (data1[0].size() != data2[0].size()) return false;
  if (data1[0].size() == 0) return true;
  auto t      = data1[0][0];
  using type1 = decltype(t);
  auto ptr1   = reinterpret_cast<type1 const*>(data1.data());
  auto ptr2   = reinterpret_cast<type1 const*>(data2.data());
  return std::equal(ptr1, ptr1 + data1.size() * data1[0].size(), ptr2);
}

//============================================================================
void subdivBaseMesh(MeshBundle&            MB,
					TriMesh<MeshType>&     subdivBase,
					uint32_t const         iterCount);

void computeVertex2triangle(uint32_t const                         vertexCount,
                            std::vector<Vector3D<uint32_t>> const& triangles,
                            AdjInfo&                               vertex2triangle);

void
tagAdjTriangles(uint32_t const idx, int8_t const tag, AdjInfo const& vertex2triangle, std::vector<int8_t>& triTags);

void computeAdjVert(uint32_t const                         idx,
                    std::vector<Vector3D<uint32_t>> const& triangles,
                    AdjInfo const&                         vertex2triangle,
                    std::vector<uint8_t>&                  visitedVert,
                    std::vector<uint32_t>&                 adjVert);

void computeAdjVert(uint32_t const                         idx,
                    std::vector<Vector3D<uint32_t>> const& triangles,
                    AdjInfo const&                         vertex2triangle,
                    std::vector<uint32_t>&                 adjVert);

void computeTriangleAdjTriangles(uint32_t const                         idx,
                                 std::vector<Vector3D<uint32_t>> const& triangles,
                                 AdjInfo const&                         vertex2triangle,
                                 std::vector<uint32_t>&                 adjTri);

uint32_t computeBoundaryVert(std::vector<Vector3D<uint32_t>> const& triangles,
                             AdjInfo const&                         vertex2triangle,
                             std::vector<uint8_t>&                  isBoundaryVertex);

//============================================================================
uint32_t floatToFixed32(float value, int fractionalBits);

MeshType fixed32ToFloat(uint32_t fixedValue, int fractionalBits);

//============================================================================
template<typename T>
bool
quantize(std::vector<T>&  data,
         const uint32_t   bitDepth,
         typename T::type& scale,
         T&               bboxMin,
         T&               bboxMax,
         bool             calScale = true) {
  typedef typename T::type type;
  if (typeid(type) != typeid(float) && typeid(type) != typeid(double)) {
    printf("Type of scale and bbox should be either float or double\n");
    return false;
  }

  if (calScale) {
    auto bbox = boundingBox(data);
    bboxMin   = bbox[0];
    bboxMax   = bbox[1];

    scale = ((1 << bitDepth) - 1) / (bboxMax - bboxMin).Linfnorm();
  }

  for (auto& vec : data) {
    T v(vec);
    v   = scale * (v - bboxMin);
    vec = v.round();
  }

  return true;
}

//============================================================================

template<typename T>
bool
dequantize(std::vector<T>& data, typename T::type scale, T& bboxMin) {
  typedef typename T::type type;
  if (typeid(type) != typeid(float) && typeid(type) != typeid(double)) {
    printf("Type of scale and bbox should be either float or double\n");
    return false;
  }

  type iscale = type(1.0) / scale;

  for (auto& vec : data) {
    T v(vec);
    vec = iscale * v + bboxMin;
  }

  return true;
}

//============================================================================

template<typename T>
T
computeTriangleNormalFitness(uint32_t const                         idx,
                             std::vector<Vector3D<uint32_t>> const& triangles,
                             std::vector<Vector3D<T>> const&        triangleNormals,
                             AdjInfo const&                         vertex2triangle) {
  std::vector<uint32_t> adjTri;
  computeTriangleAdjTriangles(idx, triangles, vertex2triangle, adjTri);
  T fitness = T(0);
  for (auto const& i : adjTri) { fitness += triangleNormals[idx] * triangleNormals[i]; }
  return fitness;
}

//============================================================================

template<typename T>
void
removeDuplicateVert(std::vector<Vector3D<T>> const&        vertices0,
                    std::vector<Vector3D<uint32_t>> const& triangles0,
                    std::vector<Vector3D<T>>&              vertices1,
                    std::vector<Vector3D<uint32_t>>&       triangles1,
                    std::vector<uint32_t>&                 uniqueOriginalIdx,
                    std::vector<uint32_t>&                 original2unique) {
  auto const vertexCount0   = vertices0.size();
  auto const triangleCount0 = triangles0.size();

  vertices1.clear();
  vertices1.reserve(vertexCount0);
  uniqueOriginalIdx.clear();
  uniqueOriginalIdx.reserve(vertexCount0);
  original2unique.resize(vertexCount0);

  std::map<Vector3D<T>, uint32_t> vertex2uniqueIdx;
  for (uint32_t idx0 = 0; idx0 < vertexCount0; idx0++) {
    auto const& vertex = vertices0[idx0];
    auto const  it     = vertex2uniqueIdx.find(vertex);
    if (it == vertex2uniqueIdx.end()) {
      auto const& uniqueIdx    = vertices1.size();
      vertex2uniqueIdx[vertex] = uniqueIdx;
      uniqueOriginalIdx.push_back(idx0);
      original2unique[idx0] = uniqueIdx;
      vertices1.push_back(vertex);
    } else {
      original2unique[idx0] = it->second;
    }
  }

  triangles1.resize(triangleCount0);
  for (uint32_t i = 0; i < triangleCount0; i++) {
    auto const& tri0 = triangles0[i];
    auto&       tri1 = triangles1[i];
    for (uint32_t j = 0; j < 3; j++) { tri1[j] = original2unique[tri0[j]]; }
  }
}

//============================================================================

template<typename T>
void
removeDuplicateTriangles(std::vector<Vector3D<uint32_t>> const& triangles0,
                         std::vector<Vector3D<T>> const&        triangleNormals,
                         AdjInfo const&                         vertex2triangle,
                         std::vector<Vector3D<uint32_t>>&       triangles1,
                         std::vector<uint32_t>&                 uniqueOriginalIdx,
                         std::vector<uint32_t>&                 original2unique,
                         bool const                             removeDegenerated = true) {
  auto const triangleCount0 = triangles0.size();
  triangles1.clear();
  triangles1.reserve(triangleCount0);
  uniqueOriginalIdx.clear();
  uniqueOriginalIdx.reserve(triangleCount0);
  original2unique.resize(triangleCount0);

  std::map<Vector3D<uint32_t>, uint32_t> triangle2uniqueIdx;
  std::map<uint32_t, T>                  triangleNormalFitness;
  for (uint32_t idx0 = 0; idx0 < triangleCount0; idx0++) {
    auto tri = triangles0[idx0];
    std::sort(tri.data(), tri.data() + 3);
    auto const it = triangle2uniqueIdx.find(tri);
    if (it == triangle2uniqueIdx.end()) {
      original2unique[idx0] = std::numeric_limits<uint32_t>::max();
      if (!removeDegenerated || !isDegenerate(tri)) {
        auto const& uniqueIdx   = triangles1.size();
        triangle2uniqueIdx[tri] = uniqueIdx;
        uniqueOriginalIdx.push_back(idx0);
        original2unique[idx0] = uniqueIdx;
        triangles1.push_back(triangles0[idx0]);
      }
    } else {
      auto const& uniqueIdx = it->second;
      if (triangleNormalFitness.find(uniqueIdx) == triangleNormalFitness.end()) {
        auto const i =
          std::distance(original2unique.begin(), std::find(original2unique.begin(), original2unique.end(), uniqueIdx));
        triangleNormalFitness[uniqueIdx] =
          computeTriangleNormalFitness(i, triangles0, triangleNormals, vertex2triangle);
      }
      auto const fitness = computeTriangleNormalFitness(idx0, triangles0, triangleNormals, vertex2triangle);
      if (fitness > triangleNormalFitness[uniqueIdx]) {
        triangleNormalFitness[uniqueIdx] = fitness;
        triangles1[uniqueIdx]            = triangles0[idx0];
        uniqueOriginalIdx[uniqueIdx]     = idx0;
      }
      original2unique[idx0] = uniqueIdx;
    }
  }
}

//============================================================================

template<typename T>
void
removeDegenerateTriangles(TriMesh<T>& mesh) {
  auto const faceCount     = mesh.FaceXYZ.size();
  auto const hasFaceUV     = mesh.FaceUV.size() == faceCount;
  auto const hasFaceNormal = mesh.FaceNormal.size() == faceCount;

  std::vector<Vector3D<uint32_t>> FaceXYZ;
  std::vector<Vector3D<uint32_t>> FaceUV;
  std::vector<Vector3D<uint32_t>> FaceNormal;

  FaceXYZ.reserve(faceCount);
  if (hasFaceUV) FaceUV.reserve(faceCount);
  if (hasFaceNormal) FaceNormal.reserve(faceCount);

  for (uint32_t i = 0; i < faceCount; i++) {
    auto const& face = mesh.FaceXYZ[i];
    if (!isDegenerate(face)) {
      FaceXYZ.push_back(face);
      if (hasFaceUV) { FaceUV.push_back(mesh.FaceUV[i]); }
      if (hasFaceNormal) { FaceNormal.push_back(mesh.FaceNormal[i]); }
    }
  }

  std::swap(mesh.FaceXYZ, FaceXYZ);
  if (hasFaceUV) std::swap(mesh.FaceUV, FaceUV);
  if (hasFaceNormal) std::swap(mesh.FaceNormal, FaceNormal);
}

//============================================================================

template<typename T>
void
subdivMidpoint1(std::vector<T>&                  vertices,
                std::vector<Vector3D<uint32_t>>& triangles,
                std::vector<uint64_t>*           midpoint2endpoints = nullptr,
                std::vector<uint32_t>*           faceIndex          = nullptr) {
  auto const vertexCount0    = vertices.size();
  auto const triangleCount0  = triangles.size();
  auto const newVertCountEst = 2 * (vertexCount0 + triangleCount0);

  vertices.reserve(vertexCount0 + newVertCountEst);
  if (midpoint2endpoints) {
    midpoint2endpoints->reserve(vertexCount0 + newVertCountEst);
    midpoint2endpoints->resize(vertexCount0);
  }
  if (faceIndex) { faceIndex->resize(4 * triangleCount0); }

  AdjInfo vertex2triangle;
  computeVertex2triangle(vertexCount0, triangles, vertex2triangle);
  std::vector<uint8_t> visitedVert(vertexCount0);

  uint32_t                               newVertIdx = vertexCount0;
  std::unordered_map<uint64_t, uint32_t> endpoints2midpoint;
  auto const                             half = vertices[0].cast(0.5);
  for (uint32_t vidx0 = 0; vidx0 < vertexCount0; vidx0++) {
    std::vector<uint32_t> adjVert;
    computeAdjVert(vidx0, triangles, vertex2triangle, visitedVert, adjVert);
    for (auto const& vidx1 : adjVert) {
      if (vidx1 < vidx0) continue;
      uint64_t const endpoints      = (uint64_t(vidx0) << 32) + vidx1;
      endpoints2midpoint[endpoints] = newVertIdx;
      vertices.push_back(half * (vertices[vidx0] + vertices[vidx1]));
      if (midpoint2endpoints) midpoint2endpoints->push_back(endpoints);
      newVertIdx++;
    }
  }

  triangles.resize(4 * triangleCount0);
  uint32_t newTriIdx = triangleCount0;
  for (uint32_t tidx0 = 0; tidx0 < triangleCount0; tidx0++) {
    auto const         tri0 = triangles[tidx0];
    Vector3D<uint32_t> tri;
    for (uint32_t i = 0; i < 3; i++) {
	  auto const     vidx0     = tri0[i];
	  auto const     vidx1     = tri0[(i + 1) % 3];
      uint64_t const endpoints = (uint64_t((std::min)(vidx0, vidx1)) << 32) + (std::max)(vidx0, vidx1);
      tri[i]                   = endpoints2midpoint[endpoints];
    }
    triangles[tidx0]         = tri;
    triangles[newTriIdx + 0] = Vector3D<uint32_t>(tri0[0], tri[0], tri[2]);
    triangles[newTriIdx + 1] = Vector3D<uint32_t>(tri[0], tri0[1], tri[1]);
    triangles[newTriIdx + 2] = Vector3D<uint32_t>(tri[2], tri[1], tri0[2]);
    if (faceIndex) {
      faceIndex->at(newTriIdx + 0) = faceIndex->at(tidx0);
      faceIndex->at(newTriIdx + 1) = faceIndex->at(tidx0);
      faceIndex->at(newTriIdx + 2) = faceIndex->at(tidx0);
    }
    newTriIdx += 3;
  }
}

//============================================================================
template<typename T>
void
TriMesh<T>::subdivMidpoint(uint32_t const          iterCount,
                           std::vector<SubdivLoD>* lod,
                           std::vector<uint64_t>*  midpoint2endpoints3D,
                           std::vector<uint64_t>*  midpoint2endpoints2D,
                           std::vector<uint32_t>*  faceIndex) {
  //printf("Subdivide Mid Point iterationCount = %d \n", iterCount);
  //fflush(stdout);

  if (lod) {
    lod->resize(iterCount + 1);
    lod->at(0).xyzCount     = XYZ.size();
    lod->at(0).uvCount      = UV.size();
    lod->at(0).faceXYZCount = FaceXYZ.size();
    lod->at(0).faceUVCount  = FaceUV.size();
  }
  uint32_t numBaseFace = FaceXYZ.size();
  //printf("numBaseFace = %u \n", numBaseFace);
  if (faceIndex) {
    auto numFacesForLod = std::pow(4, iterCount);
    faceIndex->reserve(FaceXYZ.size() * numFacesForLod);
    for (uint32_t i = 0; i < numBaseFace; i++) { faceIndex->push_back(i); }
  }

  //printf(
  //  "XYZ  = %8zu FaceXYZ = %8zu UV  = %8zu FaceXYZ = %8zu \n", XYZ.size(), FaceXYZ.size(), UV.size(), FaceUV.size());
  for (uint32_t it = 0; it < iterCount; it++) {
    subdivMidpoint1<Vector3D<T>>(XYZ, FaceXYZ, midpoint2endpoints3D, faceIndex);
    if (UV.size() && FaceUV.size()) { subdivMidpoint1<Vector2D<T>>(UV, FaceUV, midpoint2endpoints2D); }
    if (lod) {
      lod->at(it + 1).xyzCount     = XYZ.size();
      lod->at(it + 1).uvCount      = UV.size();
      lod->at(it + 1).faceXYZCount = FaceXYZ.size();
      lod->at(it + 1).faceUVCount  = FaceUV.size();
    }
  }
}

//============================================================================

template<typename T>
void
normalizeIfNeeded(Vector3D<T>& vec, std::true_type) {
  vec.normalize();
}

//============================================================================

template<typename T>
void
normalizeIfNeeded(Vector3D<T>& vec, std::false_type) {
  // Do nothing
}

//============================================================================

template<typename T, bool Normalize>
void
subdivInterpolate(std::vector<Vector3D<T>>&     signal,
                  std::vector<SubdivLoD> const& lod,
                  std::vector<uint64_t> const&  midpoint2endpoints) {
  const auto lodCount = lod.size() - 1;
  for (size_t it = 0; it < lodCount; ++it) {
    const auto xyzCount0 = lod[it].xyzCount;
    const auto xyzCount1 = lod[it + 1].xyzCount;
    for (size_t i = xyzCount0; i < xyzCount1; ++i) {
      auto const endpoint = midpoint2endpoints[i];
      auto const vidx0    = endpoint >> 32;
      auto const vidx1    = endpoint & 0xFFFFFFFF;
      signal[i]           = kSubdivWeights[0] * signal[vidx0] + kSubdivWeights[1] * signal[vidx1];
      normalizeIfNeeded(signal[i], std::integral_constant<bool, Normalize>());
    }
  }
}

//============================================================================

template<typename T>
uint32_t
extractConnectedComponents(TriMesh<T> const& mesh, std::vector<std::shared_ptr<TriMesh<T>>>* ccMeshes = nullptr) {
  if (ccMeshes) { ccMeshes->clear(); }
  if (mesh.FaceXYZ.empty()) { return 0; }

  uint32_t ccCount = 0;
  const auto pointCount    = mesh.XYZ.size();
  const auto texCoordCount = mesh.UV.size();
  const auto normalCount   = mesh.Normal.size();
  const auto triangleCount = mesh.FaceXYZ.size();

  std::vector<int32_t> partition;
  partition.assign(triangleCount, -1);

  AdjInfo vert2Face;
  computeVertex2triangle(pointCount, mesh.FaceXYZ, vert2Face);

  // clear connectedComponents
  std::vector<int32_t> lifo;
  lifo.reserve(triangleCount);
  std::vector<int32_t> posMapping;
  std::vector<int32_t> texCoordMapping;
  std::vector<int32_t> normalMapping;
  std::vector<int32_t> triangleList;

  if (ccMeshes) {
    posMapping.resize(pointCount, -1);
    if (mesh.UV.size()) { texCoordMapping.resize(texCoordCount, -1); }
    if (mesh.Normal.size()) { normalMapping.resize(normalCount, -1); }
  }

  for (uint32_t faceIdx = 0; faceIdx < triangleCount; ++faceIdx) {
    if (partition[faceIdx] == -1) {
      const auto ccIdx   = ccCount++;
      partition[faceIdx] = ccIdx;
      lifo.push_back(faceIdx);

      std::shared_ptr<TriMesh<T>> ccMesh;
      if (ccMeshes) {
        ccMesh.reset(new TriMesh<T>());
        ccMeshes->push_back(ccMesh);
        ccMesh->Mtllib = mesh.Mtllib;
        triangleList.resize(0);
      }

      int32_t ccPosCount      = 0;
      int32_t ccTexCoordCount = 0;
      int32_t ccNormalCount   = 0;

      while (!lifo.empty()) {
        const auto tIdx = lifo.back();
        lifo.pop_back();

        if (ccMeshes) {
          // add triangle to ccMesh
          triangleList.push_back(tIdx);
          const auto& face = mesh.FaceXYZ[tIdx];
          for (size_t k = 0; k < 3; ++k) {
            const auto v = face[k];
            if (posMapping[v] == -1) {
              posMapping[v] = ccPosCount++;
              ccMesh->XYZ.push_back(mesh.XYZ[v]);
            }
          }
          ccMesh->FaceXYZ.push_back(Vector3D<uint32_t>(posMapping[face[0]], posMapping[face[1]], posMapping[face[2]]));

          // add texCoord to ccMesh
          if (mesh.FaceUV.size()) {
            const auto& faceUV = mesh.FaceUV[tIdx];
            for (size_t k = 0; k < 3; ++k) {
              const auto v = faceUV[k];
              if (texCoordMapping[v] == -1) {
                texCoordMapping[v] = ccTexCoordCount++;
                ccMesh->UV.push_back(mesh.UV[v]);
              }
            }
            ccMesh->FaceUV.push_back(
              Vector3D<uint32_t>(texCoordMapping[faceUV[0]], texCoordMapping[faceUV[1]], texCoordMapping[faceUV[2]]));
          }

          // add normal to ccMesh
          if (mesh.FaceNormal.size()) {
            const auto& faceNormal = mesh.FaceNormal[tIdx];
            for (size_t k = 0; k < 3; ++k) {
              const auto v = faceNormal[k];
              if (normalMapping[v] == -1) {
                normalMapping[v] = ccNormalCount++;
                ccMesh->Normal.push_back(mesh.Normal[v]);
              }
            }
            ccMesh->FaceNormal.push_back(Vector3D<uint32_t>(
              normalMapping[faceNormal[0]], normalMapping[faceNormal[1]], normalMapping[faceNormal[2]]));
          }
        }

        const auto& face = mesh.FaceXYZ[tIdx];
        for (size_t k = 0; k < 3; ++k) {
          const auto v      = face[k];
          const auto istart = vert2Face.BeginEnd[v];
          const auto iend   = vert2Face.BeginEnd[v] + vert2Face.NeighborCount[v];

          for (size_t k = istart; k < iend; ++k) {
            const auto nIdx = vert2Face.Neighbors[k];
            if (partition[nIdx] == -1) {
              partition[nIdx] = ccIdx;
              lifo.push_back(nIdx);
            }
          }
        }
      }

      if (ccMeshes) {
        for (const auto& tIdx : triangleList) {
          const auto& face = mesh.FaceXYZ[tIdx];
          for (size_t k = 0; k < 3; ++k) {
            const auto v  = face[k];
            posMapping[v] = -1;
          }
          if (mesh.FaceUV.size()) {
            const auto& faceUV = mesh.FaceUV[tIdx];
            for (size_t k = 0; k < 3; ++k) {
              const auto v       = faceUV[k];
              texCoordMapping[v] = -1;
            }
          }
          if (mesh.FaceNormal.size()) {
            const auto& faceNormal = mesh.FaceNormal[tIdx];
            for (size_t k = 0; k < 3; ++k) {
              const auto v     = faceNormal[k];
              normalMapping[v] = -1;
            }
          }
        }
      }
    }
  }
  return ccCount;
}

