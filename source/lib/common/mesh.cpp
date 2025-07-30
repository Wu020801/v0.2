/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#include "mesh.hpp"
#include "kdtree.hpp"

//============================================================================

static std::vector<std::string>
splitString(const std::string& str, const std::string& delimiter = " ") {
  std::vector<std::string> res;
  size_t                   curr = 0, next = 0;
  while ((next = str.find(delimiter, curr)) != std::string::npos) {
    if (delimiter != " " || next > curr) res.push_back(str.substr(curr, next - curr));
    curr = next + 1;
  }
  if (curr < str.length()) res.push_back(str.substr(curr));
  return res;
}

//============================================================================

template<typename T>
bool
TriMesh<T>::readOBJ(const std::string& filePath) {
  std::ifstream ifs(filePath);
  if (ifs.is_open()) {
    clear();
    std::string line;
    while (std::getline(ifs, line)) {
      auto strs = splitString(line);
      if (strs.size() == 0) { continue; }
      if (strs[0] == "v" && strs.size() > 3) {
        Vector3D<T> v(std::stod(strs[1]), std::stod(strs[2]), std::stod(strs[3]));
        XYZ.push_back(v);
      } else if (strs[0] == "vt" && strs.size() > 2) {
        Vector2D<T> v(std::stod(strs[1]), std::stod(strs[2]));
        UV.push_back(v);
      } else if (strs[0] == "vn" && strs.size() > 3) {
        Vector3D<T> v(std::stod(strs[1]), std::stod(strs[2]), std::stod(strs[3]));
        Normal.push_back(v);
      } else if (strs[0] == "f" && strs.size() > 3) {
        std::vector<std::vector<std::string>> idx;
        for (uint32_t i = 1; i < 4; i++) { idx.push_back(splitString(strs[i], "/")); }
        uint32_t n = idx[0].size();
        if (n && idx[0][0].size()) {
          Vector3D<uint32_t> f0(std::stoi(idx[0][0]) - 1, std::stoi(idx[1][0]) - 1, std::stoi(idx[2][0]) - 1);
          FaceXYZ.push_back(f0);
        }
        if (n > 1 && idx[0][1].size()) {
          Vector3D<uint32_t> f1(std::stoi(idx[0][1]) - 1, std::stoi(idx[1][1]) - 1, std::stoi(idx[2][1]) - 1);
          FaceUV.push_back(f1);
        }
        if (n > 2 && idx[0][2].size()) {
          Vector3D<uint32_t> f2(std::stoi(idx[0][2]) - 1, std::stoi(idx[1][2]) - 1, std::stoi(idx[2][2]) - 1);
          FaceNormal.push_back(f2);
        }
      } else if (strs[0] == "mtllib" && strs.size() > 1) {
        Mtllib = strs[1];
      }
    }  // line loop
    ifs.close();
    return true;
  }
  printf("Error: cannot read %s \n", filePath.c_str());
  return false;
}

//============================================================================

template<typename T>
bool
TriMesh<T>::writeOBJ(const std::string& filePath, const std::vector<uint32_t>& mtlPos, bool isInt) const {
  if (filePath.empty()) return false;
  if (!FaceUV.empty() && FaceUV.size() != FaceXYZ.size()) {
    printf("Cannot write obj file: FaceUV.size() != FaceXYZ.size()\n");
    return false;
  }
  if (!FaceNormal.empty() && FaceNormal.size() != FaceXYZ.size()) {
    printf("Cannot write obj file: Normal.size() != XYZ.size()\n");
    return false;
  }

  std::ofstream ofs(filePath);
  if (ofs.is_open()) {
    bool const multiTexture = !mtlPos.empty();
    if (!Mtllib.empty()) ofs << "mtllib " << Mtllib << "\n";

    if (RGB.size() == XYZ.size()) {
      for (uint32_t i = 0; i < XYZ.size(); i++) {
        auto v = XYZ[i];
        auto c = RGB[i];
        if (isInt) {
		  Vector3D<int64_t> v_int;
		  Vector3D<int64_t> c_int;
          for (uint32_t j = 0; j < v.size(); ++j) {
			v_int[j] = static_cast<int>(v[j]);
            c_int[j] = static_cast<int>(c[j]);
          }
		  ofs << "v " << v_int[0] << " " << v_int[1] << " " << v_int[2] << " " << c_int[0] << " " << c_int[1] << " " << c_int[2] << "\n";
		}
		else {
			ofs << "v " << v[0] << " " << v[1] << " " << v[2] << " " << c[0] << " " << c[1] << " " << c[2] << "\n";
		}
      }
    } else {
      for (uint32_t i = 0; i < XYZ.size(); i++) {
        auto v = XYZ[i];
		if (isInt) {
			Vector3D<int64_t> v_int;
			for (uint32_t j = 0; j < v.size(); ++j) v_int[j] = static_cast<int64_t>(v[j]);
			ofs << "v " << v_int[0] << " " << v_int[1] << " " << v_int[2] << "\n";
		}
		else {
			ofs << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
		}     
      }
    }

    for (uint32_t i = 0; i < UV.size(); i++) {
      auto v = UV[i];
	  if (isInt) {
		  Vector2D<int64_t> v_int;
		  for (uint32_t j = 0; j < v.size(); ++j) v_int[j] = static_cast<int64_t>(v[j]);
		  ofs << "vt " << v_int[0] << " " << v_int[1] << "\n";
	  }
	  else {
		  ofs << "vt " << v[0] << " " << v[1] << "\n";
	  } 
    }

    for (uint32_t i = 0; i < Normal.size(); i++) {
      auto v = Normal[i];
	  if (isInt) {
		  Vector3D<int64_t> v_int;
		  for (uint32_t j = 0; j < v.size(); ++j) v_int[j] = static_cast<int>(v[j]);
		  ofs << "v " << v_int[0] << " " << v_int[1] << " " << v_int[2] << "\n";
	  }
	  else {
		  ofs << "vn " << v[0] << " " << v[1] << " " << v[2] << "\n";
	  } 
    }

    if (FaceUV.empty() && FaceNormal.empty()) {
      for (uint32_t m = 0, f = 0; f < FaceXYZ.size(); f++) {
        if (multiTexture && f == mtlPos[m]) { ofs << "usemtl material000" << m++ << "\n"; }

        const auto& face = FaceXYZ[f];
        ofs << "f " << face[0] + 1 << " " << face[1] + 1 << " " << face[2] + 1 << "\n";
      }
    } else if (!FaceUV.empty() && FaceNormal.empty()) {
      for (uint32_t m = 0, f = 0; f < FaceXYZ.size(); f++) {
        if (multiTexture && f == mtlPos[m]) { ofs << "usemtl material000" << m++ << "\n"; }

        const auto& face0 = FaceXYZ[f];
        const auto& face1 = FaceUV[f];
        ofs << "f " << face0[0] + 1 << "/" << face1[0] + 1 << " " << face0[1] + 1 << "/" << face1[1] + 1 << " "
            << face0[2] + 1 << "/" << face1[2] + 1 << "\n";
      }
    } else if (FaceUV.empty() && !FaceNormal.empty()) {
      for (uint32_t m = 0, f = 0; f < FaceXYZ.size(); f++) {
        if (multiTexture && f == mtlPos[m]) { ofs << "usemtl material000" << m++ << "\n"; }

        const auto& face0 = FaceXYZ[f];
        const auto& face1 = FaceNormal[f];
        ofs << "f " << face0[0] + 1 << "//" << face1[0] + 1 << " " << face0[1] + 1 << "//" << face1[1] + 1 << " "
            << face0[2] + 1 << "//" << face1[2] + 1 << "\n";
      }
    } else {
      for (uint32_t m = 0, f = 0; f < FaceXYZ.size(); f++) {
        if (multiTexture && f == mtlPos[m]) { ofs << "usemtl material000" << m++ << "\n"; }

        const auto& face0 = FaceXYZ[f];
        const auto& face1 = FaceUV[f];
        const auto& face2 = FaceNormal[f];
        ofs << "f " << face0[0] + 1 << "/" << face1[0] + 1 << "/" << face2[0] + 1 << " " << face0[1] + 1 << "/"
            << face1[1] + 1 << "/" << face2[1] + 1 << " " << face0[2] + 1 << "/" << face1[2] + 1 << "/" << face2[2] + 1
            << "\n";
      }
    }
    ofs.close();
    return true;
  }
  printf("Error: cannot write %s \n", filePath.c_str());
  return false;
}

//============================================================================

template<typename T>
bool
TriMesh<T>::read(const std::string& filePath) {
  auto ext = std::filesystem::path(filePath).extension().string();
  if (ext == ".obj") { return readOBJ(filePath); }
  printf("Error: reading files with format %s is not supported.\n", ext.c_str());
  return false;
}

//============================================================================

template<typename T>
bool
TriMesh<T>::write(const std::string& filePath) const {
  if (filePath.empty()) return false;
  auto ext = std::filesystem::path(filePath).extension().string();
  if (ext == ".obj") { return writeOBJ(filePath); }
  printf("Error: writing files with format %s is not supported. %s \n", ext.c_str(), filePath.c_str());
  return false;
}

//============================================================================

template<typename T>
void
TriMesh<T>::clear() {
  XYZ.clear();
  UV.clear();
  Normal.clear();
  FaceXYZ.clear();
  FaceUV.clear();
  FaceNormal.clear();
}

//============================================================================

template<typename T>
void
TriMesh<T>::print() const {
  printf("XYZ count:        %zu\n", XYZ.size());
  printf("UV count:         %zu\n", UV.size());
  printf("Normal count:     %zu\n", Normal.size());
  printf("FaceXYZ count:    %zu\n", FaceXYZ.size());
  printf("FaceUV count:     %zu\n", FaceUV.size());
  printf("FaceNormal count: %zu\n", FaceNormal.size());
  printf("Mtllib:           %s\n", Mtllib.c_str());
}

//============================================================================

template<typename T>
void
TriMesh<T>::reconGeometry(std::vector<Vector3D<T>>& disp) {
	for (uint32_t i = 0; i < XYZ.size(); ++i) {
		//std::cout << "xyz " << XYZ[i][0] << " "<< XYZ[i][1] << " " << XYZ[i][2] << std::endl;
		XYZ[i] += disp[i];
		//std::cout << "recon " << XYZ[i][0] << " " << XYZ[i][1] << " " << XYZ[i][2] << std::endl;
	}
}

//============================================================================

template<typename T>
TriMesh<T>&
TriMesh<T>::setFromFaces(TriMesh<T> const& mesh, const size_t firstFace, const size_t lastFace) {
  Mtllib     = mesh.Mtllib;
  FaceNormal = mesh.FaceNormal;
  if (mesh.FaceXYZ.size() != mesh.FaceUV.size()) {
    printf("setFromFaces: number of faces are not equal between "
           "FaceXYZ and FaceUV (%zu != %zu) \n",
           mesh.FaceXYZ.size(),
           mesh.FaceUV.size());
    exit(-1);
  }
  FaceXYZ.resize(mesh.FaceXYZ.size());
  for (uint32_t i = firstFace; i < lastFace; i++) FaceXYZ[i] = mesh.FaceXYZ[i];
  FaceUV.resize(mesh.FaceUV.size());
  for (uint32_t i = firstFace; i < lastFace; i++) FaceUV[i] = mesh.FaceUV[i];
  XYZ.resize(mesh.XYZ.size());
  for (uint32_t i = 0; i < XYZ.size(); i++) XYZ[i] = mesh.XYZ[i];
  UV.resize(mesh.UV.size());
  for (uint32_t i = 0; i < UV.size(); i++) UV[i] = mesh.UV[i];
  Normal.resize(mesh.Normal.size());
  for (uint32_t i = 0; i < Normal.size(); i++) Normal[i] = mesh.Normal[i];
  return *this;
}

//============================================================================

template<typename T>
void
TriMesh<T>::removeDegeneratedVertices() {
  std::unordered_map<uint32_t, uint32_t> vertex_map;
  std::vector<Vector3D<T>>               new_XYZ;
  uint32_t                               newIndex = 0;
  for (uint32_t faceIndex = 0; faceIndex < FaceXYZ.size(); ++faceIndex) {
    auto& face = FaceXYZ[faceIndex];
    for (uint32_t i = 0; i < 3; ++i) {
      uint32_t oldIndex = face[i];
      auto     it       = vertex_map.find(oldIndex);
      if (it == vertex_map.end()) {
        vertex_map[oldIndex] = newIndex;
        new_XYZ.push_back(XYZ[oldIndex]);
        face[i] = newIndex;
        newIndex++;
      } else {
        face[i] = it->second;
      }
    }
  }
  XYZ = std::move(new_XYZ);
}

//============================================================================

template<typename T>
void
TriMesh<T>::computeNormals(bool const normalize, uint32_t const mode) {
  Normal.resize(XYZ.size());
  std::fill(Normal.begin(), Normal.end(), T(0));
  switch (mode) {
  case 0:  // weighted by face areas
    for (auto const& face : FaceXYZ) {
      auto const i = face[0];
      auto const j = face[1];
      auto const k = face[2];
      auto const n = computeTriangleNormal(XYZ[i], XYZ[j], XYZ[k], false);
      Normal[i] += n;
      Normal[j] += n;
      Normal[k] += n;
    }
    if (normalize) {
      for (auto& n : Normal) { n.normalize(); }
    }
    break;
  }
}

//============================================================================

template<typename T>
void
TriMesh<T>::computeTriangleNormals(std::vector<Vector3D<T>>& normals, bool const normalize) const {
  normals.resize(FaceXYZ.size());
  for (uint32_t i = 0; i < FaceXYZ.size(); i++) {
    auto const& face = FaceXYZ[i];
    normals[i]       = computeTriangleNormal(XYZ[face[0]], XYZ[face[1]], XYZ[face[2]], normalize);
  }
}

//============================================================================

template<typename T>
void
TriMesh<T>::computeUV2XYZ(std::vector<uint32_t>& uv2xyz) const {
  uv2xyz.resize(UV.size());
  for (uint32_t i = 0; i < FaceXYZ.size(); i++) {
    auto const& triXYZ = FaceXYZ[i];
    auto const& triUV  = FaceUV[i];
    for (uint32_t j = 0; j < 3; j++) { uv2xyz[triUV[j]] = triXYZ[j]; }
  }
}

//============================================================================

template<typename T>
void
TriMesh<T>::computeXYZ2UV(std::vector<uint32_t>& xyz2uv) const {
  xyz2uv.resize(XYZ.size());
  for (uint32_t i = 0; i < FaceXYZ.size(); i++) {
    auto const& triXYZ = FaceXYZ[i];
    auto const& triUV  = FaceUV[i];
    for (uint32_t j = 0; j < 3; j++) { xyz2uv[triXYZ[j]] = triUV[j]; }
  }
}

//============================================================================

template<typename T>
void
TriMesh<T>::invertOrientation() {
  for (uint32_t i = 0; i < FaceXYZ.size(); i++) {
    auto& face = FaceXYZ[i];
    std::swap(face[0], face[1]);
  }
  for (uint32_t i = 0; i < FaceUV.size(); i++) {
    auto& face = FaceUV[i];
    std::swap(face[0], face[1]);
  }
  for (uint32_t i = 0; i < FaceNormal.size(); i++) {
    auto& face = FaceNormal[i];
    std::swap(face[0], face[1]);
  }
}

//============================================================================

template<typename T>
void
TriMesh<T>::append(TriMesh const& mesh) {
  auto const& xyzCount0        = XYZ.size();
  auto const& uvCount0         = UV.size();
  auto const& normalCount0     = Normal.size();
  auto const& faceXYZCount0    = FaceXYZ.size();
  auto const& faceUVCount0     = FaceUV.size();
  auto const& faceNormalCount0 = FaceNormal.size();
  auto const& faceXYZCount1    = mesh.FaceXYZ.size();
  auto const& faceUVCount1     = mesh.FaceUV.size();
  auto const& faceNormalCount1 = mesh.FaceNormal.size();
 
  XYZ.insert(XYZ.end(), mesh.XYZ.begin(), mesh.XYZ.end());
  UV.insert(UV.end(), mesh.UV.begin(), mesh.UV.end());
  Normal.insert(Normal.end(), mesh.Normal.begin(), mesh.Normal.end());
  FaceXYZ.resize(faceXYZCount0 + faceXYZCount1);
  for (uint32_t i = 0; i < faceXYZCount1; i++) { FaceXYZ[faceXYZCount0 + i] = mesh.FaceXYZ[i] + xyzCount0; }
  FaceUV.resize(faceUVCount0 + faceUVCount1);
  for (uint32_t i = 0; i < faceUVCount1; i++) { FaceUV[faceUVCount0 + i] = mesh.FaceUV[i] + uvCount0; }
  FaceNormal.resize(faceNormalCount0 + faceNormalCount1);
  for (uint32_t i = 0; i < faceNormalCount1; i++) {
    FaceNormal[faceNormalCount0 + i] = mesh.FaceNormal[i] + normalCount0;
  }
}

//============================================================================

void
subdivBaseMesh(MeshBundle&            MB,
			   TriMesh<MeshType>&	  subdivBase,
			   uint32_t const         iterCount) {
	// subdivBase = MB.base;
	subdivBase.computeNormals();
	subdivBase.subdivMidpoint(iterCount, &(MB.LoD), &(MB.midpoint2endpoints), nullptr, &MB.faceIndex);

	subdivBase.Normal.resize(subdivBase.XYZ.size());
	subdivInterpolate<MeshType, true>(subdivBase.Normal, MB.LoD, MB.midpoint2endpoints);
}

//============================================================================

void
computeVertex2triangle(uint32_t const                         vertexCount,
                       std::vector<Vector3D<uint32_t>> const& triangles,
                       AdjInfo&                               vertex2triangle) {
  vertex2triangle.resize(vertexCount);
  for (uint32_t i = 0; i < triangles.size(); i++) {
    auto const& tri = triangles[i];
    if (isDegenerate(tri)) continue;
    for (uint32_t j = 0; j < 3; j++) { vertex2triangle.NeighborCount[tri[j]]++; }
  }
  vertex2triangle.computeBeginEnd();
  std::fill(vertex2triangle.NeighborCount.begin(), vertex2triangle.NeighborCount.end(), 0);
  for (uint32_t i = 0; i < triangles.size(); i++) {
    auto const& tri = triangles[i];
    if (isDegenerate(tri)) continue;
    for (uint32_t j = 0; j < 3; j++) { vertex2triangle.addNeighbor(tri[j], i); }
  }
}

//============================================================================

void
tagAdjTriangles(uint32_t const idx, int8_t const tag, AdjInfo const& vertex2triangle, std::vector<int8_t>& triTags) {
  auto const& faceIdx = vertex2triangle.getNeighbors(idx);
  for (auto const& i : faceIdx) { triTags[i] = tag; }
}

//============================================================================

void
computeAdjVert(uint32_t const                         idx,
               std::vector<Vector3D<uint32_t>> const& faces,
               AdjInfo const&                         vertex2triangle,
               std::vector<uint32_t>&                 adjVert) {
  std::set<uint32_t> adjVertSet;
  auto const&        faceIdx = vertex2triangle.getNeighbors(idx);
  for (auto const& fIdx : faceIdx) {
    auto const& face = faces[fIdx];
    for (uint32_t j = 0; j < 3; j++) {
      if (face[j] != idx) adjVertSet.insert(face[j]);
    }
  }
  adjVert.assign(adjVertSet.begin(), adjVertSet.end());
}

//============================================================================

void
computeAdjVert(uint32_t const                         idx,
               std::vector<Vector3D<uint32_t>> const& faces,
               AdjInfo const&                         vertex2triangle,
               std::vector<uint8_t>&                  visistedVert,
               std::vector<uint32_t>&                 adjVert) {
  auto const& faceIdx = vertex2triangle.getNeighbors(idx);
  for (auto const& fIdx : faceIdx) {
    const auto& face = faces[fIdx];
    for (size_t i = 0; i < 3; ++i) { visistedVert[face[i]] = uint8_t(1); }
  }
  visistedVert[idx] = uint8_t(0);

  adjVert.resize(0);
  for (auto const& fIdx : faceIdx) {
    auto const& face = faces[fIdx];
    for (size_t j = 0; j < 3; j++) {
      const auto index = face[j];
      if (visistedVert[index] != 0) {
        adjVert.push_back(index);
        visistedVert[index] = uint8_t(0);
      }
    }
  }
}

//============================================================================

int32_t
computeEdgeAdjTriangleCount(uint32_t const       idx1,
                            uint32_t const       idx2,
                            AdjInfo const&       vert2Triangle,
                            std::vector<int8_t>& visitedTriangle) {
  tagAdjTriangles(idx1, 0, vert2Triangle, visitedTriangle);
  tagAdjTriangles(idx2, 1, vert2Triangle, visitedTriangle);
  auto const& faceIdx = vert2Triangle.getNeighbors(idx1);
  int32_t     count   = 0;

  for (const auto& fidx : faceIdx) { count += static_cast<int32_t>(visitedTriangle[fidx] != 0); }
  return count;
}

//============================================================================

void
computeTriangleAdjTriangles(uint32_t const                         idx,
                            std::vector<Vector3D<uint32_t>> const& triangles,
                            AdjInfo const&                         vertex2triangle,
                            std::vector<uint32_t>&                 adjTri) {
  auto const&        tri = triangles[idx];
  std::set<uint32_t> adjTriSet;
  for (uint32_t i = 0; i < 3; i++) {
    auto const& faceIdx = vertex2triangle.getNeighbors(tri[i]);
    for (auto const& j : faceIdx) {
      if (j != idx) adjTriSet.insert(j);
    }
  }
  adjTri.assign(adjTriSet.begin(), adjTriSet.end());
}

//============================================================================

uint32_t
computeBoundaryVert(std::vector<Vector3D<uint32_t>> const& triangles,
                    AdjInfo const&                         vertex2triangle,
                    std::vector<uint8_t>&                  isBoundaryVertex) {
  auto const vertexCount = vertex2triangle.size();
  if (vertexCount == 0) return 0;

  uint32_t boundaryVertexCount = 0;
  isBoundaryVertex.resize(vertexCount);
  std::fill(isBoundaryVertex.begin(), isBoundaryVertex.end(), 0);
  std::vector<int8_t>   visitedTriangle(triangles.size());
  std::vector<uint32_t> adjVert;

  for (uint32_t i = 0; i < vertexCount; i++) {
    computeAdjVert(i, triangles, vertex2triangle, adjVert);
    for (auto const j : adjVert) {
      if (j < i) continue;
      std::vector<uint32_t> adjTri;
      auto                  count = computeEdgeAdjTriangleCount(i, j, vertex2triangle, visitedTriangle);
      if (count == 1) {
        if (isBoundaryVertex[i] == 0) {
          isBoundaryVertex[i] = 1;
          boundaryVertexCount++;
        }
        if (isBoundaryVertex[j] == 0) {
          isBoundaryVertex[j] = 1;
          boundaryVertexCount++;
        }
      }
    }
  }
  return boundaryVertexCount;
}

uint32_t floatToFixed32(float value, int fractionalBits) {
	return static_cast<uint32_t>(std::round(value * (1 << fractionalBits)));
}

MeshType fixed32ToFloat(uint32_t fixedValue, int fractionalBits) {
	return static_cast<MeshType>(fixedValue) / (1 << fractionalBits);
}
//============================================================================

template class TriMesh<float>;
template class TriMesh<double>;
