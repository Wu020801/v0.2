
#include <chrono>
#include "misc.hpp"
#include "mesh.hpp"
#include "vector.hpp"
#include "geometryDecimate.hpp"

//============================================================================

template<typename T>
bool
GeometrySimplify<T>::removeDuplVert(TriMesh<T>& mesh) {

  // 1.1. UnifyVertices
  const auto                      xyzCount0  = mesh.XYZ.size();
  const auto                      faceCount0 = mesh.FaceXYZ.size();
  std::vector<Vector3D<T>>        xyzOut;
  std::vector<Vector3D<uint32_t>> faceOut;
  std::vector<uint32_t>                mapping;
  std::vector<uint32_t>                orgIdx;
  removeDuplicateVert(mesh.XYZ, mesh.FaceXYZ, xyzOut, faceOut, orgIdx, mapping);
  std::swap(mesh.XYZ, xyzOut);
  std::swap(mesh.FaceXYZ, faceOut);
  // 1.2. Remove degenerated triangles
  removeDegenerateTriangles(mesh);
  return true;
}

// =========================================================================

template<typename T>
bool
GeometrySimplify<T>::removeDuplFaces(TriMesh<T>& mesh) {
  AdjInfo                  vert2Face;
  const auto               noTriFace = mesh.FaceXYZ.size();
  std::vector<Vector3D<T>> triNormals;
  mesh.computeTriangleNormals(triNormals);

  computeVertex2triangle(mesh.XYZ.size(), mesh.FaceXYZ, vert2Face);

  std::vector<Vector3D<uint32_t>> faceXyzOut;
  std::vector<uint32_t>           uniqOrgIdx, org2Unique;
  removeDuplicateTriangles(mesh.FaceXYZ, triNormals, vert2Face, faceXyzOut, uniqOrgIdx, org2Unique);
  std::swap(mesh.FaceXYZ, faceXyzOut);

  return true;
}

// =========================================================================

template<typename T>
bool
GeometrySimplify<T>::removeSmallConnectedComponents(TriMesh<T>& mesh, int minCCFaceCount) {

  const auto                               xyzCount  = mesh.XYZ.size();
  const auto                               faceCount = mesh.FaceXYZ.size();
  std::vector<std::shared_ptr<TriMesh<T>>> ccMeshes;
  auto                                     ccCount0 = extractConnectedComponents(mesh, &ccMeshes);

  mesh.clear();
  int ccCount = 0;
  for (size_t i = 0; i < ccCount0; ++i) {
    const auto& ccMesh = ccMeshes[i];
    if (ccMesh->FaceXYZ.size() < minCCFaceCount) {
      // std::cout << "\t\t Skipping CC " << i << " with " << ccMesh->faceXYZ.size() << " faces\n";
      continue;
    }

    ++ccCount;
    mesh.append(*ccMesh);
  }

  return true;
}

// =========================================================================
template<typename T>
void
GeometrySimplify<T>::calAdjVert(uint32_t vindex, SmallIntVec& vadj) {
  vadj.resize(0);
  const auto& triAdj = vertices_[vindex].faces_;
  for (size_t i = 0, count = triAdj.size(); i < count; ++i) {
    const auto& triangle        = faces_[triAdj[i]].indices_;
    vertices_[triangle[0]].tag_ = 1;
    vertices_[triangle[1]].tag_ = 1;
    vertices_[triangle[2]].tag_ = 1;
  }

  vertices_[vindex].tag_ = 0;
  for (size_t i = 0, count = triAdj.size(); i < count; ++i) {
    const auto& triangle = faces_[triAdj[i]].indices_;
    for (int j = 0; j < 3; ++j) {
      const auto idx    = triangle[j];
      auto&      vertex = vertices_[idx];
      if (vertex.tag_ == 1) {
        vertex.tag_ = 0;
        vadj.push_back(idx);
      }
    }
  }
}

// =========================================================================
template<typename T>
void
GeometrySimplify<T>::tagAdjFace(uint32_t vindex, uint32_t tag) {
  const auto& tadj = vertices_[vindex].faces_;
  for (int i = 0, count = tadj.size(); i < count; ++i) { faces_[tadj[i]].tag_ = tag; }
}

// =========================================================================
template<typename T>
void
GeometrySimplify<T>::incTagAdjFace(uint32_t vindex) {
  const auto& tadj = vertices_[vindex].faces_;
  for (int i = 0, count = tadj.size(); i < count; ++i) { ++faces_[tadj[i]].tag_; }
}

// =========================================================================
template<typename T>
void
GeometrySimplify<T>::calEdgeAdjTri(uint32_t vidx0, uint32_t vidx1, SmallIntVec& triAdj) {
  triAdj.resize(0);
  tagAdjFace(vidx0, 0);
  tagAdjFace(vidx1, 1);
  const auto& triAdj0 = vertices_[vidx0].faces_;
  for (size_t i = 0, count = triAdj0.size(); i < count; ++i) {
    const auto triIdx = triAdj0[i];
    if (faces_[triIdx].tag_ != 0) { triAdj.push_back(triIdx); }
  }
}

// =========================================================================

template<typename T>
void
GeometrySimplify<T>::calAdjTriQual(int32_t vidx, const Vector3D<double>& pos, double& angleQual, double& normalQual) {
  const auto& triAdj = vertices_[vidx].faces_;
  for (int i = 0, count = triAdj.size(); i < count; ++i) {
    const auto  triIdx = triAdj[i];
    const auto& tri    = faces_[triIdx];
    if (tri.tag_ != 1 || tri.zeroArea_) { continue; }
    int j = 0, k = 0;

    if (tri.indices_[0] == vidx) {
      j = tri.indices_[1];
      k = tri.indices_[2];
    } else if (tri.indices_[1] == vidx) {
      j = tri.indices_[2];
      k = tri.indices_[0];
    } else {
      j = tri.indices_[0];
      k = tri.indices_[1];
    }

    const auto& normal0     = tri.normal_;
    const auto  e0          = vertices_[j].pos_ - pos;
    const auto  e1          = vertices_[k].pos_ - pos;
    const auto  e2          = vertices_[j].pos_ - vertices_[k].pos_;
    const auto  maxEdgeLen2 = std::max(std::max(e0.L2normSq(), e1.L2normSq()), e2.L2normSq());

    auto       normal1 = e0 ^ e1;
    const auto triArea = normal1.L2norm();
    if (triArea > 0.0) { normal1 /= triArea; }

    angleQual  = std::min(angleQual, maxEdgeLen2 > 0.0 ? triArea / maxEdgeLen2 : 0.0);
    normalQual = std::min(normalQual, normal0 * normal1);
  }
}

// =========================================================================

template<typename T>
double
GeometrySimplify<T>::calCost(const Quadratic&         q,
                             const Vector3D<double>&  pos,
                             const int32_t            vidx0,
                             const int32_t            vidx1,
                             const TriMeshDeciParams& params) {
  double     err                = q(pos);
  const bool isPenalizeTriFlip  = params.triFlipCost > 0.0;
  const bool isPenalizeElongTri = params.angleQualityThres > 0.0;

  if (isPenalizeTriFlip || isPenalizeElongTri) {
    double normalQual = 1.0;
    double angleQual  = params.angleQualityThres;

    calAdjTriQual(vidx0, pos, angleQual, normalQual);
    calAdjTriQual(vidx1, pos, angleQual, normalQual);

    if (isPenalizeElongTri) { err /= (1.0e-7 + angleQual); }

    if (isPenalizeTriFlip && normalQual < params.triFlipThres) {
      err += (params.triFlipThres - normalQual) * params.triFlipCost;
    }
  }

  return err;
}

// =========================================================================

template<typename T>
void
GeometrySimplify<T>::calEdgeCollapseCost(DecEdge& edge, const TriMeshDeciParams& params) {
  const auto  vidx0 = edge.vidx0_;
  const auto  vidx1 = edge.vidx1_;
  const auto& v0    = vertices_[vidx0];
  const auto& v1    = vertices_[vidx1];
  const auto  q     = v0.q_ + v1.q_;

  // tag the incident triangles
  tagAdjFace(vidx0, 0);
  tagAdjFace(vidx1, 1);
  incTagAdjFace(vidx0);

  // compute the best position for merged point
  auto       pos   = v0.pos_;
  auto       err   = calCost(q, v0.pos_, vidx0, vidx1, params);
  const auto errV1 = calCost(q, v1.pos_, vidx0, vidx1, params);
  if (errV1 < err) {
    err = errV1;
    pos = v1.pos_;
  }

  if (params.vertPlacement == VertPlacement::END_POINTS_OR_MID_POINT
      || params.vertPlacement == VertPlacement::OPTIMAL) {
    const auto posMidPoint(0.5 * (v0.pos_ + v1.pos_));
    const auto errMidPoint = calCost(q, posMidPoint, vidx0, vidx1, params);
    if (errMidPoint < err) {
      err = errMidPoint;
      pos = posMidPoint;
    }
  }

  Vector3D<double> posOpt{0.0, 0.0, 0.0};
  if (params.vertPlacement == VertPlacement::OPTIMAL && q.minDistPoint(posOpt) > 0.0) {
    const auto errOpt = calCost(q, posOpt, vidx0, vidx1, params);
    if (errOpt < err) {
      err = errOpt;
      pos = posOpt;
    }
  }

  edge.pos_ = pos;

  edge.setHeapKey(-err);
  if (edge.heapPosition() == -1) {
    heap_.insert(edge);
  } else {
    heap_.update(edge);
  }

}

// =========================================================================

template<typename T>
bool
GeometrySimplify<T>::init(const std::vector<Vector3D<T>>&        points,
                          const std::vector<Vector3D<uint32_t>>& triangles,
                          const TriMeshDeciParams&               params) {
 //1. init
  vertices_.reserve(xyzCount_);
  faces_.reserve(faceCount_);
  edges_.reserve(faceCount_ * 3);
  center_ = Vector3D<T>(0.0, 0.0, 0.0);

  auto createDecVertex = [](const Vector3D<T>& pos) {
    DecVertex vertex;
    vertex.pos_   = pos;
    vertex.valid_ = true;
    return vertex;
  };

  for (int32_t i = 0; i < points.size(); ++i) {
    vertices_.emplace_back(createDecVertex(points[i]));
  }

  // 2. scale mesh to improve stability
  // calculate mesh center
  for (const auto& v : vertices_) { center_ += v.pos_; }
  center_ /= vertices_.size();

  // Calculate radius
  double maxRadius = 0.0;
  for (const auto& v : vertices_) {
    const auto r = (v.pos_ - center_).L2normSq();
    if (r > maxRadius) { maxRadius = r; }
  }
  maxRadius = std::sqrt(maxRadius);

  // scale vertices
  scale_ = (maxRadius > 0.0) ? (1.0 / maxRadius) : 1.0;
  for (auto& v : vertices_) { v.pos_ = (v.pos_ - center_) * scale_; }

  // triangle initalization
  for (int t = 0, c = 0; t < faceCount_; ++t, c += 3) {
    DecTriangle triangle = triangles[t];
    const auto  i        = triangles[t][0];
    const auto  j        = triangles[t][1];
    const auto  k        = triangles[t][2];
    auto&       vi       = vertices_[i];
    auto&       vj       = vertices_[j];
    auto&       vk       = vertices_[k];
    vi.faces_.push_back(t);
    vj.faces_.push_back(t);
    vk.faces_.push_back(t);
    const Vector3D<double>& posi   = vi.pos_;
    const Vector3D<double>& posj   = vj.pos_;
    const Vector3D<double>& posk   = vk.pos_;
    Vector3D<double> normal = computeTriangleNormal(posi, posj, posk, false);

    const auto area = normal.norm();
    if (area > 0.0) {
      normal /= area;
      triangle.normal_   = normal;
      triangle.zeroArea_ = false;
    } else {
      triangle.zeroArea_ = true;
      triangle.normal_   = Vector3D<double>(0.0, 0.0, 0.0);
    }

    // compute quadratic error
	double q4 = -normal * posi;
    Quadratic q(normal[0], normal[1], normal[2], q4);
    if (params.areaWeightQuadratics) { q *= area; }

    vi.q_ += q;
    vj.q_ += q;
    vk.q_ += q;
    faces_.push_back(triangle);
  }

  // 3. Initlize edges info
  SmallIntVec verAdj;
  int         edgeCount = 0;
  for (uint32_t vidx0 = 0; vidx0 < vertices_.size(); ++vidx0) {
    calAdjVert(vidx0, verAdj);
    for (size_t i = 0, count = verAdj.size(); i < count; ++i) {
      const auto vidx1 = static_cast<uint32_t>(verAdj[i]);
      edgeCount += (vidx0 < vidx1) ? 1 : 0;
    }
  }
  edges_.resize(edgeCount);

  // 4. Triangle Update
  SmallIntVec triAdj;
  for (uint32_t vidx0 = 0, eidx = 0; vidx0 < xyzCount_; ++vidx0) {
    calAdjVert(vidx0, verAdj);
    for (size_t i = 0, count = verAdj.size(); i < count; ++i) {
      const auto vidx1 = verAdj[i];
      if (vidx0 < vidx1) {
        calEdgeAdjTri(vidx0, vidx1, triAdj);

        DecEdge& edge = edges_[eidx];
        edge.setHeapPosition(-1);
        edge.valid_ = true;
        edge.vidx0_ = vidx0;
        edge.vidx1_ = vidx1;
        edge.flag_  = vertices_[vidx0].flag_;

        vertices_[vidx0].edges_.push_back(eidx);
        vertices_[vidx1].edges_.push_back(eidx);

        auto& v0 = vertices_[vidx0];
        auto& v1 = vertices_[vidx1];
        if (triAdj.size() == 1 && params.boundWeight > 0.0) {
          const auto& pos0      = v0.pos_;
          const auto& pos1      = v1.pos_;
          const auto  delta     = pos1 - pos0;
          const auto& triNormal = faces_[triAdj[0]].normal_;
          auto        normal    = delta ^ triNormal;
          const auto  len       = normal.L2norm();
          if (len > 0.0) {
            normal /= len;
            Quadratic q(normal[0], normal[1], normal[2], -normal * pos0);
            if (false) {
              const auto weight = (params.areaWeightQuadratics ? delta.L2normSq() : 1.0) * params.boundWeight;
              q *= weight;
            } else {
              if (params.areaWeightQuadratics) { q *= delta.L2normSq(); }
              q *= params.boundWeight;
            }
            v0.q_ += q;
            v1.q_ += q;
          }
        }
        ++eidx;
      }
    }
  }

  // 5. Compute edge collapse cost
  int count = 0;
  for (auto& edge : edges_) {
    calEdgeCollapseCost(edge, params);
    count++;
  }
  return true;
}

// =========================================================================

/**
* @brief Extracts the simplified mesh data after the decimation process.
* @param points A vector of simplified mesh vertices as 3D points.
* @param triangles A vector of simplified mesh triangles, represented as indices of the corresponding vertices.
* @return true if the extraction is successful, false otherwise.
*/
template<typename T>
bool
GeometrySimplify<T>::extractSimpMesh(std::vector<Vector3D<T>>& points, std::vector<Vector3D<uint32_t>>& triangles) {
  // 1. Extract decimated vertex and vertex mapping
  std::vector<int32_t> vmap(vertices_.size(), -1);
  int                  vcount = 0;

  for (size_t i = 0; i < vertices_.size(); ++i) {
    if (vertices_[i].valid_) {
      points[vcount] = vertices_[i].pos_ / scale_ + center_;
      vmap[i] = vcount++;
    }
  }

  // 2. Clean up invalid triangle
  std::vector<Vector3D<uint32_t>> valid_triangles;
  for (const auto& triangle : faces_) {
    if (triangle.valid_) {
      //auto& tri = triangles[tcount++];
      const auto i = triangle.indices_[0];
      const auto j = triangle.indices_[1];
      const auto k = triangle.indices_[2];
      valid_triangles.push_back(Vector3D<uint32_t>(vmap[i], vmap[j], vmap[k]));
    }
  }

  if (valid_triangles.size() != triangles.size()) {
    std::cout << "Warning: invalid triangles found in the decimated mesh.\n";
    return false;
  }
  // triangles =  valid_triangles;
  std::swap(triangles, valid_triangles);

  return true;
}

// =========================================================================

/**
 * @brief Prints a progress bar to the console, representing the given progress value.
 * @param progress A value between 0.0 and 1.0, where 0.0 indicates 0% progress and
 *                 1.0 indicates 100% progress.
 */
template<typename T>
void
GeometrySimplify<T>::PrintProgress(double progress) {
  const int barWidth = 70;
  std::cout << "\t [";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) {
      std::cout << '=';
    } else if (i == pos) {
      std::cout << '>';
    } else {
      std::cout << ' ';
    }
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  std::cout << '\n';
}

// =========================================================================

/**
* @brief Checks if the vertices of the given edge are valid.
* @param edge A pointer to a DecEdge object, whose vertices' validity will be checked.
* @return Returns true if both vertices of the edge are valid, otherwise returns false.
*/
template<typename T>
bool
GeometrySimplify<T>::verticesAreValid(DecEdge* edge) {
  const auto vidx0 = edge->vidx0_;
  const auto vidx1 = edge->vidx1_;
  return vertices_[vidx0].valid_ && vertices_[vidx1].valid_;
}

// =========================================================================

/**
* @brief Extracts the adjacent triangles of a given vertex with a specific tag.
* @param vidx The index of the vertex whose adjacent triangles will be extracted.
* @param tag The specific tag value to filter the adjacent triangles.
* @param[out] adj A reference to a SmallIntVec object that will store the indices of 
*                 the extracted adjacent triangles.
*/
template<typename T>
void
GeometrySimplify<T>::extractAdjTri(int vidx, int tag, SmallIntVec& adj) {
  const auto& tadj = vertices_[vidx].faces_;
  for (int i = 0, count = tadj.size(); i < count; ++i) {
    const auto tidx = tadj[i];
    if (faces_[tidx].tag_ == tag) { adj.push_back(tidx); }
  }
}

// =========================================================================

/**
* @brief Extracts the adjacent triangles of a given vertex with two tags,
*        and stores them separately in two SmallIntVec objects.
* @param vidx The index of the vertex whose adjacent triangles will be extracted.
* @param tag0 The first tag value to filter the adjacent triangles.
* @param tag1 The second tag value to filter the adjacent triangles.
* @param[out] adj0 A reference to a SmallIntVec object that will store the indices of
*                  the extracted adjacent triangles with tag0.
* @param[out] adj1 A reference to a SmallIntVec object that will store the indices of
*                  the extracted adjacent triangles with tag1.
*/
template<typename T>
void
GeometrySimplify<T>::extractAdjTri(int vidx, int tag0, int tag1, SmallIntVec& adj0, SmallIntVec& adj1) {
  const auto& tadj = vertices_[vidx].faces_;
  for (int i = 0, count = tadj.size(); i < count; ++i) {
    const auto tidx = tadj[i];
    if (faces_[tidx].tag_ == tag0) {
      adj0.push_back(tidx);
    } else if (faces_[tidx].tag_ == tag1) {
      adj1.push_back(tidx);
    }
  }
}

// =========================================================================

/**
* @brief Updates the edges and vertices after an edge collapse operation.
* @param[in,out] eadj0 A reference to a SmallIntVec object that contains the list of
*                      edges adjacent to vidx0 before the edge collapse operation.
*                      The function will update this list with additional edges.
* @param eadj1 A SmallIntVec object that contains the list of edges adjacent to vidx1.
*              This list is used to update eadj0.
* @param vidx0 The index of the vertex that will be removed after the edge collapse operation.
* @param vidx1 The index of the vertex that will remain after the edge collapse operation.
*/
template<typename T>
void
GeometrySimplify<T>::updateEdgeAndVert(SmallIntVec& eadj0, const SmallIntVec eadj1, int vidx0, int vidx1) {
  SmallIntVec vadj;
  for (int i = 0, count = eadj0.size(); i < count; ++i) { vadj.push_back(edges_[eadj0[i]].opposite(vidx0)); }

  // update vertices
  for (int i = 0, count = eadj1.size(); i < count; ++i) {
    const auto eidx  = eadj1[i];
    auto&      cedge = edges_[eidx];
    const auto vidx  = cedge.opposite(vidx1);
    assert(vidx >= 0);

    if (vidx == vidx0) {  // remove collapsed edge
      vertices_[vidx0].edges_.erase(eidx);
      heap_.remove(cedge);
      cedge.valid_ = false;
    } else {
      const auto eidx0 = vadj.find(vidx);
      if (eidx0 != -1) {  // remove duplicated edge
        vertices_[vidx].edges_.erase(eidx);
        heap_.remove(cedge);
        cedge.valid_ = false;
      } else {  // update edge
        cedge.vidx0_ = vidx0;
        cedge.vidx1_ = vidx;
        eadj0.push_back(eidx);
      }
    }
  }
}

// =========================================================================

/**
* @brief Deletes the specified faces from the mesh.
* @param[in] deleted A SmallIntVec object that contains the indices of the faces to be deleted.
*                    The indices must be valid indices of the triangles_ vector.
*/
template<typename T>
void
GeometrySimplify<T>::deleteFaces(const SmallIntVec& deleted) {
  for (int i = 0, count = deleted.size(); i < count; ++i) {
    const auto tidx = deleted[i];
    auto&      ctri = faces_[tidx];
    ctri.valid_     = false;
    for (int j = 0; j < 3; ++j) { vertices_[ctri.indices_[j]].faces_.erase(tidx); }
  }
  faceCount_ -= deleted.size();
}

// =========================================================================

/**
* @brief Updates the modified triangles after an edge collapse operation.
* @param[in] vidx0 The index of the vertex that will be removed after the edge collapse operation.
* @param[in] vidx1 The index of the vertex that will remain after the edge collapse operation.
* @param[in] modified A SmallIntVec object that contains the indices of the triangles that
*                     are modified after the edge collapse operation.
*/
template<typename T>
void
GeometrySimplify<T>::updateModTriangles(int vidx0, int vidx1, const SmallIntVec& modified) {
  auto& tadj0 = vertices_[vidx0].faces_;
  for (int i = 0, count = modified.size(); i < count; ++i) {
    const auto tidx = modified[i];
    auto&      ctri = faces_[tidx];
    for (int j = 0; j < 3; ++j) {
      if (ctri.indices_[j] == vidx1) { ctri.indices_[j] = vidx0; }
    }
    tadj0.push_back(tidx);
  }
}

// =========================================================================

/**
* @brief Updates the position and quadratic error of the vertex after an edge collapse operation.
* @param[in] vidx0 The index of the vertex that will be removed after the edge collapse operation.
* @param[in] vidx1 The index of the vertex that will remain after the edge collapse operation.
* @param[in] edge A DecEdge object that contains information about the edge that will be collapsed.
* @param[in] params A reference to a TriMeshDeciParams object that contains the parameters for the mesh decimation.
*/
template<typename T>
void
GeometrySimplify<T>::updateVertPosAndQuadratic(int vidx0, int vidx1, DecEdge* edge, const TriMeshDeciParams& params) {
  auto& v0 = vertices_[vidx0];
  auto& v1 = vertices_[vidx1];

  v1.valid_ = false;
  --xyzCount_;

  v0.pos_ = edge->pos_;
  v0.q_ += v1.q_;
}

// =========================================================================

/**
* @brief Updates the edge collapse cost of the edges that are affected by a vertex removal operation.
* @param[in,out] eadj0 A SmallIntVec object that contains the indices of the edges that are affected
*                      by the vertex removal operation. After this function returns, the cost of the
*                      edges in the vector will be updated according to the new edge collapse information.
* @param[in] vidx0 The index of the vertex that will be removed after the edge collapse operation.
* @param[in] params A reference to a TriMeshDeciParams object that contains the parameters for the mesh decimation.
*/
template<typename T>
void
GeometrySimplify<T>::updateEdgeCollapes(SmallIntVec& eadj0,
                                        int          vidx0,
                                        const TriMeshDeciParams& params) {
  const auto tcount = vertices_[vidx0].faces_.size();
  if (tcount == 0) { return; }

  for (int i = 0, count = eadj0.size(); i < count; ++i) {
    calEdgeCollapseCost(edges_[eadj0[i]], params);
  }

  for (int i = 0; i < tcount; ++i) {
    const auto tidx = vertices_[vidx0].faces_[i];
    auto&      ctri = faces_[tidx].indices_;
    int32_t    a = 0, b = 0;
    if (ctri[0] == vidx0) {
      a = ctri[1];
      b = ctri[2];
    } else if (ctri[1] == vidx0) {
      a = ctri[0];
      b = ctri[2];
    } else {
      a = ctri[0];
      b = ctri[1];
    }

    for (int j = 0, ecount = vertices_[a].edges_.size(); j < ecount; ++j) {
      const auto eidx = vertices_[a].edges_[j];
      if ((edges_[eidx].vidx0_ == a && edges_[eidx].vidx1_ == b)
          || (edges_[eidx].vidx0_ == b && edges_[eidx].vidx1_ == a)) {
        calEdgeCollapseCost(edges_[eidx], params);
        // break;
      }
    }
  }
}

// =========================================================================

/**
* @brief Processes the given edge for mesh simplification, modifying the mesh as needed.
* @param edge A pointer to the edge to be processed.
* @param params The TriMeshDeciParams object that specifies the simplification parameters.
*/
template<typename T>
void
GeometrySimplify<T>::processEdge(DecEdge* edge, const TriMeshDeciParams& params) {
  const auto  vidx0 = edge->vidx0_;
  const auto  vidx1 = edge->vidx1_;
  SmallIntVec modified0, modified1, deleted;
  tagAdjFace(vidx1, 0);
  tagAdjFace(vidx0, 1);
  incTagAdjFace(vidx1);

  extractAdjTri(vidx0, 1, 2, modified0, deleted);
  extractAdjTri(vidx1, 1, modified1);

  // update edge and vertices
  auto&       eadj0 = vertices_[vidx0].edges_;
  const auto& eadj1 = vertices_[vidx1].edges_;
  updateEdgeAndVert(eadj0, eadj1, vidx0, vidx1);

  deleteFaces(deleted);
  updateModTriangles(vidx0, vidx1, modified1);
  updateVertPosAndQuadratic(vidx0, vidx1, edge, params);
  updateEdgeCollapes(eadj0, vidx0, params);

}

// =========================================================================

template<typename T>
bool
GeometrySimplify<T>::isSafe(const DecEdge& edge, const TriMeshDeciParams& params) {
  return true;
}

// =========================================================================

// main simplification
template<typename T>
bool
GeometrySimplify<T>::simplify(const TriMeshDeciParams& params) {

  double       lastProgress   = -1.0;
  const double triangleCount0 = faceCount_;
  double       progressStep   = 0.1;
  bool         stop           = false;

  while (!stop && faceCount_ > params.faceCount && xyzCount_ > params.xyzCount) {
    double progress = 1.0 - faceCount_ / triangleCount0;
    if (progress - lastProgress > progressStep) {
      lastProgress = progress;
    }

    if (DecEdge* const edge = heap_.extract()) {
      if (edge->heapKey_ < -params.maxErr) {
        stop = true;
        break;
      }

      if (!isSafe(*edge, params)) {
        edge->setHeapKey(std::numeric_limits<double>::lowest());
        heap_.insert(*edge);
      } else {
        if (verticesAreValid(edge)) { processEdge(edge, params); }
      }
    }
  }

  return true;
}

// =========================================================================


template<typename T>
bool
GeometrySimplify<T>::simplify(const TriMesh<T>&           imesh,
                              TriMesh<T>&                 dmesh,
                              TriMesh<T>&                 mmesh,
                              const int                   faceCount) {
  TriMeshDeciParams deciParams;
  deciParams.faceCount = faceCount;

  xyzCount_  = imesh.XYZ.size();
  faceCount_ = imesh.FaceXYZ.size();
  init(imesh.XYZ, imesh.FaceXYZ, deciParams);

  simplify(deciParams);
  dmesh.XYZ.resize(xyzCount_);
  dmesh.FaceXYZ.resize(faceCount_);
  extractSimpMesh(dmesh.XYZ, dmesh.FaceXYZ);
  mmesh.FaceXYZ = imesh.FaceXYZ;

  if (!removeDuplFaces(dmesh)) {
    std::cerr << "Error: can't remove small connected components!\n";
    return 1;
  }

  uint32_t minCCTriangleCount = 0;
  if (!removeSmallConnectedComponents(dmesh, minCCTriangleCount)) {
    std::cerr << "Error: can't remove small connected components! \n";
    return 1;
  }

  return true;
}
//============================================================================

template class GeometrySimplify<MeshType>;

