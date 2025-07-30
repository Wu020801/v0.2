
#pragma once
#include "mesh.hpp"
#include "matrix.hpp"
#include "mutablePriorityHeap.hpp"

//============================================================================

using Real = double;

//============================================================================

struct Quadratic {
  Quadratic() = default;
  Quadratic(Real nx, Real ny, Real nz, Real d) { init(nx, ny, nz, d); }
  Quadratic(const Quadratic& other) { copy(other); }
  Quadratic& operator=(const Quadratic& rhs) {
    copy(rhs);
    return *this;
  }
  ~Quadratic() = default;

  const Real& operator[](const int i) const {
    assert(i < 10);
    return data_[i];
  }

  Real& operator[](const int i) {
    assert(i < 10);
    return data_[i];
  }

  void copy(const Quadratic& other) {
    for (int i = 0; i < 10; ++i) { data_[i] = other.data_[i]; }
  }

  void zero() {
    for (double& i : data_) { i = Real(0); }
  }

  void init(Real nx, Real ny, Real nz, Real d) {
    data_[0] = nx * nx;  // a2
    data_[1] = nx * ny;  // ab
    data_[2] = nx * nz;  // ac
    data_[3] = nx * d;   // ad
    data_[4] = ny * ny;  // b2
    data_[5] = ny * nz;  // bc
    data_[6] = ny * d;   // bd
    data_[7] = nz * nz;  // c2
    data_[8] = nz * d;   // cd
    data_[9] = d * d;    // d2
  }

  friend Quadratic operator+(const Quadratic& lhs, const Quadratic& rhs) {
    Quadratic res{};
    for (int i = 0; i < 10; ++i) { res.data_[i] = lhs.data_[i] + rhs.data_[i]; }
    return res;
  }

  Quadratic& operator+=(const Quadratic& rhs) {
    for (int i = 0; i < 10; ++i) { data_[i] += rhs.data_[i]; }
    return *this;
  }

  Quadratic& operator-=(const Quadratic& rhs) {
    for (int i = 0; i < 10; ++i) { data_[i] -= rhs.data_[i]; }
    return *this;
  }

  Quadratic& operator+=(double rhs) {
    for (double& i : data_) { i += rhs; }
    return *this;
  }

  Quadratic& operator-=(double rhs) {
    for (double& i : data_) { i -= rhs; }
    return *this;
  }

  Quadratic& operator*=(double rhs) {
    for (double& i : data_) { i *= rhs; }
    return *this;
  }

  bool minDistPoint(Vector3D<double>& point) const {
    Mat3Vmc<double> M;
    Mat3Vmc<double> iM;
    M[0][0] = data_[0];
    M[0][1] = M[1][0] = data_[1];
    M[0][2] = M[2][0] = data_[2];
    M[1][1]           = data_[4];
    M[1][2] = M[2][1] = data_[5];
    M[2][2]           = data_[7];
    if (M.inverse(iM)) {
      auto point0 = iM * Vector3D<double>(-data_[3], -data_[6], -data_[8]);
      for (int i = 0; i < 3; i++) { point[i] = point0[i]; }
      return true;
    }
    return false;
  }

  double operator()(const Vector3D<double>& point) const {
    const double x   = point[0];
    const double y   = point[1];
    const double z   = point[2];
    const double xx  = x * x;
    const double yy  = y * y;
    const double zz  = z * z;
    const double x2  = 2.0 * x;
    const double y2  = 2.0 * y;
    const double z2  = 2.0 * z;
    const double xy2 = x2 * y;
    const double xz2 = x2 * z;
    const double yz2 = y2 * z;
    return data_[0] * xx + data_[1] * xy2 + data_[2] * xz2 + data_[3] * x2 + data_[4] * yy + data_[5] * yz2
           + data_[6] * y2 + data_[7] * zz + data_[8] * z2 + data_[9];
  }

  Real data_[10]{};
};

//============================================================================

template<typename T, int N>
class SmallVector {
public:
  SmallVector() {
    // std::cout << "Constructor: N = " << N << std::endl;
    static_assert(N > 0, "initial size should be strictly higher than 0");
    data_      = data0_.data();
    size_      = 0;
    allocated_ = N;
  }

  // SmallVector(const SmallVector& other) { copy(other); }
  SmallVector(const SmallVector& other) {
    data_      = data0_.data();
    size_      = 0;
    allocated_ = N;
    copy(other);
  }

  SmallVector& operator=(const SmallVector& rhs) {
    copy(rhs);
    return *this;
  }

  ~SmallVector() = default;

  void copy(const SmallVector& other) {
    resize(other.size());
    std::copy_n(other.data(), other.size(), data());
  }

  const T* data() const { return data_; }
  T*       data() { return data_; }

  const T& operator[](const int i) const {
    assert(i < size());
    return data_[i];
  }

  T& operator[](const int i) {
    assert(i < size());
    return data_[i];
  }

  int size() const { return size_; }

  void reserve(const int sz) {
    if (sz >= allocated_) {
      std::vector<T> ndata(sz);
      std::copy(data_, data_ + size_, ndata.begin());
      allocated_ = sz;
      data1_.swap(ndata);
      data_ = data1_.data();
      return;
    }
  }

  void resize(const int sz) {
    if (sz > allocated_) { reserve(sz); }
    size_ = sz;
  }

  void push_back(const T& v) {
    // std::cout << "Pushing back: current size = " << _size << ", allocated = " << _allocated << std::endl;
    if (size_ == allocated_) {
      // std::cout << " reserving new capacity: " << 2 * _allocated << std::endl;
      reserve(2 * allocated_);
    }
    data_[size_++] = v;
  }

  int find(const T& v) const {
    for (int i = 0; i < size_; ++i) {
      if (data_[i] == v) { return i; }
    }
    return -1;
  }

  void erase(const T& v) {
    const auto i0 = find(v);
    if (i0 != -1) {
      --size_;
      for (int i = i0; i < size_; ++i) { data_[i] = data_[i + 1]; }
    }
  }

private:
  int              size_{};
  int              allocated_{};
  T*               data_;
  std::array<T, N> data0_;
  std::vector<T>   data1_;
};

//============================================================================

using SmallIntVec = SmallVector<int, 8>;

//============================================================================

enum class VertPlacement {
  END_POINTS              = 0,
  END_POINTS_OR_MID_POINT = 1,
  OPTIMAL                 = 2
};

//============================================================================

struct TriMeshDeciParams {
  // Triangle flipping parameters
  double triFlipCost          = 1e+7;
  double triFlipThres         = 0.3;
  double boundWeight          = 1e+7;
  double angleQualityThres    = 1.0;
  double maxErr               = 1.0;
  bool   areaWeightQuadratics = true;
  int    faceCount            = 1;
  int    xyzCount             = 1;

  // Vertex placement mode
  VertPlacement vertPlacement = VertPlacement::OPTIMAL;

};

//============================================================================

struct DecVertex {
  Vector3D<double> pos_{};
  Quadratic        q_{};
  int              tag_{0};
  SmallIntVec      faces_;
  SmallIntVec      edges_;
  bool             valid_{false};
  bool             flag_{false};
  DecVertex() { q_.zero(); }
};

//============================================================================

struct DecTriangle {
  int                tag_{};
  Vector3D<uint32_t> indices_{};
  Vector3D<double>   normal_{};
  bool               zeroArea_{};
  bool               valid_{};

  DecTriangle(const Vector3D<uint32_t>& triangle) : valid_(true) {
    indices_[0] = triangle[0];
    indices_[1] = triangle[1];
    indices_[2] = triangle[2];
  }
};

//============================================================================

struct DecEdge {
  int    heapPosition() const { return heapPos_; }
  double heapKey() const { return heapKey_; }
  void   setHeapPosition(int hpos) { heapPos_ = hpos; }
  double setHeapKey(double hkey) { return heapKey_ = hkey; }

  int opposite(int v) const {
    assert(v == vidx0_ || v == vidx1_);
    return v == vidx0_ ? vidx1_ : vidx0_;
  }

  double           heapKey_;
  Vector3D<double> pos_;
  int32_t          vidx0_;
  int32_t          vidx1_;
  int32_t          tidx0_;
  int32_t          tid1_;
  int32_t          heapPos_;
  bool             valid_;
  bool             flag_{false};
};

//============================================================================

struct encoderParams;

//============================================================================

template<typename T>
class GeometrySimplify {
public:
  GeometrySimplify()  = default;
  ~GeometrySimplify() = default;

  bool   removeDuplVert(TriMesh<T>& mesh);
  //bool   removeDuplVert(const TriMesh<T>& mesh, TriMesh<T>& umesh);
  bool   removeDuplFaces(TriMesh<T>& mesh);
  bool   removeSmallConnectedComponents(TriMesh<T>& mesh, int minCCFaceCount);
  void   calAdjVert(uint32_t vindex, SmallIntVec& vadj);
  void   tagAdjFace(uint32_t vindex, uint32_t tag);
  void   incTagAdjFace(uint32_t vindex);
  void   calEdgeAdjTri(uint32_t vidx0, uint32_t vidx1, SmallIntVec& triAdj);
  void   calAdjTriQual(int32_t vidx, const Vector3D<double>& pos, double& angleQual, double& normalQual);
  double calCost(const Quadratic&         q,
                 const Vector3D<double>&  pos,
                 const int32_t            vidx0,
                 const int32_t            vidx1,
                 const TriMeshDeciParams& params);
  void   calEdgeCollapseCost(DecEdge& edge, const TriMeshDeciParams& params);

  bool   extractSimpMesh(std::vector<Vector3D<T>>& points, std::vector<Vector3D<uint32_t>>& triangles);
  void   PrintProgress(double progress);
  bool   verticesAreValid(DecEdge* edge);
  void   extractAdjTri(int vidx, int tag, SmallIntVec& adj);
  void   extractAdjTri(int vidx, int tag0, int tag1, SmallIntVec& adj0, SmallIntVec& adj1);
  void   updateEdgeAndVert(SmallIntVec& eadj0, const SmallIntVec eadj1, int vidx0, int vidx1);
  void   deleteFaces(const SmallIntVec& deleted);
  void   updateModTriangles(int vidx0, int vidx1, const SmallIntVec& modified);
  void   updateVertPosAndQuadratic(int vidx0, int vidx1, DecEdge* edge, const TriMeshDeciParams& params);
  void   updateEdgeCollapes(SmallIntVec& eadj0,
                            int          vidx0,
                          const TriMeshDeciParams& params);

  void processEdge(DecEdge* edge, const TriMeshDeciParams& params);
  bool isSafe(const DecEdge& edge, const TriMeshDeciParams& params);
  bool simplify(const TriMeshDeciParams& params);


  bool   init(const std::vector<Vector3D<T>>&        points,
			  const std::vector<Vector3D<uint32_t>>& triangles,
			  const TriMeshDeciParams&               params);
  bool simplify(const TriMesh<T>& mesh, TriMesh<T>& dmesh, TriMesh<T>& mmesh,
				const int faceCount);

private:
  uint32_t                     faceCount_;
  uint32_t                     xyzCount_;
  double                       scale_;
  Vector3D<double>             center_;
  std::vector<DecVertex>       vertices_;
  std::vector<DecEdge>         edges_;
  std::vector<DecTriangle>     faces_;
  MutablePriorityHeap<DecEdge> heap_;
};

