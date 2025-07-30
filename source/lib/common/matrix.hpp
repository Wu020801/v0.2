
#pragma once

#include <cassert>
#include <cstdint>
#include <ostream>
#include "vector.hpp"

//============================================================================

template<typename T>
class Matrix {
public:
  Matrix()  = default;
  ~Matrix() = default;

  Matrix(uint32_t h, uint32_t w) { resize(h, w); }

  Matrix(uint32_t h, uint32_t w, T const v) {
    resize(h, w);
    fill(v);
  }

  template<typename S>
  Matrix(Matrix<S> const& mat) {
    *this = mat;
  }

  template<typename S>
  Matrix& operator=(Matrix<S> const& rhs) {
    resize(rhs.height(), rhs.width());
    for (uint32_t r = 0; r < _height; r++) _data[r] = rhs[r];
    return *this;
  }

  VectorND<T>&        operator[](uint32_t const& r) { return _data[r]; }
  VectorND<T> const&  operator[](uint32_t const& r) const { return _data[r]; }
  inline VectorND<T>& row(uint32_t const r) { return _data[r]; }
  VectorND<T> const&  row(uint32_t const r) const { return _data[r]; }
  inline T&           at(uint32_t const r, uint32_t const c) { return _data[r][c]; }

  inline T const& at(uint32_t const r, uint32_t const c) const { return _data[r][c]; }

  template<typename S>
  Matrix& operator+=(Matrix<S> const& rhs) {
    assert(_height == rhs.height() && _width == rhs.width());
    for (uint32_t r = 0; r < _height; r++) _data[r] += rhs[r];
    return *this;
  }

  template<typename S>
  Matrix& operator-=(Matrix<S> const& rhs) {
    assert(_height == rhs.height() && _width == rhs.width());
    for (uint32_t r = 0; r < _height; r++) _data[r] -= rhs[r];
    return *this;
  }

  template<typename S>
  Matrix& operator*=(S const& rhs) {
    for (uint32_t r = 0; r < _height; r++) _data[r] *= rhs;
    return *this;
  }

  template<typename S>
  Matrix& operator/=(S const& rhs) {
    for (uint32_t r = 0; r < _height; r++) _data[r] /= rhs;
    return *this;
  }

  bool operator==(Matrix& rhs) const {
    if (dim() != rhs.dim()) return false;
    for (uint32_t r = 0; r < _height; r++)
      if (_data[r] != rhs[r]) return false;
    return true;
  }

  bool operator!=(Matrix& rhs) const {
    if (dim() != rhs.dim()) return true;
    for (uint32_t r = 0; r < _height; r++)
      if (_data[r] != rhs[r]) return true;
    return false;
  }

  Matrix operator-() const {
    Matrix mat(_height, _width);
    for (uint32_t r = 0; r < _height; r++) mat[r] = -_data[r];
    return mat;
  }

  friend Matrix operator+(Matrix const& lhs, Matrix const& rhs) {
    Matrix mat(lhs);
    mat += rhs;
    return mat;
  }

  friend Matrix operator-(Matrix const& lhs, Matrix const& rhs) {
    Matrix mat(lhs);
    mat -= rhs;
    return mat;
  }

  template<typename S>
  friend Matrix operator*(Matrix const& lhs, S const& rhs) {
    Matrix mat(lhs);
    mat *= rhs;
    return mat;
  }

  template<typename S>
  friend Matrix operator*(S const& lhs, Matrix const& rhs) {
    Matrix mat(rhs);
    mat *= lhs;
    return mat;
  }

  template<typename S>
  friend Matrix operator/(Matrix const& lhs, S const& rhs) {
    Matrix mat(lhs);
    mat /= rhs;
    return mat;
  }

  friend VectorND<T> operator*(Matrix const& lhs, VectorND<T> const& rhs) {
    assert(lhs.width() == rhs.size());
    VectorND<T> prod(lhs.height());
    for (uint32_t r = 0; r < lhs.height(); r++) prod[r] = lhs[r] * rhs;
    return prod;
  }

  friend VectorND<T> operator*(VectorND<T> const& lhs, Matrix const& rhs) {
    assert(lhs.size() == rhs.height());
    VectorND<T> prod(rhs.width());
    auto const& mat = rhs.transpose();
    for (uint32_t r = 0; r < mat.height(); r++) prod[r] = mat[r] * lhs;
    return prod;
  }

  friend Matrix operator*(Matrix const& lhs, Matrix const& rhs) {
    assert(lhs.width() == rhs.height());
    Matrix      prod(lhs.height(), rhs.width());
    auto const& mat = rhs.transpose();
    for (uint32_t r = 0; r < lhs.height(); r++)
      for (uint32_t c = 0; c < rhs.width(); c++) prod[r][c] = lhs[r] * mat[c];
    return prod;
  }

  friend std::istream& operator>>(std::istream& is, Matrix& mat) {
    for (auto& r : mat) { is >> r; }
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, Matrix const& mat) {
    for (auto const& r : mat) { os << r << "\n"; }
    return os;
  }

  Matrix transpose() const {
    Matrix mat(_width, _height);
    for (uint32_t r = 0; r < _height; r++)
      for (uint32_t c = 0; c < _width; c++) mat.at(c, r) = _data[r][c];
    return mat;
  }

  inline void resize(uint32_t h, uint32_t w) {
    _height = h;
    _width  = w;
    _data.resize(h);
    for (auto& r : _data) { r.resize(w); }
  }

  inline void fill(const T v) {
    for (uint32_t r = 0; r < _height; r++) _data[r] = v;
  }

  typename std::vector<VectorND<T>>::iterator       begin() { return _data.begin(); }
  typename std::vector<VectorND<T>>::iterator       end() { return _data.end(); }
  typename std::vector<VectorND<T>>::const_iterator begin() const { return _data.begin(); }
  typename std::vector<VectorND<T>>::const_iterator end() const { return _data.end(); }

  inline uint32_t                height() const { return _height; }
  inline uint32_t                width() const { return _width; }
  inline std::array<uint32_t, 2> dim() const { return std::array<uint32_t, 2>{_height, _width}; }

private:
  uint32_t                 _height;
  uint32_t                 _width;
  std::vector<VectorND<T>> _data;
};

//============================================================================

template<typename T>
class Mat3 : public Matrix<T> {
public:
  Mat3() : Matrix<T>(3, 3) {}
  Mat3(std::initializer_list<T> elements) : Matrix<T>(3, 3, elements) {}

  T determinant() const {
    return this->at(0, 0) * this->at(1, 1) * this->at(2, 2) + this->at(0, 1) * this->at(1, 2) * this->at(2, 0)
           + this->at(0, 2) * this->at(1, 0) * this->at(2, 1) - this->at(0, 2) * this->at(1, 1) * this->at(2, 0)
           - this->at(0, 1) * this->at(1, 0) * this->at(2, 2) - this->at(0, 0) * this->at(1, 2) * this->at(2, 1);
  }

  bool inverse(Mat3<T>& mat) const {
    T det = determinant();
    if (det == 0) return false;
    T invDet     = T(1.0) / det;
    mat.at(0, 0) = (this->at(1, 1) * this->at(2, 2) - this->at(1, 2) * this->at(2, 1)) * invDet;
    mat.at(0, 1) = (this->at(0, 2) * this->at(2, 1) - this->at(0, 1) * this->at(2, 2)) * invDet;
    mat.at(0, 2) = (this->at(0, 1) * this->at(1, 2) - this->at(0, 2) * this->at(1, 1)) * invDet;
    mat.at(1, 0) = (this->at(1, 2) * this->at(2, 0) - this->at(1, 0) * this->at(2, 2)) * invDet;
    mat.at(1, 1) = (this->at(0, 0) * this->at(2, 2) - this->at(0, 2) * this->at(2, 0)) * invDet;
    mat.at(1, 2) = (this->at(0, 2) * this->at(1, 0) - this->at(0, 0) * this->at(1, 2)) * invDet;
    mat.at(2, 0) = (this->at(1, 0) * this->at(2, 1) - this->at(1, 1) * this->at(2, 0)) * invDet;
    mat.at(2, 1) = (this->at(0, 1) * this->at(2, 0) - this->at(0, 0) * this->at(2, 1)) * invDet;
    mat.at(2, 2) = (this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0)) * invDet;

    return true;
  }

  Vector3D<T> operator*(const Vector3D<T>& rhs) const {
    Vector3D<T> prod = {0, 0, 0};
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) { prod[i] += this->at(i, j) * rhs[j]; }
    }
    return prod;
  }
};

//============================================================================

//!    3x3 Matrix
template<typename T>
class Mat3Vmc {
public:
  T operator()(int32_t rowIndex, int32_t columnIndex) const { return data[rowIndex][columnIndex]; }

  T& operator()(int32_t rowIndex, int32_t columnIndex) { return data[rowIndex][columnIndex]; }

  T* operator[](const int32_t rowIndex) {
    assert(rowIndex < 3);
    return data[rowIndex];
  }

  const T* operator[](const int32_t rowIndex) const {
    assert(rowIndex < 3);
    return data[rowIndex];
  }

  int32_t columnCount() const { return 3; }
  int32_t rowCount() const { return 3; }

  Mat3Vmc& operator=(const Mat3Vmc& rhs) {
    memcpy(data, rhs.data, sizeof(data));
    return *this;
  }

  void operator+=(const Mat3Vmc& rhs) {
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] += rhs.data[i][j]; }
    }
  }

  void operator-=(const Mat3Vmc& rhs) {
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] -= rhs.data[i][j]; }
    }
  }

  void operator-=(const T a) {
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] -= a; }
    }
  }

  void operator+=(const T a) {
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] += a; }
    }
  }

  void operator/=(const T a) {
    assert(a != 0);
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] /= a; }
    }
  }

  void operator*=(const T a) {
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { this->data[i][j] *= a; }
    }
  }

  Vector3D<T> operator*(const Vector3D<T>& rhs) const {
    Vector3D<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      res[i] = 0;
      for (int32_t j = 0; j < 3; ++j) { res[i] += this->data[i][j] * rhs[j]; }
    }
    return res;
  }

  Mat3Vmc operator*(const Mat3Vmc& rhs) const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) {
        res.data[i][j] = 0;
        for (int32_t k = 0; k < 3; ++k) { res.data[i][j] += this->data[i][k] * rhs.data[k][j]; }
      }
    }
    return res;
  }

  Mat3Vmc operator+(const Mat3Vmc& rhs) const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = this->data[i][j] + rhs.data[i][j]; }
    }
    return res;
  }

  Mat3Vmc operator-(const Mat3Vmc& rhs) const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = this->data[i][j] - rhs.data[i][j]; }
    }
    return res;
  }

  Mat3Vmc operator-() const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = -this->data[i][j]; }
    }
    return res;
  }

  Mat3Vmc operator*(T rhs) const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = this->data[i][j] * rhs; }
    }
    return res;
  }

  Mat3Vmc operator/(T rhs) const {
    assert(rhs != 0);
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = this->data[i][j] / rhs; }
    }
    return res;
  }

  Mat3Vmc transpose() const {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = this->data[j][i]; }
    }
    return res;
  }

  T determinant() const {  // compute determinant
    const auto& m = (*this);
    const T det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
                  + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    return det;
  }
  bool inverse(Mat3Vmc<T>& mat) const {
    // compute determinant
    const auto& m = (*this);
    const T det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) - m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0))
                  + m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
    const T invDet = T(1) / det;
    if (std::isinf(invDet)) {  // det == 0
      return false;
    }

    // compute inverse
    mat(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invDet;
    mat(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invDet;
    mat(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invDet;
    mat(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invDet;
    mat(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invDet;
    mat(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invDet;
    mat(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invDet;
    mat(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invDet;
    mat(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invDet;
    return true;
  }

  Mat3Vmc() = default;

  Mat3Vmc(const T a) {
    for (auto& i : data) {
      for (int32_t j = 0; j < 3; ++j) { i[j] = a; }
    }
  }

  Mat3Vmc(const Mat3Vmc& rhs) { memcpy(data, rhs.data, sizeof(data)); }

  ~Mat3Vmc() = default;

  friend inline Mat3Vmc<T> operator*(T lhs, const Mat3Vmc<T>& rhs) {
    Mat3Vmc<T> res;
    for (int32_t i = 0; i < 3; ++i) {
      for (int32_t j = 0; j < 3; ++j) { res.data[i][j] = lhs * rhs.data[i][j]; }
    }
    return res;
  }

  static void makeIdentity(Mat3Vmc<T>& mat) {
    memset(mat.data, 0, sizeof(mat.data));
    for (int32_t i = 0; i < 3; ++i) { mat[i][i] = T(1); }
  }

  static void makeScale(const T sx, const T sy, const T sz, Mat3Vmc<T>& mat) {
    makeIdentity(mat);
    mat[0][0] = sx;
    mat[1][1] = sy;
    mat[2][2] = sz;
  }

  static void makeUniformScale(const T s, Mat3Vmc<T>& mat) { makeScale(s, s, s, mat); }

  static void makeRotation(const T angle, const T ax, const T ay, const T az, Mat3Vmc<T>& mat) {
    T c       = cos(angle);
    T l_c     = 1 - c;
    T s       = sin(angle);
    mat[0][0] = ax * ax + (1 - ax * ax) * c;
    mat[0][1] = ax * ay * l_c - az * s;
    mat[0][2] = ax * az * l_c + ay * s;
    mat[1][0] = ax * ay * l_c + az * s;
    mat[1][1] = ay * ay + (1 - ay * ay) * c;
    mat[1][2] = ay * az * l_c - ax * s;
    mat[2][0] = ax * az * l_c - ay * s;
    mat[2][1] = ay * az * l_c + ax * s;
    mat[2][2] = az * az + (1 - az * az) * c;
  }

  static void makeRotationX(const T angle, Mat3Vmc<T>& mat) {
    T c       = cos(angle);
    T s       = sin(angle);
    mat[0][0] = T(1);
    mat[0][1] = T(0);
    mat[0][2] = T(0);
    mat[1][0] = T(0);
    mat[1][1] = c;
    mat[1][2] = -s;
    mat[2][0] = T(0);
    mat[2][1] = s;
    mat[2][2] = c;
  }

  static void makeRotationY(const T angle, Mat3Vmc<T>& mat) {
    T c       = cos(angle);
    T s       = sin(angle);
    mat[0][0] = c;
    mat[0][1] = T(0);
    mat[0][2] = s;
    mat[1][0] = T(0);
    mat[1][1] = T(1);
    mat[1][2] = T(0);
    mat[2][0] = -s;
    mat[2][1] = T(0);
    mat[2][2] = c;
  }

  static void makeRotationZ(const T angle, Mat3Vmc<T>& mat) {
    T c       = cos(angle);
    T s       = sin(angle);
    mat[0][0] = c;
    mat[0][1] = -s;
    mat[0][2] = T(0);
    mat[1][0] = s;
    mat[1][1] = c;
    mat[1][2] = T(0);
    mat[2][0] = T(0);
    mat[2][1] = T(0);
    mat[2][2] = T(1);
  }

  static void makeRotationEulerAngles(const T angleX, const T angleY, const T angleZ, Mat3Vmc<T>& mat) {
    Mat3<T> Rx;
    Mat3<T> Ry;
    Mat3<T> Rz;
    makeRotationX(angleX, Rx);
    makeRotationY(angleY, Ry);
    makeRotationZ(angleZ, Rz);
    mat = Rz * Ry * Rx;
  }

  friend std::ostream& operator<<(std::ostream& os, const Mat3Vmc<T>& mat) {
    os << mat[0][0] << ' ' << mat[0][1] << ' ' << mat[0][2] << '\n'
       << mat[1][0] << ' ' << mat[1][1] << ' ' << mat[1][2] << '\n'
       << mat[2][0] << ' ' << mat[2][1] << ' ' << mat[2][2] << '\n';
    return os;
  }

  friend std::istream& operator>>(std::istream& is, Mat3Vmc<T>& mat) {
    is >> mat[0][0] >> mat[0][1] >> mat[0][2] >> mat[1][0] >> mat[1][1] >> mat[1][2] >> mat[2][0] >> mat[2][1]
      >> mat[2][2];
    return is;
  }

private:
  T data[3][3];
};