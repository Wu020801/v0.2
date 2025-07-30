/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <istream>
#include <ostream>
#include <vector>
#include <array>
#if defined(WIN32) || defined(_WIN32)
#  define NOMINMAX
#  include <Windows.h>
#endif

//============================================================================

template<typename T>
class Vector2D {
public:
  Vector2D()  = default;
  ~Vector2D() = default;

  Vector2D(const T x) {
    _data[0] = x;
    _data[1] = x;
  }

  Vector2D(const T x, const T y) {
    _data[0] = x;
    _data[1] = y;
  }

  Vector2D(const Vector2D& v) {
    _data[0] = v[0];
    _data[1] = v[1];
  }

  template<typename S>
  Vector2D(const Vector2D<S>& v) {
    _data[0] = (T)v[0];
    _data[1] = (T)v[1];
  }

  Vector2D(const T* v) : Vector2D(v[0], v[1]){};

  T& operator[](const uint32_t i) {
    assert(i < 2);
    return _data[i];
  }

  const T& operator[](const uint32_t i) const {
    assert(i < 2);
    return _data[i];
  }

  Vector2D& operator=(const T rhs) {
    _data[0] = rhs;
    _data[1] = rhs;
    return *this;
  }

  template<typename S>
  Vector2D& operator=(const S rhs) {
    T x      = (T)rhs;
    _data[0] = x;
    _data[1] = x;
    return *this;
  }

  Vector2D& operator=(const Vector2D& rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    return *this;
  }

  template<typename S>
  Vector2D& operator=(const Vector2D<S>& rhs) {
    _data[0] = (T)rhs[0];
    _data[1] = (T)rhs[1];
    return *this;
  }

  Vector2D& operator=(const T* const rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    return *this;
  }

  Vector2D& operator+=(const T rhs) {
    _data[0] += rhs;
    _data[1] += rhs;
    return *this;
  }

  Vector2D& operator-=(const T rhs) {
    _data[0] -= rhs;
    _data[1] -= rhs;
    return *this;
  }

  Vector2D& operator*=(const T rhs) {
    _data[0] *= rhs;
    _data[1] *= rhs;
    return *this;
  }

  Vector2D& operator/=(const T rhs) {
    assert(rhs != 0);
    _data[0] /= rhs;
    _data[1] /= rhs;
    return *this;
  }

  Vector2D& operator+=(const Vector2D& rhs) {
    _data[0] += rhs[0];
    _data[1] += rhs[1];
    return *this;
  }

  Vector2D& operator-=(const Vector2D& rhs) {
    _data[0] -= rhs[0];
    _data[1] -= rhs[1];
    return *this;
  }

  bool operator==(const Vector2D& rhs) const { return (_data[0] == rhs[0] && _data[1] == rhs[1]); }

  bool operator!=(const Vector2D& rhs) const { return (_data[0] != rhs[0] || _data[1] != rhs[1]); }

  bool operator>(const Vector2D& rhs) const {
    if (_data[0] == rhs[0]) return _data[1] > rhs[1];
    return _data[0] > rhs[0];
  }

  bool operator<(const Vector2D& rhs) const {
    if (_data[0] == rhs[0]) return _data[1] < rhs[1];
    return _data[0] < rhs[0];
  }

  Vector2D operator-() const { return Vector2D(-_data[0], -_data[1]); }

  friend Vector2D operator+(const Vector2D& lhs, const Vector2D& rhs) {
    return Vector2D(lhs[0] + rhs[0], lhs[1] + rhs[1]);
  }

  friend Vector2D operator+(const T lhs, const Vector2D& rhs) { return Vector2D(lhs + rhs[0], lhs + rhs[1]); }

  friend Vector2D operator+(const Vector2D& lhs, const T rhs) { return Vector2D(lhs[0] + rhs, lhs[1] + rhs); }

  friend Vector2D operator-(const Vector2D& lhs, const Vector2D& rhs) {
    return Vector2D(lhs[0] - rhs[0], lhs[1] - rhs[1]);
  }

  friend Vector2D operator-(const T lhs, const Vector2D& rhs) { return Vector2D(lhs - rhs[0], lhs - rhs[1]); }

  friend Vector2D operator-(const Vector2D& lhs, const T rhs) { return Vector2D(lhs[0] - rhs, lhs[1] - rhs); }

  friend T operator^(const Vector2D& lhs, const Vector2D& rhs) { return lhs[0] * rhs[1] - lhs[1] * rhs[0]; }

  friend T operator*(const Vector2D& lhs, const Vector2D& rhs) { return lhs[0] * rhs[0] + lhs[1] * rhs[1]; }

  friend Vector2D operator*(const T lhs, const Vector2D& rhs) { return Vector2D(lhs * rhs[0], lhs * rhs[1]); }

  friend Vector2D operator*(const Vector2D& lhs, const T rhs) { return Vector2D(lhs[0] * rhs, lhs[1] * rhs); }

  friend Vector2D operator/(const Vector2D& lhs, const T rhs) {
    assert(rhs != 0);
    return Vector2D(lhs[0] / rhs, lhs[1] / rhs);
  }

  friend Vector2D operator/(const T lhs, const Vector2D& rhs) {
    assert(rhs.L0norm() == 2);
    return Vector2D(lhs / rhs[0], lhs / rhs[1]);
  }

  friend std::istream& operator>>(std::istream& s, Vector2D& v) {
    s >> v[0] >> v[1];
    return s;
  }

  friend std::ostream& operator<<(std::ostream& s, const Vector2D& v) {
    s << v[0] << ' ' << v[1];
    return s;
  }

  uint32_t L0norm() const { return uint32_t(_data[0] != T(0)) + uint32_t(_data[1] != T(0)); }
  T        L1norm() const { return std::abs(_data[0]) + std::abs(_data[1]); }
  T        L2norm() const { return std::sqrt(L2normSq()); }
  T        L2normSq() const { return (*this) * (*this); }
  T        Linfnorm() const { return std::max(std::abs(_data[0]), std::abs(_data[1])); }

  T norm() const { return L2norm(); }
  T normSq() const { return L2normSq(); }

  void normalize() {
    const T n = norm();
    if (n != T(0)) { (*this) /= n; }
  }

  Vector2D round() const {
    Vector2D v;
    v[0] = std::round(_data[0]);
    v[1] = std::round(_data[1]);
    return v;
  }

  uint32_t size() const { return 2; }

  typedef T type;

  template<typename S>
  T cast(const S x) {
    return T(x);
  }

  T*       data() { return _data; }
  const T* data() const { return _data; }

private:
  T _data[2];
};

//============================================================================

template<typename T>
class Vector3D {
public:
  Vector3D()  = default;
  ~Vector3D() = default;

  Vector3D(const T x) {
    _data[0] = x;
    _data[1] = x;
    _data[2] = x;
  }

  Vector3D(const T x, const T y, const T z) {
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
  }

  Vector3D(const Vector3D& v) {
    _data[0] = v[0];
    _data[1] = v[1];
    _data[2] = v[2];
  }

  template<typename S>
  Vector3D(const Vector3D<S>& v) {
    _data[0] = (T)v[0];
    _data[1] = (T)v[1];
    _data[2] = (T)v[2];
  }

  Vector3D(const T* v) : Vector3D(v[0], v[1], v[2]){};

  T*       begin() { return &_data[0]; }
  const T* begin() const { return &_data[0]; }

  T*       end() { return &_data[2]; }
  const T* end() const { return &_data[2]; }

  T& operator[](const uint32_t i) {
    assert(i < 3);
    return _data[i];
  }

  const T& operator[](const uint32_t i) const {
    assert(i < 3);
    return _data[i];
  }

  Vector3D& operator=(const T rhs) {
    _data[0] = rhs;
    _data[1] = rhs;
    _data[2] = rhs;
    return *this;
  }

  Vector3D& operator=(const Vector3D& rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    _data[2] = rhs[2];
    return *this;
  }

  template<typename S>
  Vector3D& operator=(const S rhs) {
    T x      = static_cast<T>(rhs);
    _data[0] = x;
    _data[1] = x;
    _data[2] = x;
    return *this;
  }

  template<typename S>
  Vector3D& operator=(const Vector3D<S>& rhs) {
    _data[0] = (T)rhs[0];
    _data[1] = (T)rhs[1];
    _data[2] = (T)rhs[2];
    return *this;
  }

  Vector3D& operator=(const T* const rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    _data[2] = rhs[2];
    return *this;
  }

  Vector3D& operator+=(const T rhs) {
    _data[0] += rhs;
    _data[1] += rhs;
    _data[2] += rhs;
    return *this;
  }

  Vector3D& operator-=(const T rhs) {
    _data[0] -= rhs;
    _data[1] -= rhs;
    _data[2] -= rhs;
    return *this;
  }

  Vector3D& operator*=(const T rhs) {
    _data[0] *= rhs;
    _data[1] *= rhs;
    _data[2] *= rhs;
    return *this;
  }

  Vector3D& operator/=(const T rhs) {
    assert(rhs != 0);
    _data[0] /= rhs;
    _data[1] /= rhs;
    _data[2] /= rhs;
    return *this;
  }

  Vector3D& operator+=(const Vector3D& rhs) {
    _data[0] += rhs[0];
    _data[1] += rhs[1];
    _data[2] += rhs[2];
    return *this;
  }

  Vector3D& operator-=(const Vector3D& rhs) {
    _data[0] -= rhs[0];
    _data[1] -= rhs[1];
    _data[2] -= rhs[2];
    return *this;
  }

  bool operator==(const Vector3D& rhs) const {
    return (_data[0] == rhs[0] && _data[1] == rhs[1] && _data[2] == rhs[2]);
  }

  bool operator!=(const Vector3D& rhs) const {
    return (_data[0] != rhs[0] || _data[1] != rhs[1] || _data[2] != rhs[2]);
  }

  bool operator>(const Vector3D& rhs) const {
    for (uint32_t i = 0; i < 3; i++) {
      if (_data[i] != rhs[i]) return _data[i] > rhs[i];
    }
    return false;
  }

  bool operator<(const Vector3D& rhs) const {
    for (uint32_t i = 0; i < 3; i++) {
      if (_data[i] != rhs[i]) return _data[i] < rhs[i];
    }
    return false;
  }

  Vector3D operator-() const { return Vector3D(-_data[0], -_data[1], -_data[2]); }

  friend Vector3D operator+(const Vector3D& lhs, const Vector3D& rhs) {
    return Vector3D(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
  }

  friend Vector3D operator+(const T lhs, const Vector3D& rhs) {
    return Vector3D(lhs + rhs[0], lhs + rhs[1], lhs + rhs[2]);
  }

  friend Vector3D operator+(const Vector3D& lhs, const T rhs) {
    return Vector3D(lhs[0] + rhs, lhs[1] + rhs, lhs[2] + rhs);
  }

  friend Vector3D operator-(const Vector3D& lhs, const Vector3D& rhs) {
    return Vector3D(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
  }

  friend Vector3D operator-(const T lhs, const Vector3D& rhs) {
    return Vector3D(lhs - rhs[0], lhs - rhs[1], lhs - rhs[2]);
  }

  friend Vector3D operator-(const Vector3D& lhs, const T rhs) {
    return Vector3D(lhs[0] - rhs, lhs[1] - rhs, lhs[2] - rhs);
  }

  friend Vector3D operator^(const Vector3D& lhs, const Vector3D& rhs) {
    return Vector3D(
      lhs[1] * rhs[2] - lhs[2] * rhs[1], lhs[2] * rhs[0] - lhs[0] * rhs[2], lhs[0] * rhs[1] - lhs[1] * rhs[0]);
  }

  friend T operator*(const Vector3D& lhs, const Vector3D& rhs) {
    return (lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]);
  }

  friend Vector3D operator*(const T lhs, const Vector3D& rhs) {
    return Vector3D(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2]);
  }

  friend Vector3D operator*(const Vector3D& lhs, const T rhs) {
    return Vector3D(lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs);
  }

  friend Vector3D operator/(const Vector3D& lhs, const T rhs) {
    assert(rhs != 0);
    return Vector3D(lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs);
  }

  friend Vector3D operator/(const T lhs, const Vector3D& rhs) {
    assert(rhs.L0norm() == 3);
    return Vector3D(lhs / rhs[0], lhs / rhs[1], lhs / rhs[2]);
  }

  friend std::istream& operator>>(std::istream& s, Vector3D& v) {
    s >> v[0] >> v[1] >> v[2];
    return s;
  }

  friend std::ostream& operator<<(std::ostream& s, const Vector3D& v) {
    s << v[0] << ' ' << v[1] << ' ' << v[2];
    return s;
  }

  uint32_t L0norm() const {
    return uint32_t(_data[0] != T(0)) + uint32_t(_data[1] != T(0)) + uint32_t(_data[2] != T(0));
  }

  T L1norm() const { return std::abs(_data[0]) + std::abs(_data[1]) + std::abs(_data[2]); }

  T L2norm() const { return std::sqrt(L2normSq()); }
  T L2normSq() const { return (*this) * (*this); }

  T Linfnorm() const { return std::max(std::max(std::abs(_data[0]), std::abs(_data[1])), std::abs(_data[2])); }

  T norm() const { return L2norm(); }
  T normSq() const { return L2normSq(); }

  void normalize1() {
    const T n = norm();
    if (n != T(0)) { (*this) /= n; }
  }

  void normalize() {
    const T n2 = normSq();
    if (n2 != T(0)) {
      T invNorm = static_cast<T>(1) / std::sqrt(n2);
      (*this) *= invNorm;
    }
  }

  Vector3D round() const {
    Vector3D v;
    v[0] = std::round(_data[0]);
    v[1] = std::round(_data[1]);
    v[2] = std::round(_data[2]);
    return v;
  }

  uint32_t size() const { return 3; }

  typedef T type;

  template<typename S>
  T cast(const S x) {
    return T(x);
  }

  T*       data() { return _data; }
  const T* data() const { return _data; }

private:
  T _data[3];
};

//============================================================================

template<typename T>
class Vector4D {
public:
  Vector4D()  = default;
  ~Vector4D() = default;

  Vector4D(const T x) {
    _data[0] = x;
    _data[1] = x;
    _data[2] = x;
    _data[3] = x;
  }

  Vector4D(const T x, const T y, const T z, const T w) {
    _data[0] = x;
    _data[1] = y;
    _data[2] = z;
    _data[3] = w;
  }

  Vector4D(const Vector4D& v) {
    _data[0] = v[0];
    _data[1] = v[1];
    _data[2] = v[2];
    _data[3] = v[3];
  }

  template<typename S>
  Vector4D(const Vector4D<S>& v) {
    _data[0] = (T)v[0];
    _data[1] = (T)v[1];
    _data[2] = (T)v[2];
    _data[3] = (T)v[3];
  }

  Vector4D(const T* v) : Vector4D(v[0], v[1], v[2], v[3]){};

  T& operator[](const uint32_t i) {
    assert(i < 4);
    return _data[i];
  }

  const T& operator[](const uint32_t i) const {
    assert(i < 4);
    return _data[i];
  }

  Vector4D& operator=(const T rhs) {
    _data[0] = rhs;
    _data[1] = rhs;
    _data[2] = rhs;
    _data[3] = rhs;
    return *this;
  }

  template<typename S>
  Vector4D& operator=(const S rhs) {
    T x      = (T)rhs;
    _data[0] = x;
    _data[1] = x;
    _data[2] = x;
    _data[3] = x;
    return *this;
  }

  Vector4D& operator=(const Vector4D& rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    _data[2] = rhs[2];
    _data[3] = rhs[3];
    return *this;
  }

  template<typename S>
  Vector4D& operator=(const Vector4D<S>& rhs) {
    _data[0] = (T)rhs[0];
    _data[1] = (T)rhs[1];
    _data[2] = (T)rhs[2];
    _data[3] = (T)rhs[3];
    return *this;
  }

  Vector4D& operator=(const T* const rhs) {
    _data[0] = rhs[0];
    _data[1] = rhs[1];
    _data[2] = rhs[2];
    _data[3] = rhs[3];
    return *this;
  }

  Vector4D& operator+=(const T rhs) {
    _data[0] += rhs;
    _data[1] += rhs;
    _data[2] += rhs;
    _data[3] += rhs;
    return *this;
  }

  Vector4D& operator-=(const T rhs) {
    _data[0] -= rhs;
    _data[1] -= rhs;
    _data[2] -= rhs;
    _data[3] -= rhs;
    return *this;
  }

  Vector4D& operator*=(const T rhs) {
    _data[0] *= rhs;
    _data[1] *= rhs;
    _data[2] *= rhs;
    _data[3] *= rhs;
    return *this;
  }

  Vector4D& operator/=(const T rhs) {
    assert(rhs != 0);
    _data[0] /= rhs;
    _data[1] /= rhs;
    _data[2] /= rhs;
    _data[3] /= rhs;
    return *this;
  }

  Vector4D& operator+=(const Vector4D& rhs) {
    _data[0] += rhs[0];
    _data[1] += rhs[1];
    _data[2] += rhs[2];
    _data[3] += rhs[3];
    return *this;
  }

  Vector4D& operator-=(const Vector4D& rhs) {
    _data[0] -= rhs[0];
    _data[1] -= rhs[1];
    _data[2] -= rhs[2];
    _data[3] -= rhs[3];
    return *this;
  }

  bool operator==(const Vector4D& rhs) const {
    return (_data[0] == rhs[0] && _data[1] == rhs[1] && _data[2] == rhs[2] && _data[3] == rhs[3]);
  }

  bool operator!=(const Vector4D& rhs) const {
    return (_data[0] != rhs[0] || _data[1] != rhs[1] || _data[2] != rhs[2] || _data[3] != rhs[3]);
  }

  bool operator>(const Vector4D& rhs) const {
    for (uint32_t i = 0; i < 4; i++) {
      if (_data[i] != rhs[i]) return _data[i] > rhs[i];
    }
    return false;
  }

  bool operator<(const Vector4D& rhs) const {
    for (uint32_t i = 0; i < 4; i++) {
      if (_data[i] != rhs[i]) return _data[i] < rhs[i];
    }
    return false;
  }

  Vector4D operator-() const { return Vector4D(-_data[0], -_data[1], -_data[2], -_data[3]); }

  friend Vector4D operator+(const Vector4D& lhs, const Vector4D& rhs) {
    return Vector4D(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]);
  }

  friend Vector4D operator+(const T lhs, const Vector4D& rhs) {
    return Vector4D(lhs + rhs[0], lhs + rhs[1], lhs + rhs[2], lhs + rhs[3]);
  }

  friend Vector4D operator+(const Vector4D& lhs, const T rhs) {
    return Vector4D(lhs[0] + rhs, lhs[1] + rhs, lhs[2] + rhs, lhs[3] + rhs);
  }

  friend Vector4D operator-(const Vector4D& lhs, const Vector4D& rhs) {
    return Vector4D(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]);
  }

  friend Vector4D operator-(const T lhs, const Vector4D& rhs) {
    return Vector4D(lhs - rhs[0], lhs - rhs[1], lhs - rhs[2], lhs - rhs[3]);
  }

  friend Vector4D operator-(const Vector4D& lhs, const T rhs) {
    return Vector4D(lhs[0] - rhs, lhs[1] - rhs, lhs[2] - rhs, lhs[3] - rhs);
  }

  friend T operator*(const Vector4D& lhs, const Vector4D& rhs) {
    return (lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2] + lhs[3] * rhs[3]);
  }

  friend Vector4D operator*(const T lhs, const Vector4D& rhs) {
    return Vector4D(lhs * rhs[0], lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]);
  }

  friend Vector4D operator*(const Vector4D& lhs, const T rhs) {
    return Vector4D(lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs, lhs[3] * rhs);
  }

  friend Vector4D operator/(const Vector4D& lhs, const T rhs) {
    assert(rhs != 0);
    return Vector4D(lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs, lhs[3] / rhs);
  }

  friend Vector4D operator/(const T lhs, const Vector4D& rhs) {
    assert(rhs.L0norm() == 4);
    return Vector4D(lhs / rhs[0], lhs / rhs[1], lhs / rhs[2], lhs / rhs[3]);
  }

  friend std::istream& operator>>(std::istream& s, Vector4D& v) {
    s >> v[0] >> v[1] >> v[2] >> v[3];
    return s;
  }

  friend std::ostream& operator<<(std::ostream& s, const Vector4D& v) {
    s << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3];
    return s;
  }

  uint32_t L0norm() const {
    return uint32_t(_data[0] != T(0)) + uint32_t(_data[1] != T(0)) + uint32_t(_data[2] != T(0))
           + uint32_t(_data[3] != T(0));
  }

  T L1norm() const { return std::abs(_data[0]) + std::abs(_data[1]) + std::abs(_data[2]) + std::abs(_data[3]); }

  T L2norm() const { return std::sqrt(L2normSq()); }
  T L2normSq() const { return (*this) * (*this); }

  T Linfnorm() const {
    return std::max(std::max(std::max(std::abs(_data[0]), std::abs(_data[1])), std::abs(_data[2])),
                    std::abs(_data[3]));
  }

  T norm() const { return L2norm(); }
  T normSq() const { return L2normSq(); }

  void normalize() {
    const T n = norm();
    if (n != T(0)) { (*this) /= n; }
  }

  Vector4D round() const {
    Vector4D v;
    v[0] = std::round(_data[0]);
    v[1] = std::round(_data[1]);
    v[2] = std::round(_data[2]);
    v[3] = std::round(_data[3]);
    return v;
  }

  uint32_t size() const { return 4; }

  typedef T type;

  template<typename S>
  T cast(const S x) {
    return T(x);
  }

  T*       data() { return _data; }
  const T* data() const { return _data; }

private:
  T _data[4];
};

//============================================================================

template<typename T>
class VectorND {
public:
  VectorND()  = default;
  ~VectorND() = default;

  VectorND(const uint32_t n) {
    N = n;
    _data.resize(N);
  }

  VectorND(const uint32_t n, const T x) {
    N = n;
    _data.resize(N);
    for (uint32_t i = 0; i < N; i++) { _data[i] = x; }
  }

  VectorND(const VectorND& v) {
    N     = v.size();
    _data = v._data;
  }

  template<typename S>
  VectorND(const VectorND<S>& v) {
    N = v.size();
    _data.resize(N);

    for (uint32_t i = 0; i < N; i++) { _data[i] = (T)v[i]; }
  }

  VectorND(const T* v, const uint32_t n) {
    N = n;
    _data.resize(N);
    memcpy(_data.data(), v, n * sizeof(T));
  }

  T& operator[](const uint32_t i) {
    assert(i < N);
    return _data[i];
  }

  const T& operator[](const uint32_t i) const {
    assert(i < N);
    return _data[i];
  }

  VectorND& operator=(const T rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] = rhs; }
    return *this;
  }

  template<typename S>
  VectorND& operator=(const S rhs) {
    T x = (T)rhs;
    for (uint32_t i = 0; i < N; i++) { _data[i] = x; }
    return *this;
  }

  VectorND& operator=(const VectorND& rhs) {
    N     = rhs.size();
    _data = rhs._data;
    return *this;
  }

  template<typename S>
  VectorND& operator=(const VectorND<S>& rhs) {
    N = rhs.size();
    _data.resize(N);

    for (uint32_t i = 0; i < N; i++) { _data[i] = (T)rhs[i]; }
    return *this;
  }

  VectorND& operator+=(const T rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] += rhs; }
    return *this;
  }

  template<typename S>
  VectorND& operator+=(const S rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] += (T)rhs; }
    return *this;
  }

  VectorND& operator-=(const T rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] -= rhs; }
    return *this;
  }

  template<typename S>
  VectorND& operator-=(const S rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] -= (T)rhs; }
    return *this;
  }

  VectorND& operator*=(const T rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] *= rhs; }
    return *this;
  }

  template<typename S>
  VectorND& operator*=(const S rhs) {
    for (uint32_t i = 0; i < N; i++) { _data[i] *= (T)rhs; }
    return *this;
  }

  VectorND& operator/=(const T rhs) {
    assert(rhs != 0);
    for (uint32_t i = 0; i < N; i++) { _data[i] /= rhs; }
    return *this;
  }

  template<typename S>
  VectorND& operator/=(const S rhs) {
    assert(rhs != 0);
    for (uint32_t i = 0; i < N; i++) { _data[i] /= (T)rhs; }
    return *this;
  }

  VectorND& operator+=(const VectorND& rhs) {
    assert(N == rhs.size());
    for (uint32_t i = 0; i < N; i++) { _data[i] += rhs[i]; }
    return *this;
  }

  template<typename S>
  VectorND& operator+=(const VectorND<S>& rhs) {
    assert(N == rhs.size());
    for (uint32_t i = 0; i < N; i++) { _data[i] += T(rhs[i]); }
    return *this;
  }

  VectorND& operator-=(const VectorND& rhs) {
    assert(N == rhs.size());
    for (uint32_t i = 0; i < N; i++) { _data[i] -= rhs[i]; }
    return *this;
  }

  template<typename S>
  VectorND& operator-=(const VectorND<S>& rhs) {
    assert(N == rhs.size());
    for (uint32_t i = 0; i < N; i++) { _data[i] -= T(rhs[i]); }
    return *this;
  }

  bool operator==(const VectorND& rhs) const {
    if (N != rhs.size()) return false;
    for (uint32_t i = 0; i < N; i++) {
      if (_data[i] != rhs[i]) return false;
    }
    return true;
  }

  bool operator!=(const VectorND& rhs) const {
    if (N != rhs.size()) return true;
    for (uint32_t i = 0; i < N; i++) {
      if (_data[i] != rhs[i]) return true;
    }
    return false;
  }

  bool operator>(const VectorND& rhs) const {
    for (uint32_t i = 0; i < N; i++) {
      if (_data[i] != rhs[i]) return _data[i] > rhs[i];
    }
    return false;
  }

  bool operator<(const VectorND& rhs) const {
    for (uint32_t i = 0; i < N; i++) {
      if (_data[i] != rhs[i]) return _data[i] < rhs[i];
    }
    return false;
  }

  VectorND operator-() const {
    VectorND v(N);
    for (uint32_t i = 0; i < N; i++) { v[i] = -_data[i]; }
    return v;
  }

  friend VectorND operator+(const VectorND& lhs, const VectorND& rhs) {
    assert(lhs.size() == rhs.size());
    VectorND v(lhs);
    v += rhs;
    return v;
  }

  friend VectorND operator+(const T lhs, const VectorND& rhs) {
    VectorND v(rhs);
    v += lhs;
    return v;
  }

  friend VectorND operator+(const VectorND& lhs, const T rhs) {
    VectorND v(lhs);
    v += rhs;
    return v;
  }

  friend VectorND operator-(const VectorND& lhs, const VectorND& rhs) {
    assert(lhs.size() == rhs.size());
    VectorND v(lhs);
    v -= rhs;
    return v;
  }

  friend VectorND operator-(const T lhs, const VectorND& rhs) {
    VectorND v(rhs.size());
    for (uint32_t i = 0; i < v.size(); i++) { v[i] = lhs - rhs[i]; }
    return v;
  }

  friend VectorND operator-(const VectorND& lhs, const T rhs) {
    VectorND v(lhs);
    v -= rhs;
    return v;
  }

  friend T operator*(const VectorND& lhs, const VectorND& rhs) {
    assert(lhs.size() == rhs.size());
    T p = 0;
    for (uint32_t i = 0; i < lhs.size(); i++) { p += lhs[i] * rhs[i]; }
    return p;
  }

  friend VectorND operator*(const T lhs, const VectorND& rhs) {
    VectorND v(rhs);
    v *= lhs;
    return v;
  }

  friend VectorND operator*(const VectorND& lhs, const T rhs) {
    VectorND v(lhs);
    v *= rhs;
    return v;
  }

  friend VectorND operator/(const VectorND& lhs, const T rhs) {
    assert(rhs != 0);
    VectorND v(lhs);
    v /= rhs;
    return v;
  }

  friend VectorND operator/(const T lhs, const VectorND& rhs) {
    uint32_t n = rhs.size();
    assert(rhs.L0norm() == n);
    VectorND v(n);
    for (uint32_t i = 0; i < n; i++) { v[i] = lhs / rhs[i]; }
    return v;
  }

  friend std::istream& operator>>(std::istream& s, VectorND& v) {
    for (uint32_t i = 0; i < v.size(); i++) { s >> v[i]; }
    return s;
  }

  friend std::ostream& operator<<(std::ostream& s, const VectorND& v) {
    for (uint32_t i = 0; i < v.size(); i++) { s << v[i] << ' '; }
    return s;
  }

  uint32_t L0norm() const {
    uint32_t nnz = 0;
    for (uint32_t i = 0; i < N; i++) { nnz += uint32_t(_data[i] != T(0)); }
    return nnz;
  }

  T L1norm() const {
    T n = 0;
    for (uint32_t i = 0; i < N; i++) { n += T(std::abs(_data[i])); }
    return n;
  }

  T L2norm() const { return std::sqrt(L2normSq()); }
  T L2normSq() const { return (*this) * (*this); }

  T Linfnorm() const {
    T n = std::abs(_data[0]);
    for (uint32_t i = 1; i < N; i++) {
      T a = std::abs(_data[i]);
      if (a > n) n = a;
    }
    return n;
  }

  T norm() const { return L2norm(); }
  T normSq() const { return L2normSq(); }

  void normalize() {
    const T n = norm();
    if (n != T(0)) { (*this) /= n; }
  }

  VectorND round() const {
    VectorND v(N);
    for (uint32_t i = 0; i < N; i++) { v[i] = std::round(_data[i]); }
    return v;
  }

  uint32_t size() const { return N; }
  void     resize(const uint32_t n) {
    N = n;
    _data.resize(N);
  }
  void clear() {
    N = 0;
    _data.clear();
  }

  typedef T type;

  template<typename S>
  T cast(const S x) {
    return T(x);
  }

  typename std::vector<T>::iterator begin() { return _data.begin(); }
  typename std::vector<T>::iterator end() { return _data.end(); }

  T*       data() { return _data.data(); }
  const T* data() const { return _data.data(); }

private:
  uint32_t       N;
  std::vector<T> _data;
};

//============================================================================

template<typename T>
T
hadamard(const T& lhs, const T& rhs) {
  uint32_t n = lhs.size();
  assert(n == rhs.size());

  T v(lhs);
  for (uint32_t i = 0; i < n; i++) { v[i] *= rhs[i]; }
  return v;
}
