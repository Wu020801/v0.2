/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include "matrix.hpp"
#include "misc.hpp"
#include "vector.hpp"
#include "sparseMatrix.hpp"

//============================================================================

template<typename T>
class Image {
public:
  Image()  = default;
  ~Image() = default;

  template<typename S>
  Image(const Image<S>& img) {
    *this = img;
  }
  Image(uint32_t h, uint32_t w, uint32_t d = 3) { resize(h, w, d); }

  inline void resize(uint32_t h, uint32_t w, uint32_t d = 3) {
    _height = h;
    _width  = w;
    _depth  = d;
    _data.resize(h);
    for (auto& row : _data) {
      row.resize(w);
      for (auto& pix : row) { pix.resize(d); }
    }
  }

  template<typename S>
  Image& operator=(const Image<S>& img) {
    resize(img.height(), img.width(), img.depth());
    for (uint32_t i = 0; i < _height; i++)
      for (uint32_t j = 0; j < _width; j++) _data[i][j] = img.get(i, j);
    return *this;
  }

  template<typename S>
  Image& operator=(const std::vector<std::vector<S>>& img) {
    resize(img.size(), img[0].size(), 1);
    for (uint32_t i = 0; i < _height; i++)
      for (uint32_t j = 0; j < _width; j++) _data[i][j][0] = img[i][j];
    return *this;
  }

  inline void fill(const T v) {
    for (uint32_t i = 0; i < _height; i++)
      for (uint32_t j = 0; j < _width; j++)
        for (uint32_t k = 0; k < _depth; k++) _data[i][j][k] = v;
  }

  inline void set(const uint32_t r, const uint32_t c, const VectorND<T>& p) { _data[r][c] = p; }

  template<typename S>
  inline void set(const uint32_t r, const uint32_t c, const VectorND<S>& p) {
    _data[r][c] = p;
  }

  inline void set(const uint32_t r, const uint32_t c, const uint32_t d, const T p) { _data[r][c][d] = p; }

  inline const VectorND<T>& get(const uint32_t r, const uint32_t c) const { return _data[r][c]; }

  inline VectorND<T>& get(const uint32_t r, const uint32_t c) { return _data[r][c]; }

  inline const T& get(const uint32_t r, const uint32_t c, const uint32_t d) const { return _data[r][c][d]; }

  inline T& get(const uint32_t r, const uint32_t c, const uint32_t d) { return _data[r][c][d]; }

  VectorND<double> interp(const double x, const double y) const;

  bool read(const std::string& filePath);
  bool write(const std::string& filePath, const uint32_t quality = 100) const;

  inline uint32_t                height() const { return _height; }
  inline uint32_t                width() const { return _width; }
  inline uint32_t                depth() const { return _depth; }
  inline std::array<uint32_t, 3> dim() const { return std::array<uint32_t, 3>{_height, _width, _depth}; }

private:
  uint32_t                              _height;
  uint32_t                              _width;
  uint32_t                              _depth = 3;
  std::vector<std::vector<VectorND<T>>> _data;
};

//============================================================================

template<typename T>
class Image1 {
public:
  Image1()  = default;
  ~Image1() = default;

  template<typename S>
  Image1(const Image1<S>& img) {
    *this = img;
  }

  Image1(uint32_t h, uint32_t w) { resize(h, w); }

  inline void resize(uint32_t h, uint32_t w) {
    _height = h;
    _width  = w;
    _data.resize(h);
    for (auto& row : _data) { row.resize(w); }
  }

  template<typename S>
  Image1& operator=(const Image1<S>& img) {
    resize(img.height(), img.width());
    for (uint32_t i = 0; i < _height; ++i)
      for (uint32_t j = 0; j < _width; ++j) _data[i][j] = img.get(i, j);
    return *this;
  }

  template<typename S>
  Image1& operator=(const std::vector<std::vector<S>>& img) {
    resize(img.size(), img[0].size());
    for (uint32_t i = 0; i < _height; ++i)
      for (uint32_t j = 0; j < _width; ++j) _data[i][j] = img[i][j];
    return *this;
  }

  inline void fill(const T v) {
    for (uint32_t i = 0; i < _height; ++i)
      for (uint32_t j = 0; j < _width; ++j) _data[i][j] = v;
  }

  inline void set(const uint32_t r, const uint32_t c, const T p) { _data[r][c] = p; }

  inline const T& get(const uint32_t r, const uint32_t c) const { return _data[r][c]; }

  inline T& get(const uint32_t r, const uint32_t c) { return _data[r][c]; }

  double interp(const double x, const double y) const;
  bool   read(const std::string& filePath);
  bool   write(const std::string& filePath, const uint32_t quality = 100) const;

  inline uint32_t height() const { return _height; }
  inline uint32_t width() const { return _width; }

private:
  uint32_t                    _height;
  uint32_t                    _width;
  std::vector<std::vector<T>> _data;
};

//============================================================================

template<typename T>
class Image3 {
public:
  Image3()  = default;
  ~Image3() = default;

  template<typename S>
  Image3(const Image3<S>& img) {
    *this = img;
  }
  Image3(uint32_t h, uint32_t w) { resize(h, w); }

  template<typename S>
  Image3& operator=(const Image3<S>& img) {
    resize(img.height(), img.width());
    for (uint32_t i = 0; i < _height; ++i)
      for (uint32_t j = 0; j < _width; ++j) _data[i][j] = img.get(i, j);
    return *this;
  }

  inline void resize(uint32_t h, uint32_t w) {
    _height = h;
    _width  = w;
    _data.resize(h);
    for (auto& row : _data) { row.resize(w); }
  }

  inline void fill(const T v) {
    for (uint32_t i = 0; i < _height; ++i)
      for (uint32_t j = 0; j < _width; ++j) _data[i][j] = v;
  }

  inline void set(const uint32_t r, const uint32_t c, const Vector3D<T>& p) { _data[r][c] = p; }

  template<typename S>
  inline void set(const uint32_t r, const uint32_t c, const Vector3D<S>& p) {
    _data[r][c] = p;
  }

  inline void set(const uint32_t r, const uint32_t c, const uint32_t d, const T p) { _data[r][c][d] = p; }

  template<typename S>
  inline void set(const uint32_t r, const uint32_t c, const uint32_t d, const S p) {
    _data[r][c][d] = p;
  }

  inline const Vector3D<T>& get(const uint32_t r, const uint32_t c) const { return _data[r][c]; }

  inline Vector3D<T>& get(const uint32_t r, const uint32_t c) { return _data[r][c]; }

  inline const T& get(const uint32_t r, const uint32_t c, const uint32_t d) const { return _data[r][c][d]; }

  inline T& get(const uint32_t r, const uint32_t c, const uint32_t d) { return _data[r][c][d]; }

  Vector3D<double> interp(const double x, const double y) const;
  bool             readImage(const std::string& filePath, uint32_t inputWidth, uint32_t inputHeight);
  bool             writeImage(const std::string& filePath, const uint32_t quality = 100) const;

  inline uint32_t height() const { return _height; }
  inline uint32_t width() const { return _width; }

  void rgb2yuv420(std::vector<unsigned char>& y_plane, std::vector<unsigned char>& u_plane, 
			      std::vector<unsigned char>& v_plane);

  void dilatePadding(std::vector<std::vector<uint8_t>>& occ, const uint32_t iter);
  void pullPushPadding(std::vector<std::vector<uint8_t>>& occ);
  void sparseLinearPadding(std::vector<std::vector<uint8_t>>& occ,
                           const double                       maxErr       = 0.001,
                           const int32_t                      maxIterCount = -1);

private:
  uint32_t                              _height;
  uint32_t                              _width;
  std::vector<std::vector<Vector3D<T>>> _data;
};

