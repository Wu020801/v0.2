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
#include <string>
#include <fstream>
#include <array>
#include <filesystem>
#include "vector.hpp"

//============================================================================

#ifdef WIN32
#  include <windows.h>
static void
mkdir(const char* pDirectory, int) {
  CreateDirectory(pDirectory, NULL);
}
#else
#  include <sys/dir.h>
#endif

//============================================================================

static inline char
separator() {
#ifndef _WIN32
  return '/';
#else
  return '\\';
#endif
}

//============================================================================

static char
separator(const std::string& string) {
  auto pos1 = string.find_last_of('/');
  auto pos2 = string.find_last_of('\\');
  auto pos  = (std::max)(pos1 != std::string::npos ? pos1 : 0, pos2 != std::string::npos ? pos2 : 0);
  return (pos != 0 ? string[pos] : separator());
}

//============================================================================

static std::string
dirname(const std::string& string) {
  auto position = string.find_last_of(separator(string));
  if (position != std::string::npos) { return string.substr(0, position + 1); }
  return string;
}

//============================================================================

static inline std::string
basename(const std::string& string) {
  auto position = string.find_last_of(separator(string));
  if (position != std::string::npos) { return string.substr(position + 1, string.length()); }
  return string;
}

//============================================================================

static inline std::string
extension(const std::string& name) {
  const auto pos = name.find_last_of('.');
  return pos == std::string::npos ? std::string("") : name.substr(pos + 1);
}

//============================================================================

static inline std::string
removeExtension(const std::string& name) {
  const auto pos = name.find_last_of('.');
  return pos == std::string::npos ? name : name.substr(0, pos);
}

//============================================================================

static int
save(const std::string& filename, std::vector<uint8_t>& buffer) {
  std::ofstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    printf("Can not save: %s \n", filename.c_str());
    return -1;
  }
  file.write(reinterpret_cast<const char*>(buffer.data()), buffer.size());
  file.close();
  return 0;
}

//============================================================================

static inline std::string
fillInNum(const std::string& str, int num) {
  char res[1000];
  snprintf(res, 1000, str.c_str(), num);
  return std::string(res);
}

//============================================================================

template<typename T>
std::array<T, 2>
boundingBox(const std::vector<T>& vectors) {
  std::array<T, 2> bbox  = {vectors[0], vectors[0]};
  auto&            bbMin = bbox[0];
  auto&            bbMax = bbox[1];

  uint32_t dim = vectors[0].size();
  for (const auto& p : vectors) {
    for (uint32_t d = 0; d < dim; d++) {
      const auto& x = p[d];
      if (x < bbMin[d]) bbMin[d] = x;
      if (x > bbMax[d]) bbMax[d] = x;
    }
  }
  return bbox;
}

//============================================================================

template<typename T>
void
computeLocalCoordinateSystem(Vector3D<T> const& n, Vector3D<T>& t, Vector3D<T>& b) {
  Vector3D<T> const e0(1, 0, 0);
  Vector3D<T> const e1(0, 1, 0);
  Vector3D<T> const e2(0, 0, 1);
  auto const        p0 = e0 * n;
  auto const        p1 = e1 * n;
  auto const        p2 = e2 * n;
  auto const        a0 = std::abs(p0);
  auto const        a1 = std::abs(p1);
  auto const        a2 = std::abs(p2);
  if (a0 <= a1 && a0 <= a2) {
    t = e0 - p0 * n;
  } else if (a1 <= a2) {
    t = e1 - p1 * n;
  } else {
    t = e2 - p2 * n;
  }
  t.normalize();
  b = n ^ t;
}

