/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include "misc.hpp"
#include "mesh.hpp"
// #include "util/verbose.hpp"
#include "vector.hpp"
// #include "vmc.hpp"

//============================================================================

struct encoderParams;

//============================================================================

class TextureParametrization {
public:
  TextureParametrization()  = default;
  ~TextureParametrization() = default;

  template<typename T>
  static bool generate(TriMesh<T> const& inputMesh, TriMesh<T>& outputMesh, encoderParams const& params);
};
