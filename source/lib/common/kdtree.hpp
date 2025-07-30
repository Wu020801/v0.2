/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once
#include <nanoflann/nanoflann.hpp>
#include <nanoflann/KDTreeVectorOfVectorsAdaptor.h>
#include <vector>
#include "vector.hpp"

//============================================================================

struct metric_L2_vvm {
  template<class T, class DataSource>
  struct traits {
    using distance_t = nanoflann::L2_Adaptor<T, DataSource, T>;
  };
};

//============================================================================

template<class T>
using KdTree = KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D<T>>, T, 3, metric_L2_vvm, uint32_t>;

//============================================================================
