/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#pragma once

#include <cstdint>
#include <vector>
#include <array>
#include "mesh.hpp"
#include "image.hpp"

//============================================================================

template<typename T>
class TriMesh;
template<typename T>
class Image3;
struct encoderParams;

//============================================================================

class TransferColor {
public:
  TransferColor()  = default;
  ~TransferColor() = default;

  void textureTransfer(TriMesh<MeshType> const&    inputMesh,
                       Image3<uint8_t> const&      inputTexture,
                       TriMesh<MeshType>&          outputMesh,
                       Image3<uint8_t>&            outputTexture,
					   encoderParams&			   params);

private:

  void textureTransferByMeshes(TriMesh<MeshType> const&           inputMesh,
                               Image3<uint8_t> const&             inputTexture,
                               TriMesh<MeshType>&                 outputMesh,
                               Image3<uint8_t>&                   outputTexture,
                               std::vector<std::vector<uint8_t>>& ocm,
							   encoderParams const&				  params);
};
