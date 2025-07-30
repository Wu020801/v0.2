
#include <chrono>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "encoder.hpp"
#include "textureParametrization.hpp"
#include <DirectXMath.h>
#include <DirectXMesh.h>
#include <UVAtlas.h>

//============================================================================

HRESULT __cdecl UVAtlasCallback(float fPercentDone) {
  // static auto prev = std::chrono::steady_clock::now();
  // const auto  tick = std::chrono::steady_clock::now();
  // if (tick - prev > std::chrono::seconds(1)) {
  //   std::cout << fPercentDone * 100. << "%   \r" << std::flush;
  //   prev = tick;
  // }
  return S_OK;
}

template<typename T>
bool
TextureParametrization::generate(TriMesh<T> const&           inputMesh,
                                 TriMesh<T>&                 outputMesh,
								 encoderParams const& params) {
  TriMesh<float> mesh(inputMesh);
  mesh.UV.clear();
  mesh.Normal.clear();
  mesh.FaceUV.clear();
  mesh.FaceNormal.clear();
  std::cout << mesh.XYZ.size() << " vertices, " << mesh.FaceXYZ.size() << " faces\n";
  if ((mesh.XYZ.size() == 0) || (mesh.FaceXYZ.size() == 0)) {
    std::cerr << "ERROR: Invalid mesh\n";
    return false;
  }

  // Prepare mesh for processing
  const float           epsilon = 0.F;
  std::vector<uint32_t> adjacency(3 * mesh.FaceXYZ.size());
  auto                  hr = DirectX::GenerateAdjacencyAndPointReps(reinterpret_cast<uint32_t*>(mesh.FaceXYZ.data()),
                                                   mesh.FaceXYZ.size(),
                                                   reinterpret_cast<DirectX::XMFLOAT3*>(mesh.XYZ.data()),
                                                   mesh.XYZ.size(),
                                                   epsilon,
                                                   nullptr,
                                                   adjacency.data());
  if (FAILED(hr)) {
    std::cerr << "ERROR: Failed generating adjacency (" << hr << ")\n";
    return false;
  }

  // Validation
  std::wstring msgs;
  DirectX::Validate(reinterpret_cast<uint32_t*>(mesh.FaceXYZ.data()),
                    mesh.FaceXYZ.size(),
                    mesh.XYZ.size(),
                    adjacency.data(),
                    DirectX::VALIDATE_BACKFACING | DirectX::VALIDATE_BOWTIES,
                    &msgs);
  if (!msgs.empty()) {
    std::cerr << "WARNING: \n";
    std::wcerr << msgs;
  }

  // Clean
  std::vector<uint32_t> dups;
  bool                  breakBowties = true;
  hr                                 = DirectX::Clean(reinterpret_cast<uint32_t*>(mesh.FaceXYZ.data()),
                      mesh.FaceXYZ.size(),
                      mesh.XYZ.size(),
                      adjacency.data(),
                      nullptr,
                      dups,
                      breakBowties);
  if (FAILED(hr)) {
    std::cerr << "ERROR: Failed mesh clean " << hr << '\n';
    return false;
  }

  if (!dups.empty()) {
    std::cout << " [" << dups.size() << " vertex dups]\n";
    mesh.XYZ.reserve(mesh.XYZ.size() + dups.size());
    for (auto dupIdx : dups) mesh.XYZ.push_back(mesh.XYZ[dupIdx]);
  }

  uint32_t textureParametrizationWidth = std::max(params.inputHeight, params.inputWidth);
  uint32_t textureParametrizationHeight = textureParametrizationWidth;
  // Perform UVAtlas isocharting
  std::cout << "Computing isochart atlas on mesh...\n";
  std::vector<DirectX::UVAtlasVertex> vb;
  std::vector<uint8_t>                ib;
  float                               outStretch = 0.F;
  size_t                              outCharts  = 0;
  std::vector<uint32_t>               facePartitioning;
  std::vector<uint32_t>               vertexRemapArray;
  const auto                          start = std::chrono::steady_clock::now();
  hr                                        = UVAtlasCreate(reinterpret_cast<DirectX::XMFLOAT3*>(mesh.XYZ.data()),
                     mesh.XYZ.size(),
                     reinterpret_cast<uint32_t*>(mesh.FaceXYZ.data()),
                     DXGI_FORMAT_R32_UINT,
                     mesh.FaceXYZ.size(),
                     params.maxCharts,
                     params.maxStretch,
                     textureParametrizationWidth,
                     textureParametrizationHeight,
                     params.gutter,
                     adjacency.data(),
                     nullptr,
                     nullptr,
                     UVAtlasCallback,
                     DirectX::UVATLAS_DEFAULT_CALLBACK_FREQUENCY,
                     params.uvOptions,
                     vb,
                     ib,
                     &facePartitioning,
                     &vertexRemapArray,
                     &outStretch,
                     &outCharts);
  if (FAILED(hr)) {
    std::cerr << "ERROR: Failed creating isocharts " << hr << '\n';
    return false;
  }

  auto stop    = std::chrono::steady_clock::now();
  auto deltams = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Processing time: " << deltams.count() << " ms\n";
  std::cout << "Output # of charts: " << outCharts << ", resulting stretching " << outStretch << ", " << vb.size()
            << " verts\n";

  assert(ib.size() == 3 * mesh.FaceXYZ.size() * sizeof(uint32_t));
  memcpy(mesh.FaceXYZ.data(), ib.data(), ib.size());

  assert(vertexRemapArray.size() == vb.size());
  std::vector<Vector3D<float>> XYZ(vertexRemapArray.size());
  hr = DirectX::UVAtlasApplyRemap(reinterpret_cast<DirectX::XMFLOAT3*>(mesh.XYZ.data()),
                                  sizeof(DirectX::XMFLOAT3),
                                  mesh.XYZ.size(),
                                  vertexRemapArray.size(),
                                  vertexRemapArray.data(),
                                  reinterpret_cast<DirectX::XMFLOAT3*>(XYZ.data()));
  if (FAILED(hr)) {
    std::cerr << "ERROR: Failed applying atlas vertex remap (" << hr << ")\n";
    return false;
  }

  mesh.XYZ.resize(XYZ.size());
  for (size_t i = 0; i < XYZ.size(); i++) mesh.XYZ[i] = XYZ[i];

  msgs.clear();
  DirectX::Validate(reinterpret_cast<uint32_t*>(mesh.FaceXYZ.data()),
                    mesh.FaceXYZ.size(),
                    mesh.XYZ.size(),
                    adjacency.data(),
                    DirectX::VALIDATE_DEFAULT,
                    &msgs);
  if (!msgs.empty()) {
    std::cerr << "WARNING: \n";
    std::wcerr << msgs;
  }

  // Copy isochart UVs into mesh
  mesh.UV.reserve(vb.size());
  std::transform(vb.begin(), vb.end(), std::back_inserter(mesh.UV), [](DirectX::UVAtlasVertex& vtx) {
    return Vector2D<float>(vtx.uv.x, vtx.uv.y);
  });

  mesh.FaceUV = mesh.FaceXYZ;
  outputMesh  = mesh;
  return true;
}

//============================================================================

template bool TextureParametrization::generate<float>(TriMesh<float> const&       inputMesh,
                                                      TriMesh<float>&             outputMesh,
													  encoderParams const& params);

template bool TextureParametrization::generate<double>(TriMesh<double> const&      inputMesh,
                                                       TriMesh<double>&            outputMesh,
													   encoderParams const& params);
