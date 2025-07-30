/***********************************************************************************
 * This software module was originally developed by Tencent America LLC, Media Lab *
 * in the course of development of dynamic mesh compression. Tencent America LLC,  *
 * Media Lab, retains full right to modify and use the code for its own purpose.   *
 * This copyright notice must be included in all copies or derivative works.       *
 * Copyright (c) Tencent 2023.                                                     *
 ***********************************************************************************/
#include <sstream>
#include "image.hpp"
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#endif
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
#endif

//============================================================================

template<typename T>
VectorND<double>
Image<T>::interp(const double x, const double y) const {
  assert(x >= 0.0 && x <= 1.0);
  assert(y >= 0.0 && y <= 1.0);

  double lx = _width - 1.0;
  double ly = _height - 1.0;
  double px = lx * x;
  double py = ly * (1.0 - y);  // y = 0 means the last row

  uint32_t x0 = std::floor(px);
  uint32_t y0 = std::floor(py);
  uint32_t x1 = std::ceil(px);
  uint32_t y1 = std::ceil(py);

  double dx0 = px - x0;
  double dy0 = py - y0;
  double dx1 = 1.0 - dx0;
  double dy1 = 1.0 - dy0;

  double w00 = dx0 * dy0;
  double w01 = dx0 * dy1;
  double w10 = dx1 * dy0;
  double w11 = dx1 * dy1;

  VectorND<double> v00(get(y0, x0));
  VectorND<double> v01(get(y1, x0));
  VectorND<double> v10(get(y0, x1));
  VectorND<double> v11(get(y1, x1));

  return w11 * v00 + w10 * v01 + w01 * v10 + w00 * v11;
}

//============================================================================

template<typename T>
bool
Image<T>::read(const std::string& filePath) {
  auto ext = std::filesystem::path(filePath).extension().string();
  if (ext == ".png" || ext == ".jpg" || ext == ".bmp" || ext == ".tga") {
    int      w, h, d;
    uint8_t* buffer = stbi_load(filePath.c_str(), &w, &h, &d, 0);
    resize(h, w, d);
    uint32_t i = 0;
    for (uint32_t r = 0; r < _height; r++) {
      for (uint32_t c = 0; c < _width; c++) {
        _data[r][c] = VectorND<uint8_t>(buffer + i, 3);
        i += 3;
      }
    }
    stbi_image_free(buffer);
  } else {
    printf("Error: reading images with format %s is not supported.\n", ext.c_str());
    return false;
  }
  return true;
}

//============================================================================

template<typename T>
bool
Image<T>::write(const std::string& filePath, const uint32_t quality) const {
  if (filePath.empty()) return false;
  auto                 ext = std::filesystem::path(filePath).extension().string();
  std::vector<uint8_t> buffer(_height * _width * _depth);
  uint32_t             count = 0;
  for (uint32_t r = 0; r < _height; r++) {
    for (uint32_t c = 0; c < _width; c++) {
      for (uint32_t d = 0; d < _depth; d++) { buffer[count++] = _data[r][c][d]; }
    }
  }
  if (ext == ".png") {
    return stbi_write_png(filePath.c_str(), _width, _height, _depth, buffer.data(), _width * _depth);
  } else if (ext == ".jpg") {
    return stbi_write_jpg(filePath.c_str(), _width, _height, _depth, buffer.data(), quality);
  } else if (ext == ".bmp") {
    return stbi_write_bmp(filePath.c_str(), _width, _height, _depth, buffer.data());
  } else if (ext == ".tga") {
    return stbi_write_tga(filePath.c_str(), _width, _height, _depth, buffer.data());
  }
  printf("Error: writing images with format %s is not supported.\n", ext.c_str());
  return false;
}

//============================================================================

template<typename T>
double
Image1<T>::interp(const double x, const double y) const {
  assert(x >= 0.0 && x <= 1.0);
  assert(y >= 0.0 && y <= 1.0);

  double lx = _width - 1.0;
  double ly = _height - 1.0;
  double px = lx * x;
  double py = ly * (1.0 - y);  // y = 0 means the last row

  uint32_t x0 = std::floor(px);
  uint32_t y0 = std::floor(py);
  uint32_t x1 = std::ceil(px);
  uint32_t y1 = std::ceil(py);

  double dx0 = px - x0;
  double dy0 = py - y0;
  double dx1 = 1.0 - dx0;
  double dy1 = 1.0 - dy0;

  double w00 = dx0 * dy0;
  double w01 = dx0 * dy1;
  double w10 = dx1 * dy0;
  double w11 = dx1 * dy1;

  double v00 = get(y0, x0);
  double v01 = get(y1, x0);
  double v10 = get(y0, x1);
  double v11 = get(y1, x1);

  return w11 * v00 + w10 * v01 + w01 * v10 + w00 * v11;
}

//============================================================================

template<typename T>
bool
Image1<T>::read(const std::string& filePath) {
  auto ext = std::filesystem::path(filePath).extension().string();
  if (ext == ".png" || ext == ".jpg" || ext == ".bmp" || ext == ".tga") {
    int      w, h, d;
    uint8_t* buffer = stbi_load(filePath.c_str(), &w, &h, &d, 0);
    if (d != 1) {
      printf("Error: Image1 cannot read images with channel number != 1.\n");
      return false;
    }
    resize(h, w);
    uint32_t i = 0;
    for (uint32_t r = 0; r < _height; ++r) {
      for (uint32_t c = 0; c < _width; ++c) {
        _data[r][c] = *(buffer + i);
        ++i;
      }
    }
    stbi_image_free(buffer);
  } else {
    printf("Error: reading images with format %s is not supported.\n", ext.c_str());
    return false;
  }
  return true;
}

//============================================================================

template<typename T>
bool
Image1<T>::write(const std::string& filePath, const uint32_t quality) const {
  if (filePath.empty()) return false;
  auto                 ext = std::filesystem::path(filePath).extension().string();
  std::vector<uint8_t> buffer(_height * _width);
  uint32_t             count = 0;
  for (uint32_t r = 0; r < _height; ++r) {
    for (uint32_t c = 0; c < _width; ++c) { buffer[count++] = _data[r][c]; }
  }
  if (ext == ".png") {
    return stbi_write_png(filePath.c_str(), _width, _height, 1, buffer.data(), _width);
  } else if (ext == ".jpg") {
    return stbi_write_jpg(filePath.c_str(), _width, _height, 1, buffer.data(), quality);
  } else if (ext == ".bmp") {
    return stbi_write_bmp(filePath.c_str(), _width, _height, 1, buffer.data());
  } else if (ext == ".tga") {
    return stbi_write_tga(filePath.c_str(), _width, _height, 1, buffer.data());
  }
  printf("Error: writing images with format %s is not supported.\n", ext.c_str());
  return false;
}

//============================================================================

template<typename T>
Vector3D<double>
Image3<T>::interp(const double x, const double y) const {
  assert(x >= 0.0 && x <= 1.0);
  assert(y >= 0.0 && y <= 1.0);

  double lx = _width - 1.0;
  double ly = _height - 1.0;
  double px = lx * x;
  double py = ly * (1.0 - y);  // y = 0 means the last row

  uint32_t x0 = std::floor(px);
  uint32_t y0 = std::floor(py);
  uint32_t x1 = std::ceil(px);
  uint32_t y1 = std::ceil(py);

  double dx0 = px - x0;
  double dy0 = py - y0;
  double dx1 = 1.0 - dx0;
  double dy1 = 1.0 - dy0;

  double w00 = dx0 * dy0;
  double w01 = dx0 * dy1;
  double w10 = dx1 * dy0;
  double w11 = dx1 * dy1;

  Vector3D<double> v00(get(y0, x0));
  Vector3D<double> v01(get(y1, x0));
  Vector3D<double> v10(get(y0, x1));
  Vector3D<double> v11(get(y1, x1));

  return w11 * v00 + w10 * v01 + w01 * v10 + w00 * v11;
}

//============================================================================

template<typename T>
bool
Image3<T>::readImage(const std::string& filePath, uint32_t inputWidth, uint32_t inputHeight) {
  auto ext = std::filesystem::path(filePath).extension().string();
  if (ext == ".png" || ext == ".jpg" || ext == ".bmp" || ext == ".tga") {
    int      w, h, d;
    uint8_t* buffer = stbi_load(filePath.c_str(), &w, &h, &d, 0);
    //if (d != 3) {
    //  printf("Error: Image3 cannot read images with channel number != 3.\n");
    //  return false;
    //}
	assert((w == inputWidth) && (h == inputHeight));
    resize(h, w);
    uint32_t i = 0;
    for (uint32_t r = 0; r < _height; ++r) {
      for (uint32_t c = 0; c < _width; ++c) {
        _data[r][c] = Vector3D<uint8_t>(buffer + i);
        i += d;
      }
    }
    stbi_image_free(buffer);
  } else {
    printf("Error: reading images with format %s is not supported.\n", ext.c_str());
    return false;
  }
  return true;
}

//============================================================================

template<typename T>
bool
Image3<T>::writeImage(const std::string& filePath, const uint32_t quality) const {
  if (filePath.empty()) return false;
  auto                 ext = std::filesystem::path(filePath).extension().string();
  std::vector<uint8_t> buffer(_height * _width * 3);
  uint32_t             count = 0;
  for (uint32_t r = 0; r < _height; ++r) {
    for (uint32_t c = 0; c < _width; ++c) {
      for (uint32_t d = 0; d < 3; ++d) { buffer[count++] = _data[r][c][d]; }
    }
  }
  if (ext == ".png") {
    return stbi_write_png(filePath.c_str(), _width, _height, 3, buffer.data(), _width * 3);
  } else if (ext == ".jpg") {
    return stbi_write_jpg(filePath.c_str(), _width, _height, 3, buffer.data(), quality);
  } else if (ext == ".bmp") {
    return stbi_write_bmp(filePath.c_str(), _width, _height, 3, buffer.data());
  } else if (ext == ".tga") {
    return stbi_write_tga(filePath.c_str(), _width, _height, 3, buffer.data());
  }
  printf("Error: writing images with format %s is not supported.\n", ext.c_str());
  return false;
}

//============================================================================
template<typename T> 
void 
Image3<T>::rgb2yuv420(std::vector<unsigned char>& y_plane, 
					  std::vector<unsigned char>& u_plane,
					  std::vector<unsigned char>& v_plane){

	const float R2Y = 0.2126f;
	const float G2Y = 0.7152f;
	const float B2Y = 0.0722f;
	const float R2U = -0.1146f;
	const float G2U = -0.3854f;
	const float B2U = 0.5000f;
	const float R2V = 0.5000f;
	const float G2V = -0.4542f;
	const float B2V = -0.0458f;

	for (uint32_t y = 0; y < _height; ++y) {
		for (uint32_t x = 0; x < _width; ++x) {
			float r = _data[y][x][0];
			float g = _data[y][x][1];
			float b = _data[y][x][2];
			unsigned char Y = (unsigned char)((R2Y * r + G2Y * g + B2Y * b));
			unsigned char U = (unsigned char)(((R2U * r + G2U * g + B2U * b) + 128));
			unsigned char V = (unsigned char)(((R2V * r + G2V * g + B2V * b) + 128));

			y_plane[y * _width + x] = Y;
			if (x % 2 == 0 && y % 2 == 0) {
				int index_uv = (y / 2 * (_width / 2) + (x / 2));
				u_plane[index_uv] = U;
				v_plane[index_uv] = V;
			}
		}
	}
}

//============================================================================

template<typename T>
void
Image3<T>::dilatePadding(std::vector<std::vector<uint8_t>>& occ, const uint32_t iter) {
  auto&          img    = *this;
  const uint32_t height = img.height();
  const uint32_t width  = img.width();

  assert(occ.size() == height && occ[0].size() == width);

  for (uint32_t it = 0; it < iter; ++it) {
    auto occ2 = occ;
    for (uint32_t i = 0; i < height; ++i) {
      for (uint32_t j = 0; j < width; ++j) {
        if (!occ[i][j]) {
          uint32_t         count = 0;
          Vector3D<double> avg(0.0);
          if (i > 0 && occ[i - 1][j]) {
            avg += img.get(i - 1, j);
            count++;
          }
          if (j > 0 && occ[i][j - 1]) {
            avg += img.get(i, j - 1);
            count++;
          }
          if (i < height - 1 && occ[i + 1][j]) {
            avg += img.get(i + 1, j);
            count++;
          }
          if (j < width - 1 && occ[i][j + 1]) {
            avg += img.get(i, j + 1);
            count++;
          }
          if (count) {
            avg /= count;
            img.set(i, j, avg.round());
            occ2[i][j] = 255;
          }
        }
      }
    }
    occ = occ2;
  }
}

//============================================================================

template class Image3<uint8_t>;
template class Image1<uint8_t>;
template class Image<uint8_t>;
