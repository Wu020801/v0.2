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
#include "vector.hpp"

//============================================================================

template<class T>
class SparseMatrix {
public:
  SparseMatrix()  = default;
  ~SparseMatrix() = default;

  SparseMatrix(uint32_t const rowCount, uint32_t const colCount) {
    _rowCount = rowCount;
    _colCount = colCount;
    _pointers.resize(_rowCount + 1, 0);
  }

  void addElementInOrder(uint32_t const r, uint32_t const c, T const v) {
    assert(r < _rowCount);
    assert(c < _colCount);
    ++_pointers[r + 1];
    _indexes.push_back(c);
    _data.push_back(v);
  }
  void initialize(const uint32_t rowCount, const uint32_t columnCount) {
    _colCount = columnCount;
    _rowCount = rowCount;
    _data.resize(0);
    _pointers.resize(0);
    _indexes.resize(0);
    _pointers.resize(_rowCount + 1, 0);
  }

  void updatePointers() {
    for (uint32_t i = 1; i <= _rowCount; ++i) { _pointers[i] += _pointers[i - 1]; }
  }

  uint32_t rowBegin(uint32_t const r) const {
    assert(r < _rowCount);
    return _pointers[r];
  }

  uint32_t rowEnd(uint32_t const r) const {
    assert(r < _rowCount);
    return _pointers[r + 1];
  }

  SparseMatrix transpose() const {
    SparseMatrix res(_colCount, _rowCount);
    auto&        rpointers = res._pointers;
    for (auto const i : _indexes) { ++rpointers[i + 1]; }
    res.updatePointers();
    auto const nnz      = rpointers.back();
    auto&      rindexes = res._indexes;
    auto&      rbuffer  = res._data;
    rindexes.resize(nnz);
    rbuffer.resize(nnz);
    for (uint32_t r = 0; r < _rowCount; ++r) {
      for (uint32_t i = rowBegin(r), end = rowEnd(r); i < end; ++i) {
        auto const p = rpointers[_indexes[i]]++;
        rindexes[p]  = r;
        rbuffer[p]   = _data[i];
      }
    }
    for (int64_t c = _colCount - 1; c >= 0; --c) { rpointers[c + 1] = rpointers[c]; }
    rpointers[0] = 0;
    return res;
  }

  friend VectorND<T> operator*(SparseMatrix const& lhs, VectorND<T> const& rhs) {
    assert(lhs.colCount() == rhs.size());
    VectorND<T> res(lhs._rowCount);
    for (uint32_t r = 0; r < lhs._rowCount; ++r) {
      T p = T(0);
      for (uint32_t i = lhs.rowBegin(r), end = lhs.rowEnd(r); i < end; ++i) {
        p += lhs._data[i] * rhs[lhs._indexes[i]];
      }
      res[r] = p;
    }
    return res;
  }

  friend VectorND<T> operator*(SparseMatrix const& lhs, const T* rhs) {
    assert(rhs);
    VectorND<T> res(lhs._rowCount);
    for (uint32_t r = 0; r < lhs._rowCount; ++r) {
      T p = T(0);
      for (uint32_t i = lhs.rowBegin(r), end = lhs.rowEnd(r); i < end; ++i) {
        p += lhs._data[i] * rhs[lhs._indexes[i]];
      }
      res[r] = p;
    }
    return res;
  }

  friend SparseMatrix operator*(SparseMatrix const& lhs, SparseMatrix const& rhs) {
    assert(lhs.colCount() == rhs.rowCount());
    auto const   rhsT = rhs.transpose();
    SparseMatrix res(lhs._rowCount, rhs._colCount);
    VectorND<T>  tmp(rhs._rowCount);
    tmp = T(0);
    for (uint32_t r = 0; r < res._rowCount; ++r) {
      for (uint32_t c = 0; c < res._colCount; ++c) {
        for (uint32_t i = rhsT.rowBegin(c), end = rhsT.rowEnd(c); i < end; ++i) {
          tmp[rhsT._indexes[i]] = rhsT._data[i];
        }

        T p = T(0);
        for (uint32_t i = lhs.rowBegin(r), end = lhs.rowEnd(r); i < end; ++i) {
          p += lhs._data[i] * tmp[lhs._indexes[i]];
        }

        if (p != T(0)) { res.addElementInOrder(r, c, p); }

        for (uint32_t i = rhsT.rowBegin(c), end = rhsT.rowEnd(c); i < end; ++i) { tmp[rhsT._indexes[i]] = T(0); }
      }
    }

    res.updatePointers();
    return res;
  }

  uint32_t rowCount() const { return _rowCount; }
  uint32_t colCount() const { return _colCount; }

private:
  std::vector<T>        _data;
  std::vector<uint32_t> _pointers;
  std::vector<uint32_t> _indexes;
  uint32_t              _colCount = 0;
  uint32_t              _rowCount = 0;
};