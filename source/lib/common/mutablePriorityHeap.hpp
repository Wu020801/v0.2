
#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

//============================================================================

// T should define
//  heapPosition() returns the index in the heap or -1 if not in heap
//  heapKey() returns the key
//  setHeapPosition() sets the element position
//  setHeapKey() sets the key
template<typename T>
class MutablePriorityHeap {
public:
  MutablePriorityHeap()                                      = default;
  MutablePriorityHeap(const MutablePriorityHeap&)            = default;
  MutablePriorityHeap& operator=(const MutablePriorityHeap&) = default;
  ~MutablePriorityHeap()                                     = default;
  int32_t  size() const { return static_cast<int32_t>(elements_.size()); }
  const T* element(const int32_t index) {
    assert(index < size());
    return elements_[index];
  }
  bool empty() const { return elements_.empty(); }
  void reserve(const int32_t sz) { elements_.reserve(sz); }
  void clear() { elements_.resize(0); }
  T*   top() { return elements_.empty() ? nullptr : elements_[0]; }
  T*   remove(T& e) {
    const auto i = e.heapPosition();
    if (i < 0) { return nullptr; }

    assert(i >= 0 && i < size());

    swap(i, size() - 1);
    elements_.pop_back();
    e.setHeapPosition(-1);  // not in heap

    if (i == size()) upheap(i);
    else if (elements_[i]->heapKey() < e.heapKey()) {
      downheap(i);
    } else {
      upheap(i);
    }
    return &e;
  }
  T* extract() {
    if (size() == 0) { return nullptr; }
    swap(0, size() - 1);
    T* dead = elements_.back();  //[size() - 1];
    dead->setHeapPosition(-1);
    elements_.pop_back();
    downheap(0);
    return dead;
  }
  void insert(T& e) {
    const auto i = size();
    e.setHeapPosition(i);
    elements_.push_back(&e);
    upheap(i);
  }
  void update(T& e) {
    const auto i = e.heapPosition();
    if (i < 0) { return; }
    assert(i < size());
    if (i > 0 && e.heapKey() > elements_[parent(i)]->heapKey()) {
      upheap(i);
    } else {
      downheap(i);
    }
  }

private:
  void place(T& e, const int32_t position) {
    elements_[position] = &e;
    e.setHeapPosition(position);
  }

  void swap(const int32_t position1, const int32_t position2) {
    std::swap(elements_[position1], elements_[position2]);
    elements_[position1]->setHeapPosition(position1);
    elements_[position2]->setHeapPosition(position2);
  }

  int32_t parent(const int32_t i) const {
    assert(i >= 0 && i < size());
    return (i - 1) / 2;
  }

  int32_t left(const int32_t i) const {
    assert(i >= 0 && i < size());
    return 2 * i + 1;
  }

  int32_t right(const int32_t i) const {
    assert(i >= 0 && i < size());
    return 2 * i + 2;
  }

  void upheap(const int32_t i) {
    if (!size() || i == size()) { return; }

    assert(i >= 0 && i < size());
    T*   moving = elements_[i];
    auto index  = i;
    auto p      = parent(i);

    while (index > 0 && moving->heapKey() > elements_[p]->heapKey()) {
      place(*elements_[p], index);
      index = p;
      p     = parent(p);
    }

    if (index != i) { place(*moving, index); }
  }

  void downheap(const int32_t i) {
    if (!size()) { return; }

    assert(i >= 0 && i < size());
    T*      moving  = elements_[i];
    auto    index   = i;
    int32_t largest = 0;

    const auto elementCount = size();
    while (true) {
      auto l = left(index);
      auto r = right(index);

      if (l >= elementCount) { break; }

      if (r < elementCount && elements_[l]->heapKey() < elements_[r]->heapKey()) {
        largest = r;
      } else {
        largest = l;
      }

      if (moving->heapKey() >= elements_[largest]->heapKey()) { break; }

      place(*elements_[largest], index);
      index = largest;
    }

    if (index != i) { place(*moving, index); }
  }

  std::vector<T*> elements_;
};
