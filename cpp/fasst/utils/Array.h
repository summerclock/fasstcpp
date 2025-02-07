#pragma once

#include <vector>
#include <array>
#include <cassert>
#include <algorithm>

namespace fasst {
namespace utils {

template<typename T, int M, int N, bool Fixed = (M != -1 && N != -1)>
class Array;

// Specialization for fixed-size array (both dimensions are fixed)
template<typename T, int M, int N>
class Array<T, M, N, true>
{
public:

  static Array Zero()
  {
    Array mat{};
    mat.fill(T(0));
    return mat;
  }

  void fill(const T& val)
  {
    data.fill(val);
  }

  int rows() const { return M; }
  int cols() const { return N; }
  int size() const { return M * N; }

  T& operator()(int i)
  {
    assert(i >= 0 && i < M * N);
    return data[i];
  }
  const T& operator()(int i) const
  {
    assert(i >= 0 && i < M * N);
    return data[i];
  }
  T& operator()(int i, int j)
  {
    assert(i >= 0 && i < M && j >= 0 && j < N);
    return data[i * N + j];
  }
  const T& operator()(int i, int j) const
  {
    assert(i >= 0 && i < M && j >= 0 && j < N);
    return data[i * N + j];
  }

  // data now public for aggregate initialization.
  std::array<T, M * N> data;
};

// Specialization for array with dynamic dimensions (if one dimension is dynamic).
template<typename T, int M, int N>
class Array<T, M, N, false>
{
public:
  Array() :
      rows_(M != -1 ? M : 0), cols_(N != -1 ? N : 0) {}

  Array(int rows, int cols = 1) :
      rows_(rows), cols_(cols)
  {
    if (M != -1) {
      assert(rows == M);
    }
    if (N != -1) {
      assert(cols == N);
    }
    data_.resize(rows * cols);
  }

  // New: Initializer-list constructor to support initializer-list initialization.
  Array(int rows, int cols, std::initializer_list<T> init)
      : rows_(rows), cols_(cols), data_(init)
  {
      assert(init.size() == static_cast<size_t>(rows * cols));
  }

  // Static function to create and return a zero array.
  // For matrices with a fixed dimension on one side, the respective parameter must match.
  static Array Zero(int rows, int cols)
  {
    if (M != -1) {
      assert(rows == M);
    }
    if (N != -1) {
      assert(cols == N);
    }
    Array mat(rows, cols);
    mat.fill(T(0));
    return mat;
  }

  void resize(int rows, int cols)
  {
    if (M != -1) {
      assert(rows == M);
    }
    if (N != -1) {
      assert(cols == N);
    }
    rows_ = rows;
    cols_ = cols;
    data_.resize(rows * cols);
  }
  // only for lists
  void resize(int size)
  {
    if (M != -1) {
      assert(size == M);
    }
    rows_ = size;
    cols_ = 1;
    data_.resize(size);
  }

  void fill(const T& val)
  {
    std::fill(data_.begin(), data_.end(), val);
  }

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  int size() const { return rows_ * cols_; }

  T& operator()(int i)
  {
    assert(i >= 0 && i < rows_ * cols_);
    return data_[i];
  }
  const T& operator()(int i) const
  {
    assert(i >= 0 && i < rows_ * cols_);
    return data_[i];
  }
  T& operator()(int i, int j)
  {
    assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
    return data_[i * cols_ + j];
  }
  const T& operator()(int i, int j) const
  {
    assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
    return data_[i * cols_ + j];
  }

private:
  int rows_, cols_;
  std::vector<T> data_;
};

using ArrayXXd = Array<double, -1, -1>;
using ArrayXXi = Array<int, -1, -1>;

// Define a list as a special case of a array with one dimension fixed to 1.
template<typename T, int M>
using ArrayX = Array<T, M, 1>;

using ArrayXd = ArrayX<double, -1>;
using ArrayXi = ArrayX<int, -1>;

} // namespace utils
} // namespace fasst