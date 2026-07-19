// SPDX-License-Identifier: Apache-2.0

#include <type_traits>
#include <utility>

#include "hundun/runtime/field_view.hpp"

namespace {

struct CellCallback {
  template <class View>
  void operator()(View) const;
};

struct FaceCallback {
  template <class View>
  void operator()(View) const;
};

template <class T, class = void>
struct HasAdlCellHelper : std::false_type {};

template <class T>
struct HasAdlCellHelper<
    T, std::void_t<decltype(with_kernel_cell_view(
           std::declval<const hundun::runtime::FieldView<T> &>(),
           CellCallback{}))>> : std::true_type {};

template <class T, class = void>
struct HasAdlFaceHelper : std::false_type {};

template <class T>
struct HasAdlFaceHelper<
    T, std::void_t<decltype(with_kernel_face_view(
           std::declval<const hundun::runtime::FaceFieldView<T> &>(),
           FaceCallback{}))>> : std::true_type {};

static_assert(!HasAdlCellHelper<double>::value,
              "field_view.hpp must not advertise the cell kernel helper");
static_assert(!HasAdlFaceHelper<double>::value,
              "field_view.hpp must not advertise the face kernel helper");

}  // namespace

int main() { return 0; }
