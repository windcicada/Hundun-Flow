// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_error.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string_view>
#include <type_traits>
#include <utility>

namespace hundun::execution {

using BackendIdentity = std::uint64_t;
using AllocationIdentity = std::uint64_t;

struct AllocationCounters final {
  std::uint64_t allocation_events{};
  std::uint64_t allocated_bytes{};
  std::uint64_t deallocation_events{};
  std::uint64_t deallocated_bytes{};
  std::uint64_t live_bytes{};
  std::uint64_t peak_live_bytes{};
};

[[nodiscard]] AllocationCounters allocation_counters() noexcept;

enum class ExecutionSpace { host, device };
enum class ScalarFormat { float64 };
enum class ExecutionCapability {
  buffer_allocation,
  host_access,
  transfer,
  asynchronous_event
};

class ExecutionContext {
 public:
  virtual ~ExecutionContext() noexcept = default;
  virtual std::string_view backend_name() const noexcept = 0;
  virtual BackendIdentity backend_identity() const noexcept = 0;
  virtual ExecutionSpace space() const noexcept = 0;
  virtual bool ordered() const noexcept = 0;
  virtual bool supports(ExecutionCapability capability) const noexcept = 0;
};

class CpuReferenceContext final : public ExecutionContext {
 public:
  CpuReferenceContext() noexcept = default;
  std::string_view backend_name() const noexcept override;
  BackendIdentity backend_identity() const noexcept override;
  ExecutionSpace space() const noexcept override;
  bool ordered() const noexcept override;
  bool supports(ExecutionCapability capability) const noexcept override;
};

namespace detail {
struct AllocationControl;
struct EventState;

double* checked_view_element(
    const std::weak_ptr<AllocationControl>& control,
    AllocationIdentity allocation_identity, std::uint64_t epoch,
    BackendIdentity backend_identity, ExecutionSpace space,
    std::size_t offset_bytes, std::size_t element_count,
    std::size_t stride_elements, bool mutable_access_allowed,
    std::size_t index, bool require_index);
}  // namespace detail

namespace test {
class ExecutionTestAccess;
}

template <class T>
class VectorView;
class Buffer;
class ExecutionEvent;
ExecutionEvent transfer(VectorView<const double> source,
                        VectorView<double> destination,
                        ExecutionContext& context);

template <class T>
class VectorView final {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, const double>,
                "VectorView supports only double");

 public:
  VectorView(const VectorView&) noexcept = default;
  VectorView& operator=(const VectorView&) noexcept = default;
  VectorView(VectorView&&) noexcept = default;
  VectorView& operator=(VectorView&&) noexcept = default;

  template <class U,
            std::enable_if_t<std::is_same_v<T, const double> &&
                                 std::is_same_v<U, double>,
                             int> = 0>
  VectorView(const VectorView<U>& other) noexcept
      : control_(other.control_),
        offset_bytes_(other.offset_bytes_),
        element_count_(other.element_count_),
        stride_elements_(other.stride_elements_),
        allocation_identity_(other.allocation_identity_),
        epoch_(other.epoch_),
        backend_identity_(other.backend_identity_),
        space_(other.space_),
        writable_(false) {}

  std::size_t size() const noexcept { return element_count_; }
  std::size_t stride() const noexcept { return stride_elements_; }
  ScalarFormat scalar_format() const noexcept { return ScalarFormat::float64; }
  AllocationIdentity allocation_identity() const noexcept {
    return allocation_identity_;
  }
  std::size_t offset_bytes() const noexcept { return offset_bytes_; }
  std::uint64_t epoch() const noexcept { return epoch_; }
  BackendIdentity backend_identity() const noexcept {
    return backend_identity_;
  }
  ExecutionSpace space() const noexcept { return space_; }
  bool writable() const noexcept {
    return writable_ && std::is_same_v<T, double>;
  }

  T* data() const {
    double* pointer = detail::checked_view_element(
        control_, allocation_identity_, epoch_, backend_identity_, space_,
        offset_bytes_, element_count_, stride_elements_,
        !std::is_same_v<T, double> || writable_, 0, false);
    return pointer;
  }

  T& operator[](std::size_t index) const {
    double* pointer = detail::checked_view_element(
        control_, allocation_identity_, epoch_, backend_identity_, space_,
        offset_bytes_, element_count_, stride_elements_,
        !std::is_same_v<T, double> || writable_, index, true);
    return *pointer;
  }

 private:
  template <class>
  friend class VectorView;
  friend class Buffer;
  friend class test::ExecutionTestAccess;
  friend class ExecutionEvent;
  friend ExecutionEvent transfer(VectorView<const double>, VectorView<double>,
                                 ExecutionContext&);

  VectorView(std::weak_ptr<detail::AllocationControl> control,
             std::size_t offset_bytes, std::size_t element_count,
             std::size_t stride_elements,
             AllocationIdentity allocation_identity, std::uint64_t epoch,
             BackendIdentity backend_identity, ExecutionSpace space,
             bool writable) noexcept
      : control_(std::move(control)),
        offset_bytes_(offset_bytes),
        element_count_(element_count),
        stride_elements_(stride_elements),
        allocation_identity_(allocation_identity),
        epoch_(epoch),
        backend_identity_(backend_identity),
        space_(space),
        writable_(writable) {}

  std::weak_ptr<detail::AllocationControl> control_;
  std::size_t offset_bytes_{0};
  std::size_t element_count_{0};
  std::size_t stride_elements_{1};
  AllocationIdentity allocation_identity_{0};
  std::uint64_t epoch_{0};
  BackendIdentity backend_identity_{0};
  ExecutionSpace space_{ExecutionSpace::host};
  bool writable_{false};
};

class Buffer final {
 public:
  Buffer(ExecutionContext& context, std::size_t byte_size);
  ~Buffer();
  Buffer(Buffer&& other) noexcept;
  Buffer& operator=(Buffer&& other) noexcept;
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;

  std::size_t byte_size() const noexcept;
  AllocationIdentity allocation_identity() const noexcept;
  std::uint64_t epoch() const noexcept;
  BackendIdentity backend_identity() const noexcept;
  ExecutionSpace space() const noexcept;

  void reallocate(std::size_t byte_size);

  VectorView<double> view(std::size_t offset_elements,
                          std::size_t element_count,
                          std::size_t stride_elements = 1);
  VectorView<const double> view(std::size_t offset_elements,
                                std::size_t element_count,
                                std::size_t stride_elements = 1) const;

 private:
  friend class test::ExecutionTestAccess;

  void invalidate() noexcept;
  void require_live() const;
  void validate_view_range(std::size_t offset_elements,
                           std::size_t element_count,
                           std::size_t stride_elements) const;

  std::shared_ptr<detail::AllocationControl> control_;
};

class ExecutionEvent final {
 public:
  static ExecutionEvent completed();

  ExecutionEvent(const ExecutionEvent&) noexcept = default;
  ExecutionEvent& operator=(const ExecutionEvent&) noexcept = default;
  ExecutionEvent(ExecutionEvent&&) noexcept = default;
  ExecutionEvent& operator=(ExecutionEvent&&) noexcept = default;
  ~ExecutionEvent() = default;

  bool ready() const noexcept;
  void wait() const;

 private:
  friend class test::ExecutionTestAccess;
  friend ExecutionEvent transfer(VectorView<const double>, VectorView<double>,
                                 ExecutionContext&);

  explicit ExecutionEvent(std::shared_ptr<detail::EventState> state) noexcept
      : state_(std::move(state)) {}

  std::shared_ptr<detail::EventState> state_;
};

ExecutionEvent transfer(VectorView<const double> source,
                        VectorView<double> destination,
                        ExecutionContext& context);

}  // namespace hundun::execution
