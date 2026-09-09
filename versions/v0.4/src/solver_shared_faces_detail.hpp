// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "solver_cartesian_detail.hpp"
#include <array>

namespace hundun::v04::detail {
// The caller validates the box, activity and nonaliasing immutable face inputs.
// A cell consumes faces in x-/x+/y-/y+/z-/z+ order. No cache survives this call.
template<class Face, class Cell>
Status shared_cell_faces(Int3 cells, KernelBox box,
    Span<const std::uint8_t> activity, Face face_value, Cell consume) noexcept {
  const Int3 end{box.begin.x+box.cells.x,box.begin.y+box.cells.y,box.begin.z+box.cells.z};
  // Eight cells per edge bounds scratch to below 2 KiB. A face shared
  // by active cells in this tile is evaluated once; tile/solid boundaries
  // invalidate the cache. The face callback reads only immutable inputs.
  constexpr int width=8;
  struct Saved { double value{}; bool valid{}; };
  for(int bz=box.begin.z;bz<end.z;bz+=width)
    for(int by=box.begin.y;by<end.y;by+=width)
      for(int bx=box.begin.x;bx<end.x;bx+=width) {
        std::array<Saved,width*width> z_faces{};
        for(int z=bz;z<std::min(bz+width,end.z);++z) {
          std::array<Saved,width> y_faces{};
          for(int y=by;y<std::min(by+width,end.y);++y) {
            Saved x_face{};
            for(int x=bx;x<std::min(bx+width,end.x);++x) {
              const Int3 c{x,y,z};
              const auto index=(static_cast<std::size_t>(z)*cells.y+y)*cells.x+x;
              Saved* saved[]{&x_face,&y_faces[x-bx],&z_faces[(y-by)*width+x-bx]};
              if(activity.size && activity.data[index]==0U) {
                for(auto* entry:saved) entry->valid=false;
                continue;
              }
              std::array<double,6U> values{};
              for(unsigned a=0;a<3;++a) {
                auto& entry=*saved[a];
                if(entry.valid) values[2*a]=entry.value;
                else {
                  auto status=face_value(static_cast<CartesianAxis>(a),c,values[2*a]);
                  if(!status) return status;
                }
                auto plus=c; ++(a==0 ? plus.x : a==1 ? plus.y : plus.z);
                auto status=face_value(static_cast<CartesianAxis>(a),plus,values[2*a+1]);
                if(!status) return status;
                entry={values[2*a+1],true};
              }
              auto status=consume(c,values);
              if(!status) return status;
            }
          }
        }
      }
  return {};
}
} // namespace hundun::v04::detail
