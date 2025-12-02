#!/bin/bash
cd /home/ghosh/Codes/hypar

# Replace all usages of _NavierStokes3D_stride_ with s__NavierStokes3D_stride_
find . -type f \( -name "*.c" -o -name "*.cu" -o -name "*.cpp" -o -name "*.h" \) \
  -exec sed -i 's/_NavierStokes3D_stride_/s__NavierStokes3D_stride_/g' {} \;

echo "Fixed _NavierStokes3D_stride_ variable usages"
