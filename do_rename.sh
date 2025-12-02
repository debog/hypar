#!/bin/bash
cd /home/ghosh/Codes/hypar/src

# Rename static variables (not functions) by adding s_ prefix
# Using word boundaries to avoid partial matches

find . -type f \( -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" \) -exec sed -i \
  -e 's/\bstatic double one_twelve\b/static double s_one_twelve/g' \
  -e 's/\bstatic int count\b/static int s_count/g' \
  -e 's/\bstatic const double tolerance\b/static const double s_tolerance/g' \
  -e 's/\bstatic const double one_by_thirty\b/static const double s_one_by_thirty/g' \
  -e 's/\bstatic const double one_half\b/static const double s_one_half/g' \
  -e 's/\bstatic const double one_third\b/static const double s_one_third/g' \
  -e 's/\bstatic const double one_sixth\b/static const double s_one_sixth/g' \
  -e 's/\bstatic const double three_by_ten\b/static const double s_three_by_ten/g' \
  -e 's/\bstatic const double c1\b/static const double s_c1/g' \
  -e 's/\bstatic const double c2\b/static const double s_c2/g' \
  -e 's/\bstatic const double thirteen_by_twelve\b/static const double s_thirteen_by_twelve/g' \
  -e 's/\bstatic const double one_fourth\b/static const double s_one_fourth/g' \
  -e 's/\bstatic const int ndims\b/static const int s_ndims/g' \
  -e 's/\bstatic const int nvars\b/static const int s_nvars/g' \
  -e 's/\bstatic const int JacSize\b/static const int s_JacSize/g' \
  -e 's/\bstatic double two_third\b/static double s_two_third/g' \
  -e 's/\bstatic const int dummy\b/static const int s_dummy/g' \
  {} \;

# Now replace usages (not in declarations)
find . -type f \( -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" \) -exec sed -i \
  -e 's/\([^_]\)one_twelve\b/\1s_one_twelve/g' \
  -e 's/^\(\s*\)one_twelve\b/\1s_one_twelve/g' \
  -e 's/\*one_twelve\b/*s_one_twelve/g' \
  -e 's/\bone_by_thirty\b/s_one_by_thirty/g' \
  -e 's/\bone_half\b/s_one_half/g' \
  -e 's/\bone_third\b/s_one_third/g' \
  -e 's/\bone_sixth\b/s_one_sixth/g' \
  -e 's/\bthree_by_ten\b/s_three_by_ten/g' \
  -e 's/\bc1\b/s_c1/g' \
  -e 's/\bc2\b/s_c2/g' \
  -e 's/\bthirteen_by_twelve\b/s_thirteen_by_twelve/g' \
  -e 's/\bone_fourth\b/s_one_fourth/g' \
  -e 's/\btwo_third\b/s_two_third/g' \
  {} \;

echo "Static variables renamed successfully!"
