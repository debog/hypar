#!/bin/bash
cd /home/ghosh/Codes/hypar/src

# Fix usages of static variables that were already renamed in declarations
# This handles cases where the variable is used without the s_ prefix

find . -type f \( -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" \) -exec sed -i \
  -e 's/\([^s_]\)one_twelve/\1s_one_twelve/g' \
  -e 's/^one_twelve/s_one_twelve/g' \
  -e 's/(one_twelve/(s_one_twelve/g' \
  -e 's/\*one_twelve/*s_one_twelve/g' \
  -e 's/ one_twelve/ s_one_twelve/g' \
  -e 's/\[one_twelve/[s_one_twelve/g' \
  -e 's/\([^s_]\)count\b/\1s_count/g' \
  -e 's/^count\b/s_count/g' \
  -e 's/(count\b/(s_count/g' \
  -e 's/ count\b/ s_count/g' \
  -e 's/++count/++s_count/g' \
  -e 's/count++/s_count++/g' \
  -e 's/\([^s_]\)tolerance/\1s_tolerance/g' \
  -e 's/^tolerance/s_tolerance/g' \
  -e 's/(tolerance/(s_tolerance/g' \
  -e 's/ tolerance/ s_tolerance/g' \
  -e 's/\([^s_]\)one_by_thirty/\1s_one_by_thirty/g' \
  -e 's/^one_by_thirty/s_one_by_thirty/g' \
  -e 's/\*one_by_thirty/*s_one_by_thirty/g' \
  -e 's/ one_by_thirty/ s_one_by_thirty/g' \
  -e 's/\([^s_]\)one_half/\1s_one_half/g' \
  -e 's/^one_half/s_one_half/g' \
  -e 's/\*one_half/*s_one_half/g' \
  -e 's/ one_half/ s_one_half/g' \
  -e 's/\([^s_]\)one_third/\1s_one_third/g' \
  -e 's/^one_third/s_one_third/g' \
  -e 's/\*one_third/*s_one_third/g' \
  -e 's/ one_third/ s_one_third/g' \
  -e 's/\([^s_]\)one_sixth/\1s_one_sixth/g' \
  -e 's/^one_sixth/s_one_sixth/g' \
  -e 's/\*one_sixth/*s_one_sixth/g' \
  -e 's/ one_sixth/ s_one_sixth/g' \
  -e 's/\([^s_]\)three_by_ten/\1s_three_by_ten/g' \
  -e 's/^three_by_ten/s_three_by_ten/g' \
  -e 's/\*three_by_ten/*s_three_by_ten/g' \
  -e 's/ three_by_ten/ s_three_by_ten/g' \
  -e 's/\([^s_]\)thirteen_by_twelve/\1s_thirteen_by_twelve/g' \
  -e 's/^thirteen_by_twelve/s_thirteen_by_twelve/g' \
  -e 's/\*thirteen_by_twelve/*s_thirteen_by_twelve/g' \
  -e 's/ thirteen_by_twelve/ s_thirteen_by_twelve/g' \
  -e 's/\([^s_]\)one_fourth/\1s_one_fourth/g' \
  -e 's/^one_fourth/s_one_fourth/g' \
  -e 's/\*one_fourth/*s_one_fourth/g' \
  -e 's/ one_fourth/ s_one_fourth/g' \
  -e 's/\([^s_]\)two_third/\1s_two_third/g' \
  -e 's/^two_third/s_two_third/g' \
  -e 's/\*two_third/*s_two_third/g' \
  -e 's/ two_third/ s_two_third/g' \
  {} \;

# Fix c1 and c2 carefully (avoiding matching ac1, c10, etc.)
find . -type f \( -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" \) -exec sed -i \
  -e 's/\([^a-zA-Z0-9_]\)c1\([^a-zA-Z0-9_]\)/\1s_c1\2/g' \
  -e 's/^c1\([^a-zA-Z0-9_]\)/s_c1\1/g' \
  -e 's/\([^a-zA-Z0-9_]\)c1$/\1s_c1/g' \
  -e 's/\([^a-zA-Z0-9_]\)c2\([^a-zA-Z0-9_]\)/\1s_c2\2/g' \
  -e 's/^c2\([^a-zA-Z0-9_]\)/s_c2\1/g' \
  -e 's/\([^a-zA-Z0-9_]\)c2$/\1s_c2/g' \
  {} \;

# Fix model constants (ndims, nvars, JacSize) - but be careful as these might be used in other contexts
find . -type f \( -name "*.c" -o -name "*.cpp" -o -name "*.cu" -o -name "*.h" \) -exec sed -i \
  -e 's/\([^m_]\)ndims\b/\1s_ndims/g' \
  -e 's/\([^m_]\)nvars\b/\1s_nvars/g' \
  -e 's/\bJacSize\b/s_JacSize/g' \
  -e 's/\bdummy\b/s_dummy/g' \
  {} \;

echo "Static variable usages fixed!"
