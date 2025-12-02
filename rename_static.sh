#!/bin/bash

# Navigate to source directory
cd /home/ghosh/Codes/hypar/src

# List of static variables to rename (found through grep)
# Format: "filename:old_name:new_name"

declare -a renames=(
    "SecondDerivative/SecondDerivativeFourthOrder.c:one_twelve:s_one_twelve"
    "IOFunctions/WriteArray.c:count:s_count"
    "HyParFunctions/CalculateROMDiff.c:tolerance:s_tolerance"
    "HyParFunctions/CalculateError.c:tolerance:s_tolerance"
    "TimeIntegration/TimeError.c:tolerance:s_tolerance"
    "InterpolationFunctions/Interp1PrimFifthOrderUpwindChar.c:one_by_thirty:s_one_by_thirty"
    "InterpolationFunctions/Interp1PrimFifthOrderHCWENOChar.c:one_half:s_one_half"
    "InterpolationFunctions/Interp1PrimFifthOrderHCWENOChar.c:one_third:s_one_third"
    "InterpolationFunctions/Interp1PrimFifthOrderHCWENOChar.c:one_sixth:s_one_sixth"
    "InterpolationFunctions/Interp1PrimFifthOrderCRWENOChar.c:one_third:s_one_third"
    "InterpolationFunctions/Interp1PrimFifthOrderCRWENOChar.c:one_sixth:s_one_sixth"
    "InterpolationFunctions/Interp1PrimFifthOrderCRWENO.c:one_third:s_one_third"
    "InterpolationFunctions/Interp1PrimFifthOrderCRWENO.c:one_sixth:s_one_sixth"
    "InterpolationFunctions/Interp1PrimFourthOrderCentral.c:c1:s_c1"
    "InterpolationFunctions/Interp1PrimFourthOrderCentral.c:c2:s_c2"
    "InterpolationFunctions/Interp1PrimFourthOrderCentralChar.c:c1:s_c1"
    "InterpolationFunctions/Interp1PrimFourthOrderCentralChar.c:c2:s_c2"
    "InterpolationFunctions/Interp1PrimFifthOrderUpwind.c:one_by_thirty:s_one_by_thirty"
    "InterpolationFunctions/Interp1PrimFifthOrderCompactUpwind.c:one_third:s_one_third"
    "InterpolationFunctions/Interp1PrimFifthOrderWENO.c:one_sixth:s_one_sixth"
    "InterpolationFunctions/WENOFifthOrderCalculateWeights.c:thirteen_by_twelve:s_thirteen_by_twelve"
    "InterpolationFunctions/WENOFifthOrderCalculateWeights.c:one_fourth:s_one_fourth"
)

for entry in "${renames[@]}"; do
    IFS=':' read -r file old new <<< "$entry"
    if [ -f "$file" ]; then
        echo "Processing $file: $old -> $new"
        sed -i "s/\b$old\b/$new/g" "$file"
    fi
done

echo "Done renaming static variables"
