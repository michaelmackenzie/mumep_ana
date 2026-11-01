#! /bin/bash

for TARGET in "nominal" "st37plus6up6dn" "st49hole1to5" "st49hole9to41"; do
    echo ">>> Target: ${TARGET}"
    root.exe -q -b "optimize_selection_fast.C(true, 81, false, \"${TARGET}\", 1e3)"
done
