#!/bin/bash

set -e

echo "`date` -------------"  >> stdout 2>&1

make

# Run Strip (rather than deconvolution) on raw data.
# echo "" >> stdout 2>&1
# ./0_SubtractData.out >> stdout 2>&1


# Run synthetics models.
# echo "" >> stdout 2>&1
# ./RunModels.out >> stdout 2>&1


# Run 1D synthetics (for reflectivity)
# echo "" >> stdout 2>&1
# ./Run1D.out >> stdout 2>&1


# Run 2D synthetics (for SHAXI)
# echo "" >> stdout 2>&1
# ./Run2.5D.out >> stdout 2>&1


# Modeling ...
echo "" >> stdout 2>&1
# ./2_StackStripPREM.out >> stdout 2>&1
./4_Figure2.out >> stdout 2>&1

# or Modeling the direct subtraction
# echo "" >> stdout 2>&1
# ./2_SubtractBinStack.out >> stdout 2>&1




# echo "" >> stdout 2>&1
# ../plot_data_all >> stdout 2>&1
# mv 3_Figure2.pdf C9dataStackAll.pdf
# ../plot_data >> stdout 2>&1
# mv 3_Figure2.pdf C9dataStack.pdf
# ../plot_model_all >> stdout 2>&1
# mv 3_Figure2.pdf C9modelStackAll.pdf
# ../plot_model >> stdout 2>&1
# mv 3_Figure2.pdf C9modelStack.pdf

echo "" >> stdout 2>&1
echo "------------- `date`"  >> stdout 2>&1

exit 0
