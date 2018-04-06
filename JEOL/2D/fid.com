#!/bin/csh

# X-Axis Size after DF Adjustment: 960 Complex Points.

delta2pipe -in ./indanone_hsqc_std-1.jdf \
  \
  -nodf  \
  -xN              1924  -yN               256  \
  -xT               962  -yT               128  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW         6410.256  -ySW        20105.555  \
  -xORIG       -803.422  -yORIG       -927.535  \
  -xOBS         399.756  -yOBS         100.519  \
  -xCAR           6.000  -yCAR          90.000  \
  -xFT             Time  -yFT             Time  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
