# wxQuantizeCrash

The MSVC++ optimizer changes the output of this program.  This could
be the optimizer's fault.  Changing the "curN" variables in
quantize.cpp's pass2_fs_dither() to "volatile" makes the output
from the optimized code match the non-optimized output.

This code is an attempt to isolate the problem referenced in
[wxWidgets Issue #17764](http://trac.wxwidgets.org/ticket/17764).