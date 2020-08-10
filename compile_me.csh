# Switches probably OK on a Mac, but had to be skipped on Linux:
# -bind_at_load -Wl,-stack_size  -Wl,-stack_size -Wl,4000000 -L/usr/local/lib -lgsl -L/usr/X11/lib -lX11 -lXt

#--- Non-debug version:

## g++ -std=c++11 -Wall -O2 -Df2cFortran -DSI_PIXEL_TEMPLATE_USE_BOOST -DSI_PIXEL_TEMPLATE_STANDALONE PixelResolutionHistograms.cc $1.cxx `root-config --cflags --libs` -I/usr/local/include  -o $1

#--- Debug version:
g++ -std=c++11 -Wall -g -O0 -Df2cFortran -DSI_PIXEL_TEMPLATE_USE_BOOST -DSI_PIXEL_TEMPLATE_STANDALONE PixelResolutionHistograms.cc $1.cxx `root-config --cflags --libs` -I/usr/local/include  -o ${1}
