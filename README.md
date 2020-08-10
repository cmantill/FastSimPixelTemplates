# Legacy version of PixelAV Templates for FastSim

The `templatehitreco` will produce position resolution histograms for a reconstruction based on PIXELAV pixel templates for FastSim.

`.compile_me.csh` should contain the command to compile the main executable `templatehitreco.cxx`, e.g.:
```
source compile_me.csh templatehitreco
```

This executable takes as input a configuration file e.g. `phase2_barrel_100x25x150_58360_292.txt` which would point to the PIXELAV template number.

In addition you should have the PIXELAV templates *.out, e.g. `template_events_d58360.out`.
