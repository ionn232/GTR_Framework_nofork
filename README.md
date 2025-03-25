# GTR
Implementation of a forward and deferred graphics pipeline, using C++, GLSL and the provided framework, with the following features:

Deferred:
- Ambient, directional, spotlights and pointlights using PBR
- Multiple lights supported with optimizations (not using quads)
- Dynamic shadowmap atlas for all spotlights and directional light sorces
- SSAO, SSAO+
- HDR and tonemapper
- irradiance using probes
- reflection using probes
- volumetric fog
- various post-processing effects: chromatic aberration, DOF, bloom, lens distortion...

Forward:
- Alpha sorting
- Ambient, directional, spotlights and pointlights using blinn-phong model
- Multiple light support (single pass and multi pass)
- Dynamic shadowmap atlas for all spotlights and directional light sorces

# Framework
OpenGL C++ Framework used for teach the Real-time Grapchics course at Universitat Pompeu Fabra.

Course material can be found here: https://tamats.com/learn/realtime-graphics/
