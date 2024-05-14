Simón Gasión
u186412
240126

Gbuffers:
FragColor (xyz) --> texture color
NormalColor (xyz) --> normals
MaterialProperties (x) --> occlusion
MaterialProperties (y) --> metalness
MaterialProperties (z) --> shininess (inverse roughness)
FragColor (w) --> emissive red component
NormalColor (w) --> emissive green component
MaterialProperties (w) --> emissive blue component

Done:
- Deferred:
Gbuffers
Render illumination
Forward for alpha nodes
Geometry instead of quads (point&spotlight)

- PBR:
TODO

- SSAO:
SSAO
SSAO+
blurring (tried using previous results, final version uses neighboring pixels)
normal shading component (if normal is facing down occlusion is more intense - proposed in practice class)

- linear space:
implemented for all pipelines (deferred+forward+flat) (operations in linear space then converts to gamma for output)
converts to linear when reading, FBOs formatted for higher quality

- HDR and tonemapper:
Basic tonemapper from slides implemented. Configure from UI.


TODO:
refactor shadowmap atlas
report