Simón Gasión
u186412
240126

Gbuffers:
FragColor (xyz) --> texture color
NormalColor (xyz) --> normals
MaterialProperties (x) --> occlusion
MaterialProperties (y) --> metalness
MaterialProperties (z) --> roughness
FragColor (w) --> emissive red component
NormalColor (w) --> emissive green component
MaterialProperties (w) --> emissive blue component

Done:
- Deferred:
Gbuffers
Render illumination
Forward for alpha nodes
Geometry instead of quads

- PBR:
todo

- SSAO:
SSAO
SSAO+
blurring (using previous frame)
range check (TODO)
normal weights (TODO)

- linear space:
todo

- HDR and tonemapper:
todo

TODO: duda ajenjo: 
visualizar gbuffers y emissive / 
dithering cuenta nota si hemos hecho forward también?