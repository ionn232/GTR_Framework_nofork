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
Geometry instead of quads (point&spotlight)

- PBR:
todo

- SSAO:
SSAO
SSAO+
blurring (tried using previous results, final version uses neighboring pixels)
normal component (if normal is facing down occlusion is more intense)

- linear space:
implemented for all pipelines (deferred+forward+flat) (operations in linear space then converts to gamma for output)
linear transformation for: texture color, emissive color, light colors, ambient color, OCCLUSION?¿?¿?

- HDR and tonemapper:
todo

ajenjo seminar: 
linear en color textura y luz, occlusion factor?
linearizar valor bruto de la textura o (textura * u_color)? (igual para emissive)
contar pixeles de fondo para blur? (enseñar arbol ssao)
metallic/roughness usa metal_factor y rough_factor como emissive_factor? o solo textura?
