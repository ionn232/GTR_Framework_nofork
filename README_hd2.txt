Simón Gasión
u186412
240126

Gbuffers:
FragColor (xyz) --> texture color
NormalColor (xyz) --> normals
MaterialProperties (x) --> occlusion
MaterialProperties (y) --> metalness
MaterialProperties (z) --> shininess
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
Enabled in deferred rendering, forward rendering keeps the previous implementation.

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



- Changes to previous submissions:
Reworked lightning in forward rendering (PBR is implemented in deferred) to correctly display shadows and spotlights.

- Known issues:
Because PBR is only implemented for deferred, semitransparent nodes (especially the car windows) look out of place when under high-intensity lights as they use forward rendering even when deferred mode is selected.



- new render options:
Render mode:		Select rendering pipeline (flat, forward, deferred)
Display channel:	Visualize something other than the final render in the deferred pipeline (Gbuffer values, SSAO, depth)
SSAO radius:		Sphere radius for generating the random points used in SSAO.
Use tonemapper:		Enable tonemapper. Default values don't change the render in respect to just gamma pass.
Tonemapper properties:	Sliders for each of the tonemapper components.