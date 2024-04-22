Simón Gasión
simon.gasion01@estudiant.upf.edu
u186412
240126


Done:
- Alpha sorting:
renderables (didnt use structure seen in the lesson - done using the alredy existing Node class with slight modifications since it alredy contains the data needed)
alpha elements rendered last and sorted by distance

- Phong: 
ambient light + emissive + occlusion (occlusion only implemented for red channel of roughness texture - not act)
point light
directional light
spot light
ranged attenuation
specular light (using blinn-phong model as suggested)

- Several lights: 
single pass
multi pass (some inconsistencies in directional lights)
frustum culling 
further optimization (TODO, test if sphere overlaps with bounding boxes and camera BoundingBoxSphereOverlap no implementado asi testear)

- Shadows:
Spotlight
Directional light
Shadow atlas
Single pass + multi pass
shadowmap resizing
optimizations: direccional se mueve, interpolacion lineal, multi-pass (TO-DO)


(new) Render options:
Added a toggle for each of the lab's tasks for debugging and building the report.
- Render with lights		
Activate lab 1 implementation.
- Disable lights
Ignore light entities. (for debug purposes)
- Multipass lights		
Render lights using single pass or multi pass implementation.
- Enable normals
Account for normalmaps when rendering
- Enable occlusion
Account for occlusion when rendering
- Enable emissive
Accouns for emissive materials when rendering
- Enable specular
Accounts for specular lighting when rendering
- Enable shadowmaps
Directional and spot lights cast shadows
