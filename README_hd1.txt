Simón Gasión
simon.gasion01@estudiant.upf.edu
u186412
240126


Done:
- Alpha sorting:
renderables (didnt use structure seen in the lesson - done using the alredy existing Node class with slight modifications since it alredy contains the data needed)
alpha elements rendered last and sorted by distance

- Phong: 
ambient light + emissive + occlusion (occlusion only implemented for red channel of roughness texture - not occlusion texture)
point light
directional light
spot light
normalmaps
ranged attenuation
specular light (using blinn-phong model)

- Several lights: 
single pass
multi pass
frustum culling for lights

- Shadows:
Spotlight
Directional light
Shadow atlas
Single pass + multi pass
other: shadowmap (+atlas) resizing, directional light follows camera

(new) Render options:
Added a toggle for each of the lab's tasks under a dropbox for debugging and building the report.
- Render with lights		
Activate lab 1 implementation.
- Disable lights
Ignore light entities. (for debug purposes)
- Multipass lights		
Render lights using single pass (unchecked) or multi pass (checked) approach.
- use normals
Account for normalmaps when rendering
- use occlusion
Account for occlusion when rendering
- use emissive
Accouns for emissive materials when rendering
- use specular
Accounts for specular lighting when rendering
- use shadowmaps
Directional and spot lights cast shadows
- shadowmap resolution
Set the exact shadowmap resolution (10-5000, initially 1024)
