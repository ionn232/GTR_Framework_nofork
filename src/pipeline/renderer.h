#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"

#define MAX_LIGHTS_SP 10

//forward declarations
class Camera;
class Skeleton;
namespace GFX {
	class Shader;
	class Mesh;
	class FBO;
}

namespace SCN {

	class Prefab;
	class Material;

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
	public:
		bool render_wireframe;
		bool render_boundaries;
		bool use_multipass;
		bool render_lights;
		bool disable_lights;
		bool gui_use_normalmaps = true;
		bool gui_use_emissive = true;
		bool gui_use_occlusion = true;
		bool gui_use_specular = true;

		GFX::Texture* skybox_cubemap;

		GFX::Texture* shadowmapAtlas;

		SCN::Scene* scene;

		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene();

		//add here your functions
		//...

		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera);

		//sorts node and children nodes to their respective container
		void categorizeNodes(SCN::Node* node, Camera* camera);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		//lab1
		void renderMeshWithMaterialLights(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void renderShadowmap(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void showUI();

		void generateShadowmaps();

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader
		void lightToShaderSP(GFX::Shader* shader); //send light uniforms to shader for single-pass rendering
		void lightToShaderMP(LightEntity* light, GFX::Shader* shader); //send light uniforms to shader for multi-pass rendering (one light)
		void baseRenderMP(GFX::Mesh* mesh, GFX::Shader* shader); //draws first render of multi-pass using only ambien light (blends others on top)
	};

};