#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"

#define MAX_LIGHTS_SP 15

//forward declarations
class Camera;
class Skeleton;
namespace GFX {
	class Shader;
	class Mesh;
	class FBO;
}

enum eRenderTypes {
	FLAT,
	FORWARD,
	DEFERRED
};

enum eDeferredDisplay {
	DEFAULT,
	COLOR,
	NORMALS,
	MATERIAL_PROPERTIES,
	DEPTH, //IDEA: traducir de logaritmico a lineal, adaptar shader depth.fs
	EMISSIVE,
	SSAO_result,
	COUNT
};

enum eSSAO {
	TEXTURE,
	SSAO,
	SSAOplus,
};

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
		bool gui_use_multipass = false;
		bool disable_lights= false;
		bool gui_use_normalmaps = true;
		bool gui_use_emissive = true;
		bool gui_use_occlusion = true;
		bool gui_use_specular = true;
		bool gui_use_shadowmaps = true;
		int gui_shadowmap_res = 1024;
		bool gui_use_tonemapper = false;
		eRenderTypes render_mode = eRenderTypes::DEFERRED;
		eDeferredDisplay deferred_display = eDeferredDisplay::DEFAULT;
		eSSAO occlusion_mode = eSSAO::SSAOplus;

		vec2 prevScreenSize = vec2(0.0, 0.0); //to detect if fbos must be rebuilt

		GFX::Texture* skybox_cubemap;

		GFX::FBO* shadowmapAtlas;
		int prevAtlasSize = 0;
		int numShadowmaps; //to avoid building an unnecesarily large atlas

		GFX::FBO* gBuffersFBO;
		GFX::FBO* deferred_fbo;

		GFX::FBO* ssao_fbo; //current frame's occlusion
		GFX::FBO* blurred_ssao; //occlusion results
		float ssao_radius = 10.0f;
		std::vector<vec3> ssao_positions;
		Matrix44 prevViewProj; //for temporal reprojection in SSAO

		SCN::Scene* scene;

		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene();

		void categorizeNodes(Camera* camera);

		//add here your functions
		//...

		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);
		void renderSceneFlat(SCN::Scene* scene, Camera* camera);
		void renderSceneForward(SCN::Scene* scene, Camera* camera);
		void renderSceneDeferred(SCN::Scene* scene, Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera, eRenderTypes render_type);

		//sorts node and children nodes to their respective container
		void categorizeNodes(SCN::Node* node, Camera* camera);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		//lab1
		void renderMeshWithMaterialLights(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, bool use_multipass);

		void renderMeshWithMaterialGBuffers(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void renderShadowmap(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void showUI();

		void generateShadowmaps(Camera* main_camera);

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader
		void lightToShaderSP(GFX::Shader* shader); //send light uniforms to shader for single-pass rendering
		void lightToShaderMP(LightEntity* light, GFX::Shader* shader); //send light uniforms to shader for multi-pass rendering (one light)
		void baseRenderMP(GFX::Mesh* mesh, GFX::Shader* shader); //draws first render of multi-pass using only ambien light (blends others on top)
	};

};