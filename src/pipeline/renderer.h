#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"
#include "../gfx/sphericalharmonics.h"

#define MAX_LIGHTS_SP 25

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
	IRRADIANCE,
	VOLUMETRIC,
	COUNT
};

enum eSSAO {
	TEXTURE,
	SSAO,
	SSAOplus,
};

//data structure for irradiance probes
struct sProbe {
	vec3 pos; //where is located
	vec3 local; //its ijk pos in the matrix
	int index; //its index in the linear array
	SphericalHarmonics sh; //coeffs
};

//struct to store grid info
struct sIrradianceInfo {
	vec3 start;
	vec3 end;
	vec3 dim;
	vec3 delta;
	int num_probes;
};

//struct to store reflection probes info
struct sReflectionProbe {
	vec3 pos;
	GFX::Texture* cubemap = NULL;
};



namespace SCN {

	class Prefab;
	class Material;

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
	public:
		//render options
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
		bool gui_use_tonemapper = true;
		bool use_irradiance = true;
		bool show_irr_probes = false;
		float irr_probe_size = 5.0f;
		float irr_factor = 1.0f;
		float irr_normal_distance = 7.0f;
		bool show_ref_probes = false;
		bool volumetric_light = true;
		float air_density = 0.06f;
		bool blur_volumetric = true;
		bool blur_render = false;
		bool circular_blur = true;
		bool use_motion_blur = true;
		float motion_blur_intensity = 0.1f;
		vec2 fx_blur = vec2(5.0, 5.0);

		//tonemapper parameters
		float tmp_scale = 1.0f;
		float tmp_avg_lum = 1.0f;
		float tmp_lumwhite = 1.0f;

		//program configurations
		eRenderTypes render_mode = eRenderTypes::DEFERRED;
		eDeferredDisplay deferred_display = eDeferredDisplay::DEFAULT;
		eSSAO occlusion_mode = eSSAO::SSAOplus;

		GFX::Texture* skybox_cubemap;

		//FBOs and configuration
		GFX::FBO* shadowmapAtlas;
		int prevAtlasSize = 0;
		int numShadowmaps; //to avoid building an unnecesarily large atlas

		GFX::FBO* gBuffersFBO;

		GFX::FBO* linear_fbo;

		GFX::FBO* ssao_fbo; //current frame's occlusion
		GFX::FBO* fx_fbo; //fbo to render post-processing effect such as blurring
		float ssao_radius = 7.0f;
		std::vector<vec3> ssao_positions;
		Matrix44 prevViewProj; //for temporal reprojection in SSAO (unused)

		GFX::FBO* volumFBO;

		GFX::FBO* multi_probes_fbo;

		//scene container
		SCN::Scene* scene;

		//window size
		vec2 size;
		vec2 prevScreenSize = vec2(0.0, 0.0); //to detect if fbos must be rebuilt

		//variables for irradiance computation
		std::vector<sProbe> probes;
		sIrradianceInfo probes_info;
		Camera* probeCam;
		GFX::Texture* probes_texture;

		//reflection probes
		std::vector<sReflectionProbe*> reflection_probes;
		Camera* reflectionCam;

		GFX::FBO* prev_motionblur;

		bool partial_render; //stop rendering at some intermediate point (to visualize gbuffers, ssao, etc.)

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

		//forward scene render adapted for irradiance probes (no shadowmaps, renders to whatever context is active previous to all)
		void renderProbeFaces(SCN::Scene*, Camera* camera, bool render_skybox);
	
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

		//irradiance methods
		void renderAllProbes(float size);
		void renderProbe(vec3 pos, float size, float* coeffs);
		void captureAllProbes();
		void captureProbe(sProbe& p);
		void renderIrradianceTexture();

		//reflection methods
		void renderReflectionProbe(sReflectionProbe& s, float scale);
		void captureReflectionProbe(sReflectionProbe& s);
	};

};