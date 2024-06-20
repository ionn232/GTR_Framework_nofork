#include "renderer.h"

#include <algorithm> //sort
#include <format>

#include "camera.h"
#include "../gfx/gfx.h"
#include "../gfx/shader.h"
#include "../gfx/mesh.h"
#include "../gfx/texture.h"
#include "../gfx/fbo.h"
#include "../pipeline/prefab.h"
#include "../pipeline/material.h"
#include "../pipeline/animation.h"
#include "../utils/utils.h"
#include "../extra/hdre.h"
#include "../core/ui.h"

#include "scene.h"

using namespace SCN;

//some globals
GFX::Mesh sphere;

std::vector<SCN::Node*> opaque_objects;
std::vector<SCN::Node*> semitransparent_objects;
std::vector<LightEntity*> lights;

bool compareDist(Node* s1, Node* s2) { 
	return s1->distance_to_camera > s2->distance_to_camera;
}

Renderer::Renderer(const char* shader_atlas_filename)
{
	//general variables initialization
	render_wireframe = false;
	render_boundaries = false;
	scene = nullptr;
	skybox_cubemap = nullptr;
	shadowmapAtlas = nullptr;
	gBuffersFBO = nullptr;
	linear_fbo = nullptr;
	ssao_fbo = nullptr;
	fx_fbo = nullptr;
	multi_probes_fbo = nullptr;
	volumFBO = nullptr;
	prevViewProj = Matrix44::IDENTITY;
	prev_motionblur = nullptr;
	ssao_positions = generateSpherePoints(128, 1, false); //IDEA: cambiar num puntos UI

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
	sphere.uploadToVRAM();

	//irradiance + reflection fbo
	multi_probes_fbo = new GFX::FBO();
	multi_probes_fbo->create(64, 64, 1, GL_RGB, GL_FLOAT, true);
	multi_probes_fbo->color_textures[0]->setName("PROBE FBO");

	//irradiance variables initialization
	probeCam = new Camera();
	probeCam->setPerspective(90, 1, 0.1, 1000);
	probes_texture = nullptr;

	FILE* f = fopen("irradiance.bin", "rb");
	if (!f) {
		std::cout << "previous irradiance data not found, please use 'Capture Irradiance' option" << std::endl;
		//define and fill irradiance probe grid
		//define bounding of the grid and num probes
		probes_info.start.set(-400, 0, -500);
		probes_info.end.set(600, 300, 500);
		probes_info.dim.set(25, 20, 25); //TODO: ajustar en UI

		//compute the vector from one corner to the other
		vec3 delta = (probes_info.end - probes_info.start);
		//compute delta from one probe to the next one
		delta.x /= (probes_info.dim.x - 1);
		delta.y /= (probes_info.dim.y - 1);
		delta.z /= (probes_info.dim.z - 1);
		probes_info.delta = delta; //store

		//lets compute the centers
		//pay attention at the order at which we add them
		for (int z = 0; z < probes_info.dim.z; ++z) {
			for (int y = 0; y < probes_info.dim.y; ++y) {
				for (int x = 0; x < probes_info.dim.x; ++x) {
					sProbe p;
					p.local.set(x, y, z);

					//index in the linear array
					p.index = x + y * probes_info.dim.x + z *
						probes_info.dim.x * probes_info.dim.y;

					//and its position
					p.pos = probes_info.start +
						probes_info.delta * Vector3f(x, y, z);
					probes.push_back(p);
				}
			}
		}
		probes_info.num_probes = probes.size();
	}
	else {
		std::cout << "Loading previous irradiance data..." << std::endl;
		//read header
		fread(&probes_info, sizeof(sIrradianceInfo), 1, f);

		//allocate space for the probes
		probes.resize(probes_info.num_probes);

		//read from disk directly to our probes container in memory
		fread(&probes[0], sizeof(sProbe), probes.size(), f);
		fclose(f);

		//render irradiance texture using loaded data.
		renderIrradianceTexture();
	}

	//reflection probes
	//hardcoded for the scene
	reflection_probes.push_back(new sReflectionProbe({ vec3(-80,80,-100), nullptr}));
	reflection_probes.push_back(new sReflectionProbe({ vec3(-50,80,233), nullptr }));
	reflection_probes.push_back(new sReflectionProbe({ vec3(350,100,-150), nullptr }));
	reflection_probes.push_back(new sReflectionProbe({ vec3(339,100,265), nullptr }));
}

void Renderer::setupScene()
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;

	//clear node and light containers
	lights.clear();
	semitransparent_objects.clear();
	opaque_objects.clear();
	numShadowmaps = 0;

	size = CORE::getWindowSize();
	partial_render = true; //set to false when necessary on render scene methods (flat, forward, deferred)

	if (!linear_fbo || prevScreenSize.distance(size) > 0.0) {
		delete linear_fbo;
		linear_fbo = new GFX::FBO();
		linear_fbo->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, true);
		linear_fbo->color_textures[0]->setName("LINEAR RENDER");
	}

}

void Renderer::categorizeNodes(Camera* camera) {
	//store entities in containers
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible)
			continue;

		//categorize all objects into containers
		if (ent->getType() == eEntityType::PREFAB) //prefabs
		{
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
			if (pent->prefab)
				categorizeNodes(&pent->root, camera);
		}
		else if (ent->getType() == eEntityType::LIGHT && !disable_lights) { //light objects
			//downcast to EntityLight and store in light array if it affects objects in the scene
			LightEntity* light = (SCN::LightEntity*)ent;
			if (light->light_type == eLightType::DIRECTIONAL || camera->testSphereInFrustum(light->root.model.getTranslation(), light->max_distance) == CLIP_INSIDE) { //simple frustum culling (for point and spotlight)
				lights.push_back(light);
				if (light->cast_shadows) {
					light->shadowmap_index = numShadowmaps;
					numShadowmaps += 1;
				}
			}
		}
	}
}

int getSMapdDimensions(int numlights) {
	return(ceil(sqrt((float)numlights)));
}

void Renderer::generateShadowmaps(Camera* main_camera) {
	//IDEA: move directional lights to 'player' relative distance to better fit shadowmaps
	int shadowmap_size = gui_shadowmap_res;
	int atlas_dimensions = getSMapdDimensions(numShadowmaps);
	int shadowatlas_size = atlas_dimensions * shadowmap_size;
	if ((shadowmapAtlas == nullptr || (prevAtlasSize != shadowatlas_size && prevAtlasSize > 0)) && (shadowatlas_size > 0)) {
		delete shadowmapAtlas;
		shadowmapAtlas = new GFX::FBO();
		shadowmapAtlas->setDepthOnly(shadowatlas_size, shadowatlas_size);
		shadowmapAtlas->depth_texture->setName("SHADOWMAP ATLAS");
	}

	shadowmapAtlas->bind();
	Camera camera;
	int light_index_i = 0;
	int light_index_j = 0;
	glClear(GL_DEPTH_BUFFER_BIT);

	if (numShadowmaps == 0) { //if shadowmaps dont exist just render a 
		GFX::Texture* texture = GFX::Texture::getWhiteTexture();
		texture->toViewport();
		shadowmapAtlas->unbind();
		prevAtlasSize = shadowatlas_size;
		return;
	}

	for (LightEntity* currentLight : lights) {
		if (currentLight->cast_shadows) {
			glViewport((light_index_i * shadowmap_size), (light_index_j * shadowmap_size), shadowmap_size, shadowmap_size);
			vec3 position = currentLight->root.model.getTranslation();
			if (currentLight->light_type == eLightType::DIRECTIONAL) { //IDEA: re-do main light system (there is just one directional light anyways)
				//set camera
				camera.setOrthographic(currentLight->area * -0.5, currentLight->area * 0.5, currentLight->area * -0.5, currentLight->area * 0.5, currentLight->near_distance, currentLight->max_distance);
				//move with player
				position = main_camera->eye;
			}
			else if (currentLight->light_type == eLightType::SPOT) {
				//set camera PREGUNTAR
				camera.setPerspective(currentLight->cone_info.y * 2.0f, 1.0f, currentLight->near_distance, currentLight->max_distance); 
			}
			camera.lookAt(position, position + currentLight->root.model.frontVector().normalize() * -1.0f, vec3(0, 1, 0));
			camera.enable();

			//smooth directional position
			float grid = currentLight->area / (float)shadowmap_size;
			camera.view_matrix.M[3][0] = round(camera.view_matrix.M[3][0] / grid) * grid;
			camera.view_matrix.M[3][1] = round(camera.view_matrix.M[3][1] / grid) * grid;
			camera.view_matrix.M[3][2] = round(camera.view_matrix.M[3][2] / grid) * grid;
			camera.viewprojection_matrix = camera.view_matrix * camera.projection_matrix;

			//render all opaque nodes
			for (Node* currentObj : opaque_objects) {
				BoundingBox world_bounding = transformBoundingBox(currentObj->getGlobalMatrix(), currentObj->getBoundingBox());
				//BoundingBox world_bounding = transformBoundingBox(node_model, node->mesh->box);
				//frustum culling for shadowmaps
				if (camera.testBoxInFrustum(world_bounding.center, world_bounding.halfsize)) {
					renderShadowmap(currentObj->getGlobalMatrix(), currentObj->mesh, currentObj->material); //IDEA: just use renderMeshWithMaterial
				}
			}
			currentLight->shadowmap_viewprojection = camera.viewprojection_matrix;
			//update shadowmap indexes
			light_index_i += 1;
			if (light_index_i >= atlas_dimensions) { light_index_i = 0; light_index_j += 1; }
		}
	}
	shadowmapAtlas->unbind();
	prevAtlasSize = shadowatlas_size;
}

void Renderer::renderShadowmap(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material) {
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	Camera* camera = Camera::current;

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* texture = NULL;

	texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	if (texture == NULL)
		texture = GFX::Texture::getWhiteTexture(); //a 1x1 white texture

	//select the blending
	glDisable(GL_BLEND);
	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);


	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("texture");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	if (texture)
		shader->setUniform("u_texture", texture, 0);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();
}

void Renderer::renderScene(SCN::Scene* scene, Camera* camera) {
	GFX::Mesh* quad = GFX::Mesh::getQuad();
	this->scene = scene;
	setupScene();
	categorizeNodes(camera);
	std::sort(std::begin(semitransparent_objects), std::end(semitransparent_objects), compareDist);	//sort blending vector - sorts nodes by distance in descending order

	//render to linear_fbo
	switch (render_mode) {
	case eRenderTypes::DEFERRED: //deferred rendering
		renderSceneDeferred(scene, camera);
		break;
	case eRenderTypes::FORWARD: //forward rendering
		renderSceneForward(scene, camera);
		break;
	default: //flat rendering
		renderSceneFlat(scene, camera);
	}

	//skip further computations if an intermediate process has been selected as deferred type (viewport alredy filled)
	if (!partial_render) {
		//apply PostFX transformations to the render in linear color space
		//each of these should be a function probably

		if (!fx_fbo || prevScreenSize.distance(size) > 0.0) {
			delete fx_fbo;
			fx_fbo = new GFX::FBO();
			fx_fbo->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, false);
			fx_fbo->color_textures[0]->setName("FX FBO");
		}

		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
		glClearColor(0.0, 0.0, 0.0, 1.0);
		vec2 invRes = vec2(1.0 / linear_fbo->color_textures[0]->width, 1.0 / linear_fbo->color_textures[0]->height);

		//FX0: color banding
		if (color_banding) {
			//apply color banding
			GFX::Shader* limit_shader = GFX::Shader::Get("color_banding");
			limit_shader->enable();
			fx_fbo->bind();
			vec2 invRes = vec2(1.0 / linear_fbo->color_textures[0]->width, 1.0 / linear_fbo->color_textures[0]->height);
			limit_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			limit_shader->setUniform("u_invRes", invRes);
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			limit_shader->disable();
			fx_fbo->unbind();

			linear_fbo->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			fx_fbo->color_textures[0]->toViewport();
			linear_fbo->unbind();
		}

		//FX1: DoF
		if (use_dof) {
			//get blurred image
			GFX::Shader* blur_shader = circular_blur ? GFX::Shader::Get("blur_circular") : GFX::Shader::Get("blur_neighbors");
			blur_shader->enable();
			fx_fbo->bind();
			blur_shader->setUniform("u_raw", linear_fbo->color_textures[0], 0);
			blur_shader->setUniform("u_depth_texture", linear_fbo->depth_texture, 1);
			blur_shader->setUniform("u_invRes", invRes);
			blur_shader->setUniform("u_blur_far", true);
			blur_shader->setUniform("u_blur_dimensions", vec2(8,8));
			blur_shader->setUniform("u_blur_distance", vec2(2,2));
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			fx_fbo->unbind();
			blur_shader->disable();

			//perform DoF
			GFX::Shader* dof_shader = GFX::Shader::Get("depth_of_field");
			dof_shader->enable();
			linear_fbo->bind();
			dof_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			dof_shader->setUniform("u_blurred_render", fx_fbo->color_textures[0], 1);
			dof_shader->setUniform("u_depth_texture", linear_fbo->depth_texture, 2);
			dof_shader->setUniform("u_min_dist", dof_min);
			dof_shader->setUniform("u_max_dist", dof_max);
			dof_shader->setUniform("u_invRes", invRes);
			dof_shader->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
			cameraToShader(camera, dof_shader);
			quad->render(GL_TRIANGLES);
			linear_fbo->unbind();
			dof_shader->disable();
		}

		//FX2: bloom
		//this can easily be made with less FBOs
		//and better performance
		//but I have an exam tomorrow
		if (use_bloom) {
			//(re)create as many FBOS as needed
			if (bloom_samples.size() != bloom_iterations || prevScreenSize.distance(size) > 0.0) {
				//free memory of previous FBOs
				for (std::vector<GFX::FBO*>::iterator i = bloom_samples.begin(); i != bloom_samples.end(); ++i) {
					delete* i;
				}
				bloom_samples.clear();

				//create FBO for each layer
				GFX::FBO* bloom_layer;
				for (int i = 0; i < bloom_iterations; i++) {
					//create fbo
					bloom_layer = new GFX::FBO();
					bloom_layer->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, false);
					//set name
					char buffer[50];
					std::snprintf(buffer, sizeof(buffer), "BLOOM LAYER %d", i+1);
					char* name = new char[std::strlen(buffer) + 1];
					std::strcpy(name, buffer);
					bloom_layer->color_textures[0]->setName(name);
					bloom_samples.push_back(bloom_layer);
				}
			}

			//isolate bright pixels
			GFX::Shader* bloom_shader = GFX::Shader::Get("bloom_pass");
			bloom_shader->enable();
			fx_fbo->bind();
			bloom_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			bloom_shader->setUniform("u_threshold", bloom_threshold);
			bloom_shader->setUniform("u_invRes", invRes);
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			bloom_shader->disable();
			fx_fbo->unbind();

			//downsample iteratively
			GFX::Shader* blur_shader = circular_blur ? GFX::Shader::Get("blur_circular") : GFX::Shader::Get("blur_neighbors");
			blur_shader->enable();
			for (int i = 0; i < bloom_iterations; i++) {
				bloom_samples.at(i)->bind();
				blur_shader->setUniform("u_raw", i==0 ? fx_fbo->color_textures[0] : bloom_samples.at(i-1)->color_textures[0], 0);
				blur_shader->setUniform("u_invRes", invRes);
				blur_shader->setUniform("u_blur_far", true);
				blur_shader->setUniform("u_blur_dimensions", vec2(9, 9));
				blur_shader->setUniform("u_blur_distance", vec2(4, 4));
				blur_shader->setUniform("u_set_alpha", 1.0f/float(bloom_iterations));
				glClear(GL_COLOR_BUFFER_BIT);
				quad->render(GL_TRIANGLES);
				bloom_samples.at(i)->unbind();
			}
			blur_shader->disable();

			//blend to render
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
			linear_fbo->bind();
			for (int i = 0; i < bloom_iterations; i++) {
				bloom_samples.at(i)->color_textures[0]->toViewport();
			}
			linear_fbo->unbind();
			glDisable(GL_BLEND);
		}

		//FX3: grain
		if (use_grain) {
			GFX::Shader* grain_shader = GFX::Shader::Get("grain");
			grain_shader->enable();
			fx_fbo->bind();
			grain_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			grain_shader->setUniform("u_intensity", grain_intensity);
			grain_shader->setUniform("u_grain_size", 10.0f - grain_size); //10 is the max size, this should not be hardcoded
			grain_shader->setUniform("u_invRes", invRes);
			grain_shader->setUniform("u_rand_seed", vec3(std::rand(),std::rand(),std::rand()));
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			grain_shader->disable();
			fx_fbo->unbind();

			linear_fbo->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			fx_fbo->color_textures[0]->toViewport();
			linear_fbo->unbind();
		}

		//FX4: image blur
		if (blur_render) {
			GFX::Shader* blur_shader = circular_blur ? GFX::Shader::Get("blur_circular") : GFX::Shader::Get("blur_neighbors");
			blur_shader->enable();
			fx_fbo->bind();
			blur_shader->setUniform("u_raw", linear_fbo->color_textures[0], 0);
			blur_shader->setUniform("u_depth_texture", linear_fbo->depth_texture, 1);
			blur_shader->setUniform("u_invRes", invRes);
			blur_shader->setUniform("u_blur_far", true);
			blur_shader->setUniform("u_blur_dimensions", fx_blur_res);
			blur_shader->setUniform("u_blur_distance", fx_blur_dist);
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			fx_fbo->unbind();
			blur_shader->disable();

			linear_fbo->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			fx_fbo->color_textures[0]->toViewport();
			linear_fbo->unbind();
		}

		//FX5: motion blur
		if (use_motion_blur) {
			glDisable(GL_DEPTH_TEST);
			glDisable(GL_BLEND);
			if (!prev_motionblur || prevScreenSize.distance(size) > 0.0) {
				prev_motionblur = new GFX::FBO();
				prev_motionblur->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, false);
				prev_motionblur->color_textures[0]->setName("MOTION BLUR RESULTS");

				prev_motionblur->bind();
				linear_fbo->color_textures[0]->toViewport();
				prev_motionblur->unbind();
			}
			linear_fbo->bind();
			GFX::Shader* motion_blur = GFX::Shader::Get("motion_blur");
			motion_blur->enable();
			motion_blur->setUniform("u_render", linear_fbo->color_textures[0], 0);
			motion_blur->setUniform("u_last_results", prev_motionblur->color_textures[0], 1);
			motion_blur->setUniform("u_intensity", motion_blur_intensity);
			motion_blur->setUniform("u_burn_in", frame_smoothing_burnin);
			quad->render(GL_TRIANGLES);
			motion_blur->disable();
			linear_fbo->unbind();

			//save previous results
			prev_motionblur->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			linear_fbo->color_textures[0]->toViewport();
			prev_motionblur->unbind();
		}

		//FX6: chromatic aberration
		if (use_chromatic_aberration) {
			GFX::Shader* CA_shader = GFX::Shader::Get("chromatic_aberration");
			CA_shader->enable();
			fx_fbo->bind();
			CA_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			CA_shader->setUniform("u_intensity", chromatic_aberration_factor);
			CA_shader->setUniform("u_invRes", invRes);
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			CA_shader->disable();
			fx_fbo->unbind();

			linear_fbo->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			fx_fbo->color_textures[0]->toViewport();
			linear_fbo->unbind();
		}


		//FX7: lens distortion
		if (use_lens_distortion) {
			GFX::Shader* LD_shader = GFX::Shader::Get("lens_distortion");
			LD_shader->enable();
			fx_fbo->bind();
			LD_shader->setUniform("u_render", linear_fbo->color_textures[0], 0);
			LD_shader->setUniform("u_intensity", lens_distortion_intensity);
			LD_shader->setUniform("u_dist_mode", lens_distortion_mode);
			LD_shader->setUniform("u_invRes", invRes);
			glClear(GL_COLOR_BUFFER_BIT);
			quad->render(GL_TRIANGLES);
			LD_shader->disable();
			fx_fbo->unbind();

			linear_fbo->bind();
			glClear(GL_COLOR_BUFFER_BIT);
			fx_fbo->color_textures[0]->toViewport();
			linear_fbo->unbind();
		}

		//FX8: tonemapper / gamma
		if (gui_use_tonemapper) {
			GFX::Shader* tonemapper_shader = GFX::Shader::Get("tonemapper");
			tonemapper_shader->enable();
			tonemapper_shader->setUniform("u_texture", linear_fbo->color_textures[0], 0);
			tonemapper_shader->setUniform("u_scale", tmp_scale);
			tonemapper_shader->setUniform("u_average_lum", tmp_avg_lum);
			tonemapper_shader->setUniform("u_lumwhite2", (float)pow(tmp_lumwhite, 2.0f));
			tonemapper_shader->setUniform("u_inv_gamma", 1.0f / 2.2f);
			quad->render(GL_TRIANGLES);
			tonemapper_shader->disable();
		}
		else {
			GFX::Shader* gamma_shader = GFX::Shader::Get("gamma");
			gamma_shader->enable();
			gamma_shader->setUniform("u_texture", linear_fbo->color_textures[0], 0);
			gamma_shader->setUniform("u_inv_gamma", 1.0f / 2.2f);
			quad->render(GL_TRIANGLES);
			gamma_shader->disable();
		}
	}
	//log size to re-create FBOs on window resize
	prevScreenSize = size;
}

void Renderer::renderSceneFlat(SCN::Scene* scene, Camera* camera) {
	
	linear_fbo->bind();

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if (skybox_cubemap)
		renderSkybox(skybox_cubemap);

	//render entities
	for (int i = 0; i < opaque_objects.size(); i++)
	{
		renderNode(opaque_objects[i], camera, eRenderTypes::FLAT);
	}
	//render semitransparent entities
	for (int i = 0; i < semitransparent_objects.size(); i++)
	{
		renderNode(semitransparent_objects[i], camera, eRenderTypes::FLAT);
	}

	linear_fbo->unbind();
	partial_render = false;
}

void Renderer::renderSceneForward(SCN::Scene* scene, Camera* camera)
{
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	generateShadowmaps(camera); //render shadowmap atlas
	camera->enable(); //reactivate scene camera

	linear_fbo->bind();

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if (skybox_cubemap)
		renderSkybox(skybox_cubemap);

	//render entities
	for (int i = 0; i < opaque_objects.size(); i++)
	{
		renderNode(opaque_objects[i], camera, eRenderTypes::FORWARD);
	}
	//render semitransparent entities
	for (int i = 0; i < semitransparent_objects.size(); i++)
	{
		renderNode(semitransparent_objects[i], camera, eRenderTypes::FORWARD);
	}

	linear_fbo->unbind();
	partial_render = false;
}

void Renderer::renderSceneDeferred(SCN::Scene* scene, Camera* camera) {
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	if (!gBuffersFBO || prevScreenSize.distance(size) > 0.0) {
		delete gBuffersFBO;
		gBuffersFBO = new GFX::FBO();
		gBuffersFBO->create(size.x, size.y, 3, GL_RGBA, GL_FLOAT, true);
	}
	gBuffersFBO->bind();

	//set the clear color
	glClearColor(0.0, 0.0, 0.0, 1.0);
	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render entities
	for (int i = 0; i < opaque_objects.size(); i++)
	{
		renderNode(opaque_objects[i], camera, eRenderTypes::DEFERRED);
	}

	gBuffersFBO->unbind();

	GFX::Mesh* quad = GFX::Mesh::getQuad(); //create quad to draw to screen (shit laptop does not support VBOs. Sorry!)

	switch (deferred_display) { //if special display is selected, show respective buffer and abort further rendering
	case eDeferredDisplay::COLOR: gBuffersFBO->color_textures[0]->toViewport();
		return;
	case eDeferredDisplay::NORMALS: gBuffersFBO->color_textures[1]->toViewport();
		return;
	case eDeferredDisplay::MATERIAL_PROPERTIES: gBuffersFBO->color_textures[2]->toViewport();
		return;
	case eDeferredDisplay::DEPTH: gBuffersFBO->depth_texture->toViewport();
		return;
	case eDeferredDisplay::EMISSIVE:
		GFX::Shader* emissive = GFX::Shader::Get("view_emissive");
		emissive->enable();
		emissive->setUniform("u_color_texture", gBuffersFBO->color_textures[0], 0);
		emissive->setUniform("u_normal_texture", gBuffersFBO->color_textures[1], 1);
		emissive->setUniform("u_mat_properties_texture", gBuffersFBO->color_textures[2], 2);
		emissive->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 3);
		quad->render(GL_TRIANGLES);
		emissive->disable();
		return;
	}

	//prepare shadowmaps
	generateShadowmaps(camera); //render shadowmap atlas
	camera->enable(); //reactivate scene camera

	//update reflection probes
	for (int i = 0; i < reflection_probes.size(); i++) {
		captureReflectionProbe(*reflection_probes.at(i));
	}
	camera->enable(); //reactivate scene camera

	//prepare SSAO
	if (occlusion_mode != eSSAO::TEXTURE) {
		//bool firstIteration = false;
		//create fbo
		if (!ssao_fbo || prevScreenSize.distance(size) > 0.0) {
			delete ssao_fbo;
			delete fx_fbo;
			ssao_fbo = new GFX::FBO();
			fx_fbo = new GFX::FBO();
			ssao_fbo->create(size.x/2, size.y/2, 1, GL_RGB, GL_UNSIGNED_BYTE, false);
			ssao_fbo->color_textures[0]->setName("RAW SSAO");
			fx_fbo->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, false);
			fx_fbo->color_textures[0]->setName("FX FBO");
			//firstIteration = true;
		}
		ssao_fbo->bind();
		glClearColor(1.0, 1.0, 1.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT);
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_BLEND);
		//upload necessary shader uniforms
		//IDEA: range check (diapos)
		GFX::Shader* ssao_shader = GFX::Shader::Get("ssao");
		ssao_shader->enable();
		ssao_shader->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 0);
		ssao_shader->setUniform("u_normal_texture", gBuffersFBO->color_textures[1], 1);
		ssao_shader->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
		ssao_shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
		ssao_shader->setUniform("u_invRes", vec2(1.0 / ssao_fbo->color_textures[0]->width, 1.0 / ssao_fbo->color_textures[0]->height));
		ssao_shader->setUniform("u_radius", ssao_radius);
		ssao_shader->setUniform("u_use_plus", occlusion_mode - 1); //SSAO --> 0, SSAO+ --> 1
		ssao_shader->setUniform3Array("u_random_pos", (float*) &ssao_positions[0], ssao_positions.size());

		//render SSAO
		quad->render(GL_TRIANGLES);
		//configure linear interpolation
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		ssao_shader->disable();
		ssao_fbo->unbind();

		//blur image using reprojection (results look like shit, no time to fix) (points were recomputed each frame, not anymore)
		/*
		if (!firstIteration) {
			GFX::Shader* ssao_blur_shader = GFX::Shader::Get("blur_reprojection");
			blurred_ssao->bind();
			ssao_blur_shader->enable();
			//ssao occlusion textures
			ssao_blur_shader->setUniform("u_current_frame", ssao_fbo->color_textures[0], 0);
			ssao_blur_shader->setUniform("u_last_results", blurred_ssao->color_textures[0], 1);
			//data needed for temporal reprojection
			ssao_blur_shader->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 2);
			ssao_shader->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
			ssao_shader->setUniform("u_previous_viewprojection", prevViewProj);
			ssao_shader->setUniform("u_invRes", vec2(1.0 / ssao_fbo->color_textures[0]->width, 1.0 / ssao_fbo->color_textures[0]->height));

			quad->render(GL_TRIANGLES);
			ssao_blur_shader->disable();
			blurred_ssao->unbind();

		}
		else {
			ssao_fbo->color_textures[0]->copyTo(blurred_ssao->color_textures[0]);
		} */
		
		//blur raw SSAO using neighbor pixels
		GFX::Shader* ssao_blur_shader = circular_blur ? GFX::Shader::Get("blur_circular") : GFX::Shader::Get("blur_neighbors");
		ssao_blur_shader->enable();
		fx_fbo->bind();
		vec2 invRes = vec2(1.0 / fx_fbo->color_textures[0]->width, 1.0 / fx_fbo->color_textures[0]->height);
		ssao_blur_shader->setUniform("u_raw", ssao_fbo->color_textures[0], 0);
		ssao_blur_shader->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 1);
		ssao_blur_shader->setUniform("u_invRes", invRes);
		ssao_blur_shader->setUniform("u_blur_far", false);
		ssao_blur_shader->setUniform("u_blur_dimensions", vec2(4.0, 4.0));
		ssao_blur_shader->setUniform("u_blur_distance", vec2(2.0, 2.0));
		ssao_blur_shader->setUniform("u_set_alpha", (float)0.0);
		glClear(GL_COLOR_BUFFER_BIT);
		quad->render(GL_TRIANGLES);
		fx_fbo->unbind();
		ssao_blur_shader->disable();
	}
	//SSAO to viewport if necessary
	if (deferred_display == eDeferredDisplay::SSAO_result) {
		fx_fbo->color_textures[0]->toViewport();
		return;
	}

	//deferred render to FBO
	linear_fbo->bind();

	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

	//render skybox
	if (skybox_cubemap && !(deferred_display == eDeferredDisplay::IRRADIANCE)) { //irradiance reuses the linear fbo so skip drawing skybox to display consistent results
		renderSkybox(skybox_cubemap);
	}

	GFX::Shader* deferred_global = GFX::Shader::Get("deferred_global");
	deferred_global->enable();

	//upload uniforms to draw emissive quad
	deferred_global->setUniform("u_is_quad", 1);

	//upload uniforms to shader
	cameraToShader(camera, deferred_global);

	deferred_global->setUniform("u_color_texture", gBuffersFBO->color_textures[0], 0);
	deferred_global->setUniform("u_normal_texture", gBuffersFBO->color_textures[1], 1);
	deferred_global->setUniform("u_mat_properties_texture", gBuffersFBO->color_textures[2], 2);
	deferred_global->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 3);

	deferred_global->setUniform("u_occlusion_type", occlusion_mode);
	deferred_global->setUniform("u_ssao_map", (occlusion_mode!=eSSAO::TEXTURE ? fx_fbo->color_textures[0] : GFX::Texture::getWhiteTexture()), 4);

	deferred_global->setUniform("u_invRes", vec2(1.0 / size.x, 1.0 / size.y));
	deferred_global->setUniform("u_ambient_light", scene->ambient_light);
	deferred_global->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);

	deferred_global->setUniform("u_shadowmap", shadowmapAtlas->depth_texture, 8);
	deferred_global->setUniform("u_shadowmap_dimensions", getSMapdDimensions(numShadowmaps));

	//initial pass (color + ambient + emissive + occlusion)
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_ALWAYS);
	baseRenderMP(quad, deferred_global); //this forces the deferred depthmap into the linear fbo
	deferred_global->disable();


	//render reflections if necessary
	if (render_reflections) {
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		glDepthMask(false);
		GFX::Shader* reflections = GFX::Shader::Get("reflection");
		reflections->enable();
		cameraToShader(camera, reflections);
		reflections->setUniform("u_invRes", vec2(1.0 / size.x, 1.0 / size.y));
		reflections->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
		reflections->setUniform("u_render", gBuffersFBO->color_textures[0], 0);
		reflections->setUniform("u_normals_gbuffer", gBuffersFBO->color_textures[1], 1);
		reflections->setUniform("u_mat_properties_gbuffer", gBuffersFBO->color_textures[2], 2);
		reflections->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 3);

		//pass probe data to gpu
		Vector3f probe_positions[NUM_REFL_PROBES];
		for (int i = 0; i < reflection_probes.size(); i++) {
			sReflectionProbe* current_probe = reflection_probes.at(i);
			probe_positions[i] = current_probe->pos;
			GFX::Texture* reflection_texture = current_probe == nullptr ? skybox_cubemap : current_probe->cubemap;

			char buffer[50];  // Temporary buffer to hold the concatenated string
			std::snprintf(buffer, sizeof(buffer), "u_reflection_%d", i); // Safely format the string
			// Allocate memory for the char* name
			char* name = new char[std::strlen(buffer) + 1];
			// Copy the formatted string into name
			std::strcpy(name, buffer);

			reflections->setUniform(name, reflection_texture, 4 + i);
		}
		reflections->setUniform3Array("u_probe_pos", (float*)&probe_positions, NUM_REFL_PROBES);

		quad->render(GL_TRIANGLES);
		reflections->disable();
		glDepthMask(true);
	}


	//variables and flags for illumination pass
	deferred_global->enable();
	GFX::Mesh light_sphere;
	Vector3f light_pos;
	Matrix44 sphere_model;
	float radius = 0.0f;
	light_sphere.createSphere(1);
	glEnable(GL_BLEND);
	glEnable(GL_CULL_FACE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);


	//illumination pass
	for (int i = 0; i < lights.size(); i++) {
		//upload uniforms
		lightToShaderMP(lights[i], deferred_global); //IDEA: trim unnecesary uniform uploads with a new, more specific function
		//point and spotlight (rendering geometry + frustum culling)
		if ((lights[i]->light_type == eLightType::POINT || lights[i]->light_type == eLightType::SPOT) && camera->testSphereInFrustum(lights[i]->root.model.getTranslation(), lights[i]->max_distance) == CLIP_INSIDE) {
			//update sphere to current light 
			//hack: setscale and settranslation break the program - update matrix manually
			sphere_model.setIdentity();
			//translation
			light_pos = lights[i]->root.model.getTranslation();
			sphere_model.m[12] = light_pos.x;
			sphere_model.m[13] = light_pos.y;
			sphere_model.m[14] = light_pos.z;
			//scale
			radius = lights[i]->max_distance;
			sphere_model.m[0] = radius;
			sphere_model.m[5] = radius;
			sphere_model.m[10] = radius;
			//check if camera is inside light sphere
			if (camera->eye.distance(light_pos) < lights[i]->max_distance) {
				glFrontFace(GL_CW);
			}
			//upload uniforms and render
			deferred_global->setUniform("u_is_quad", 0);
			deferred_global->setUniform("u_model", sphere_model);
			light_sphere.render(GL_TRIANGLES);
			//reset cull face method (prevent rendering mistakes)
			glFrontFace(GL_CCW); //IDEA: only do this if inside current light to optimize
		}
		//directional light
		else {
			deferred_global->setUniform("u_is_quad", 1);
			quad->render(GL_TRIANGLES);
		}
	}
	deferred_global->disable();


	//render probes if necessary (important to do before rendering alpha nodes)
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDisable(GL_BLEND);
	if (show_irr_probes) {
		renderAllProbes(irr_probe_size);
	}
	if (show_ref_probes) {
		for (int i = 0; i < reflection_probes.size(); i++) { //TODO: create method to iterate. Or not.
			renderReflectionProbe(*reflection_probes.at(i), 10.0);
		}
	}

	//render semitransparent objects using forward rendering
	glDepthMask(false);
	for (int i = 0; i < semitransparent_objects.size(); i++)
	{
		renderNode(semitransparent_objects[i], camera, eRenderTypes::FORWARD); //this shouldnt override the depthmask
	}
	glDepthMask(true);

	//irradiance pass
	if (probes_texture && use_irradiance) {
		deferred_display == eDeferredDisplay::IRRADIANCE ? glDisable(GL_BLEND) : glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glDepthMask(false);

		GFX::Shader* irradiance_shader = GFX::Shader::Get("irradiance");
		irradiance_shader->enable();

		// we send every data necessary
		cameraToShader(camera, irradiance_shader);

		irradiance_shader->setUniform("u_irr_start", probes_info.start);
		irradiance_shader->setUniform("u_irr_end", probes_info.end);
		irradiance_shader->setUniform("u_irr_dims", probes_info.dim);
		irradiance_shader->setUniform("u_irr_delta", probes_info.delta);
		irradiance_shader->setUniform("u_num_probes", probes_info.num_probes);
		irradiance_shader->setUniform("u_irr_normal_distance", irr_normal_distance);

		irradiance_shader->setUniform("u_color_texture", gBuffersFBO->color_textures[0], 0);
		irradiance_shader->setUniform("u_normal_texture", gBuffersFBO->color_textures[1], 1);
		irradiance_shader->setUniform("u_mat_properties_texture", gBuffersFBO->color_textures[2], 2);
		irradiance_shader->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 3);

		irradiance_shader->setUniform("u_probes_texture", probes_texture, 4);

		irradiance_shader->setUniform("u_probe_depth", linear_fbo->depth_texture, 5); //depth comparisons using gl_FragCoord are driving me insane

		irradiance_shader->setUniform("u_irr_factor", irr_factor);

		irradiance_shader->setUniform("u_invRes", vec2(1.0 / size.x, 1.0 / size.y));
		irradiance_shader->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
		quad->render(GL_TRIANGLES);

		glDisable(GL_BLEND);
		glDepthMask(true);

		//show results if necessary
		if (deferred_display == eDeferredDisplay::IRRADIANCE) {
			linear_fbo->unbind();
			linear_fbo->color_textures[0]->toViewport();
			return;
		}
	}

	linear_fbo->unbind();

	if (volumetric_light) {
		if (!volumFBO || prevScreenSize.distance(size) > 0.0) {
			delete volumFBO;
			delete fx_fbo;
			volumFBO = new GFX::FBO();
			volumFBO->create(size.x/2, size.y/2, 1, GL_RGBA, GL_FLOAT, false);
			volumFBO->color_textures[0]->setName("Volumetric light render");
			fx_fbo = new GFX::FBO();
			fx_fbo->create(size.x, size.y, 1, GL_RGBA, GL_FLOAT, false);
			fx_fbo->color_textures[0]->setName("FX FBO");
		}
		volumFBO->bind();
		glDisable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);

		//upload necessary shader uniforms
		GFX::Shader* vol_shader = GFX::Shader::Get("volumetric");
		vol_shader->enable();
		vol_shader->setUniform("u_depth_texture", linear_fbo->depth_texture, 0);
		vol_shader->setUniform("u_normal_texture", gBuffersFBO->color_textures[1], 1);
		vol_shader->setUniform("u_shadowmap", shadowmapAtlas->depth_texture, 2);
		vol_shader->setUniform("u_shadowmap_dimensions", getSMapdDimensions(numShadowmaps));
		vol_shader->setUniform("u_inverse_viewprojection", camera->inverse_viewprojection_matrix);
		cameraToShader(camera, vol_shader);
		vol_shader->setUniform("u_invRes", vec2(1.0 / volumFBO->width, 1.0 / volumFBO->height));
		vol_shader->setUniform("u_ambient_light", scene->ambient_light);
		vol_shader->setUniform("u_air_density", air_density);
		vol_shader->setUniform("u_time", (float)getTime() * 0.001f);
		lightToShaderSP(vol_shader);

		//render
		quad->render(GL_TRIANGLES);

		//configure linear interpolation
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		vol_shader->disable();

		volumFBO->unbind();

		//upscale volumetric render
		volumFBO->color_textures[0]->bind();
		volumFBO->color_textures[0]->toViewport();
		volumFBO->color_textures[0]->unbind();

		// blur raw volumetric render using neighbor pixels
		GFX::Shader* vol_blur_shader = circular_blur ? GFX::Shader::Get("blur_circular") : GFX::Shader::Get("blur_neighbors");
		vol_blur_shader->enable();
		fx_fbo->bind();
		vec2 invRes = vec2(1.0 / fx_fbo->color_textures[0]->width, 1.0 / fx_fbo->color_textures[0]->height);
		vol_blur_shader->setUniform("u_raw", volumFBO->color_textures[0], 0);
		vol_blur_shader->setUniform("u_depth_texture", gBuffersFBO->depth_texture, 1);
		vol_blur_shader->setUniform("u_invRes", invRes);
		vol_blur_shader->setUniform("u_blur_far", true);
		vol_blur_shader->setUniform("u_blur_dimensions", vec2(3.0f, 3.0f));
		vol_blur_shader->setUniform("u_blur_distance", vec2(2.0f, 2.0f));
		vol_blur_shader->setUniform("u_set_alpha", (float)0.0);
		glClear(GL_COLOR_BUFFER_BIT);
		quad->render(GL_TRIANGLES);
		fx_fbo->unbind();
		vol_blur_shader->disable();
		
		//show if necessary
		if (deferred_display == eDeferredDisplay::VOLUMETRIC) {
			blur_volumetric ? fx_fbo->color_textures[0]->toViewport() : volumFBO->color_textures[0]->toViewport();
			return;
		}

		//blend results to main render
		linear_fbo->bind();
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDepthMask(false);
		blur_volumetric ? fx_fbo->color_textures[0]->toViewport() : volumFBO->color_textures[0]->toViewport();
		glDepthMask(true);
		glDisable(GL_BLEND);
		linear_fbo->unbind();
	}

	partial_render = false;
}

void Renderer::renderSkybox(GFX::Texture* cubemap, float intensity)
{
	Camera* camera = Camera::current;

	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	GFX::Shader* shader = GFX::Shader::Get("skybox");
	if (!shader)
		return;
	shader->enable();

	Matrix44 m;
	m.setTranslation(camera->eye.x, camera->eye.y, camera->eye.z);
	m.scale(10, 10, 10);
	shader->setUniform("u_model", m);
	cameraToShader(camera, shader);
	shader->setUniform("u_texture", cubemap, 0);
	shader->setUniform("u_intensity", intensity);
	sphere.render(GL_TRIANGLES);
	shader->disable();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_DEPTH_TEST);
}

void Renderer::renderProbeFaces(SCN::Scene*, Camera* camera, bool render_skybox) {
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable(GL_CULL_FACE);
	GFX::checkGLErrors();

	if (render_skybox) {
		if (skybox_cubemap)
			renderSkybox(skybox_cubemap, 0.001f); //this is before color banding so value is really small to counter the colors (hardcoded to current scene TODO adjust)
	}

	//render entities
	for (int i = 0; i < opaque_objects.size(); i++)
	{
		renderNode(opaque_objects[i], camera, eRenderTypes::FORWARD);
	}
	//render semitransparent entities
	for (int i = 0; i < semitransparent_objects.size(); i++)
	{
		renderNode(semitransparent_objects[i], camera, eRenderTypes::FORWARD);
	}
}

//renders a node
void Renderer::renderNode(SCN::Node* node, Camera* camera, eRenderTypes render_type)
{
	if (!node->visible)
		return;

	//compute global matrix
	Matrix44 node_model = node->getGlobalMatrix(true);

	//does this node have a mesh? then we must render it
	if (node->mesh && node->material)
	{
		//compute the bounding box of the object in world space (by using the mesh bounding box transformed to world space)
		BoundingBox world_bounding = transformBoundingBox(node_model,node->mesh->box);
		
		//if bounding box is inside the camera frustum then the object is probably visible
		if (camera->testBoxInFrustum(world_bounding.center, world_bounding.halfsize) )
		{
			if(render_boundaries)
				node->mesh->renderBounding(node_model, true);
			switch (render_type) {
				case eRenderTypes::FLAT: renderMeshWithMaterial(node_model, node->mesh, node->material);
					break;
				case eRenderTypes::FORWARD: renderMeshWithMaterialLights(node_model, node->mesh, node->material, render_mode == eRenderTypes::DEFERRED ? false : gui_use_multipass); //little hack to optimize deferred performance using single pass
					break;
				case eRenderTypes::DEFERRED: renderMeshWithMaterialGBuffers(node_model, node->mesh, node->material);
					break;
			}
		}
	}
}

void Renderer::categorizeNodes(SCN::Node* node, Camera* camera) { //adds node and children nodes to their respective container

	if (node->material && node->material->alpha_mode == SCN::eAlphaMode::BLEND) { //objects with transparency
		Matrix44 global_model = node->getGlobalMatrix(); //'fast' option does not work - global_model not set correctly?
		float dist = std::sqrt(std::pow(global_model.getTranslation().x - camera->eye.x, 2) + std::pow(global_model.getTranslation().y - camera->eye.y, 2) + std::pow(global_model.getTranslation().z - camera->eye.z, 2));
		node->distance_to_camera = dist;
		semitransparent_objects.push_back(node);
	}
	else //other objects
	{
		//distance is only really used for semitransparent nodes, ignore here to optimize resources
		/*Matrix44 global_model = node->getGlobalMatrix();
		float dist = std::sqrt(std::pow(global_model.getTranslation().x - camera->eye.x, 2) + std::pow(global_model.getTranslation().y - camera->eye.y, 2) + std::pow(global_model.getTranslation().z - camera->eye.z, 2));
		node->distance_to_camera = dist; */
		opaque_objects.push_back(node);
	}
	//iterate recursively with children
	for (int i = 0; i < node->children.size(); ++i) {
		categorizeNodes(node->children[i], camera);
	}
}

//renders a mesh given its transform and material
void Renderer::renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material )
		return;
    assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;
	
	texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	if (texture == NULL)
		texture = GFX::Texture::getWhiteTexture(); //a 1x1 white texture

	//select the blending
	//lab1: send to global renderable vector, postpone rendering
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else {
		glDisable(GL_BLEND);
	}
	//select if render both sides of the triangles
	if(material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);


    assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("texture");

    assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t );

	shader->setUniform("u_color", material->color);
	if(texture)
		shader->setUniform("u_texture", texture, 0);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	if (render_wireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

void Renderer::renderMeshWithMaterialLights(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, bool use_multipass)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* colorTexture = NULL;
	GFX::Texture* normalMap = NULL;
	GFX::Texture* emissive = NULL;
	GFX::Texture* occlusion = NULL;
	GFX::Texture* metal_roughness = NULL; //contains occlusion in red channel
	Camera* camera = Camera::current;

	colorTexture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	normalMap = material->textures[SCN::eTextureChannel::NORMALMAP].texture; int useNormalmap = 1;
	emissive = material->textures[SCN::eTextureChannel::EMISSIVE].texture; int useEmissive = 1;
	occlusion = material->textures[SCN::eTextureChannel::OCCLUSION].texture; int useOcclusion = 1;
	metal_roughness = material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture; int useSpecular = 1;
	

	//get dummy textures if anything is missing
	if (colorTexture == NULL) 
		colorTexture = GFX::Texture::getWhiteTexture(); //a 1x1 white texture
	if (normalMap == NULL) {
		normalMap = GFX::Texture::getWhiteTexture();
		useNormalmap = 0;
	}
	if (emissive == NULL) {
		emissive = GFX::Texture::getWhiteTexture();
		if (material->emissive_factor.length() == 0.0) { useEmissive = 0; } //some emissive objects dont have emissive textures
	}
	if (occlusion == NULL) {
		occlusion = GFX::Texture::getWhiteTexture();
		if (metal_roughness == NULL) { //occlusions can be here too
			useOcclusion = 0;
		}
	}
	if (metal_roughness == NULL) {
		metal_roughness = GFX::Texture::getWhiteTexture();
		if (material->metallic_factor == 0 && material->roughness_factor == 1) 
		{ useSpecular = 0; }
	}

	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else {
		glDisable(GL_BLEND);
	}
	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);


	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = use_multipass ? GFX::Shader::Get("lightMP") : GFX::Shader::Get("lightSP"); 

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);
	shader->setUniform("u_ambient_light", scene->ambient_light);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_metal_factor", material->metallic_factor);
	shader->setUniform("u_rough_factor", material->roughness_factor);

	shader->setUniform("u_color", material->color);

	shader->setUniform("u_texture", colorTexture, 0);
	shader->setUniform("u_normalmap", normalMap, 1);
	useNormalmap = gui_use_normalmaps ? useNormalmap : 0;
	shader->setUniform("u_use_normalmap", useNormalmap);
	shader->setUniform("u_emissive", emissive, 2);
	useEmissive = gui_use_emissive ? useEmissive : 0;
	shader->setUniform("u_use_emissive", useEmissive);
	shader->setUniform("u_occlusion", occlusion, 3);
	shader->setUniform("u_metal_roughness", metal_roughness, 4);
	useOcclusion = gui_use_occlusion ? useOcclusion : 0;
	useSpecular = gui_use_specular ? useSpecular : 0;
	shader->setUniform("u_use_occlusion", useOcclusion);
	shader->setUniform("u_use_specular", useSpecular);
	shader->setUniform("u_shadowmap", shadowmapAtlas->depth_texture, 8);
	shader->setUniform("u_shadowmap_dimensions", getSMapdDimensions(numShadowmaps));

	shader->setUniform("u_use_shadowmaps", gui_use_shadowmaps);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	if (use_multipass) {
		if (material->alpha_mode != SCN::eAlphaMode::BLEND) {
			glDisable(GL_BLEND);
		}
		glDepthFunc(GL_LEQUAL);
		baseRenderMP(mesh, shader);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		for (int i = 0; i < lights.size(); i++) {
			lightToShaderMP(lights[i], shader);
			mesh->render(GL_TRIANGLES);
		}
		glDepthFunc(GL_LESS);
	}
	else {
		lightToShaderSP(shader);
		mesh->render(GL_TRIANGLES);	//do the draw call that renders the mesh into the screen
	}

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Renderer::renderMeshWithMaterialGBuffers(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	GFX::Texture* texture = NULL;
	GFX::Texture* normalMap = NULL;
	GFX::Texture* metal_roughness = NULL;
	GFX::Texture* emissive = NULL;


	Camera* camera = Camera::current;

	int useNormalmap = 1; //to avoid visual bugs when an object has no normalmap
	texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	normalMap = material->textures[SCN::eTextureChannel::NORMALMAP].texture;
	emissive = material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	metal_roughness = material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	if (texture == NULL) {
		texture = GFX::Texture::getWhiteTexture(); //a 1x1 white texture
	}
	if (normalMap == NULL) {
		normalMap = GFX::Texture::getWhiteTexture();
		useNormalmap = 0;
	}
	if (emissive == NULL) {
		emissive = GFX::Texture::getWhiteTexture();
	}
	if (metal_roughness == NULL) {
		metal_roughness = GFX::Texture::getWhiteTexture();
	}

	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else {
		glDisable(GL_BLEND);
	}
	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);


	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("gbuffers");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_metal_factor", material->metallic_factor);
	shader->setUniform("u_rough_factor", material->roughness_factor);

	shader->setUniform("u_texture", texture, 0);
	shader->setUniform("u_normalmap", normalMap, 1);
	shader->setUniform("u_emissive", emissive, 2);
	shader->setUniform("u_metal_roughness", metal_roughness, 3);

	shader->setUniform("u_use_normalmap", useNormalmap);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void SCN::Renderer::cameraToShader(Camera* camera, GFX::Shader* shader)
{
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix );
	shader->setUniform("u_camera_position", camera->eye);
}

void SCN::Renderer::lightToShaderSP(GFX::Shader* shader) {
	Vector3f light_positions[MAX_LIGHTS_SP];
	Vector3f light_fronts[MAX_LIGHTS_SP];
	Vector3f light_colors[MAX_LIGHTS_SP];
	Vector2f cones_info[MAX_LIGHTS_SP];
	float max_distances[MAX_LIGHTS_SP];
	int light_types[MAX_LIGHTS_SP];
	int lights_cast_shadows[MAX_LIGHTS_SP];
	mat4 shadow_viewprojections[MAX_LIGHTS_SP];
	float shadow_biases[MAX_LIGHTS_SP];
	int shadowmap_indexes[MAX_LIGHTS_SP];
	int num_lights = lights.size();
	for (int i = 0; i < num_lights; i++) {
		light_positions[i] = lights[i]->root.model.getTranslation();
		light_fronts[i] = lights[i]->root.model.frontVector().normalize();
		light_colors[i] = lights[i]->color * lights[i]->intensity; //IDEA: this should be done in gamma space, send to gpu
		Vector2f currentCone = Vector2f(cos(lights[i]->cone_info.x * (PI/180.0)), cos(lights[i]->cone_info.y * (PI / 180.0))); //convert angles to cosinus
		cones_info[i] = currentCone;
		max_distances[i] = lights[i]->max_distance;
		light_types[i] = (int)lights[i]->light_type;
		lights_cast_shadows[i] = (int)lights[i]->cast_shadows;
		shadow_viewprojections[i] = lights[i]->shadowmap_viewprojection;
		shadow_biases[i] = lights[i]->shadow_bias;
		shadowmap_indexes[i] = lights[i]->shadowmap_index;
	}
	shader->setUniform("u_num_lights", num_lights);
	shader->setUniform3Array("u_light_pos", (float*)&light_positions, MAX_LIGHTS_SP);
	shader->setUniform3Array("u_light_front", (float*)&light_fronts, MAX_LIGHTS_SP);
	shader->setUniform3Array("u_light_col", (float*)&light_colors, MAX_LIGHTS_SP);
	shader->setUniform2Array("u_cone_info", (float*)&cones_info, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_max_distance", (float*)&max_distances, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_light_type", (int*)&light_types, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_light_cast_shadows", (int*)&lights_cast_shadows, MAX_LIGHTS_SP);
	shader->setMatrix44Array("u_shadow_viewproj", (Matrix44*)&shadow_viewprojections, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_shadow_bias", (float*)&shadow_biases, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_shadowmap_index", (int*)&shadowmap_indexes, MAX_LIGHTS_SP);
}

void SCN::Renderer::lightToShaderMP(LightEntity* light, GFX::Shader* shader) {
	Vector3f light_position = light->root.model.getTranslation();
	Vector3f light_front = light->root.model.frontVector().normalize();
	Vector3f light_color = light->color * light->intensity;
	Vector2f cone_info = Vector2f(cos(light->cone_info.x * (PI / 180.0)), cos(light->cone_info.y * (PI / 180.0))); //convert angle to cosinus
	float max_distance = light->max_distance;
	int light_type = (int)light->light_type;
	int light_cast_shadows = light->cast_shadows;
	Matrix44 shadow_viewprojection = light->shadowmap_viewprojection;
	float shadow_bias = light->shadow_bias;
	int shadowmap_index = light->shadowmap_index;

	shader->setUniform("u_light_pos", light_position);
	shader->setUniform("u_light_front", light_front);
	shader->setUniform("u_light_col", light_color);
	shader->setUniform("u_cone_info", cone_info);
	shader->setUniform("u_max_distance", max_distance);
	shader->setUniform("u_light_type", light_type);
	shader->setUniform("u_light_cast_shadows", light_cast_shadows);
	shader->setUniform("u_shadow_viewproj", shadow_viewprojection);
	shader->setUniform("u_shadow_bias", shadow_bias);
	shader->setUniform("u_shadowmap_index", shadowmap_index);
}

void SCN::Renderer::baseRenderMP(GFX::Mesh* mesh, GFX::Shader* shader) {
	int light_type = 4; //defined as ambient light (u_ambient_light alredy passed to shader)
	shader->setUniform("u_light_type", light_type);
	mesh->render(GL_TRIANGLES);
}

void SCN::Renderer::renderAllProbes(float size) {
	for (int iP = 0; iP < probes_info.num_probes; ++iP)
	{
		renderProbe(probes.at(iP).pos, size, probes.at(iP).sh.coeffs->v);
	}
}

void SCN::Renderer::renderProbe(vec3 pos, float size, float* coeffs)
{
	Camera* camera = Camera::current;
	GFX::Shader* shader = GFX::Shader::Get("probe");

	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	Matrix44 model;
	model.setTranslation(pos.x, pos.y, pos.z);
	model.scale(size, size, size);

	shader->enable();
	cameraToShader(camera, shader);
	shader->setUniform("u_model", model);
	shader->setUniform3Array("u_coeffs", coeffs, 9);

	sphere.render(GL_TRIANGLES);
}

void SCN::Renderer::captureAllProbes() {
	//set coefficients
	for (int iP = 0; iP < probes_info.num_probes; ++iP)
	{
		captureProbe(probes.at(iP));
	}
	renderIrradianceTexture();
}

void SCN::Renderer::renderIrradianceTexture() {

	//create the texture to store the probes
	if (!probes_texture) {
		probes_texture = new GFX::Texture(
			9, //9 coefficients per probe
			probes.size(), //as many rows as probes
			GL_RGB, //3 channels per coefficient
			GL_FLOAT); //they require a high range
	}

	//we must create the color information for the texture. because every SH are 27 floats in the RGB,RGB,... order, we can create an array of SphericalHarmonics and use it as pixels of the texture
	SphericalHarmonics* sh_data = NULL;
	sh_data = new SphericalHarmonics[probes_info.dim.x * probes_info.dim.y * probes_info.dim.z];

	//here we fill the data of the array with our probes in x,y,z order
	for (int i = 0; i < probes.size(); ++i)
		sh_data[i] = probes[i].sh;

	//now upload the data to the GPU as a texture
	probes_texture->upload(GL_RGB, GL_FLOAT, false, (uint8*)sh_data);

	//disable any texture filtering when reading
	probes_texture->bind();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	//always free memory after allocating it!!!
	delete[] sh_data;
}

void SCN::Renderer::captureProbe(sProbe& p) {
	FloatImage images[6]; //here we will store the six views
	GFX::Texture temp;
	temp.create(64, 64, GL_RGB, GL_FLOAT);

	for (int i = 0; i < 6; ++i) //for every cubemap face
	{
		multi_probes_fbo->setTexture(&temp);
		//compute camera orientation using defined vectors
		vec3 eye = p.pos;
		vec3 front = cubemapFaceNormals[i][2];
		vec3 center = p.pos + front;
		vec3 up = cubemapFaceNormals[i][1];
		probeCam->lookAt(eye, center, up);
		probeCam->enable();
		//render the scene from this point of view
		multi_probes_fbo->bind();
		renderProbeFaces(scene, probeCam, true);
		multi_probes_fbo->unbind();

		//read the pixels back and store in a FloatImage
		images[i].fromTexture(multi_probes_fbo->color_textures[0]);
	}
	//compute the coefficients given the six images
	p.sh = computeSH(images);
}

void SCN::Renderer::renderReflectionProbe(sReflectionProbe& s, float scale) {
	Camera* camera = Camera::current;
	GFX::Shader* shader = GFX::Shader::Get("reflection_probe");

	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	GFX::Texture* reflection_texture = s.cubemap == nullptr ? skybox_cubemap : s.cubemap;


	Matrix44 model;
	model.setTranslation(s.pos.x, s.pos.y, s.pos.z);
	model.scale(scale, scale, scale);

	shader->enable();
	cameraToShader(camera, shader);
	shader->setUniform("u_model", model);
	shader->setUniform("u_environment_texture", reflection_texture, 0);
	shader->setUniform("u_linearize_colors", (int)(s.cubemap == nullptr));

	sphere.render(GL_TRIANGLES);
}

void SCN::Renderer::captureReflectionProbe(sReflectionProbe& s) {
	if (!s.cubemap) {
		s.cubemap = new GFX::Texture();
		s.cubemap->createCubemap(356,356, nullptr, GL_RGB, GL_FLOAT, true);
	}

	//render the view from every side
	for (int i = 0; i < 6; ++i)
	{
		//assign cubemap face to FBO
		multi_probes_fbo->setTexture(s.cubemap, i);

		//bind FBO
		multi_probes_fbo->bind();

		//render view
		vec3 eye = s.pos;
		vec3 center = s.pos + cubemapFaceNormals[i][2];
		vec3 up = cubemapFaceNormals[i][1];
		probeCam->lookAt(eye, center, up);
		probeCam->enable();
		renderProbeFaces(scene, probeCam, true);
		multi_probes_fbo->unbind();
	}
	//generate the mipmaps
	glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	s.cubemap->generateMipmaps();

}


#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ImGui::Checkbox("Wireframe", &render_wireframe);
	ImGui::Checkbox("Boundaries", &render_boundaries);
	ImGui::Checkbox("Disable lights", &disable_lights);
	ImGui::SliderInt("10-5000", &gui_shadowmap_res, 10, 5000, "shadowmap resolution");
	ImGui::Checkbox("Circular blurring", &circular_blur);
	if (render_mode == eRenderTypes::FLAT) {
		if (ImGui::BeginCombo("Render mode", "Flat")) {
			ImGui::RadioButton("Flat", (int*)&render_mode, (int)eRenderTypes::FLAT);
			ImGui::RadioButton("Forward", (int*)&render_mode, (int)eRenderTypes::FORWARD);
			ImGui::RadioButton("Deferred", (int*)&render_mode, (int)eRenderTypes::DEFERRED);
			ImGui::EndCombo();
		}
	}
	else if (render_mode == eRenderTypes::FORWARD) {
		if (ImGui::BeginCombo("Render mode", "Forward")) {
			ImGui::RadioButton("Flat", (int*) & render_mode, (int)eRenderTypes::FLAT);
			ImGui::RadioButton("Forward", (int*)&render_mode, (int)eRenderTypes::FORWARD);
			ImGui::RadioButton("Deferred", (int*)&render_mode, (int)eRenderTypes::DEFERRED);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Lab1", "Show options")) {
			ImGui::Checkbox("Multipass lights", &gui_use_multipass);
			ImGui::Checkbox("use normalmaps", &gui_use_normalmaps);
			ImGui::Checkbox("use emissive", &gui_use_emissive);
			ImGui::Checkbox("use occlusion", &gui_use_occlusion);
			ImGui::Checkbox("use specular", &gui_use_specular);
			ImGui::Checkbox("use shadowmaps", &gui_use_shadowmaps);
			ImGui::EndCombo();
		}
	}
	else if (render_mode == eRenderTypes::DEFERRED) {
		if (ImGui::BeginCombo("Render mode", "Deferred")) {
			ImGui::RadioButton("Flat", (int*)&render_mode, (int)eRenderTypes::FLAT);
			ImGui::RadioButton("Forward", (int*)&render_mode, (int)eRenderTypes::FORWARD);
			ImGui::RadioButton("Deferred", (int*)&render_mode, (int)eRenderTypes::DEFERRED);
			ImGui::EndCombo();
		}
		ImGui::Combo("Display channel", (int*)&deferred_display, "DEFAULT\0COLOR\0NORMALS\0MATERIAL_PROPERTIES\0DEPTH\0EMISSIVE\0SSAO_result\0IRRADIANCE\0VOLUMETRIC\0");
		if (ImGui::BeginCombo("Occlusion", "Show options")) {
			ImGui::Combo("Occlusion mode", (int*)&occlusion_mode, "TEXTURE\0SSAO\0SSAOplus\0");
			ImGui::DragFloat("SSAO radius", &ssao_radius, 0.01f, 0.0f, 20.0f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Irradiance", "Show options")) {
			ImGui::Checkbox("Use irradience", &use_irradiance);
			if (ImGui::Button("Capture Irradiance")) {
				//compute probe coefficients and save to disk
				captureAllProbes();
				FILE* f = fopen("irradiance.bin", "wb");
				fwrite(&probes_info, sizeof(sIrradianceInfo), 1, f);
				fwrite(&(probes[0]), sizeof(sProbe), probes.size(), f);
				fclose(f);
			}
			if (ImGui::Button("Delete irradiance data")) {
				//compute probe coefficients and save to disk
				if (remove("irradiance.bin") != 0) {
					std::cout << "irradiance data not found" << std::endl;
				}
				else {
					std::cout << "succesfully removed irradiance data" << std::endl;
				}
			}
			ImGui::Checkbox("Show irradiance probes", &show_irr_probes);
			ImGui::DragFloat("Probe size", &irr_probe_size, 0.01f, 1.0f, 10.0f);
			ImGui::DragFloat("Irradiance weight", &irr_factor, 0.01f, 0.01f, 10.0f);
			ImGui::DragFloat("Irradiance normal factor", &irr_normal_distance, 0.01f, 0.01f, 15.0f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Reflections", "Show options")) {
			ImGui::Checkbox("Render reflections", &render_reflections);
			ImGui::Checkbox("Show reflection probes", &show_ref_probes);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Volumetric", "Show options")) {
			ImGui::Checkbox("Volumetric lights", &volumetric_light);
			ImGui::DragFloat("Air density", &air_density, 0.001f, 0.0f, 10.0f);
			ImGui::Checkbox("Blurring for volumetric render", &blur_volumetric);
			ImGui::EndCombo();
		}
	}
	if (ImGui::BeginCombo("PostFX", "Show options")) {
		ImGui::Checkbox("Color banding", &color_banding);
		if (ImGui::BeginCombo("Chromatic Aberration", "Show options")) {
			ImGui::Checkbox("Chromatic aberration", &use_chromatic_aberration);
			ImGui::DragFloat("CA intensity", &chromatic_aberration_factor, 0.01f, 0.0f, 1.0f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Depth of Field", "Show options")) {
			ImGui::Checkbox("Depth of Field", &use_dof);
			ImGui::DragFloat("DoF - min distance", &dof_min, 1.0f, 0.0f, dof_max - 1.0);
			ImGui::DragFloat("DoF - max distance", &dof_max, 1.0f, dof_min, 5000.00f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Blur", "Show options")) {
			ImGui::Checkbox("Blur render", &blur_render);
			ImGui::DragFloat("Blur (X)", &fx_blur_res.x, 1.00f, 1.0f, 100.0f);
			ImGui::DragFloat("Blur (Y)", &fx_blur_res.y, 1.00f, 1.0f, 100.0f);
			ImGui::DragFloat("Blur dist (X)", &fx_blur_dist.x, 1.00f, 1.0f, 10.0f);
			ImGui::DragFloat("Blur dist (Y)", &fx_blur_dist.y, 1.00f, 1.0f, 10.0f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Frame smoothing", "Show options")) {
			ImGui::Checkbox("Frame smoothing", &use_motion_blur);
			ImGui::Checkbox("Frame smoothing burn-in", &frame_smoothing_burnin);
			ImGui::DragFloat("Frame smoothing intensity", &motion_blur_intensity, 0.01f, 0.01f, 0.99f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Bloom", "Show options")) {
			ImGui::Checkbox("Activate bloom", &use_bloom);
			ImGui::DragFloat("Bloom threshold", &bloom_threshold, 0.01f, 0.01f, 2.0f);
			ImGui::DragInt("Bloom iterations", &bloom_iterations, 1, 1, 10);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Grain", "Show options")) {
			ImGui::Checkbox("Grain filter", &use_grain);
			ImGui::DragFloat("Grain intensity", &grain_intensity, 0.001f, 0.00f, 1.0f);
			ImGui::DragFloat("Grain size", &grain_size, 0.1f, 0.00f, 10.0f);
			ImGui::EndCombo();
		}
		if (ImGui::BeginCombo("Lens Distortion", "Show options")) {
			ImGui::Checkbox("Lens distortion", &use_lens_distortion);
			if (ImGui::BeginCombo("select mode", "Lens distortion")) {
				if (ImGui::Button("Pincushion")) {
					lens_distortion_mode = 0;
				}
				else if (ImGui::Button("Barrel")) {
					lens_distortion_mode = 1;
				}
				else if (ImGui::Button("Fisheye")) {
					lens_distortion_mode = 2;
				}
				else if (ImGui::Button("CRT")) {
					lens_distortion_mode = 3;
				}
				ImGui::EndCombo();
			}
			ImGui::DragFloat("LD intensity", &lens_distortion_intensity, 0.01f, 0.0f, 1.0f);
			ImGui::EndCombo();
		}
		ImGui::EndCombo();
	}
	ImGui::Checkbox("Use tonemapper", &gui_use_tonemapper);
	if (gui_use_tonemapper) {
		if (ImGui::BeginCombo(" ", "Tonemapper properties")) {
			ImGui::DragFloat("color scale", &tmp_scale, 0.005f, 0.0f, 2.0f);
			ImGui::DragFloat("average lum", &tmp_avg_lum, 0.005f, 0.0f, 2.0f);
			ImGui::DragFloat("lumwhite", &tmp_lumwhite, 0.005f, 0.0f, 2.0f);
			ImGui::EndCombo();
		}
	}
}

#else
void Renderer::showUI() {}
#endif