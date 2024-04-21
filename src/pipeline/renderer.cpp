#include "renderer.h"

#include <algorithm> //sort

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

//struct sRenderable {
//	Material* material;
//	GFX::Mesh* mesh;
//	mat4 model;
//	float distanceToCamera;
//};

std::vector<SCN::Node*> opaque_objects;
std::vector<SCN::Node*> semitransparent_objects;
std::vector<LightEntity*> lights;
LightEntity* mainLight; //sunlight or equivalent (directional) light

bool compareDist(Node* s1, Node* s2) { 
	return s1->distance_to_camera > s2->distance_to_camera;
}

Renderer::Renderer(const char* shader_atlas_filename)
{
	render_wireframe = false;
	render_boundaries = false;
	scene = nullptr;
	skybox_cubemap = nullptr;
	use_multipass = false;
	render_lights = true;
	disable_lights = false;
	shadowmapAtlas = nullptr; //TODO remove shadowmaps in lights

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
	sphere.uploadToVRAM();
}

void Renderer::setupScene()
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;

	//clear node and light containers
	mainLight = nullptr;
	lights.clear();
	semitransparent_objects.clear();
	opaque_objects.clear();

}

int getSMapdDimensions(int numlights) {
	return(ceil(sqrt((float)numlights)));
}

void Renderer::generateShadowmaps() {
	//IDEA: move directional lights to 'player' relative distance to better fit shadowmaps
	LightEntity* currentLight = mainLight;
	int shadowmap_size = 1024; //TODO: customizar tama�o en variable
	int shadowatlas_size = getSMapdDimensions(1) * shadowmap_size; //TODO 1 -> lights.size()
	if (currentLight != nullptr) {
		if (currentLight->cast_shadows && (shadowmapAtlas == nullptr || prevAtlasSize != shadowatlas_size)) {
			delete shadowmapAtlas;
			shadowmapAtlas = new GFX::FBO();
			shadowmapAtlas->setDepthOnly(shadowatlas_size, shadowatlas_size); 
		}
		shadowmapAtlas->bind();
		glClear(GL_DEPTH_BUFFER_BIT);
		Camera camera;
		vec3 position = currentLight->root.model.getTranslation();
		camera.lookAt(position, position + currentLight->root.model.frontVector() * -1.0f, vec3(0, 1, 0));
		//TODO: directional ortografic, spot perspective
		camera.setOrthographic(currentLight->area * -0.5, currentLight->area * 0.5, currentLight->area * -0.5, currentLight->area * 0.5, currentLight->near_distance, currentLight->max_distance);
		camera.enable();
		//render all opaque nodes
		for (Node* currentObj : opaque_objects) {
			BoundingBox world_bounding = transformBoundingBox(currentObj->model, currentObj->getBoundingBox());
			//frustum culling for shadowmaps
			if (camera.testBoxInFrustum(world_bounding.center, world_bounding.halfsize)) {
				renderShadowmap(currentObj->getGlobalMatrix(), currentObj->mesh, currentObj->material);
			}
		}
		shadowmapAtlas->unbind();
		currentLight->shadowmap_viewprojection = camera.viewprojection_matrix;
		prevAtlasSize = shadowatlas_size;
	}
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

void Renderer::renderScene(SCN::Scene* scene, Camera* camera)
{
	this->scene = scene;
	setupScene();

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if(skybox_cubemap)
		renderSkybox(skybox_cubemap);

	//pass 1: store entities in containers
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible )
			continue;

		//categorize all objects into containers
		if (ent->getType() == eEntityType::PREFAB) //prefabs
		{
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
			if (pent->prefab)
				categorizeNodes(&pent->root, camera);
		}
		else if (ent->getType() == eEntityType::LIGHT && !disable_lights) { //light objects
			//IDEA: test spheres to bounding boxes and cull invisible lights
			//downcast to EntityLight and store in light array if it affects objects in the scene
			LightEntity* light = (SCN::LightEntity*)ent; 
			if (light->light_type == eLightType::DIRECTIONAL || camera->testSphereInFrustum(light->root.model.getTranslation(), light->max_distance) == CLIP_INSIDE) { //simple frustum culling
				lights.push_back(light);
			}
			//assign first encountered directional light to main light
			if (!mainLight && light->light_type == eLightType::DIRECTIONAL) {
				mainLight = light;
			}
		}
	}
	generateShadowmaps(); //uses light container
	camera->enable(); //reactivate scene camera
	//pass 2: render entities
	for (int i = 0; i < opaque_objects.size(); i++)
	{
		renderNode(opaque_objects[i], camera);
	}
	//render semitransparent entities
	//sort blending vector - sorts nodes by distance in descending order
	std::sort(std::begin(semitransparent_objects), std::end(semitransparent_objects), compareDist);
	for (int i = 0; i < semitransparent_objects.size(); i++)
	{
		renderNode(semitransparent_objects[i], camera);
	}
}


void Renderer::renderSkybox(GFX::Texture* cubemap)
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
	sphere.render(GL_TRIANGLES);
	shader->disable();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_DEPTH_TEST);
}

//renders a node of the prefab and its children
void Renderer::renderNode(SCN::Node* node, Camera* camera)
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
			render_lights ? renderMeshWithMaterialLights(node_model, node->mesh, node->material) : renderMeshWithMaterial(node_model, node->mesh, node->material);
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
	//texture = material->emissive_texture;
	//texture = material->metallic_roughness_texture;
	//texture = material->normal_texture;
	//texture = material->occlusion_texture;
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


void Renderer::renderMeshWithMaterialLights(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
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
	int num_lights = lights.size();
	for (int i = 0; i < num_lights; i++) {
		light_positions[i] = lights[i]->root.model.getTranslation();
		light_fronts[i] = lights[i]->root.model.frontVector().normalize();
		light_colors[i] = lights[i]->color * lights[i]->intensity;
		Vector2f currentCone = Vector2f(cos(lights[i]->cone_info.x * (PI/180.0)), cos(lights[i]->cone_info.y * (PI / 180.0))); //convert angles to cosinus
		cones_info[i] = currentCone;
		max_distances[i] = lights[i]->max_distance;
		light_types[i] = (int)lights[i]->light_type;
		lights_cast_shadows[i] = (int)lights[i]->cast_shadows;
		shadow_viewprojections[i] = lights[i]->shadowmap_viewprojection;
		shadow_biases[i] = lights[i]->shadow_bias;

	}
	shader->setUniform("u_num_lights", num_lights);
	shader->setUniform3Array("u_light_pos", (float*)&light_positions, MAX_LIGHTS_SP);
	shader->setUniform3Array("u_light_front", (float*)&light_fronts, MAX_LIGHTS_SP);
	shader->setUniform3Array("u_light_col", (float*)&light_colors, MAX_LIGHTS_SP);
	shader->setUniform2Array("u_cone_info", (float*)&cones_info, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_max_distance", (float*)&max_distances, MAX_LIGHTS_SP);
	shader->setUniform1Array("u_light_type", (int*)&light_types, MAX_LIGHTS_SP);
	//shader->setUniform1Array("u_light_cast_shadows", (int*)&lights_cast_shadows, MAX_LIGHTS_SP);
	//shader->setMatrix44Array("u_shadow_viewproj", (Matrix44*)&shadow_viewprojections, MAX_LIGHTS_SP);
	//shader->setUniform1Array("u_shadow_bias", (float*)&shadow_biases, MAX_LIGHTS_SP);
	shader->setUniform("u_shadow_viewproj", mainLight->shadowmap_viewprojection);
	shader->setUniform("u_shadow_bias", mainLight->shadow_bias);
	shader->setUniform("u_light_cast_shadows", 1);
}

void SCN::Renderer::lightToShaderMP(LightEntity* light, GFX::Shader* shader) {
	Vector3f light_position = light->root.model.getTranslation();
	Vector3f light_front = light->root.model.frontVector().normalize();
	Vector3f light_color = light->color * light->intensity;
	Vector2f cone_info = Vector2f(cos(light->cone_info.x * (PI / 180.0)), cos(light->cone_info.y * (PI / 180.0))); //convert angle to cosinus
	float max_distance = light->max_distance;
	int light_type = (int)light->light_type;

	shader->setUniform("u_light_pos", light_position);
	shader->setUniform("u_light_front", light_front);
	shader->setUniform("u_light_col", light_color);
	shader->setUniform("u_cone_info", cone_info);
	shader->setUniform("u_max_distance", max_distance);
	shader->setUniform("u_light_type", light_type);
}

void SCN::Renderer::baseRenderMP(GFX::Mesh* mesh, GFX::Shader* shader) {
	int light_type = 4; //defined as ambient light (u_ambient_light alredy passed to shader)
	shader->setUniform("u_light_type", light_type);
	mesh->render(GL_TRIANGLES);
}

#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ImGui::Checkbox("Wireframe", &render_wireframe);
	ImGui::Checkbox("Boundaries", &render_boundaries);
	ImGui::Checkbox("Render with lights", &render_lights);
	ImGui::Checkbox("Disable lights", &disable_lights);
	ImGui::Checkbox("Multipass lights", &use_multipass);
	ImGui::Checkbox("use normalmaps", &gui_use_normalmaps);
	ImGui::Checkbox("use emissive", &gui_use_emissive);
	ImGui::Checkbox("use occlusion", &gui_use_occlusion);
	ImGui::Checkbox("use specular", &gui_use_specular);



	//add here your stuff
	//...
}

#else
void Renderer::showUI() {}
#endif