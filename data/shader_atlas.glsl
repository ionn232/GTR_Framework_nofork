//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
lightSP basic.vs lightSP.fs
lightMP basic.vs lightMP.fs
gbuffers basic.vs gbuffers.fs
deferred_global deferred.vs deferred_global.fs
view_emissive quad.vs emissive.fs
ssao quad.vs ssao.fs
blur_reprojection quad.vs blur_reprojection.fs
blur_neighbors quad.vs blur_neighbors.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs
gamma quad.vs gamma.fs
tonemapper quad.vs tonemapper.fs

\basic.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;
in vec4 a_color;

uniform vec3 u_camera_pos;

uniform mat4 u_model;
uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;
out vec4 v_color;

uniform float u_time;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( v_position, 1.0) ).xyz;
	
	//store the color in the varying var to use it from the pixel shader
	v_color = a_color;

	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\quad.vs

#version 330 core

in vec3 a_vertex;
in vec2 a_coord;
out vec2 v_uv;

void main()
{	
	v_uv = a_coord;
	gl_Position = vec4( a_vertex, 1.0 );
}

\deferred.vs

#version 330 core

uniform int u_is_quad;
uniform mat4 u_model;
uniform mat4 u_viewprojection;

in vec3 a_vertex;
in vec2 a_coord;
out vec2 v_uv;

void main()
{	

	//draw light passes (directional and ambient)
	if (u_is_quad == 1) {
		v_uv = a_coord;
		gl_Position = vec4( a_vertex, 1.0 );
	}
	//draw light passes (point and spot)
	else{
		//calcule the vertex in object space
		vec3 position = a_vertex;
		vec3 world_position = (u_model * vec4( position, 1.0) ).xyz;

		//calcule the position of the vertex using the matrices
		gl_Position = u_viewprojection * vec4( world_position, 1.0 );

		//store the texture coordinates
		v_uv = a_coord;
	}
}

\flat.fs

#version 330 core

uniform vec4 u_color;

out vec4 FragColor;

void main()
{
	FragColor = u_color;
}


\texture.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );
	color.xyz = pow(color.xyz, vec3(2.2));

	if(color.a < u_alpha_cutoff)
		discard;

	FragColor = color;
}

\lightSP.fs

#version 330 core

const int MAX_LIGHTS = 15;
in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

uniform vec3 u_camera_position;

uniform sampler2D u_normalmap;
uniform int u_use_normalmap;

uniform sampler2D u_emissive;
uniform vec3 u_emissive_factor;
uniform int u_use_emissive;

uniform sampler2D u_occlusion;
uniform sampler2D u_metal_roughness;
uniform int u_use_occlusion;
uniform int u_use_specular;
uniform float u_metal_factor;
uniform float u_rough_factor;

uniform vec3 u_ambient_light;

uniform vec3 u_light_pos[MAX_LIGHTS];
uniform vec3 u_light_front[MAX_LIGHTS];
uniform vec3 u_light_col[MAX_LIGHTS];
uniform vec2 u_cone_info[MAX_LIGHTS];
uniform float u_max_distance[MAX_LIGHTS];
uniform int u_light_type[MAX_LIGHTS];
		//NO_LIGHT = 0,
		//POINT = 1,
		//SPOT = 2,
		//DIRECTIONAL = 3

uniform int u_num_lights;

uniform int u_use_shadowmaps;
uniform int u_light_cast_shadows[MAX_LIGHTS];
uniform sampler2D u_shadowmap;
uniform mat4 u_shadow_viewproj[MAX_LIGHTS];
uniform float u_shadow_bias[MAX_LIGHTS];
uniform int u_shadowmap_index[MAX_LIGHTS];
uniform int u_shadowmap_dimensions;

out vec4 FragColor;

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx( p );
	vec3 dp2 = dFdy( p );
	vec2 duv1 = dFdx( uv );
	vec2 duv2 = dFdy( uv );
	
	// solve the linear system
	vec3 dp2perp = cross( dp2, N );
	vec3 dp1perp = cross( N, dp1 );
	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;
 
	// construct a scale-invariant frame 
	float invmax = inversesqrt( max( dot(T,T), dot(B,B) ) );
	return mat3( T * invmax, B * invmax, N );
}

// assume N, the interpolated vertex normal and 
// WP the world position
//vec3 normal_pixel = texture2D( normalmap, uv ).xyz; 
vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255./127. - 128./127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}

float computeShadowSP(vec3 wp, int i){
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj[i] * vec4(wp,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = (proj_pos.xy / proj_pos.w);

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);
	//it is outside on the sides
	if( shadow_uv.x < 0.0 || shadow_uv.x > 1.0 || shadow_uv.y < 0.0 || shadow_uv.y > 1.0 ) {
		return 1.0;
	}

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_bias[i]) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;
	//it is before near or behind far plane
	if(real_depth < 0.0 || real_depth > 1.0) {
		return 1.0;
	}

	//read depth from depth buffer in [0..+1] non-linear
	//accounts for offset using shadowmap dimensions and id
	float shadow_depth = texture( u_shadowmap, vec2(shadow_uv.x*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*(u_shadowmap_index[i]%u_shadowmap_dimensions),	 shadow_uv.y*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*floor(u_shadowmap_index[i]/u_shadowmap_dimensions))).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if( shadow_depth < real_depth ) {
		shadow_factor = 0.0;
	}
	return shadow_factor;
}

void main()
{
	
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );
	color.xyz =  pow(color.xyz, vec3(2.2));

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 light = vec3(0.0, 0.0, 0.0);
	
	float occlusion = texture( u_metal_roughness, v_uv).r;
	vec3 ambient_light = pow(u_ambient_light, vec3(2.2));
	if (u_use_occlusion == 1) {
		light += ambient_light * occlusion;
	}
	else{
		light += ambient_light;
	}

	vec3 V = normalize(u_camera_position - v_world_position);

	vec3 N = normalize(v_normal);
	vec3 normal_pixel = texture( u_normalmap, v_uv ).xyz;
	if (u_use_normalmap == 1) {
		N = perturbNormal(N,v_world_position, v_uv , normal_pixel);
	}

	vec3 emissive_pixel = texture( u_emissive, v_uv ).xyz;

	float spec_ks = texture( u_metal_roughness, v_uv).g;
	float spec_a =  texture( u_metal_roughness, v_uv).b;

	for (int i=0; i<MAX_LIGHTS; i++) {
		if (i<u_num_lights) {
			//linearize light color
			vec3 light_color =  pow(u_light_col[i], vec3(2.2));
			if (u_light_type[i] == 1) { 		//point lights
				//diffuse value
				vec3 L = u_light_pos[i] - v_world_position;
				L= normalize(L);
				float NdotL = clamp(dot(N, L), 0.0, 1.0);
				
				float lightDist = distance(u_light_pos[i], v_world_position);
				float att_factor = u_max_distance[i] - lightDist;
				att_factor = att_factor/u_max_distance[i];
				att_factor = max(att_factor, 0.0);

				light += (NdotL * light_color) * att_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;
					if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * light_color * att_factor; }
				}

			}
			else if (u_light_type[i] == 2) { 		//spot lights
				//diffuse value
				vec3 L = u_light_pos[i] - v_world_position;
				L= normalize(L);
				float NdotL = clamp(dot(N, L), 0.0, 1.0);
				
				float lightDist = distance(u_light_pos[i], v_world_position);
				float att_factor = u_max_distance[i] - lightDist;
				att_factor = att_factor/u_max_distance[i];
				att_factor = max(att_factor, 0.0);

				float cos_angle = dot(u_light_front[i], L);
				if (cos_angle < u_cone_info[i].y) {
					NdotL = 0.0;
				}
				else if (cos_angle < u_cone_info[i].x) {
					NdotL *= (cos_angle - u_cone_info[i].y) / (u_cone_info[i].x - u_cone_info[i].y);
				}

				//shadow value
				float shadow_factor = 1.0;
				if (u_light_cast_shadows[i] == 1 && u_use_shadowmaps == 1) {
					shadow_factor = computeShadowSP(v_world_position, i);
				}

				light += NdotL * light_color * att_factor * shadow_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;

					if (final_a != 0) {light += NdotL * final_ks * pow(NdotH, final_a) * light_color * att_factor * shadow_factor; }
				}
			}
			else if (u_light_type[i] == 3) {		//directional lights
				//diffuse value
				vec3 L = u_light_front[i];
				L = normalize(L);
				float NdotL = clamp(dot(N, L), 0.0, 1.0);

				//shadow value
				float shadow_factor = 1.0;
				if (u_light_cast_shadows[i] == 1 && u_use_shadowmaps == 1) {
					shadow_factor = computeShadowSP(v_world_position, i);
				}
				light += NdotL * light_color * shadow_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;
					if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * light_color * shadow_factor; }
				}
			}
		}
	}
	FragColor.xyz = color.xyz * light;
	if (u_use_emissive == 1) {FragColor.xyz += pow(emissive_pixel * u_emissive_factor, vec3(2.2));}
	FragColor.a = color.a;
}


\lightMP.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

uniform vec3 u_camera_position;

uniform sampler2D u_normalmap;
uniform int u_use_normalmap;

uniform sampler2D u_emissive;
uniform vec3 u_emissive_factor;
uniform int u_use_emissive;

uniform sampler2D u_occlusion;
uniform sampler2D u_metal_roughness;
uniform int u_use_occlusion;
uniform int u_use_specular;
uniform float u_metal_factor;
uniform float u_rough_factor;

uniform vec3 u_ambient_light;

uniform vec3 u_light_pos;
uniform vec3 u_light_front;
uniform vec3 u_light_col;
uniform vec2 u_cone_info;
uniform float u_max_distance;
uniform int u_light_type;
		//NO_LIGHT = 0,
		//POINT = 1,
		//SPOT = 2,
		//DIRECTIONAL = 3,
		//AMBIENT = 4

uniform int u_use_shadowmaps;
uniform int u_light_cast_shadows;
uniform sampler2D u_shadowmap;
uniform mat4 u_shadow_viewproj;
uniform float u_shadow_bias;
uniform int u_shadowmap_index;
uniform int u_shadowmap_dimensions;

out vec4 FragColor;

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx( p );
	vec3 dp2 = dFdy( p );
	vec2 duv1 = dFdx( uv );
	vec2 duv2 = dFdy( uv );
	
	// solve the linear system
	vec3 dp2perp = cross( dp2, N );
	vec3 dp1perp = cross( N, dp1 );
	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;
 
	// construct a scale-invariant frame 
	float invmax = inversesqrt( max( dot(T,T), dot(B,B) ) );
	return mat3( T * invmax, B * invmax, N );
}

// assume N, the interpolated vertex normal and 
// WP the world position
//vec3 normal_pixel = texture2D( normalmap, uv ).xyz; 
vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255./127. - 128./127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}

float computeShadowMP(vec3 wp){
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj * vec4(wp,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = (proj_pos.xy / proj_pos.w);

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);
	//it is outside on the sides
	if( shadow_uv.x < 0.0 || shadow_uv.x > 1.0 || shadow_uv.y < 0.0 || shadow_uv.y > 1.0 ) {
		return 1.0;
	}

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_bias) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;
	//it is before near or behind far plane
	if(real_depth < 0.0 || real_depth > 1.0) {
		return 1.0;
	}

	//read depth from depth buffer in [0..+1] non-linear
	//accounts for offset using shadowmap dimensions and id
	float shadow_depth = texture( u_shadowmap, vec2(shadow_uv.x*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*(u_shadowmap_index%u_shadowmap_dimensions),	 shadow_uv.y*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*floor(u_shadowmap_index/u_shadowmap_dimensions))).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if( shadow_depth < real_depth ) {
		shadow_factor = 0.0;
	}
	return shadow_factor;
}

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );
	color.xyz =  pow(color.xyz, vec3(2.2));

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 light = vec3(0.0, 0.0, 0.0);

	vec3 V = normalize(u_camera_position - v_world_position);
	vec3 N = normalize(v_normal);
	vec3 normal_pixel = texture( u_normalmap, v_uv ).xyz;
	if (u_use_normalmap == 1) {
		N = perturbNormal(N,v_world_position, v_uv , normal_pixel);
	}
	float spec_ks = texture( u_metal_roughness, v_uv).g;
	float spec_a =  texture( u_metal_roughness, v_uv).b;

	vec3 emissive_pixel = texture( u_emissive, v_uv ).xyz;

	vec3 light_color =  pow(u_light_col, vec3(2.2));

	if (u_light_type == 1) { 		//point lights
		//diffuse value
		vec3 L = u_light_pos - v_world_position;
		L= normalize(L);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);
		
		float lightDist = distance(u_light_pos, v_world_position);
		float att_factor = u_max_distance - lightDist;
		att_factor = att_factor/u_max_distance;
		att_factor = max(att_factor, 0.0);

		light += (NdotL * light_color) * att_factor;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor; 
			if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * light_color * att_factor; }
		}
	}
	else if (u_light_type == 2) { 		//spot lights
		//diffuse value
		vec3 L = u_light_pos - v_world_position;
		L= normalize(L);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);
		
		float lightDist = distance(u_light_pos, v_world_position);
		float att_factor = u_max_distance - lightDist;
		att_factor = att_factor/u_max_distance;
		att_factor = max(att_factor, 0.0);

		float cos_angle = dot(u_light_front, L);
		if (cos_angle < u_cone_info.y) {
			NdotL = 0.0;
		}
		else if (cos_angle < u_cone_info.x) {
			NdotL *= (cos_angle - u_cone_info.y) / (u_cone_info.x - u_cone_info.y);
		}

		//shadow value
		float shadow_factor = 1.0;
		if (u_light_cast_shadows == 1 && u_use_shadowmaps == 1) {
			shadow_factor = computeShadowMP(v_world_position);
		}

		light += (NdotL * light_color) * att_factor * shadow_factor;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor;

			if (final_a != 0) {light += NdotL * final_ks * pow(NdotH, final_a) * light_color * att_factor * shadow_factor; }
		}
	}
	else if (u_light_type == 3) {		//directional lights
		//diffuse value
		vec3 L = u_light_front;
		L = normalize(L);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);

		//shadow value
		float shadow_factor = 1.0;
		if (u_light_cast_shadows == 1 && u_use_shadowmaps == 1) {
			shadow_factor = computeShadowMP(v_world_position);
		}

		light += NdotL * light_color * shadow_factor;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor;
			if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * light_color * shadow_factor; }
		}
	}
	else if (u_light_type == 4) {		//ambient light (first pass)
		float occlusion = texture( u_metal_roughness, v_uv).r;
		vec3 ambient_light = pow(u_ambient_light, vec3(2.2));
		if (u_use_occlusion == 1) {
			light += ambient_light * occlusion;
		}
		else{
			light += ambient_light;
		}
	}
	FragColor.xyz = color.xyz * light;
	if (u_light_type == 4 && u_use_emissive == 1) {FragColor.xyz += pow(emissive_pixel * u_emissive_factor, vec3(2.2));}
	FragColor.a = color.a;
}

\gbuffers.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

uniform sampler2D u_metal_roughness;
uniform sampler2D u_normalmap;
uniform int u_use_normalmap; //avoids visual bugs when perturbing a normal when there is no normalmap
uniform sampler2D u_emissive;

uniform vec3 u_emissive_factor;
uniform float u_metal_factor;
uniform float u_rough_factor;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;
layout(location = 2) out vec4 MaterialProperties;

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx( p );
	vec3 dp2 = dFdy( p );
	vec2 duv1 = dFdx( uv );
	vec2 duv2 = dFdy( uv );
	
	// solve the linear system
	vec3 dp2perp = cross( dp2, N );
	vec3 dp1perp = cross( N, dp1 );
	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;
 
	// construct a scale-invariant frame 
	float invmax = inversesqrt( max( dot(T,T), dot(B,B) ) );
	return mat3( T * invmax, B * invmax, N );
}

// assume N, the interpolated vertex normal and 
// WP the world position
vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255./127. - 128./127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color * texture( u_texture, uv ); //TODO linearize both components before operating (no effective result as u_color is always (1,1,1,1) but just in case another scene is used)
	color.xyz = pow(color.xyz, vec3(2.2)); 
	vec3 colorTexture = color.xyz;

	float occlusion = texture(u_metal_roughness, uv).x;
	float shininess =  1.0 - (texture(u_metal_roughness, uv).b * u_rough_factor);
	float metalness = texture(u_metal_roughness, uv).g * u_metal_factor;

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 emissive_pixel = texture( u_emissive, uv ).xyz * u_emissive_factor;
	emissive_pixel = pow(emissive_pixel, vec3(2.2));

	vec3 N = normalize(v_normal);
	vec3 normal_pixel = texture( u_normalmap, uv ).xyz;

	if (u_use_normalmap == 1) {N = perturbNormal(N,v_world_position, uv , normal_pixel); }

	FragColor.xyz = colorTexture;
	FragColor.w = emissive_pixel.r;

	NormalColor.xyz = N * 0.5 + vec3(0.5);
	NormalColor.w = emissive_pixel.g;

	MaterialProperties.x = occlusion;
	MaterialProperties.y = metalness;
	MaterialProperties.z = shininess;
	MaterialProperties.w = emissive_pixel.b;
}

\deferred_global.fs

#version 330 core

#define RECIPROCAL_PI 0.3183098861837697
#define PI 3.1416

in vec3 v_position;
in vec2 v_uv;

uniform sampler2D u_color_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_mat_properties_texture;
uniform sampler2D u_depth_texture;

uniform vec2 u_invRes;
uniform mat4 u_inverse_viewprojection;

uniform vec3 u_ambient_light;

uniform vec3 u_camera_position;

uniform vec3 u_light_pos;
uniform vec3 u_light_front;
uniform vec3 u_light_col;
uniform vec2 u_cone_info;
uniform float u_max_distance;
uniform int u_light_type;
		//NO_LIGHT = 0,
		//POINT = 1,
		//SPOT = 2,
		//DIRECTIONAL = 3
		//AMBIENT = 4

uniform int u_light_cast_shadows;
uniform sampler2D u_shadowmap;
uniform mat4 u_shadow_viewproj;
uniform float u_shadow_bias;
uniform int u_shadowmap_index;
uniform int u_shadowmap_dimensions;

uniform sampler2D u_ssao_map;
uniform int u_occlusion_type;
		//TEXTURE = 0,
		//SSAO = 1,
		//SSAO+ = 2

out vec4 FragColor;

// Normal Distribution Function using GGX Distribution
float D_GGX (const in float NoH, const in float linearRoughness )
{
	float a2 = linearRoughness * linearRoughness;
	float f = (NoH * NoH) * (a2 - 1.0) + 1.0;
	return a2 / (PI * f * f);
}

// Fresnel term with scalar optimization(f90=1)
float F_Schlick( const in float VoH, const in float f0)
{
	float f = pow(1.0 - VoH, 5.0);
	return f0 + (1.0 - f0) * f;
}

// Fresnel term with colorized fresnel
vec3 F_Schlick( const in float VoH, const in vec3 f0)
{
	float f = pow(1.0 - VoH, 5.0);
	return f0 + (vec3(1.0) - f0) * f;
}

// Geometry Term: Geometry masking/shadowing due to microfacets
float GGX(float NdotV, float k){
	return NdotV / (NdotV * (1.0 - k) + k);
}
	
float G_Smith( float NdotV, float NdotL, float roughness)
{
	float k = pow(roughness + 1.0, 2.0) / 8.0;
	return GGX(NdotL, k) * GGX(NdotV, k);
}

// Diffuse Reflections: Disney BRDF using retro-reflections using F term, this is much more complex!!
// might be doing something wrong w/linear roughness but this looks like shit
float Fd_Burley ( const in float NoV, const in float NoL, const in float LoH, const in float linearRoughness)
{
        float f90 = 0.5 + 2.0 * linearRoughness * LoH * LoH;
        float lightScatter = F_Schlick(NoL, f90);
        float viewScatter  = F_Schlick(NoV, f90);
        return lightScatter * viewScatter * RECIPROCAL_PI;
}

vec3 specularBRDF( float roughness, vec3 f0, float NoH, float NoV, float NoL, float LoH )
{
	float a = roughness * roughness;

	// Normal Distribution Function
	float D = D_GGX( NoH, a );

	// Fresnel Function
	vec3 F = F_Schlick( LoH, f0 );

	// Visibility Function (shadowing/masking)
	float G = G_Smith( NoV, NoL, roughness );
		
	// Norm factor
	vec3 spec = D * G * F;
	spec /= (4.0 * NoL * NoV + 1e-6);

	return spec;
}

float computeShadowMP(vec3 wp){
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj * vec4(wp,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = (proj_pos.xy / proj_pos.w);

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);
	//it is outside on the sides
	if( shadow_uv.x < 0.0 || shadow_uv.x > 1.0 || shadow_uv.y < 0.0 || shadow_uv.y > 1.0 ) {
		return 1.0;
	}

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_bias) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;
	//it is before near or behind far plane
	if(real_depth < 0.0 || real_depth > 1.0) {
		return 1.0;
	}

	//read depth from depth buffer in [0..+1] non-linear
	//accounts for offset using shadowmap dimensions and id
	float shadow_depth = texture( u_shadowmap, vec2(shadow_uv.x*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*(u_shadowmap_index%u_shadowmap_dimensions),	 shadow_uv.y*(1.0/u_shadowmap_dimensions)+(1.0/u_shadowmap_dimensions)*floor(u_shadowmap_index/u_shadowmap_dimensions))).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if( shadow_depth < real_depth ) {
		shadow_factor = 0.0;
	}
	return shadow_factor;
}

void main()
{
	vec2 uv = gl_FragCoord.xy * u_invRes.xy;
	vec3 light = vec3(0.0, 0.0, 0.0);
	vec3 color = texture( u_color_texture, uv ).xyz;
	float depth = texture(u_depth_texture, uv).x;
	vec3 N = texture(u_normal_texture, uv).xyz * 2.0 - vec3(1.0);
	N = normalize(N);
	float metalness = texture( u_mat_properties_texture, uv).g;
	float shininess =  texture( u_mat_properties_texture, uv).b;

	if(depth == 1.0) {
		discard;
	}

	//reconstruct world position from depth and inv. viewproj
	vec4 screen_pos = vec4(uv.x*2.0-1.0, uv.y*2.0-1.0, depth*2.0-1.0, 1.0);
	vec4 proj_worldpos = u_inverse_viewprojection * screen_pos;
	vec3 worldpos = proj_worldpos.xyz / proj_worldpos.w;

	vec3 V = normalize(u_camera_position - worldpos);

	//linearize light (+ambient) color
	vec3 light_color =  pow(u_light_col, vec3(2.2));
	vec3 ambient_light = pow(u_ambient_light, vec3(2.2));

	if (u_light_type == 1 || u_light_type == 2) { 		//point lights and spot lights
		//shadow factor initialization
		float shadow_factor = 1.0;
		
		//get vectors
		vec3 L = u_light_pos - worldpos;
		L= normalize(L);
		vec3 H = normalize(L + V);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);

		//correction for spotlights
		float spotlight_attenuation = 1.0;
		if (u_light_type == 2) {
			float cos_angle = dot(u_light_front, L);
			if (cos_angle < u_cone_info.y) {
				spotlight_attenuation = 0.0;
			}
			else if (cos_angle < u_cone_info.x) {
				spotlight_attenuation = (cos_angle - u_cone_info.y) / (u_cone_info.x - u_cone_info.y);
			}
		//shadow factor
			if (u_light_cast_shadows == 1) {
				shadow_factor = computeShadowMP(worldpos);
			}
		}
		
		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix( vec3(0.5), color, metalness );

		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness) * color;

		//specular value
		vec3 Fr_d = specularBRDF((1.0-shininess), f0, NdotH, NdotV, NdotL, LdotH);
		
		//diffuse value
		vec3 Fd_d = diffuseColor * NdotL;
		
		//add diffuse and specular reflection
		vec3 direct = Fd_d + Fr_d;

		//attenuation factor
		float lightDist = distance(u_light_pos, worldpos);
		float att_factor = u_max_distance - lightDist;
		att_factor = att_factor/u_max_distance;
		att_factor = max(att_factor, 0.0);

		//final light computation
		light += direct * att_factor * shadow_factor * light_color * spotlight_attenuation;
	}
	else if (u_light_type == 3) {		//directional lights
		//get vectors
		vec3 L = u_light_front;
		L = normalize(L);
		vec3 H = normalize(L + V);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);
		float NdotH = dot(N, H);
		float NdotV = dot(N, V);
		float LdotH = dot(L, H);

		float shadow_factor = 1.0;
		if (u_light_cast_shadows == 1) {
			shadow_factor = computeShadowMP(worldpos);
		}

		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix( vec3(0.5), color, metalness );

		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness) * color;

		//specular value
		vec3 Fr_d = specularBRDF((1.0-shininess), f0, NdotH, NdotV, NdotL, LdotH);
		
		//diffuse value
		vec3 Fd_d = diffuseColor * NdotL;
		
		//add diffuse and specular reflection
		vec3 direct = Fr_d + Fd_d;

		//final light computation
		light += direct * shadow_factor * light_color;
	}
	else if (u_light_type == 4) {	//ambient light
		//occlusion from texture
		float occlusion;
		if (u_occlusion_type == 0) {
			occlusion = texture( u_mat_properties_texture, uv).r;
		}
		else { //occlusion from SSAO map (no blending because some prefabs have no occlusion texture and results are inconsistent)
			occlusion = texture( u_ssao_map, uv).r;
		}
		light += ambient_light * occlusion;
	}

	FragColor.xyz = color * light;
	
	//emissive light on first draw
	if (u_light_type == 4) {
		//TODO: optimize to avoid multiple reads on same texture (save previous as vec4)
		FragColor.r += texture(u_color_texture, uv).w;
		FragColor.g += texture(u_normal_texture, uv).w;
		FragColor.b += texture(u_mat_properties_texture, uv).w;
	}
	FragColor.a = 1.0;
	gl_FragDepth = depth;
}


\emissive.fs

#version 330 core

in vec3 v_position;
in vec2 v_uv;

uniform sampler2D u_color_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_mat_properties_texture;
uniform sampler2D u_depth_texture;

out vec4 FragColor;

void main() {
	vec2 uv = v_uv;

	float depth = texture(u_depth_texture, uv).x;
	if(depth == 1.0)
		discard;

	FragColor.r = texture(u_color_texture, uv).w;
	FragColor.g = texture(u_normal_texture, uv).w;
	FragColor.b = texture(u_mat_properties_texture, uv).w;
	FragColor.a = 1.0;
}

\ssao.fs

#version 330 core

const int RANDOM_POINTS = 128;

in vec3 v_position;
in vec2 v_uv;

uniform sampler2D u_depth_texture;
uniform sampler2D u_normal_texture;

uniform mat4 u_inverse_viewprojection;
uniform mat4 u_viewprojection;
uniform vec2 u_invRes; 

uniform float u_radius;
uniform float u_dist_threshold; //TODO range check

uniform vec3 u_random_pos[RANDOM_POINTS];

uniform int u_use_plus;

out vec4 FragColor;

void main() {
	vec2 uv = gl_FragCoord.xy * u_invRes.xy;
	float depth = texture(u_depth_texture, uv).x;
	vec3 N = texture(u_normal_texture, uv).xyz * 2.0 - vec3(1.0);
	N = normalize(N);

	if(depth == 1.0)
		discard;

	//reconstruct world position from depth and inv. viewproj
	vec4 screen_pos = vec4(uv.x*2.0-1.0, uv.y*2.0-1.0, depth*2.0-1.0, 1.0);
	vec4 proj_worldpos = u_inverse_viewprojection * screen_pos;
	vec3 worldpos = proj_worldpos.xyz / proj_worldpos.w;

	int num = RANDOM_POINTS;
	for (int i=0; i<RANDOM_POINTS; i++) {
		//get point's world position
		vec3 p = worldpos + u_random_pos[i] * u_radius;
		if (u_use_plus == 1) {
			//project points to the correct hemisphere if they're not
				vec3 point_vector = u_random_pos[i] * u_radius;
				if( dot(N, point_vector) < 0.0) {
					p -= 2.0 * point_vector;
				}
		}

		//convert to texture space
		vec4 p_proj = u_viewprojection * vec4(p, 1.0);
		p_proj.xy /= p_proj.w; //convert to clipspace from homogeneous
		
		//apply a tiny bias to its z before converting to clip-space
		p_proj.z = (p_proj.z - 0.005) / p_proj.w;
		p_proj.xyz = p_proj.xyz * 0.5 + vec3(0.5); //to [0..1]
		
		//read p true depth
		float pdepth = texture( u_depth_texture, p_proj.xy ).x;
		
		//compare true depth with its depth
		if( pdepth < p_proj.z) //if true depth smaller, is inside
			num--; //remove this point from the list of visible
	}

	float occlusion_factor = float(num)/RANDOM_POINTS;

	//shade more if normal is pointing down
	vec3 down_vector = vec3(0.0,-1.0,0.0);
	float extra_shading_factor = dot(N, down_vector);
	if (extra_shading_factor > 0.0) {
		occlusion_factor *= (1.0 - extra_shading_factor * 0.95);
	}

	FragColor.xyz = vec3(occlusion_factor);
	FragColor.a = 1.0;
}

\blur_reprojection.fs

//NOT WORKING CORRECTLY

#version 330 core

in vec3 v_position;
in vec2 v_uv;

uniform sampler2D u_current_frame;
uniform sampler2D u_last_results;
uniform sampler2D u_depth_texture;

uniform mat4 u_inverse_viewprojection;
uniform mat4 u_previous_viewprojection;
uniform vec2 u_invRes;

out vec4 FragColor;

void main() {
	vec2 uv = gl_FragCoord.xy * u_invRes.xy;
	float current_depth = texture(u_depth_texture, uv).x;

	float current_occlusion = texture(u_current_frame, uv).x;
	float previous_occlusion = 1.0;

	//reconstruct world position from depth and inv. viewproj
	vec4 screen_pos = vec4(uv.x*2.0-1.0, uv.y*2.0-1.0, current_depth*2.0-1.0, 1.0);
	vec4 proj_worldpos = u_inverse_viewprojection * screen_pos;
	vec3 worldpos = proj_worldpos.xyz / proj_worldpos.w;

	//project world position to the previous results
	vec4 proj_pos = u_previous_viewprojection * vec4(worldpos,1.0);

	//from homogeneus space to clip space
	vec2 prev_uv = (proj_pos.xy / proj_pos.w);

	//from clip space to uv space
	prev_uv = prev_uv * 0.5 + vec2(0.5);

	//if pixel was not seen in previous frame consider current results only
	if( prev_uv.x < 0.0 || prev_uv.x > 1.0 || prev_uv.y < 0.0 || prev_uv.y > 1.0 || current_depth == 1.0) {
		previous_occlusion = current_occlusion;
	}
	//get value from previous results
	else {
		previous_occlusion = texture(u_last_results, prev_uv).x;
	}

	//get final occlusion value
	float occlusion_factor = 0.1 * current_occlusion + 0.9 * previous_occlusion;

	//for some reason the depth is inconsistent and results incorrect :(
	float prev_depth = texture(u_depth_texture, prev_uv).x;
	if (prev_depth == 1.0) {occlusion_factor = 1.0;}

	FragColor.xyz = vec3(occlusion_factor);
	FragColor.a = 1.0;
}

\blur_neighbors.fs

#version 330 core

in vec3 v_position;
in vec2 v_uv;

uniform sampler2D u_ssao_map;
uniform sampler2D u_depth_texture;

uniform vec2 u_invRes;

out vec4 FragColor;

void main() {
	vec2 uv = gl_FragCoord.xy * u_invRes.xy;

	float depth = texture(u_depth_texture, uv).x;
	if (depth==1.0) {
		discard;
	}

	float occlusion_factor = 0.0;
	float valid_pixels = pow(3.0, 2.0);

	float current_depth;
	//blur by averaging the pixels on a 3*3 square centered on the current pixel.
	vec2 current_uv;
	vec2 offset = vec2(-1.0, -1.0) * u_invRes;
	for (int i=0; i < 3; i++) {
		for (int j=0; j < 3; j++) {
				current_uv = vec2(float(i),float(j)) * u_invRes;
				current_uv += uv;
				current_uv += offset;
				current_depth = texture(u_depth_texture, uv).x; //to avoid counting skybox pixels as valid occlusion data
				if (current_uv.x < 0.0 || current_uv.x > 1.0 || current_uv.y < 0.0 || current_uv.y > 1.0 || current_depth == 1.0) {
					valid_pixels -= 1.0;
				}
				else {
					occlusion_factor += texture(u_ssao_map, current_uv).x;
				}
		}
	}

	FragColor.xyz = vec3(occlusion_factor/ valid_pixels);
	FragColor.a = 1.0;
}

\skybox.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;

uniform samplerCube u_texture;
uniform vec3 u_camera_position;
out vec4 FragColor;

void main()
{
	vec3 E = v_world_position - u_camera_position;
	vec4 color = texture( u_texture, E );
	color.xyz = pow(color.xyz, vec3(2.2));
	FragColor = color;


}


\multi.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	FragColor = color;
	NormalColor = vec4(N,1.0);
}


\depth.fs

#version 330 core

uniform vec2 u_camera_nearfar;
uniform sampler2D u_texture; //depth map
in vec2 v_uv;
out vec4 FragColor;

void main()
{
	float n = u_camera_nearfar.x;
	float f = u_camera_nearfar.y;
	float z = texture2D(u_texture,v_uv).x;
	if( n == 0.0 && f == 1.0 )
		FragColor = vec4(z);
	else
		FragColor = vec4( n * (z + 1.0) / (f + n - z * (f - n)) );
}


\instanced.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;

in mat4 u_model;

uniform vec3 u_camera_pos;

uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( a_vertex, 1.0) ).xyz;
	
	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\gamma.fs
#version 330 core

in vec2 v_uv;

uniform sampler2D u_texture;
uniform float u_inv_gamma;


out vec4 FragColor;

void main() {
	vec3 color = texture2D( u_texture, v_uv ).xyz;
	color = pow( color, vec3( u_inv_gamma ) );
	FragColor = vec4( color, 1.0 );
}

\tonemapper.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_texture;
uniform float u_scale; //color scale before tonemapper
uniform float u_average_lum; 
uniform float u_lumwhite2;
uniform float u_inv_gamma; //inverse gamma


void main() {
	vec4 color = texture2D( u_texture, v_uv );
	vec3 rgb = color.xyz;

	float lum = dot(rgb, vec3(0.2126, 0.7152, 0.0722));
	float L = (u_scale / u_average_lum) * lum;
	float Ld = (L * (1.0 + L / u_lumwhite2)) / (1.0 + L);

	rgb = (rgb / lum) * Ld;
	rgb = max(rgb,vec3(0.001));

	//linear to gamma
	rgb = pow( rgb, vec3( u_inv_gamma ) );
	gl_FragColor = vec4( rgb, color.a );
}

