//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
lightSP basic.vs lightSP.fs
lightMP basic.vs lightMP.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs

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

	if(color.a < u_alpha_cutoff)
		discard;

	FragColor = color;
}

\lightSP.fs

#version 330 core

const int MAX_LIGHTS = 10;
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

float computeShadow(vec3 wp, int i){
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj[i] * vec4(wp,1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = (proj_pos.xy / proj_pos.w);

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_bias[i]) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;

	//read depth from depth buffer in [0..+1] non-linear
	//float shadow_depth = texture( u_shadowmap, shadow_uv).x;
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

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 light = vec3(0.0, 0.0, 0.0);
	
	float occlusion = texture( u_metal_roughness, v_uv).r;
	if (u_use_occlusion == 1) {
		light += u_ambient_light * occlusion;
	}
	else{
		light += u_ambient_light;
	}

	vec3 V = normalize(u_camera_position - v_world_position);

	vec3 N = normalize(v_normal);
	vec3 normal_pixel = texture( u_normalmap, v_uv ).xyz;
	if (u_use_normalmap == 1) {
		N = perturbNormal(N,v_world_position, v_uv , normal_pixel);
	}

	vec3 emissive_pixel = texture( u_emissive, v_uv ).xyz;
	if (u_use_emissive == 1) {
		light += emissive_pixel * u_emissive_factor;
	}

	float spec_ks = texture( u_metal_roughness, v_uv).g;
	float spec_a =  texture( u_metal_roughness, v_uv).b;

	for (int i=0; i<MAX_LIGHTS; i++) {
		if (i<u_num_lights) {
			if (u_light_type[i] == 1) { 		//point lights
				//diffuse value
				vec3 L = u_light_pos[i] - v_world_position;
				L= normalize(L);
				float NdotL = clamp(dot(N, L), 0.0, 1.0);
				
				float lightDist = distance(u_light_pos[i], v_world_position);
				float att_factor = u_max_distance[i] - lightDist;
				att_factor = att_factor/u_max_distance[i];
				att_factor = max(att_factor, 0.0);

				light += (NdotL * u_light_col[i]) * att_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;
					if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col[i] * att_factor; }
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
					shadow_factor = computeShadow(v_world_position, i);
				}

				light += (NdotL * u_light_col[i]) * att_factor * shadow_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;

					float cos_angle = dot(u_light_front[i], L);
					if (cos_angle < u_cone_info[i].y) {
						NdotH = 0;
					}
					else if (cos_angle < u_cone_info[i].x) {
						NdotH *= (cos_angle - u_cone_info[i].y) / (u_cone_info[i].x - u_cone_info[i].y);
					}

					if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col[i] * att_factor; }
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
					shadow_factor = computeShadow(v_world_position, i);
				}
				light += NdotL * u_light_col[i] * shadow_factor;

				//specular value (blinn-phong)
				if (u_use_specular == 1 && NdotL > 0.0) {
					vec3 H = normalize(L + V);
					float NdotH = clamp(dot(N, H), 0.0, 1.0);
					float final_a = 1-(spec_a * u_rough_factor);
					float final_ks = spec_ks * u_metal_factor;
					if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col[i]; }
				}
			}
		}
	}
	FragColor.xyz = color.xyz * light;
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
		//DIRECTIONAL = 3

uniform int u_num_lights;

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

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );

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

	if (u_light_type == 1) { 		//point lights
		//diffuse value
		vec3 L = u_light_pos - v_world_position;
		L= normalize(L);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);
		
		float lightDist = distance(u_light_pos, v_world_position);
		float att_factor = u_max_distance - lightDist;
		att_factor = att_factor/u_max_distance;
		att_factor = max(att_factor, 0.0);

		light += (NdotL * u_light_col) * att_factor;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor;
			if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col * att_factor; }
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

		light += (NdotL * u_light_col) * att_factor;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor;

			float cos_angle = dot(u_light_front, L);
			if (cos_angle < u_cone_info.y) {
				NdotH = 0;
			}
			else if (cos_angle < u_cone_info.x) {
				NdotH *= (cos_angle - u_cone_info.y) / (u_cone_info.x - u_cone_info.y);
			}

			if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col * att_factor; }
		}
	}
	else if (u_light_type == 3) {		//directional lights
		//diffuse value
		vec3 L = u_light_front;
		L = normalize(L);
		float NdotL = clamp(dot(N, L), 0.0, 1.0);

		light += NdotL * u_light_col;

		//specular value (blinn-phong)
		if (u_use_specular == 1 && NdotL > 0.0) {
			vec3 H = normalize(L + V);
			float NdotH = clamp(dot(N, H), 0.0, 1.0);
			float final_a = 1-(spec_a * u_rough_factor);
			float final_ks = spec_ks * u_metal_factor;
			if (final_a != 0) {light += final_ks * pow(NdotH, final_a) * u_light_col; }
		}
	}
	else if (u_light_type == 4) {		//ambient light (first pass)
		float occlusion = texture( u_metal_roughness, v_uv).r;
		if (u_use_occlusion == 1) {
			light += u_ambient_light * occlusion;
		}
		else{
			light += u_ambient_light;
		}
		
		vec3 emissive_pixel = texture( u_emissive, v_uv ).xyz;
		if (u_use_emissive == 1) {
			light += emissive_pixel * u_emissive_factor;
		}
	}
	FragColor.xyz = color.xyz * light;
	FragColor.a = color.a;
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