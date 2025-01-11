/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/shadow:
  Render shadow map

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

out vec2 uv;

flat out uint material_mask;
flat out vec3 tint;

#ifdef WATER_CAUSTICS
out vec3 scene_pos;

#if defined (PHYSICS_MOD_OCEAN) && defined (PHYSICS_OCEAN)
	#include "/include/misc/oceans.glsl"
#endif
#endif

// --------------
//   Attributes
// --------------

attribute vec3 at_midBlock;
attribute vec4 at_tangent;
attribute vec3 mc_Entity;
attribute vec2 mc_midTexCoord;

// ------------
//   Uniforms
// ------------

uniform sampler2D tex;
uniform sampler2D noisetex;

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferProjectionInverse;

uniform mat4 shadowModelView;
uniform mat4 shadowModelViewInverse;

uniform vec3 cameraPosition;

uniform float near;
uniform float far;

uniform float frameTimeCounter;
uniform float rainStrength;

uniform vec2 taa_offset;
uniform vec3 light_dir;

#ifdef COLORED_LIGHTS
writeonly uniform uimage3D voxel_img;

uniform int renderStage;
#endif

// ------------
//   Includes
// ------------

#include "/include/lighting/distortion.glsl"
#include "/include/vertex/displacement.glsl"

#ifdef COLORED_LIGHTS
#include "/include/lighting/lpv/voxelization.glsl"
#endif

void main() {
	uv            = gl_MultiTexCoord0.xy;
	material_mask = uint(mc_Entity.x - 10000.0);
	tint          = gl_Color.rgb;

#ifdef COLORED_LIGHTS
	update_voxel_map(material_mask);
#endif

	bool is_top_vertex = uv.y < mc_midTexCoord.y;

	vec4 vert = gl_Vertex;

#if defined (WATER_CAUSTICS) && defined (PHYSICS_MOD_OCEAN) && defined (PHYSICS_OCEAN)
	if(physics_iterationsNormal >= 1.0 && material_mask == 1) {
		// basic texture to determine how shallow/far away from the shore the water is
		physics_localWaviness = texelFetch(physics_waviness, ivec2(gl_Vertex.xz) - physics_textureOffset, 0).r;
		// transform gl_Vertex (since it is the raw mesh, i.e. not transformed yet)
		vec4 finalPosition = vec4(gl_Vertex.x, gl_Vertex.y + physics_waveHeight(gl_Vertex.xz, PHYSICS_ITERATIONS_OFFSET, physics_localWaviness, physics_gameTime), gl_Vertex.z, gl_Vertex.w);
		// pass this to the fragment shader to fetch the texture there for per fragment normals
		physics_localPosition = finalPosition.xyz;

		// now use finalPosition instead of gl_Vertex
		vert.xyz = finalPosition.xyz;
	}
#endif

	vec3 pos = transform(gl_ModelViewMatrix, vert.xyz);
	     pos = transform(shadowModelViewInverse, pos);
	     pos = pos + cameraPosition;
	     pos = animate_vertex(pos, is_top_vertex, clamp01(rcp(240.0) * gl_MultiTexCoord1.y), material_mask);
		 pos = pos - cameraPosition;
		 pos = world_curvature(pos);

#ifdef WATER_CAUSTICS
	scene_pos = pos;
#endif

	vec3 shadow_view_pos = transform(shadowModelView, pos);
	vec3 shadow_clip_pos = project_ortho(gl_ProjectionMatrix, shadow_view_pos);
	     shadow_clip_pos = distort_shadow_space(shadow_clip_pos);

	gl_Position = vec4(shadow_clip_pos, 1.0);
}

