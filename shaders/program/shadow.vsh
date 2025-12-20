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
#ifdef COLORED_LIGHTS_ENTITIES
uniform int blockEntityId;
uniform int entityId;
#ifdef IS_IRIS
uniform int currentRenderedItemId;
#endif
#endif
#endif

// ------------
//   Includes
// ------------

#include "/include/lighting/shadows/distortion.glsl"
#include "/include/vertex/displacement.glsl"

#ifdef COLORED_LIGHTS
#include "/include/lighting/lpv/voxelization.glsl"
#endif

void main() {
	uv            = gl_MultiTexCoord0.xy;
	material_mask = uint(mc_Entity.x - 10000.0);
	tint          = gl_Color.rgb;

#ifdef COLORED_LIGHTS

	#ifdef COLORED_LIGHTS_ENTITIES

	if ((renderStage == MC_RENDER_STAGE_ENTITIES || renderStage == MC_RENDER_STAGE_BLOCK_ENTITIES)) {

		#ifdef IS_IRIS
		if (currentRenderedItemId > 0 && entityId != 10103) material_mask = uint(currentRenderedItemId - 10000);
		else
		#endif

		if (entityId > 0) material_mask = uint(entityId - 10000);
		else if (blockEntityId > 0) material_mask = uint(blockEntityId - 10000);
		else material_mask = 181u; // Transparent

		if (100u <= material_mask && material_mask < 164u ) {
			if (material_mask == 102u) material_mask = 80u; // Lightning
			else if (material_mask == 103u) material_mask = 181u; // Player
			else if (material_mask == 104u) material_mask = 48u; // Drowned
		#ifdef WORLD_END
			else if (material_mask == 105u) material_mask = 62u; // Dragon breath
		#endif
		}
		//else if (264u <= material_mask && material_mask < 280u) material_mask += 32u; // Candle Items
	} else if (renderStage == MC_RENDER_STAGE_PARTICLES) {
		// Make enderman/nether portal particles glow
		if (gl_Color.r > gl_Color.g && gl_Color.g < 0.6 && gl_Color.b > 0.4) material_mask = 47u;
		else material_mask = 27u;
	}
	#endif

	if (renderStage != MC_RENDER_STAGE_HAND_SOLID && renderStage != MC_RENDER_STAGE_HAND_TRANSLUCENT && material_mask != 181u) {
		update_voxel_map(material_mask);
	}
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

