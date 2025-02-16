/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/c1_blend_layers
  Apply volumetric fog

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

layout (location = 0) out vec3 fragment_color;

/* RENDERTARGETS: 0 */

#ifdef BLOOMY_FOG
layout (location = 1) out float bloomy_fog;

/* RENDERTARGETS: 0,3 */
#endif

in vec2 uv;

flat in vec3 ambient_color;
flat in vec3 light_color;

#ifdef WORLD_OVERWORLD
#include "/include/fog/overworld/coeff_struct.glsl"
flat in AirFogCoefficients air_fog_coeff;
#endif

// ------------
//   Uniforms
// ------------

uniform sampler2D noisetex;

uniform sampler2D colortex0;  // scene color
uniform sampler2D colortex1;  // gbuffer 0
uniform sampler2D colortex2;  // gbuffer 1
uniform sampler2D colortex3;  // refraction data
uniform sampler2D colortex4;  // sky map
uniform sampler2D colortex5;  // scene history
uniform sampler2D colortex6;  // volumetric fog scattering
uniform sampler2D colortex7;  // volumetric fog transmittance
uniform sampler2D colortex11; // clouds history
uniform sampler2D colortex12; // clouds data
uniform sampler2D colortex13; // rendered translucent layer

#ifdef SHADOW

#ifdef AIR_FOG_COLORED_LIGHT_SHAFTS
uniform sampler2D shadowcolor0;
uniform sampler2D shadowtex0;
#endif
uniform sampler2D shadowtex1;
#endif

uniform sampler2D depthtex0;
uniform sampler2D depthtex1;

uniform mat4 gbufferModelView;
uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferProjection;
uniform mat4 gbufferProjectionInverse;

uniform mat4 gbufferPreviousModelView;
uniform mat4 gbufferPreviousProjection;

uniform mat4 shadowModelView;
uniform mat4 shadowProjection;

uniform vec3 cameraPosition;
uniform vec3 previousCameraPosition;

uniform float near;
uniform float far;

uniform float frameTimeCounter;
uniform float sunAngle;
uniform float rainStrength;
uniform float wetness;

uniform int worldTime;
uniform int moonPhase;
uniform int frameCounter;

uniform int isEyeInWater;
uniform float eyeAltitude;
uniform float blindness;
uniform float nightVision;
uniform float darknessFactor;

uniform vec3 light_dir;
uniform vec3 sun_dir;
uniform vec3 moon_dir;

uniform vec2 view_res;
uniform vec2 view_pixel_size;
uniform vec2 taa_offset;

uniform float eye_skylight;

uniform float biome_cave;
uniform float biome_may_rain;
uniform float biome_may_snow;

uniform float time_sunrise;
uniform float time_noon;
uniform float time_sunset;
uniform float time_midnight;

// ------------
//   Includes
// ------------

#define TEMPORAL_REPROJECTION

#include "/include/fog/simple_fog.glsl"
#include "/include/misc/distant_horizons.glsl"
#include "/include/utility/color.glsl"
#include "/include/utility/encoding.glsl"
#include "/include/utility/fast_math.glsl"
#include "/include/utility/space_conversion.glsl"

#ifdef WORLD_OVERWORLD
#include "/include/fog/overworld/analytic.glsl"
#endif

#ifdef DISTANT_HORIZONS
#include "/include/misc/distant_water.glsl"
#endif

vec3 blend_layers_with_fog(
	vec3 background_color,
	vec4 translucent_color,
	vec3 front_position_world,
	vec3 back_position_world,
	bool is_translucent,
	bool is_sky
) {
	// Apply analytic fog behind translucents

#if defined WORLD_OVERWORLD
	if (is_translucent) {
		mat2x3 analytic_fog = air_fog_analytic(
			front_position_world,
			back_position_world,
			is_sky,
			eye_skylight
		);

		background_color = background_color * analytic_fog[1] + analytic_fog[0];
	}
#endif

	return background_color * (1.0 - translucent_color.a) + translucent_color.rgb;
}

#if defined (PHYSICS_MOD_OCEAN) && defined (PHYSICS_OCEAN)
#include "/include/misc/oceans.glsl"
#endif

// https://iquilezles.org/www/articles/texture/texture.htm
vec4 smooth_filter(sampler2D sampler, vec2 coord) {
	vec2 res = vec2(textureSize(sampler, 0));

	coord = coord * res + 0.5;

	vec2 i, f = modf(coord, i);
	f = f * f * f * (f * (f * 6.0 - 15.0) + 10.0);
	coord = i + f;

	coord = (coord - 0.5) / res;
	return texture(sampler, coord);
}

vec4 read_clouds(out float apparent_distance) {
#if defined WORLD_OVERWORLD
	// Soften clouds for new pixels
	float pixel_age = texelFetch(colortex12, ivec2(gl_FragCoord.xy), 0).y;
	int ld = int(3.0 * dampen(max0(1.0 - 0.1 * pixel_age)));

	apparent_distance = min_of(textureGather(colortex12, uv * taau_render_scale, 0));

	return bicubic_filter_lod(colortex11, uv * taau_render_scale, ld);
#else
	return vec4(0.0, 0.0, 0.0, 1.0);
#endif
}

void main() {
	ivec2 texel = ivec2(gl_FragCoord.xy);

	// Sample textures

	float front_depth      = texelFetch(depthtex0, texel, 0).x;
	float back_depth       = texelFetch(depthtex1, texel, 0).x;
	
	vec4 refraction_data   = texelFetch(colortex3, texel, 0);
	vec4 translucent_color = texelFetch(colortex13, texel, 0);

#ifdef VL
	vec3 fog_transmittance = smooth_filter(colortex6, uv).rgb;
	vec3 fog_scattering    = smooth_filter(colortex7, uv).rgb;
#endif

	// Distant Horizons support

#ifdef DISTANT_HORIZONS
    float front_depth_dh   = texelFetch(dhDepthTex, texel, 0).x;
    float back_depth_dh    = texelFetch(dhDepthTex1, texel, 0).x;

    bool front_is_dh_terrain = is_distant_horizons_terrain(front_depth, front_depth_dh);
    bool back_is_dh_terrain = is_distant_horizons_terrain(back_depth, back_depth_dh);
#else
	#define front_depth_dh      front_depth
	#define back_depth_dh       back_depth
	#define front_is_dh_terrain false
	#define back_is_dh_terrain  false
#endif

	bool is_translucent = front_depth != back_depth;
	bool is_sky = back_depth == 1.0 && back_depth_dh == 1.0;

	// Space conversions

	vec3 front_position_screen = vec3(uv, front_is_dh_terrain ? front_depth_dh : front_depth);
	vec3 front_position_view   = screen_to_view_space(front_position_screen, true, front_is_dh_terrain);
	vec3 front_position_scene  = view_to_scene_space(front_position_view);
	vec3 front_position_world  = front_position_scene + cameraPosition;

	vec3 back_position_screen  = vec3(uv, back_is_dh_terrain ? back_depth_dh : back_depth);
	vec3 back_position_view    = screen_to_view_space(back_position_screen, true, back_is_dh_terrain);
	vec3 back_position_world   = view_to_scene_space(back_position_view) + cameraPosition;

	vec3 direction_world; float view_distance;
	length_normalize(front_position_scene - gbufferModelViewInverse[3].xyz, direction_world, view_distance);

	// Refraction

	vec2 refracted_uv = uv;
	float layer_dist = abs(view_distance - length(back_position_view));

#if REFRACTION != REFRACTION_OFF
	if (is_translucent && refraction_data != vec4(0.0)) {
		vec2 normal_tangent = vec2(
			unsplit_2x8(refraction_data.xy) * 2.0 - 1.0,
			unsplit_2x8(refraction_data.zw) * 2.0 - 1.0
		);

		refracted_uv = uv + normal_tangent.xy * rcp(max(view_distance, 1.0)) * min(layer_dist, 8.0) * (0.1 * REFRACTION_INTENSITY);

		// Make sure the refracted fragment is behind the fragment position
		float depth_refracted = texture(depthtex1, refracted_uv).x;
		refracted_uv = mix(refracted_uv, uv, float(depth_refracted < front_depth));
	}
#endif

	fragment_color = texture(colortex0, refracted_uv * taau_render_scale).rgb;

	// Draw DH water

#ifdef DISTANT_HORIZONS
	if (front_depth_dh != back_depth_dh) {
		// if there is a layer of DH water behind the translucent layer, these 
		// will store the position of that layer
		vec3 dh_position_screen = front_position_screen;
		vec3 dh_position_view   = front_position_view;
		vec3 dh_position_world  = front_position_world;

		// detect whether translucent DH terrain may be behind the translucent layer
		float z_mc = screen_to_view_space_depth(gbufferProjectionInverse, front_depth);
		float z_dh = screen_to_view_space_depth(dhProjectionInverse, front_depth_dh);

		const float error_margin = 1.0;
		bool dh_behind_translucent = z_dh > z_mc + error_margin && back_depth == 1.0;

		if (front_is_dh_terrain || dh_behind_translucent) {
			if (dh_behind_translucent) {
				dh_position_screen = vec3(uv, front_depth_dh);
				dh_position_view = screen_to_view_space(dhProjectionInverse, dh_position_screen, true);
				dh_position_world = view_to_scene_space(dh_position_view) + cameraPosition;
			}

			// Unpack gbuffer data

			vec4 gbuffer_data = texelFetch(colortex1, texel, 0);

			mat4x2 data = mat4x2(
				unpack_unorm_2x8(gbuffer_data.x),
				unpack_unorm_2x8(gbuffer_data.y),
				unpack_unorm_2x8(gbuffer_data.z),
				unpack_unorm_2x8(gbuffer_data.w)
			);

			vec3 tint          = vec3(data[0], data[1].x);
			uint material_mask = uint(255.0 * data[1].y);
			vec3 flat_normal   = decode_unit_vector(data[2]);
			vec2 light_levels  = data[3];
			if (material_mask == 1) { // Water
				vec4 water_color = draw_distant_water(
				dh_position_screen,
				dh_position_view,
				dh_position_world,
				direction_world,
				flat_normal,
				tint,
				light_levels,
				length_knowing_direction(cameraPosition - dh_position_world, direction_world),
				length_knowing_direction(dh_position_world - back_position_world, direction_world)
				);

				fragment_color = fragment_color * (1.0 - water_color.a) + water_color.rgb;
			}

			back_position_world = dh_behind_translucent
			? dh_position_world
			: back_position_world;
		}
	}
#endif

	// Blend layers

	fragment_color = blend_layers_with_fog(
		fragment_color,
		translucent_color,
		front_position_world,
		back_position_world,
		is_translucent,
		is_sky
	);

	// Blend clouds in front of translucents
	
	if (is_translucent) {
		float clouds_dist;
		vec4 clouds = read_clouds(clouds_dist);

		if (clouds_dist < view_distance) {
			fragment_color = fragment_color * clouds.w + clouds.xyz;
		}
	}

	// Blend fog

#if (defined WORLD_OVERWORLD || defined WORLD_END) && defined VL
	// Volumetric fog

	fragment_color = fragment_color * fog_transmittance + fog_scattering;
	
	#ifdef BLOOMY_FOG
	bloomy_fog = clamp01(dot(fog_transmittance, vec3(luminance_weights_rec2020)));
	bloomy_fog = isEyeInWater == 1.0 ? sqrt(bloomy_fog) : bloomy_fog;
	#endif
#else
	// Analytic fog

	if (isEyeInWater == 1) {
		// water fog
		float LoV = dot(direction_world, light_dir);
		mat2x3 analytic_fog = water_fog_simple(
			light_color,
			ambient_color,
			water_absorption_coeff,
			vec2(0.0, eye_skylight),
			view_distance,
			LoV,
			15.0 * eye_skylight
		);

		fragment_color *= analytic_fog[1];
		fragment_color += analytic_fog[0];

	#ifdef BLOOMY_FOG
		bloomy_fog = sqrt(clamp01(dot(analytic_fog[1], vec3(0.33))));
	#endif
	} else {
		// air fog

	#if defined WORLD_OVERWORLD
		mat2x3 analytic_fog = air_fog_analytic(
			cameraPosition,
			front_position_world,
			is_sky,
			eye_skylight
		);

		fragment_color *= analytic_fog[1];
		fragment_color += analytic_fog[0];

		#ifdef BLOOMY_FOG
		bloomy_fog = clamp01(dot(fog_transmittance, vec3(luminance_weights_rec2020)));
		bloomy_fog = isEyeInWater == 1.0 ? sqrt(bloomy_fog) : bloomy_fog;
		#endif
	#else 
		#ifdef BLOOMY_FOG
		bloomy_fog = 1.0;
		#endif
	#endif
	}
#endif

#ifdef BLOOMY_FOG
	#if   defined WORLD_NETHER
	bloomy_fog = spherical_fog(view_distance, nether_fog_start, nether_bloomy_fog_density) * 0.33 + 0.67;
	#elif defined WORLD_END
	bloomy_fog = bloomy_fog * 0.5 + 0.5;
	#endif
#endif
}