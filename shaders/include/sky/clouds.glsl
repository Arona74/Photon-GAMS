#if !defined INCLUDE_SKY_CLOUDS
#define INCLUDE_SKY_CLOUDS

#include "clouds/altocumulus.glsl"
#include "clouds/cumulus.glsl"
#include "clouds/cumulus_congestus.glsl"
#include "clouds/cumulonimbus.glsl"
#include "clouds/cirrus.glsl"
#include "clouds/noctilucent.glsl"

#if !defined CLOUDS_CUMULUS && !defined CLOUDS_CUMULUS_CONGESTUS && !defined CLOUDS_CUMULONIMBUS && !defined CLOUDS_ALTOCUMULUS && !defined CLOUDS_CIRRUS
CloudsResult draw_clouds(
	vec3 air_viewer_pos,
	vec3 ray_dir,
	vec3 clear_sky,
	float distance_to_terrain,
	float dither
) { return clouds_not_hit; }
#else

// cursed sort optimization
const uint defined_clouds =
#ifdef CLOUDS_CUMULUS
	1u |
#endif
#ifdef CLOUDS_CUMULUS_CONGESTUS
	2u |
#endif
#ifdef CLOUDS_CUMULONIMBUS
	4u |
#endif
#ifdef CLOUDS_ALTOCUMULUS
	8u |
#endif
#ifdef CLOUDS_CIRRUS
	16u
#else
	0u
#endif
;

const uint cloud_layers = max(
    ((defined_clouds      ) & 1u) +
    ((defined_clouds >> 1u) & 1u) +
    ((defined_clouds >> 2u) & 1u) +
    ((defined_clouds >> 3u) & 1u) +
    ((defined_clouds >> 4u) & 1u), 1u
);

// Insertion sort from https://github.com/OpenGLInsights/OpenGLInsightsCode/blob/master/Chapter%2020%20Efficient%20Layered%20Fragment%20Buffer%20Techniques/sorting.glsl
vec2[cloud_layers] sort_clouds(vec2[cloud_layers] cloud_types) {
	if (cloud_layers < 2) return cloud_types;

	else if (cloud_layers == 2) {
		if (cloud_types[0].y > cloud_types[cloud_layers - 1].y) {
			vec2 temp = cloud_types[0];
			cloud_types[0] = cloud_types[cloud_layers - 1];
			cloud_types[cloud_layers - 1] = temp;
		}
		return cloud_types;
	}

	for (int j = 1; j < cloud_layers; ++j) {
		vec2 key = cloud_types[j];
		int i = j - 1;
		while (i >= 0 && cloud_types[i].y > key.y) {
			cloud_types[i + 1] = cloud_types[i];
			--i;
		}
		cloud_types[i + 1] = key;
	}
	return cloud_types;
}

float abs_min(float a, float b) {
	if(abs(b) < abs(a)) return b;
	return a;
}

CloudsResult draw_clouds(
	vec3 air_viewer_pos,
	vec3 ray_dir,
	vec3 clear_sky,
	float distance_to_terrain,
	float dither
) {

	CloudsResult result = clouds_not_hit;
	if (defined_clouds == 0u) return result;
	float r = length(air_viewer_pos);

	int ct_index = 0;
	vec2 cloud_types[cloud_layers];

	#ifdef CLOUDS_CUMULUS_CONGESTUS
		cloud_types[ct_index++] = vec2(0., abs(max0(clouds_cumulus_congestus_radius) - r));
	#endif
	#ifdef CLOUDS_CUMULUS
		cloud_types[ct_index++] = vec2(1., abs(clouds_cumulus_radius - r));
	#endif
	#ifdef CLOUDS_CUMULONIMBUS
		cloud_types[ct_index++] = vec2(2., mix(abs(clouds_cumulonimbus_radius - r), 0.0, rainStrength));
	#endif
	#ifdef CLOUDS_ALTOCUMULUS
		cloud_types[ct_index++] = vec2(3., abs(clouds_altocumulus_radius - r));
	#endif
	#ifdef CLOUDS_CIRRUS
		cloud_types[ct_index++] = vec2(4., abs(clouds_cirrus_radius - r));
	#endif

	cloud_types = sort_clouds(cloud_types);

	for(uint i = 0; i < cloud_types.length(); ++i) {

		int cloud_type = int(cloud_types[i].x + 0.5);
		switch (cloud_type) {
#ifdef CLOUDS_CUMULUS_CONGESTUS
		case 0:
			if(clouds_cumulus_congestus_amount >= 1e-3) {
				result = blend_layers(result, draw_cumulus_congestus_clouds(air_viewer_pos, ray_dir, clear_sky, distance_to_terrain, dither), i);
			}
			break;
#endif
#ifdef CLOUDS_CUMULUS
		case 1:
			if(max(clouds_cumulus_coverage.x, clouds_cumulus_coverage.y) >= 1e-3) {
				result = blend_layers(result, draw_cumulus_clouds(air_viewer_pos, ray_dir, clear_sky, distance_to_terrain, dither), i);
			}
			break;
#endif
#ifdef CLOUDS_CUMULONIMBUS
		case 2:
			if(clouds_cumulonimbus_amount >= 1e-3) {
				result = blend_layers(result, draw_cumulonimbus_clouds(air_viewer_pos, ray_dir, clear_sky, distance_to_terrain, dither), i);
			}
			break;
#endif
#ifdef CLOUDS_ALTOCUMULUS
		case 3:
			if(max(clouds_altocumulus_coverage.x, clouds_altocumulus_coverage.y) >= 1e-3) {
				result = blend_layers(result, draw_altocumulus_clouds(air_viewer_pos, ray_dir, clear_sky, distance_to_terrain, dither), i);
			}
			break;
#endif
#ifdef CLOUDS_CIRRUS
		case 4:
			if(max(clouds_cirrus_coverage.x, clouds_cirrus_coverage.y) >= 1e-3) {
				result = blend_layers(result, draw_cirrus_clouds(air_viewer_pos, ray_dir, clear_sky, distance_to_terrain, dither), i);
			}
			break;
#endif
		default: break;
		}

		if (result.transmittance < 1e-3 /*&& (cloud_type != 0 && cloud_type != 2)*/) return result;
	}
#ifdef CLOUDS_NOCTILUCENT	
	vec4 result_nlc = draw_noctilucent_clouds(air_viewer_pos, ray_dir, clear_sky);
	result.scattering += result_nlc.xyz * result.transmittance;
	result.transmittance *= result_nlc.w;
#endif
	return result;
}

#endif // HAS CLOUDS

#endif // INCLUDE_SKY_CLOUDS
