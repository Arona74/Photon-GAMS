#ifndef INCLUDE_LIGHTING_LPV_FLOODFILL
#define INCLUDE_LIGHTING_LPV_FLOODFILL

#include "voxelization.glsl"
#include "/include/lighting/lpv/light_colors.glsl"

bool is_emitter(uint block_id) {
	return 32u <= block_id && block_id < 64u || block_id > 184u && block_id < 192u || block_id > 191u && block_id < 240u ;
}

bool is_translucent(uint block_id) {
	return 164u <= block_id && block_id < 180u;
}

// Workaround for emitter ids 61 and >=64 not working in compute - TODO
bool is_custom(uint block_id) {
	return 64u <= block_id && block_id < 96u;
}
bool is_candle(uint block_id) {
	return 264u <= block_id && block_id < 332u;
}

float get_candle_intensity(uint level) {
	//return level > 0 ? (level > 1 ? (level > 2 ? 10.0 : 8.0) : 6.0) : 3.0;
	return sqr(float(level + 1u)) + 1u;
}

vec3 get_emitted_light(uint block_id) {
	if (is_emitter(block_id)) {
		if (32u <= block_id && block_id < 64u) {
		return texelFetch(light_data_sampler, ivec2(int(block_id) - 32, 0), 0).rgb;
		} else if (block_id > 184u && block_id < 192u) {
			#ifdef GLOWING_ORE
			if (block_id == 185) {
		return vec3(0.85, 0.59, 0.45);
			}
			if (block_id == 186) {
		return light_color[36u];
			}
			if (block_id == 187) {
		return light_color[9u];
			}
			if (block_id == 188) {
		return light_color[38u];
			}
			if (block_id == 189) {
		return vec3(0.76, 0.86, 0.41);
			}
			if (block_id == 190) {
		return light_color[40u];
			}
			if (block_id == 191) {
		return light_color[6u] / 1.5;
			}
			#endif
		} else if (block_id > 191u && block_id < 240u) {
			if (block_id == 192) {
		return light_color[32u] * 3;
			}
			if (block_id == 193) {
		return light_color[33u] * 3;
			}
			if (block_id == 194) {
		return light_color[34u] * 3;
			}
			if (block_id == 195) {
		return light_color[35u] * 3;
			}
			if (block_id == 196) {
		return light_color[36u] * 3;
			}
			if (block_id == 197) {
		return light_color[37u] * 3;
			}
			if (block_id == 198) {
		return light_color[38u] * 3;
			}
			if (block_id == 199) {
		return light_color[39u] * 3;
			}
			if (block_id == 200) {
		return light_color[40u] * 3;
			}
			if (block_id == 201) {
		return light_color[41u] * 3;
			}
			if (block_id == 202) {
		return light_color[42u] * 3;
			}
			if (block_id == 203) {
		return light_color[43u] * 3;
			}
			if (block_id == 204) {
		return light_color[44u] * 3;
			}
			if (block_id == 205) {
		return light_color[45u] * 3;
			}
			if (block_id == 206) {
		return light_color[46u] * 3;
			}
			if (block_id == 207) {
		return light_color[47u] * 3;
			}


			if (block_id == 208) {
		return light_color[32u] * 6;
			}
			if (block_id == 209) {
		return light_color[33u] * 6;
			}
			if (block_id == 210) {
		return light_color[34u] * 6;
			}
			if (block_id == 211) {
		return light_color[35u] * 6;
			}
			if (block_id == 212) {
		return light_color[36u] * 6;
			}
			if (block_id == 213) {
		return light_color[37u] * 6;
			}
			if (block_id == 214) {
		return light_color[38u] * 6;
			}
			if (block_id == 215) {
		return light_color[39u] * 6;
			}
			if (block_id == 216) {
		return light_color[40u] * 6;
			}
			if (block_id == 217) {
		return light_color[41u] * 6;
			}
			if (block_id == 218) {
		return light_color[42u] * 6;
			}
			if (block_id == 219) {
		return light_color[43u] * 6;
			}
			if (block_id == 220) {
		return light_color[44u] * 6;
			}
			if (block_id == 221) {
		return light_color[45u] * 6;
			}
			if (block_id == 222) {
		return light_color[46u] * 6;
			}
			if (block_id == 223) {
		return light_color[47u] * 6;
			}


			if (block_id == 224) {
		return light_color[32u] * 9;
			}
			if (block_id == 225) {
		return light_color[33u] * 9;
			}
			if (block_id == 226) {
		return light_color[34u] * 9;
			}
			if (block_id == 227) {
		return light_color[35u] * 9;
			}
			if (block_id == 228) {
		return light_color[36u] * 9;
			}
			if (block_id == 229) {
		return light_color[37u] * 9;
			}
			if (block_id == 230) {
		return light_color[38u] * 9;
			}
			if (block_id == 231) {
		return light_color[39u] * 9;
			}
			if (block_id == 232) {
		return light_color[40u] * 9;
			}
			if (block_id == 233) {
		return light_color[41u] * 9;
			}
			if (block_id == 234) {
		return light_color[42u] * 9;
			}
			if (block_id == 235) {
		return light_color[43u] * 9;
			}
			if (block_id == 236) {
		return light_color[44u] * 9;
			}
			if (block_id == 237) {
		return light_color[45u] * 9;
			}
			if (block_id == 238) {
		return light_color[46u] * 9;
			}
			if (block_id == 239) {
		return light_color[47u] * 9;
			}
		}
	} else if (is_custom(block_id)) {
		return light_color[block_id - 32u];
	} else if (is_candle(block_id)) {
		if(block_id > 327) { // Uncolored Candle
			return light_color[4u] / 8.0 * get_candle_intensity(block_id - 328u);
		}
		block_id -= 264u;
		uint level = uint(floor(float(block_id) / 16.0));
		float intensity = get_candle_intensity(level);
		return tint_color[block_id - level * 16u] * intensity;
	} else {
		return vec3(0.0);
	}
}

vec3 get_tint(uint block_id, bool is_transparent) {
	if (is_translucent(block_id)) {
		return texelFetch(light_data_sampler, ivec2(int(block_id) - 164, 1), 0).rgb;
	} else {
		return vec3(is_transparent);
	}
}

ivec3 clamp_to_voxel_volume(ivec3 pos) {
	return clamp(pos, ivec3(0), voxel_volume_size - 1);
}

vec3 gather_light(sampler3D light_sampler, ivec3 pos) {
	const ivec3[6] face_offsets = ivec3[6](
		ivec3( 1,  0,  0),
		ivec3( 0,  1,  0),
		ivec3( 0,  0,  1),
		ivec3(-1,  0,  0),
		ivec3( 0, -1,  0),
		ivec3( 0,  0, -1)
	);

	return texelFetch(light_sampler, pos, 0).rgb +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[0]), 0).xyz +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[1]), 0).xyz +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[2]), 0).xyz +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[3]), 0).xyz +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[4]), 0).xyz +
	       texelFetch(light_sampler, clamp_to_voxel_volume(pos + face_offsets[5]), 0).xyz;
}

void update_lpv(writeonly image3D light_img, sampler3D light_sampler) {
	ivec3 pos = ivec3(gl_GlobalInvocationID);
	ivec3 previous_pos = ivec3(vec3(pos) - floor(previousCameraPosition) + floor(cameraPosition));

	uint block_id      = texelFetch(voxel_sampler, pos, 0).x;
	bool transparent   = block_id == 0u || block_id >= 1024u || block_id == 50u || block_id == 44u || block_id == 45u || block_id == 46u || block_id == 49u || block_id == 36u;
	     block_id      = block_id & 1023;
	vec3 light_avg     = gather_light(light_sampler, previous_pos) * rcp(7.0);
	vec3 emitted_light = get_emitted_light(block_id);
	     emitted_light = sqr(pow(emitted_light, vec3(COLOREDLIGHT_SATURATION))) * (1,3 + (1-pow(COLOREDLIGHT_SATURATION, 3.0))) * sign(emitted_light);
	vec3 tint          = sqr(get_tint(block_id, transparent));

	vec3 light = emitted_light + light_avg * tint;

	imageStore(light_img, pos, vec4(light, 0.0));
}

#endif // INCLUDE_LIGHTING_LPV_FLOODFILL
