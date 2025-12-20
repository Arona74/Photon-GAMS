#if !defined INCLUDE_LIGHTING_LPV_BLOCKLIGHT
#define INCLUDE_LIGHTING_LPV_BLOCKLIGHT

#include "voxelization.glsl"
#include "/include/surface/material.glsl"

#ifdef DIRECTIONAL_LIGHTMAPS
//#include "/include/lighting/specular_lighting.glsl"
#include "/include/lighting/directional_lightmaps.glsl"
#endif

vec3 read_lpv_linear(vec3 pos) {
	if ((frameCounter & 1) == 0) {
		return texture(light_sampler_a, pos).rgb;
	} else {
		return texture(light_sampler_b, pos).rgb;
	}
}

float lpv_distance_fade(vec3 scene_pos) {
	float distance_fade  = 2.0 * max_of(abs(scene_pos / vec3(voxel_volume_size)));
		  distance_fade  = linear_step(0.75, 1.0, distance_fade);

	return distance_fade;
}

vec3 get_lpv_linear_voxel(vec3 voxel_pos, vec3 scene_pos) {
	vec3 sample_pos = clamp01(voxel_pos / vec3(voxel_volume_size));
	vec3 lpv_linear = read_lpv_linear(sample_pos);

	float distance_fade = lpv_distance_fade(scene_pos);

	return lpv_linear = mix(lpv_linear, vec3(0.0), distance_fade);
}

vec3 get_lpv_linear_voxel(vec3 voxel_pos) {
	return get_lpv_linear_voxel(voxel_pos, voxel_to_scene_space(voxel_pos));
}

vec3 get_lpv_linear(vec3 scene_pos) {
	vec3 voxel_pos = scene_to_voxel_space(scene_pos);

	if (is_inside_voxel_volume(voxel_pos)) {
		return get_lpv_linear_voxel(voxel_pos, scene_pos);
	} else {
		return vec3(0.0);
	}
}

vec3 get_lpv_basic_voxel(vec3 voxel_pos, vec3 scene_pos) {
	vec3 sample_pos = clamp01(voxel_pos / vec3(voxel_volume_size));
	vec3 lpv_blocklight = sqrt(read_lpv_linear(sample_pos));

	float distance_fade = lpv_distance_fade(scene_pos);

	return lpv_blocklight = mix(lpv_blocklight, vec3(0.0), distance_fade);
}

vec3 get_lpv_basic(vec3 scene_pos) {
	vec3 voxel_pos = scene_to_voxel_space(scene_pos);

	if (is_inside_voxel_volume(voxel_pos)) {
		return get_lpv_basic_voxel(voxel_pos, scene_pos);
	} else {
		return vec3(0.0);
	}
}

vec3 lpv_fog_curve(vec3 lpv_blocklight) {
	float lpv_max = max_of(lpv_blocklight) + 1e-3;
	lpv_blocklight *= (pow(1.08, 4.0 * lpv_max) - 1.0) / lpv_max;
	return lpv_blocklight * COLORED_LIGHT_I;
}

vec3 get_lpv_fog(vec3 scene_pos) { return lpv_fog_curve(get_lpv_basic(scene_pos)); }

vec3 get_lpv_direction(vec3 scene_pos) {
	const vec3 epsilon = vec3(0.1);

	vec3 lpv_light1 = get_lpv_basic(scene_pos - epsilon);
	vec3 lpv_light2 = get_lpv_basic(scene_pos + epsilon);
	vec3 light_dir = ((lpv_light2.r + lpv_light2.g + lpv_light2.b) - (lpv_light1.r + lpv_light1.g + lpv_light1.b)) / (epsilon * 2.0);

	return normalize(light_dir);
}

#define average3(v) (((v).x + (v).y + (v).z) / 3.0)

mat3 get_lpv_gradient_rgb(vec3 scene_pos, float epsilon/*, vec3 normal, vec3 lpv_blocklight*/) {
	//float epsilon = 0.25;//exp2(-length(lpv_blocklight));
	//scene_pos += normal * (2.0 * epsilon + eps);
	//vec3 light_dir;
	float epl = epsilon;

	mat3 lpv_light_1 = mat3(
		get_lpv_linear(scene_pos - vec3(epl, 0.0, 0.0)),
		get_lpv_linear(scene_pos - vec3(0.0, epl, 0.0)),
		get_lpv_linear(scene_pos - vec3(0.0, 0.0, epl))
	);

	mat3 lpv_light_2 = mat3(
		get_lpv_linear(scene_pos + vec3(epl, 0.0, 0.0)),
		get_lpv_linear(scene_pos + vec3(0.0, epl, 0.0)),
		get_lpv_linear(scene_pos + vec3(0.0, 0.0, epl))
	);

	/*mat3 lpv_light_3 = mat3(
		get_lpv_linear(scene_pos - vec3(epl, 0.0, 0.0)),
		get_lpv_linear(scene_pos - vec3(epl, epl, 0.0)),
		get_lpv_linear(scene_pos - vec3(0.0, 0.0, epl))
	)*/

	/*mat3 gradient = mat3(
		(lpv_light_2[0] - lpv_light_1[0]) / (2.0 * epsilon),
		(lpv_light_2[1] - lpv_light_1[1]) / (2.0 * epsilon),
		(lpv_light_2[2] - lpv_light_1[2]) / (2.0 * epsilon)
	);*/
	mat3 gradient = (lpv_light_2 - lpv_light_1) / (2.0 * epsilon);

	/*vec3 lpv_light1 = get_lpv_linear(scene_pos - epsilon);
	vec3 lpv_light2 = get_lpv_linear(scene_pos + epsilon);
	vec3 dir = ((lpv_light2.r + lpv_light2.g + lpv_light2.b) - (lpv_light1.r + lpv_light1.g + lpv_light1.b)) / vec3(epsilon * 2.0);
	mat3 gradient = mat3(dir, dir, dir);*/

	//light_dir = vec3(length(gradient[0]))

	//light_dir = vec3(average3(gradient[0]), average3(gradient[1]), average3(gradient[2]));

	/*vec3 lpv_light1 = get_lpv_basic(scene_pos - epsilon);
	vec3 lpv_light2 = get_lpv_basic(scene_pos + epsilon);
	float light_diff = length(lpv_light2 - lpv_light1);

	light_dir = light_diff / (epsilon * 2.0);*/
	//light_dir[1] = (lpv_light2.g - lpv_light1.g) / (epsilon * 2.0);
	//light_dir[2] = (lpv_light2.b - lpv_light1.b) / (epsilon * 2.0);
	return gradient;
	//return light_dir;
}

vec3 get_directional_lpv_mult(vec3 normal, vec3 scene_pos, vec3 lpv_blocklight) {
	vec3 blocklight_mul = vec3(1.0);

	mat3 blocklight_dir = transpose(get_lpv_gradient_rgb(scene_pos, 0.25/*, normal, lpv_blocklight*/));

	//if(length_squared(blocklight_dir) > 1e-12) {
	//if (max_of(vec3(length_squared(blocklight_dir[0]), length_squared(blocklight_dir[1]), length_squared(blocklight_dir[2]))) > eps) {
		vec3 NoL = vec3(dot(normalize_safe(blocklight_dir[0]), normal), dot(normalize_safe(blocklight_dir[1]), normal), dot(normalize_safe(blocklight_dir[2]), normal));
		//NoL = vec3(1.0, 0.0, 0.5);
		blocklight_mul = (clamp(NoL, 0.0, 1.0) * DIRECTIONAL_LIGHTMAPS_INTENSITY + (1.0 - DIRECTIONAL_LIGHTMAPS_INTENSITY)) * inversesqrt(sqrt(lpv_blocklight) + eps);
	//}
	return blocklight_mul;  //clamp01(normalize(normal * blocklight_dir) + 0.8)
}

vec3 get_lpv_direction_fog(vec3 scene_pos) {
	//dither = clamp(dither, 0.1, 1.0);
	/*const vec3 epsilon = vec3(0.1);

	//vec3 lpv_light1 = get_lpv_basic;clamp01(voxel_pos / vec3(voxel_volume_size));
	//vec3 lpv_light2 = get_lpv_basic(scene_pos + epsilon);
	vec3 lpv_light1 = get_lpv_linear(scene_pos - epsilon);
	vec3 lpv_light2 = get_lpv_linear(scene_pos + epsilon);
	vec3 light_dir = ((lpv_light2.x + lpv_light2.y + lpv_light2.z) - (lpv_light1.x + lpv_light1.y + lpv_light1.z)) / (epsilon * 2.0);*/
	mat3 light_dir_rgb = get_lpv_gradient_rgb(scene_pos, 0.1);
	light_dir_rgb = transpose(light_dir_rgb);
	vec3 light_dir = light_dir_rgb[0] + light_dir_rgb[1] + light_dir_rgb[2];

	return light_dir;
}

vec3 get_lpv_direction_clouds(vec3 scene_pos, vec2 dither) {
	vec3 voxel_pos = scene_to_voxel_space(scene_pos);
	//dither = clamp(dither, 0.1, 1.0);
	const float epsilon = 1.0;

	vec3 offset;
	float sin_x = sin(dither.x);
	offset.y = cos(dither.x);
	offset.x = sin_x * cos(dither.y);
	offset.z = sin_x * sin(dither.y);
	//offset = normalize(offset) * epsilon;
	offset *= epsilon;
	//vec3 lpv_light1 = get_lpv_basic;clamp01(voxel_pos / vec3(voxel_volume_size));
	//vec3 lpv_light2 = get_lpv_basic(scene_pos + offset);
	vec3 lpv_light1 = get_lpv_linear_voxel(voxel_pos - offset);
	vec3 lpv_light2 = get_lpv_linear_voxel(voxel_pos + offset);
	vec3 light_dir = ((lpv_light2.x + lpv_light2.y + lpv_light2.z) - (lpv_light1.x + lpv_light1.y + lpv_light1.z)) / (offset * 2.0);

	return light_dir;
}

#if !defined(PROGRAM_COMPOSITE0) && !defined(PROGRAM_DEFERRED0)

vec3 get_lpv_blocklight(vec3 scene_pos, vec3 normal, vec3 flat_normal, vec3 mc_blocklight, float ao, Material material) {
	vec3 voxel_pos = scene_to_voxel_space(scene_pos);
	//vec3 lpv_specular = vec3(0.0);
	float normal_diff_angle = 0.0;
	vec3 lpv_normal = vec3(0.0);

	if (is_inside_voxel_volume(voxel_pos)) {

#ifndef NO_NORMAL
	vec3 lpv_specular = vec3(1.0);
	voxel_pos += flat_normal * 0.5;

	#ifdef DIRECTIONAL_LIGHTMAPS
	vec3 normal_diff = normal - flat_normal;
	//normal_diff = sign(normal_diff) * sqrt(abs(normal_diff));
	// 0 at angle == 0 rad, 1 at angle >= PI rad
	normal_diff_angle = 1.0 - clamp01(dot(normalize(normal_diff), normalize(flat_normal)));
	vec3 normal_boost = normal_diff + sqrt(normal_diff_angle) * normal_diff * 1.0;
	lpv_normal += normal_boost /* 2.5*/;
	#endif
	
#else
	ao = 1.0;
#endif

		vec3 sample_pos = clamp01((voxel_pos + lpv_normal) / vec3(voxel_volume_size));
		vec3 lpv_blocklight = read_lpv_linear(sample_pos);

/*#if !defined(NO_NORMAL) && defined(DIRECTIONAL_LIGHTMAPS)
		// Fix directional lightmap occlusion
		const float min_lpv = 0.05;
		const float max_lpv_diff = 5.0;
		sample_pos = clamp01((voxel_pos) / vec3(voxel_volume_size));
		vec3 lpv_blocklight2 = read_lpv_linear(sample_pos);

		lpv_blocklight = min(lpv_blocklight, lpv_blocklight2 * max_lpv_diff);
		lpv_blocklight = max(lpv_blocklight, lpv_blocklight2 * min_lpv);
#endif*/

		lpv_blocklight = sqrt(lpv_blocklight) * ao;

/*#ifdef DIRECTIONAL_LIGHTMAPS
		//lpv_blocklight *= 1.0 + normal_diff_angle;
		//lpv_blocklight *= mix(vec3(1.0), vec3(8.0 * normal_diff_angle), ao);
		lpv_blocklight *= material.f0 * (0.5 + normal_diff_angle * 8.0) + (1.0 - material.f0);
#endif*/

#if !defined(NO_NORMAL) && defined(DIRECTIONAL_LIGHTMAPS) && DIRECTIONAL_LIGHTMAPS_INTENSITY > 0.0
		//lpv_blocklight *= get_directional_lpv(normal, scene_pos, lpv_blocklight);
		vec3 sample_pos_normal = clamp01((voxel_pos + lpv_normal) / vec3(voxel_volume_size));
		vec3 lpv_blocklight_normal = read_lpv_linear(sample_pos_normal);
		lpv_blocklight_normal = sqrt(lpv_blocklight_normal) * ao;
		vec3 lpv_blocklight_directional = lpv_blocklight * get_directional_lpv_mult(normal - flat_normal, scene_pos, lpv_blocklight);

		vec3 lpv_blocklight_x_directional = lpv_blocklight * lpv_blocklight_directional;
		float dot_nd = dot(lpv_blocklight_normal, lpv_blocklight_x_directional);
		if (abs(max_of(lpv_blocklight_normal) - max_of(lpv_blocklight * lpv_blocklight_directional)) < 0.2)
			lpv_blocklight += lpv_blocklight_normal * lpv_blocklight_directional;
		else 
			lpv_blocklight += lpv_blocklight * lpv_blocklight_directional;
#endif

/*#ifdef DIRECTIONAL_LIGHTMAPS
	vec3 world_dir = normalize(scene_pos - gbufferModelViewInverse[3].xyz);
	float light_length = length_squared(lpv_blocklight);
	vec2 light_gradient = vec2(dFdx(light_length), dFdy(light_length));

	if (length_squared(light_gradient) > 1e-12) {
		mat2x3 pos_gradient = mat2x3(dFdx(scene_pos), dFdy(scene_pos));
		vec3 light_dir = pos_gradient * light_gradient;

		float NoL = dot(normal, light_dir);
		float NoV = clamp01(dot(normal, -world_dir));
		float LoV = dot(light_dir, -world_dir);
		float halfway_norm = inversesqrt(2.0 * LoV + 2.0);
		float NoH = (NoL + NoV) * halfway_norm;
		float LoH = LoV * halfway_norm + halfway_norm;

		vec3 lpv_specular = get_specular_highlight(material, NoL, NoV, NoH, LoV, LoH);
		lpv_blocklight += lpv_specular * lpv_blocklight;
	}
#endif*/

#ifdef COLORED_LIGHTS_VANILLA_LIGHTMAP_CONTRIBUTION
		float vanilla_lightmap_contribution = exp2(-4.0 * dot(lpv_blocklight, luminance_weights_rec2020));
		lpv_blocklight += mc_blocklight * vanilla_lightmap_contribution;
#endif

        // Darkness effect
        float darkness_factor =
            mix(1.0,
                dampen(abs(cos(2.0 * frameTimeCounter))) * 0.67 + 0.2,
                darknessFactor) *
                0.75 +
            0.25;
        lpv_blocklight *= darkness_factor;

		float distance_fade = lpv_distance_fade(scene_pos);

		return mix(lpv_blocklight * COLORED_LIGHT_I, mc_blocklight, distance_fade);
	} else {
		return mc_blocklight;
	}
}

#endif // PROGRAMS

#endif // INCLUDE_LIGHTING_LPV_BLOCKLIGHT
