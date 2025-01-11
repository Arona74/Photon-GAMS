#if !defined INCLUDE_SKY_CLOUDS_CUMULONIMBUS
#define INCLUDE_SKY_CLOUDS_CUMULONIMBUS

// Towering layer: cumulonimbus clouds

#include "common.glsl"

const float clouds_cumulonimbus_radius           = planet_radius + CLOUDS_CUMULONIMBUS_ALTITUDE;
const float clouds_cumulonimbus_thickness        = /*CLOUDS_CUMULONIMBUS_ALTITUDE*/ 1146 * CLOUDS_CUMULONIMBUS_THICKNESS;
const float clouds_cumulonimbus_top_radius       = clouds_cumulonimbus_radius + clouds_cumulonimbus_thickness * CLOUDS_CUMULONIMBUS_END_DISTANCE * 0.00004;
float clouds_cumulonimbus_distance               = mix(CLOUDS_CUMULONIMBUS_DISTANCE, 0.0, rainStrength);
const float clouds_cumulonimbus_end_distance     = CLOUDS_CUMULONIMBUS_END_DISTANCE;
const float clouds_cumulonimbus_blend_distance   = 6.0;
const float clouds_cumulonimbus_extinction_coeff = 0.05 * CLOUDS_CUMULONIMBUS_DENSITY;
float clouds_cumulonimbus_scattering_coeff       = clouds_cumulonimbus_extinction_coeff * (1.0 - 0.33 * rainStrength);

// from https://www.shadertoy.com/view/csSfRK
// x: Coverage signal [0, 1]
// y: Height signal [0, 1]
// shapeParams.x: Horizontal shape [0, 1]
// shapeParams.y: Vertical shape [0, 1]
float CloudShape(float x, float y, vec2 shapeParams) {
	shapeParams.x *= shapeParams.y;
	shapeParams.y = 1.0 / (1.0 - shapeParams.y);

	float anvil = 1.0 - sqr(abs(y - 0.5) * 2.0);
	return clamp01(x - anvil * shapeParams.x - pow(y, shapeParams.y));
}

// altitude_fraction := 0 at the bottom of the cloud layer and 1 at the top
float clouds_cumulonimbus_altitude_shaping(float noise, float altitude_fraction, vec2 blend_distance) {
	float boundary = 1.0;
	vec2 shapeParams = vec2(0.25, 0.90);  // CUMULUS
	// Override shape params
	//shapeParams = vec2(0.00, 0.00); // LINEAR_RAMP
	shapeParams = vec2(0.45, 1.00);  // CUMULONIMBUS

	// World cloud coverage signal. Usually represented by a 2D texture/noise function
	float coverage = noise;
	// Linear height within the cloud layer
	float height = altitude_fraction <= boundary ? linear_step(0.0, boundary, altitude_fraction) : linear_step(boundary, 1.0, altitude_fraction);
	//float shape = CloudShape(1.0, height, altitude_fraction <= boundary ? vec2(0.45, 1.00) : vec2(0.00, 0.50));

	float density = 1.2 * linear_step(0.53, 1.0, noise) * linear_step(0.5, 0.75, daily_weather_variation.clouds_cumulonimbus_amount);

	//density -= 1.0 - CloudShape(1.0, height, shapeParams) /* * 0.5 */;

	// Carve egg shape
	//density -= smoothstep(0.75, 1.0, linear_step(boundary, 1.0, altitude_fraction)) * 0.6;
	float egg = pow1d5(smoothstep(0.2, 1.1, altitude_fraction)) * 0.5;
	//float anvil = 1.0 - CloudShape(1.0, height, shapeParams);
	//float blend = linear_step(0.0, 1.0, clamp01((blend_distance.y - rcp(blend_distance.x)) * blend_distance.x));
	//blend = mix(blend, 1.0, rainStrength);
	//blend = 1.0;
	density -= egg;// * (1.0 - blend);
	//density -= anvil * blend;

	/*if (altitude_fraction < 0.3) density = density;
	else if (altitude_fraction < 0.6) density -= linear_step(0.3, 1.5, altitude_fraction);
	else if (altitude_fraction < 0.8) density -= linear_step(0.2, 1.0, 1.0 - altitude_fraction);*/

	// Reduce density at the top and bottom of the cloud
	density *= smoothstep(0.0, 0.1, altitude_fraction);
	density *= smoothstep(0.0, 0.0001, 1.0 - altitude_fraction);

	return density;
}

float clouds_cumulonimbus_thickness2(float dist2) {
	return clouds_cumulonimbus_thickness * pow1d5(max1(dist2 * 0.2));
}

float clouds_cumulonimbus_dist2(vec3 pos) {
	return max0((length(pos.xz) - clouds_cumulonimbus_distance) * rcp(min(clouds_cumulonimbus_end_distance, 230000)) * 5.75); // [0.0 - 5.75]
}

float clouds_cumulonimbus_density(vec3 pos) {
	const float wind_angle = CLOUDS_CUMULONIMBUS_WIND_ANGLE * degree;
	const vec2 wind_velocity = CLOUDS_CUMULONIMBUS_WIND_SPEED * vec2(cos(wind_angle), sin(wind_angle));

	float r = length(pos);
	if (r < clouds_cumulonimbus_radius || r > clouds_cumulonimbus_top_radius) return 0.0;

	float distance_fraction = linear_step(clouds_cumulonimbus_distance, clouds_cumulonimbus_end_distance, length(pos.xz));
	float distance_fraction_start = linear_step(clouds_cumulonimbus_distance, min(clouds_cumulonimbus_distance * 4.0, clouds_cumulonimbus_end_distance), length(pos.xz));
	float dist2  = clouds_cumulonimbus_dist2(pos);
	float altitude_fraction = (r - clouds_cumulonimbus_radius) * rcp(clouds_cumulonimbus_thickness2(dist2));

	pos.xz += cameraPosition.xz * CLOUDS_SCALE + wind_velocity * (world_age + 50.0 * sqr(altitude_fraction));

	// 2D noise for base shape and coverage
	vec4 noise = texture(noisetex, (0.000004 / CLOUDS_CUMULONIMBUS_SIZE) * pos.xz);

	float density  = clouds_cumulonimbus_altitude_shaping((max1(daily_weather_variation.clouds_cumulonimbus_amount + max0(sqrt(min(dist2, 2.0) / 1.8) - 0.65)) * noise.y - 0.1737 * noise.w), altitude_fraction, vec2(clouds_cumulonimbus_blend_distance, distance_fraction)/*vec2(length(pos.xz), clouds_cumulonimbus_end_distance * 0.25)*/); //.4 * noise.y + .07 * noise.z + .3 * noise.x + .23 * noise.w
	density *= 4.0 * distance_fraction_start * (1.0 - distance_fraction);

	if (density < eps) return 0.0;

	#ifndef PROGRAM_PREPARE
		// Curl noise used to warp the 3D noise into swirling shapes
	vec3 curl = (0.181 * CLOUDS_CUMULONIMBUS_CURL_STRENGTH) * texture(colortex7, 0.0002 * pos).xyz * smoothstep(0.4, 1.0, 1.0 - altitude_fraction);
	vec3 wind = vec3(wind_velocity * world_age, 0.0).xzy;

	// 3D worley noise for detail
	float worley_0 = texture(colortex6, (pos + 0.2 * wind) * 0.00016 + curl * 1.0).x;
	float worley_1 = texture(colortex6, (pos + 0.4 * wind) * 0.0010 + curl * 3.0).x;
	#else
	const float worley_0 = 0.5;
	const float worley_1 = 0.5;
	#endif

		float detail_fade = 0.20 * smoothstep(0.85, 1.0, 1.0 - altitude_fraction)
	- 0.35 * smoothstep(0.05, 0.5, altitude_fraction) + 0.6;
	
	//float detail_strength = CLOUDS_CUMULONIMBUS_DETAIL_STRENGTH * (1.00 - clamp(dist2 - clouds_cumulonimbus_dist2(vec3(30000.0, 0., 0.)), 0.1, 1.0));
	density -= (0.25 * CLOUDS_CUMULONIMBUS_DETAIL_STRENGTH) * sqr(worley_0) * dampen(clamp01(1.0 - density));
	density -= (0.20 * CLOUDS_CUMULONIMBUS_DETAIL_STRENGTH) * sqr(worley_1) * dampen(clamp01(1.0 - density)) * detail_fade;

	// Adjust density so that the clouds are wispy at the bottom and hard at the top
	density  = max0(density);
	density  = 1.0 - pow(max0(1.0 - density), mix(2.0, 7.0, altitude_fraction));
	density *= sqr(linear_step(0.0, 0.5, altitude_fraction));

	return density;
}

float clouds_cumulonimbus_optical_depth(
	vec3 ray_origin,
	vec3 ray_dir,
	float dither,
	const uint step_count
) {
	const float step_growth = 2.0;

	float step_length = 0.1 * clouds_cumulonimbus_thickness / float(step_count); // m

	vec3 ray_pos = ray_origin;
	vec4 ray_step = vec4(ray_dir, 1.0) * step_length;

	float optical_depth = 0.0;

	for (uint i = 0u; i < step_count; ++i, ray_pos += ray_step.xyz) {
		ray_step *= step_growth;
		optical_depth += clouds_cumulonimbus_density(ray_pos + ray_step.xyz * dither) * ray_step.w;
	}

	return optical_depth;
}

vec2 clouds_cumulonimbus_scattering(
	float density,
	float light_optical_depth,
	float sky_optical_depth,
	float ground_optical_depth,
	float step_transmittance,
	float cos_theta,
	float bounced_light
) {
	vec2 scattering = vec2(0.0);

	float scatter_amount = clouds_cumulonimbus_scattering_coeff;
	float extinct_amount = clouds_cumulonimbus_extinction_coeff;

	float scattering_integral_times_density = (1.0 - step_transmittance) / clouds_cumulonimbus_extinction_coeff;

	float powder_effect = clouds_powder_effect(4.0 * density, cos_theta);

	float phase = clouds_phase_single(cos_theta);
	vec3 phase_g = pow(vec3(0.6, 0.9, 0.3), vec3(1.0 + light_optical_depth));

	for (uint i = 0u; i < 8u; ++i) {
		scattering.x += scatter_amount * exp(-extinct_amount *  light_optical_depth) * phase * 1.33;
		scattering.x += scatter_amount * exp(-extinct_amount * ground_optical_depth) * isotropic_phase * bounced_light;
		scattering.y += scatter_amount * exp(-extinct_amount *    sky_optical_depth) * isotropic_phase;

		scatter_amount *= 0.55 * mix(lift(clamp01(clouds_cumulonimbus_scattering_coeff / 0.1), 0.33), 1.0, cos_theta * 0.5 + 0.5) * powder_effect;
		extinct_amount *= 0.5;
		phase_g *= 0.5;

		powder_effect = mix(powder_effect, sqrt(powder_effect), 0.5);

		phase = clouds_phase_multi(cos_theta, phase_g);
	}

	return scattering * scattering_integral_times_density;
}

CloudsResult draw_cumulonimbus_clouds(
	vec3 air_viewer_pos,
	vec3 ray_dir,
	vec3 clear_sky,
	float distance_to_terrain,
	float dither
) {
	// ---------------------
	//   Raymarching Setup
	// ---------------------

	const uint  primary_steps     = CLOUDS_CUMULONIMBUS_PRIMARY_STEPS;
	const uint  lighting_steps    = CLOUDS_CUMULONIMBUS_LIGHTING_STEPS;
	const uint  ambient_steps     = CLOUDS_CUMULONIMBUS_AMBIENT_STEPS;

	const float min_transmittance = 0.075;

	const float planet_albedo     = 0.4;
	const vec3  sky_dir           = vec3(0.0, 1.0, 0.0);

	float r = length(air_viewer_pos);

	vec2 sphere_dists   = intersect_spherical_shell(air_viewer_pos, ray_dir, clouds_cumulonimbus_radius, clouds_cumulonimbus_top_radius);
	vec2 cylinder_dists = intersect_cylindrical_shell(air_viewer_pos, ray_dir, clouds_cumulonimbus_distance, clouds_cumulonimbus_end_distance);
	vec2 dists          = vec2(max(sphere_dists.x, cylinder_dists.x), min(sphere_dists.y, cylinder_dists.y));

	bool planet_intersected = intersect_sphere(air_viewer_pos, ray_dir, min(r - 10.0, planet_radius)).y >= 0.0;
	bool terrain_intersected = distance_to_terrain >= 0.0 && r < clouds_cumulonimbus_radius && distance_to_terrain * CLOUDS_SCALE < dists.y;

	if (dists.y < 0.0                                        // volume not intersected
	 || planet_intersected && r < clouds_cumulonimbus_radius // planet blocking clouds
	 || terrain_intersected                                  // terrain blocking clouds
	) { return clouds_not_hit; }

	float ray_length = (distance_to_terrain >= 0.0) ? distance_to_terrain : dists.y;
	      ray_length = max0(ray_length - dists.x);
	float step_length = ray_length * rcp(float(primary_steps));

	vec3 ray_step = ray_dir * step_length;

	vec3 ray_origin = air_viewer_pos + ray_dir * (dists.x + step_length * dither);

	vec2 scattering = vec2(0.0); // x: direct light, y: skylight
	float transmittance = 1.0;

	float distance_sum = 0.0;
	float distance_weight_sum = 0.0;

	// ------------------
	//   Lighting Setup
	// ------------------

	bool moonlit = sun_dir.y < -0.04;
	vec3 light_dir = moonlit ? moon_dir : sun_dir;
	float cos_theta = dot(ray_dir, light_dir);
	float bounced_light = planet_albedo * light_dir.y * rcp_pi;
	
	// --------------------
	//   Raymarching Loop
	// --------------------

	for (uint i = 0u; i < primary_steps; ++i) {
		if (transmittance < min_transmittance) break;

		vec3 ray_pos = ray_origin + ray_step * i;

		float altitude_fraction = (length(ray_pos) - clouds_cumulonimbus_radius) * rcp(clouds_cumulonimbus_thickness2(clouds_cumulonimbus_dist2(ray_pos)));

		float density = clouds_cumulonimbus_density(ray_pos);

		if (density < eps) continue;

		float distance_to_sample = distance(ray_origin, ray_pos);

		float step_optical_depth = density * clouds_cumulonimbus_extinction_coeff * step_length;
		float step_transmittance = exp(-step_optical_depth);

#if defined PROGRAM_DEFERRED0
		vec2 hash = vec2(0.0);
#else
		vec2 hash = hash2(fract(ray_pos)); // used to dither the light rays
#endif

		float light_optical_depth  = clouds_cumulonimbus_optical_depth(ray_pos, light_dir, hash.x, lighting_steps);
		float sky_optical_depth    = clouds_cumulonimbus_optical_depth(ray_pos, sky_dir, hash.y, ambient_steps);
		float ground_optical_depth = mix(density, 1.0, clamp01(altitude_fraction * 2.0 - 1.0)) * altitude_fraction * clouds_cumulonimbus_thickness2(clouds_cumulonimbus_dist2(ray_pos)); // guess optical depth to the ground using altitude fraction and density from this sample

		scattering += clouds_cumulonimbus_scattering(
		density,
		light_optical_depth,
		sky_optical_depth,
		ground_optical_depth,
		step_transmittance,
		cos_theta,
		bounced_light
		) * transmittance;

		transmittance *= step_transmittance;

		// Update distance to cloud
		distance_sum += distance_to_sample * density;
		distance_weight_sum += density;
	}

	// Get main light color for this layer
	vec3 light_color  = sunlight_color * atmosphere_transmittance(ray_origin, light_dir);
		 light_color  = atmosphere_post_processing(light_color);
	     light_color *= moonlit ? moon_color : sun_color;
		 light_color *= 1.0 - rainStrength;

	// Remap the transmittance so that min_transmittance is 0
	float clouds_transmittance = linear_step(min_transmittance, 1.0, transmittance);

	// Aerial perspective
	vec3 clouds_scattering = scattering.x * light_color + scattering.y * sky_color;
	if (distance_to_terrain < 0.0) clouds_scattering = clouds_aerial_perspective(clouds_scattering, clouds_transmittance, air_viewer_pos, ray_origin, ray_dir, clear_sky);

	float apparent_distance = (distance_weight_sum == 0.0)
		? 1e6
		: (distance_sum / distance_weight_sum) + distance(air_viewer_pos, ray_origin);

	return CloudsResult(
		clouds_scattering,
		clouds_transmittance,
		apparent_distance
	);

}
#endif
