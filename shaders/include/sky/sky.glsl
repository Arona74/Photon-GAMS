#if !defined INCLUDE_SKY_SKY
#define INCLUDE_SKY_SKY

#include "/include/utility/color.glsl"
#include "/include/utility/dithering.glsl"
#include "/include/utility/fast_math.glsl"
#include "/include/utility/random.glsl"

// Stars based on https://www.shadertoy.com/view/Md2SR3

vec3 unstable_star_field(vec2 coord, float star_threshold) {
	const float min_temp = 3500.0;
	const float max_temp = 9500.0;

	vec4 noise = hash4(coord);

	float star = linear_step(star_threshold, 1.0, noise.x);
	      star = pow4(star) * STARS_INTENSITY;

	float temp = mix(min_temp, max_temp, noise.y);
	vec3 color = blackbody(temp);

	const float twinkle_speed = 2.0;
	float twinkle_amount = noise.z;
	float twinkle_offset = tau * noise.w;
	star *= 1.0 - twinkle_amount * cos(frameTimeCounter * twinkle_speed + twinkle_offset);

	return star * color;
}

// Stabilizes the star field by sampling at the four neighboring integer coordinates and
// interpolating
vec3 stable_star_field(vec2 coord, float star_threshold) {
	coord = abs(coord) + 33.3 * step(0.0, coord);
	vec2 i, f = modf(coord, i);

	f.x = cubic_smooth(f.x);
	f.y = cubic_smooth(f.y);

	return unstable_star_field(i + vec2(0.0, 0.0), star_threshold) * (1.0 - f.x) * (1.0 - f.y)
	     + unstable_star_field(i + vec2(1.0, 0.0), star_threshold) * f.x * (1.0 - f.y)
	     + unstable_star_field(i + vec2(0.0, 1.0), star_threshold) * f.y * (1.0 - f.x)
	     + unstable_star_field(i + vec2(1.0, 1.0), star_threshold) * f.x * f.y;
}

vec3 draw_stars(vec3 ray_dir, float galaxy_luminance) {
	// Adjust star threshold so that brightest stars appear first
#if defined WORLD_OVERWORLD
	float star_threshold = 1.0 - 0.008 * STARS_COVERAGE * smoothstep(-0.2, 0.5, -sun_dir.y) - 0.5 * cube(galaxy_luminance);
#else
	float star_threshold = 1.0 - 0.008 * STARS_COVERAGE;
#endif

	// Project ray direction onto the plane
	vec2 coord  = ray_dir.xy * rcp(abs(ray_dir.z) + length(ray_dir.xy)) + 41.21 * sign(ray_dir.z);
	     coord *= 600.0;

	return stable_star_field(coord, star_threshold);
}

//----------------------------------------------------------------------------//
#if defined WORLD_OVERWORLD

#include "/include/lighting/colors/light_color.glsl"
#include "/include/lighting/colors/weather_color.glsl"
#include "/include/lighting/bsdf.glsl"
#include "/include/sky/atmosphere.glsl"
#include "/include/sky/projection.glsl"
#include "/include/utility/geometry.glsl"
#include "/include/sky/shooting_stars.glsl"
#include "/include/sky/nebula.glsl"

const float sun_luminance  = SUN_LUMINANCE * SUN_DISK_INTENSITY; // luminance of sun disk
const float moon_luminance = MOON_LUMINANCE * MOON_DISK_INTENSITY; // luminance of moon disk

vec3 draw_sun(vec3 ray_dir) {
	float nu = dot(ray_dir, sun_dir);

	// Limb darkening model from http://www.physics.hmc.edu/faculty/esin/a101/limbdarkening.pdf
	const vec3 alpha = vec3(0.429, 0.522, 0.614);
	float center_to_edge = max0(sun_angular_radius - fast_acos(nu));
	vec3 limb_darkening = pow(vec3(1.0 - sqr(1.0 - center_to_edge)), 0.5 * alpha);

	return sun_luminance * sun_color * step(0.0, center_to_edge) * limb_darkening;
}

#if defined GALAXY

#if defined GALAXY_GAMS

// Galaxy from old Photon-GAMS

vec3 draw_galaxy(vec3 ray_dir, out float galaxy_luminance) {
	//const float galaxy_intensity = GALAXY_INTENSITY;
	const vec3 galaxy_tint = vec3(GALAXY_TINT_R, GALAXY_TINT_G, GALAXY_TINT_B) * GALAXY_INTENSITY;
	// Check if it's night time
	//if (sun_dir.y > -0.05) return vec3(0.0); // Return black if it's not night
	// Convert ray direction to spherical coordinates
	float phi = atan(ray_dir.y, ray_dir.x);
	float theta = acos(ray_dir.z);
	// Map spherical coordinates to UV coordinates
	vec2 uv = vec2(phi / (2.0 * pi) + 0.5, theta / pi);

	vec3 galaxy = from_srgb(texture(colortex13, uv).rgb);

	// Fade in/out at twilight
	float night_factor = smoothstep(0.01, 0.5, -sun_dir.y);

	return galaxy * galaxy_tint * night_factor;
}

#else

// GALAXY from Photon

vec3 draw_galaxy(vec3 ray_dir, out float galaxy_luminance) {
	const vec3 galaxy_tint = vec3(GALAXY_TINT_R, GALAXY_TINT_G, GALAXY_TINT_B) * GALAXY_INTENSITY;

	float galaxy_intensity = 0.05 + 1.0 * linear_step(-0.1, 0.25, -sun_dir.y);

	float lon = atan(ray_dir.x, ray_dir.z);
	float lat = fast_acos(-ray_dir.y);

	vec3 galaxy = texture(
		galaxy_sampler,
		vec2(lon * rcp(tau) + 0.5, lat * rcp(pi))
	).rgb;

	galaxy = srgb_eotf_inv(galaxy) * rec709_to_working_color;

	galaxy *= 2 * galaxy_intensity * galaxy_tint;

	galaxy_luminance = dot(galaxy, luminance_weights_rec709);

	galaxy = mix(
		vec3(galaxy_luminance),
		galaxy,
		2.0
	);

	return max0(galaxy);
}
#endif

#endif

vec3 adjust_night_atmosphere(vec3 atmosphere, vec3 ray_dir) {
	#ifdef BLACK_NIGHT_SKY
	float night_factor = smoothstep(0.1, -0.1, sun_dir.y);
	float height_fade = smoothstep(-0.1, 0.3, ray_dir.y);

	float blue_hour = linear_step(0.05, 1.0, exp(-190.0 * sqr(sun_dir.y + 0.09604)));
	vec3 blue_hour_tint = vec3(0.95, 0.80, 1.0);
	vec3 blue_hour_sky = mix(atmosphere, atmosphere * blue_hour_tint, blue_hour);

	vec3 night_sky = mix(atmosphere * 0.1, vec3(0.0), height_fade);
	vec3 blended_sky = mix(blue_hour_sky, night_sky, night_factor);

	return mix(atmosphere, blended_sky, smoothstep(0.2, -0.2, sun_dir.y));
	#else
	return atmosphere;
	#endif
}

vec4 get_clouds_and_aurora(vec3 ray_dir, vec3 clear_sky) {
#if defined PROGRAM_DEFERRED0
	ivec2 texel   = ivec2(gl_FragCoord.xy);
	      texel.x = texel.x % (sky_map_res.x - 4);

	float dither = interleaved_gradient_noise(vec2(texel));

	// Render clouds
	#ifndef BLOCKY_CLOUDS
	const vec3 air_viewer_pos = vec3(0.0, planet_radius, 0.0);
	CloudsResult result = draw_clouds(air_viewer_pos, ray_dir, clear_sky, -1.0, dither);
	#else
	CloudsResult result = clouds_not_hit;
	#endif

	// Render aurora
	vec3 aurora = draw_aurora(ray_dir, dither);

	return vec4(
		result.scattering + aurora * result.transmittance,
		result.transmittance
	);
#else
	return vec4(0.0, 0.0, 0.0, 1.0);
#endif
}

vec3 draw_sky(vec3 ray_dir, vec3 atmosphere) {
	vec3 sky = vec3(0.0);

#if defined SHADOW
	// Trick to make stars rotate with sun and moon
	mat3 rot = (sunAngle < 0.5)
		? mat3(shadowModelViewInverse)
		: mat3(-shadowModelViewInverse[0].xyz, shadowModelViewInverse[1].xyz, -shadowModelViewInverse[2].xyz);

	vec3 celestial_dir = ray_dir * rot;
#else
	vec3 celestial_dir = ray_dir;
#endif

	// Galaxy
#ifdef GALAXY
	float galaxy_luminance;
	#ifdef GALAXY_ROTATION
	sky += draw_galaxy(celestial_dir, galaxy_luminance);
	#else
	sky += draw_galaxy(ray_dir, galaxy_luminance);
	#endif
	
#else
	const float galaxy_luminance = 0.0;
#endif

	// Sun, moon and stars

#if defined PROGRAM_DEFERRED4
	// Output of skytextured
	sky += texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0).rgb;

#ifdef STARS
	// Stars
	#ifdef STARS_ROTATION
	sky += draw_stars(celestial_dir, galaxy_luminance);
	#else
	sky += draw_stars(ray_dir, galaxy_luminance);
	#endif
#endif

#ifndef VANILLA_SUN
	// Sun
	sky += draw_sun(ray_dir);
#endif
	
	// DEBUG Try to fix dimmed moon ############################################
#if MOON_TYPE == MOON_VANILLA
	//vec4 vanilla_sky = texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0);
	vec3 vanilla_sky_rgb = texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0).rgb;
	vec3 vanilla_sky_color = from_srgb(vanilla_sky_rgb);
	//uint vanilla_sky_id = uint(255.0 * vanilla_sky.a);
	
	/*if (max_of(vanilla_sky_color) > 0.2) {
		//const vec3 brightness_scale = sunlight_color * moon_luminance;
		//if(dot(vanilla_sky_color, vec3(1.0)) > 1e-3) sky *= 0.0; // Hide stars behind moon
		sky *= 0.0;
	}*/
	
	sky += vanilla_sky_color * (sunlight_color * moon_luminance);
	
#elif MOON_TYPE == MOON_PHOTON
	//vec4 vanilla_sky = texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0);
	vec3 vanilla_sky_rgb = texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0).rgb;
	vec3 vanilla_sky_color = from_srgb(vanilla_sky_rgb);
	//uint vanilla_sky_id = uint(255.0 * vanilla_sky.a);
	
	if (max_of(vanilla_sky_color) > 0.00001) {
		//const vec3 brightness_scale = sunlight_color * moon_luminance;
		//if(dot(vanilla_sky_color, vec3(1.0)) > 1e-3) sky *= 0.0; // Hide stars behind moon
		//sky += vanilla_sky_color * brightness_scale;
		sky *= 0.0;
	}
	
	sky += vanilla_sky_color * (sunlight_color * moon_luminance);
	
#endif
	// END OF DEBUG Try to fix dimmed moon ############################################

#endif

	// Sky gradient
	atmosphere = adjust_night_atmosphere(atmosphere, ray_dir);
	sky *= atmosphere_transmittance(ray_dir.y, planet_radius) * (1.0 - rainStrength);
	sky += atmosphere;

	// Rain
	vec3 rain_sky = get_weather_color() * (1.0 - exp2(-0.8 / clamp01(ray_dir.y)));
	sky = mix(sky, rain_sky, rainStrength * mix(1.0, 0.9, time_sunrise + time_sunset));

	// Clouds
	vec4 clouds = get_clouds_and_aurora(ray_dir, sky);
	sky *= clouds.a;   // transmittance
	sky += clouds.rgb; // scattering

	// Shooting stars
#if defined SHOOTING_STARS && !defined PROGRAM_DEFERRED0
	sky = DrawShootingStars(sky, ray_dir);
#endif
	
	// Nebula
	sky = draw_nebula(ray_dir, sky);

	// Fade lower part of sky into cave fog color when underground so that the sky isn't visible
	// beyond the render distance
	float underground_sky_fade = biome_cave * smoothstep(-0.1, 0.1, 0.4 - ray_dir.y);
	sky = mix(sky, vec3(0.0), underground_sky_fade);

	return sky;
}

vec3 draw_sky(vec3 ray_dir) {
	
	vec3 atmosphere = atmosphere_scattering(ray_dir, sun_color, sun_dir, moon_color, moon_dir, true);
	return draw_sky(ray_dir, atmosphere);
}

//----------------------------------------------------------------------------//
#elif defined WORLD_NETHER

vec3 draw_sky(vec3 ray_dir) {
	return ambient_color;
}

//----------------------------------------------------------------------------//
#elif defined WORLD_END

#include "/include/misc/end_lighting_fix.glsl"
#include "/include/sky/atmosphere.glsl"

const float sun_solid_angle = cone_angle_to_solid_angle(sun_angular_radius);
const vec3 end_sun_color = vec3(END_SOLAR_FLARE_COLOR_R, END_SOLAR_FLARE_COLOR_G, END_SOLAR_FLARE_COLOR_B);

vec3 draw_sun(vec3 ray_dir) {
	float nu = dot(ray_dir, sun_dir);
	float r = fast_acos(nu);

	// Sun disk

	const vec3 alpha = vec3(0.6, 0.5, 0.4);
	float center_to_edge = max0(sun_angular_radius - r);
	vec3 limb_darkening = pow(vec3(1.0 - sqr(1.0 - center_to_edge)), 0.5 * alpha);
	vec3 sun_disk = vec3(r < sun_angular_radius);

    // Solar flare effect

#ifdef END_SOLAR_FLARE_ENABLED

	// Transform the coordinate space such that z is parallel to sun_dir
    vec3 tangent = sun_dir.y == 11.0 ? vec3(1.0, 0.0, 0.0) : normalize(cross(vec3(1.0, 0.0, 0.0), sun_dir));
    vec3 bitangent = normalize(cross(tangent, sun_dir));
    mat3 rot = mat3(tangent, bitangent, sun_dir);

	// Vector from ray dir to sun dir
    vec2 q = ((ray_dir - sun_dir) * rot).xy;

    float theta = fract(linear_step(-pi, pi, atan(q.y, q.x)) + 0.015 * frameTimeCounter - 0.33 * r);

    float flare = texture(noisetex, vec2(theta, r - END_SOLAR_FLARE_SPEED * frameTimeCounter)).x;
          flare = pow5(flare) * exp(END_SOLAR_FLARE_FALLOFF * (r - sun_angular_radius));
          flare = r < sun_angular_radius ? 0.0 : flare;

    return end_sun_color * rcp(sun_solid_angle) * max0(sun_disk + END_SOLAR_FLARE_INTENSITY * flare);
#else
    return end_sun_color * rcp(sun_solid_angle) * sun_disk;
#endif
}

vec3 draw_sky(vec3 ray_dir) {

	// Sky gradient

	float up_gradient = linear_step(0.0, 0.4, ray_dir.y) + linear_step(0.1, 0.8, -ray_dir.y);
	vec3 sky = ambient_color * mix(0.1, 0.04, up_gradient);
	float mie_phase = cornette_shanks_phase(dot(ray_dir, sun_dir), 0.6);
	sky += 0.1 * (ambient_color + 0.5 * end_sun_color) * mie_phase;

#if defined PROGRAM_DEFERRED4
	// Sun

	sky += draw_sun(ray_dir);

	// Stars

	vec3 stars_fade = exp2(-0.1 * max0(1.0 - ray_dir.y) / max(ambient_color, eps)) * linear_step(-0.2, 0.0, ray_dir.y);
	sky += draw_stars(ray_dir, 0.0).xzy * stars_fade;
#endif

	return sky;
}

//----------------------------------------------------------------------------//
#elif defined WORLD_SPACE

vec3 draw_space_moon(vec3 ray_dir, vec3 color) {
	float nu = dot(ray_dir, moon_dir);

	// Limb darkening model from http://www.physics.hmc.edu/faculty/esin/a101/limbdarkening.pdf
	float center_to_edge = max0(sun_angular_radius - fast_acos(nu));
	float limb_darkening = pow(1.0 - sqr(1.0 - center_to_edge), 0.25);

	return color * step(0.0, center_to_edge) * limb_darkening;
}

vec3 draw_sky(vec3 ray_dir) {
	vec3 sky = vec3(0.0);

	// Sun and stars

#if defined PROGRAM_DEFERRED4
	// Output of skytextured
	sky += texelFetch(colortex0, ivec2(gl_FragCoord.xy), 0).rgb;

#ifdef STARS
	// Stars
	sky += draw_stars(ray_dir);
#endif

#ifndef VANILLA_SUN
	// Sun
	sky += draw_sun(ray_dir);
#endif

	return sky;
}

#endif

#endif // INCLUDE_SKY_SKY
