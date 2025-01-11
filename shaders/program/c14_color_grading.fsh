/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/c14_color_grading:
  Apply bloom, color grading and tone mapping then convert to rec. 709

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"


//----------------------------------------------------------------------------//
#ifdef vsh

out vec2 uv;

void main() {
	uv = gl_MultiTexCoord0.xy;

	gl_Position = vec4(gl_Vertex.xy * 2.0 - 1.0, 0.0, 1.0);
}

#endif
//----------------------------------------------------------------------------//



//----------------------------------------------------------------------------//
#ifdef fsh

layout (location = 0) out vec3 scene_color;

/* RENDERTARGETS: 0 */

in vec2 uv;

// ------------
//   Uniforms
// ------------

uniform sampler2D colortex0; // bloom tiles
uniform sampler2D colortex3; // fog transmittance
uniform sampler2D colortex5; // scene color

uniform float aspectRatio;
uniform float blindness;
uniform float darknessFactor;
uniform float frameTimeCounter;

uniform float biome_cave;
uniform float time_noon;
uniform float eye_skylight;

uniform vec2 view_pixel_size;

#include "/include/tonemapping/aces/aces.glsl"
#include "/include/tonemapping/agx.glsl"
#include "/include/tonemapping/bottosson.glsl"
#include "/include/tonemapping/zcam_justjohn.glsl"
//#include "/include/tonemapping/hatchling_oklab.glsl"

#if (tonemap == tonemap_opendrt) || (defined(TONEMAP_COMPARISON) && (tonemap_left == tonemap_opendrt || tonemap_right == tonemap_opendrt))
#undef tonemap_opendrt
//#include "/include/tonemapping/opendt/opendrt.glsl"
#endif

#if (tonemap == tonemap_jzdt) || (defined(TONEMAP_COMPARISON) && (tonemap_left == tonemap_jzdt || tonemap_right == tonemap_jzdt))
#undef tonemap_jzdt
//#include "/include/tonemapping/opendt/jzdt.glsl"
#endif

#if (tonemap == tonemap_rgbdrt) || (defined(TONEMAP_COMPARISON) && (tonemap_left == tonemap_rgbdrt || tonemap_right == tonemap_rgbdrt))
#undef tonemap_rgbdrt
#include "/include/tonemapping/rgbdrt.glsl"
#endif

#if (tonemap == tonemap_zcam_drt) || (defined(TONEMAP_COMPARISON) && (tonemap_left == tonemap_zcam_drt || tonemap_right == tonemap_zcam_drt))
#undef tonemap_zcam_drt
#include "/include/tonemapping/zcam_drt.glsl"
#endif

#if (tonemap == tonemap_fidelityfx_lpm) || (defined(TONEMAP_COMPARISON) && (tonemap_left == tonemap_fidelityfx_lpm || tonemap_right == tonemap_fidelityfx_lpm))
#undef tonemap_fidelityfx_lpm
//#include "/include/tonemapping/fidelityfx_lpm.glsl"
#endif

#ifdef TONEMAP_HUE_FADE
#include "/include/tonemapping/hue_fade.glsl"
#endif

#include "/include/utility/bicubic.glsl"
#include "/include/utility/color.glsl"

// Bloom

vec3 get_bloom(out vec3 fog_bloom) {
	const int tile_count = 6;
	const float radius  = 1.0;

	vec3 tile_sum = vec3(0.0);

	float weight = 1.0;
	float weight_sum = 0.0;

#if defined (BLOOMY_FOG) || defined (BLOOMY_RAIN)
	const float fog_bloom_radius = 1.5;

	fog_bloom = vec3(0.0); // large-scale bloom for bloomy fog
	float fog_bloom_weight = 1.0;
	float fog_bloom_weight_sum = 0.0;
#endif

	for (int i = 0; i < tile_count; ++i) {
		float a = exp2(float(-i));

		float tile_scale = 0.5 * a;
		vec2 tile_offset = vec2(1.0 - a, float(i & 1) * (1.0 - 0.5 * a));

		vec2 tile_coord = uv * tile_scale + tile_offset;

		vec3 tile = bicubic_filter(colortex0, tile_coord).rgb;

		tile_sum += tile * weight;
		weight_sum += weight;

		weight *= radius;

#if defined (BLOOMY_FOG) || defined (BLOOMY_RAIN)
		fog_bloom += tile * fog_bloom_weight;

		fog_bloom_weight_sum += fog_bloom_weight;
		fog_bloom_weight *= fog_bloom_radius;
#endif
	}

#if defined (BLOOMY_FOG) || defined (BLOOMY_RAIN)
	fog_bloom /= fog_bloom_weight_sum;
#endif

	return tile_sum / weight_sum;
}

// Color grading

vec3 gain(vec3 x, float k) {
    vec3 a = 0.5 * pow(2.0 * mix(x, 1.0 - x, step(0.5, x)), vec3(k));
    return mix(a, 1.0 - a, step(0.5, x));
}

// Color grading applied before tone mapping
// rgb := color in acescg [0, inf]
vec3 grade_input(vec3 rgb) {
	const float brightness = 0.83 * GRADE_BRIGHTNESS;
	const float contrast   = 1.00 * GRADE_CONTRAST;
	const float saturation = 0.98 * GRADE_SATURATION;
	const float vibrance   = GRADE_VIBRANCE;

	// Brightness
	rgb *= brightness;

	// Contrast
	const float log_midpoint = log2(0.18);
	rgb = log2(rgb + eps);
	rgb = contrast * (rgb - log_midpoint) + log_midpoint;
	rgb = max0(exp2(rgb) - eps);

	// Saturation
	float lum = dot(rgb, luminance_weights);
	rgb = max0(mix(vec3(lum), rgb, saturation));

	// Vibrance
	float max_c = max_of(rgb);
	float min_c = min_of(rgb);
	float sat = (max_c - min_c);
	rgb = max0(mix(rgb, mix(vec3(lum), rgb, vibrance), clamp01(1.0 - sat)));

	//vec3 hsl = rgb_to_hsl(rgb);
	//hsl.y = mix(hsl.y, hsl.y * vibrance, 1.0 - hsl.y);
	////hsl.y *= 1.0 + (vibrance - 1.0) * (1.0 - hsl.y);
	//rgb = hsl_to_rgb(hsl);


#if GRADE_WHITE_BALANCE != 6500
	// White balance (slow)
	vec3 src_xyz = blackbody(float(GRADE_WHITE_BALANCE)) * rec2020_to_xyz;
	vec3 dst_xyz = blackbody(                    6500.0) * rec2020_to_xyz;
	mat3 cat = get_chromatic_adaptation_matrix(src_xyz, dst_xyz);

	rgb = rgb * rec2020_to_xyz;
	rgb = rgb * cat;
	rgb = rgb * xyz_to_rec2020;

	rgb = max0(rgb);
#endif

	return rgb;
}

// Color grading applied after tone mapping
// rgb := color in linear rec.709 [0, 1]
vec3 grade_output(vec3 rgb) {
	// Convert to roughly perceptual RGB for color grading
	rgb = sqrt(rgb);

	// HSL color grading inspired by Tech's color grading setup in Lux Shaders

	const float orange_sat_boost = GRADE_ORANGE_SAT_BOOST;
	const float teal_sat_boost   = GRADE_TEAL_SAT_BOOST;
	const float green_sat_boost  = GRADE_GREEN_SAT_BOOST;
	const float green_hue_shift  = GRADE_GREEN_HUE_SHIFT / 360.0;

	vec3 hsl = rgb_to_hsl(rgb);

	// Oranges
	float orange = isolate_hue(hsl, 30.0, 20.0); //isolate_hue(hsl, 30.0, 20.0) // custom : 20.0, 30.0
	hsl.y *= 1.0 + orange_sat_boost * orange;

	// Teals
	float teal = isolate_hue(hsl, 210.0, 20.0);
	hsl.y *= 1.0 + teal_sat_boost * teal;

	// Greens
	float green = isolate_hue(hsl, 90.0, 53.0); //isolate_hue(hsl, 90.0, 44.0) // custom : 90.0, 53.0
	hsl.x += green_hue_shift * green;
	hsl.y *= 1.0 + green_sat_boost * green;

	rgb = hsl_to_rgb(hsl);

	rgb = gain(rgb, 1.05 * GRADE_GAIN);

	return sqr(rgb);
}

// Tonemapping operators

vec3 tonemap_none(vec3 rgb) { return rgb; }

// ACES RRT and ODT
vec3 academy_rrt(vec3 rgb) {
	rgb *= 1.6; // Match the exposure to the RRT

	rgb = rgb * rec709_to_ap0;

	rgb = aces_lmt(rgb);
	rgb = aces_rrt(rgb);
	rgb = aces_odt(rgb);

	return rgb * ap1_to_rec709;
}

// ACES RRT and ODT approximation
vec3 academy_fit(vec3 rgb) {
	rgb *= 1.6; // Match the exposure to the RRT

	rgb = rgb * rec709_to_ap0;

	rgb = aces_lmt(rgb);
	rgb = rrt_sweeteners(rgb);
	rgb = rrt_and_odt_fit(rgb);

	// Global desaturation
	vec3 grayscale = vec3(dot(rgb, luminance_weights));
	rgb = mix(grayscale, rgb, odt_sat_factor);

	return rgb * ap1_to_rec709;
}

vec3 tonemap_hejl_2015(vec3 rgb) {
	const float white_point = 5.0;

	vec4 vh = vec4(rgb, white_point);
	vec4 va = (1.425 * vh) + 0.05; // eval filmic curve
	vec4 vf = ((vh * va + 0.004) / ((vh * (va + 0.55) + 0.0491))) - 0.0821;

	return vf.rgb / vf.www; // white point correction
}

// Filmic tonemapping operator made by Jim Hejl and Richard Burgess
// Modified by Tech to not lose color information below 0.004
vec3 tonemap_hejl_burgess(vec3 rgb) {
	rgb = rgb * min(vec3(1.0), 1.0 - 0.8 * exp(rcp(-0.004) * rgb));
	rgb = (rgb * (6.2 * rgb + 0.5)) / (rgb * (6.2 * rgb + 1.7) + 0.06);
	return srgb_eotf_inv(rgb); // Revert built-in sRGB conversion
}

// Timothy Lottes 2016, "Advanced Techniques and Optimization of HDR Color Pipelines"
// https://gpuopen.com/wp-content/uploads/2016/03/GdcVdrLottes.pdf
vec3 tonemap_lottes(vec3 rgb) {
	const vec3 a       = vec3(LOTTES_CONTRAST); // Contrast
	const vec3 d       = vec3(LOTTES_SHOULDER); // Shoulder contrast
	const vec3 hdr_max = vec3(8.0);  // White point
	const vec3 mid_in  = vec3(0.26); // Fixed midpoint x
	const vec3 mid_out = vec3(0.32); // Fixed midput y

	const vec3 b =
		(-pow(mid_in, a) + pow(hdr_max, a) * mid_out) /
		((pow(hdr_max, a * d) - pow(mid_in, a * d)) * mid_out);
	const vec3 c =
		(pow(hdr_max, a * d) * pow(mid_in, a) - pow(hdr_max, a) * pow(mid_in, a * d) * mid_out) /
		((pow(hdr_max, a * d) - pow(mid_in, a * d)) * mid_out);

	return pow(rgb, a) / (pow(rgb, a * d) * b + c);
}

// Filmic tonemapping operator made by John Hable for Uncharted 2
vec3 tonemap_uncharted_2_partial(vec3 rgb) {
	const float a = UNCHARTED_2_SHOULDER_STRENGTH; // Shoulder strength
	const float b = UNCHARTED_2_LINEAR_STRENGTH; // Linear strength
	const float c = UNCHARTED_2_LINEAR_ANGLE; // Linear angle
	const float d = UNCHARTED_2_TOE_STRENGTH; // Toe strength
	const float e = UNCHARTED_2_TOE_NUMERATOR; // Toe numerator
	const float f = UNCHARTED_2_TOE_DENOMINATOR; // Toe denominator

	return ((rgb * (a * rgb + (c * b)) + (d * e)) / (rgb * (a * rgb + b) + d * f)) - e / f;
}

vec3 tonemap_uncharted_2_filmic(vec3 rgb) {
	const float exposure_bias = 2.0;
	const vec3 w = vec3(11.2);

	vec3 curr = tonemap_uncharted_2_partial(rgb * exposure_bias);
	vec3 white_scale = vec3(1.0) / tonemap_uncharted_2_partial(w);
	return curr * white_scale;
}

vec3 tonemap_uncharted_2(vec3 rgb) {
#ifdef UNCHARTED_2_PARTIAL
	rgb *= 3.0;
	return tonemap_uncharted_2_partial(rgb);
#else
	return tonemap_uncharted_2_filmic(rgb);
#endif
}

// Tone mapping operator made by Tech for his shader pack Lux
vec3 tonemap_tech(vec3 rgb) {
	vec3 a = rgb * min(vec3(1.0), 1.0 - exp(-1.0 / 0.038 * rgb));
	a = mix(a, rgb, rgb * rgb);
	return a / (a + 0.6);
}

// Tonemapping operator made by Zombye for his old shader pack Ozius
// It was given to me by Jessie
vec3 tonemap_ozius(vec3 rgb) {
    const vec3 a = vec3(0.46, 0.46, 0.46);
    const vec3 b = vec3(0.60, 0.60, 0.60);

	rgb *= 1.6;

    vec3 cr = mix(vec3(dot(rgb, luminance_weights_ap1)), rgb, 0.5) + 1.0;

    rgb = pow(rgb / (1.0 + rgb), a);
    return pow(rgb * rgb * (-2.0 * rgb + 3.0), cr / b);
}

vec3 tonemap_reinhard(vec3 rgb) {
	return rgb / (rgb + 1.0);
}

vec3 tonemap_reinhard_jodie(vec3 rgb) {
	vec3 reinhard = rgb / (rgb + 1.0);
	return mix(rgb / (dot(rgb, luminance_weights) + 1.0), reinhard, reinhard);
}

// Uchimura 2017, "HDR theory and practice"
// Math: https://www.desmos.com/calculator/gslcdxvipg
// Source: https://www.slideshare.net/nikuque/hdr-theory-and-practicce-jp
vec3 tonemap_uchimura(vec3 rgb) {
	const float P = UCHIMURA_MAX_BRIGHTNESS;  // max display brightness
	const float a = UCHIMURA_CONTRAST;  // contrast
	const float m = UCHIMURA_LINEAR_SECTION_START; // linear section start
	const float l = UCHIMURA_LINEAR_SECTION_LENGTH;  // linear section length
	const float c = UCHIMURA_BLACK_TIGHTNESS; // black
	const float b = UCHIMURA_BLACK_PEDESTAL;  // pedestal

	float l0 = ((P - m) * l) / a;
	float L0 = m - m / a;
	float L1 = m + (1.0 - m) / a;
	float S0 = m + l0;
	float S1 = m + a * l0;
	float C2 = (a * P) / (P - S1);
	float CP = -C2 / P;

	vec3 w0 = vec3(1.0 - smoothstep(0.0, m, rgb));
	vec3 w2 = vec3(step(m + l0, rgb));
	vec3 w1 = vec3(1.0 - w0 - w2);

	vec3 T = vec3(m * pow(rgb / m, vec3(c)) + b);
	vec3 S = vec3(P - (P - S1) * exp(CP * (rgb - S0)));
	vec3 L = vec3(m + a * (rgb - m));

	return T * w0 + L * w1 + S * w2;
}

vec3 tonemap_justjohn(vec3 rgb) {
	rgb *= 1.6;
#ifdef JJS_ZCAM_REC2020
	rgb = zcam_tonemap_rec2020(rgb);
#else
	vec3 sRGB = rgb * working_to_display_color;
	sRGB = zcam_tonemap(sRGB);
	//sRGB = zcam_gamma_correct(sRGB);
	rgb = sRGB * display_to_working_color;
#endif
	return rgb;
}

// Minimal implementation of Troy Sobotka's AgX display transform by bwrensch
// Source: https://www.shadertoy.com/view/cd3XWr
//         https://iolite-engine.com/blog_posts/minimal_agx_implementation
// Original: https://github.com/sobotka/AgX
vec3 tonemap_agx(vec3 rgb) {
	//rgb = srgb_eotf(rgb);

	rgb = agx_pre(rgb);

	// Apply sigmoid function approximation
	rgb = agx_default_contrast_approx(rgb);
#if AGX_LOOK != 0
	rgb = agx_look(rgb);
#endif
#ifdef AGX_EOTF
	rgb = agx_eotf(rgb);
#endif

	return srgb_eotf_inv(rgb);
}

// Source: https://modelviewer.dev/examples/tone-mapping#commerce
vec3 tonemap_khronos_neutral(vec3 rgb) {
	const float start_compression = KHRONOS_NEUTRAL_SCMP - 0.04; // 0.8
	const float desaturation = KHRONOS_NEUTRAL_DESAT; // 0.15

	rgb *= 1.2;

	float x = min_of(rgb);
	float offset = x < 0.08 ? x - 6.25 * sqr(x) : 0.04;
	rgb -= offset;

	float peak = max_of(rgb);
	if(peak < start_compression) return rgb;

	float d = 1.0 - start_compression;
	float new_peak = 1.0 - sqr(d) / (peak + d - start_compression);
	rgb *= new_peak / peak;

	float g = 1.0 - rcp(desaturation * (peak - new_peak) + 1.0);
	return mix(rgb, vec3(1.0), g);
}

vec3 tonemap_bottosson(vec3 rgb) {
	rgb *= 0.8;
	vec3 sRGB = rgb * working_to_display_color;
	sRGB = btsn_tonemap_hue_preserving(sRGB);
#ifdef BOTTOSSON_TONEMAP_CURVE
	sRGB = btsn_tonemap_per_channel(sRGB);
#endif
#ifdef BOTTOSSON_SOFT_CLIP
	sRGB = btsn_soft_clip_color(sRGB);
#endif
	rgb = sRGB * display_to_working_color;
	return rgb;
}

// Cinematic Tonemapper designed to fix issues with ACES
// Source: https://github.com/ltmx/Melon-Tonemapper/tree/main
vec3 melon_hueshift(vec3 col) {
	float A = max(col.r, col.g);
	return vec3(A, max(A, col.b), col.b);
}

vec3 tonemap_melon(vec3 rgb) {
	rgb *= 1.8;

	// remaps the colors to [0-1] range
	// tested to be as close to ACES contrast levels as possible
	rgb = pow(rgb, vec3(1.56));
	rgb = rgb / (rgb + 0.84);

	// governs the transition to white for high color intensities
	float factor = max_of(rgb) * MELON_TONEMAP_SHIFT_STRENGTH; // multiply by 0.15 to get a similar look to ACES
	factor = factor / (factor + 1); // remaps the factor to [0-1] range
	factor *= factor; // smooths the transition to white

#ifdef MELON_TONEMAP_HUESHIFT
	// shift the hue for high intensities (for a more pleasing look).
	rgb = mix(rgb, melon_hueshift(rgb), factor); // can be removed for more neutral colors
#endif
#ifdef MELON_TONEMAP_WHITESHIFT
	rgb = mix(rgb, vec3(1.0), factor); // shift to white for high intensities
#endif

	return clamp01(rgb);
}

#ifndef tonemap_opendrt
vec3 tonemap_opendrt(vec3 rgb) {
	rgb *= OPENDRT_EXPOSURE;
	//rgb = opendrtransform(rgb);
	return rgb;
}
#endif

#ifndef tonemap_jzdt
vec3 tonemap_jzdt(vec3 rgb) {
	rgb *= 1.2;
	//rgb = jzdtransform(rgb);
	return rgb;
}
#endif

#ifndef tonemap_rgbdrt
vec3 tonemap_rgbdrt(vec3 rgb) {
	rgb *= 1.4;
	rgb = rgbdrtransform(rgb);
	return rgb;
}
#endif

#ifndef tonemap_zcam_drt
vec3 tonemap_zcam_drt(vec3 rgb) {
	rgb *= 1.6;
	rgb = zcamdrtransform(rgb);
	return rgb;
}
#endif

#ifndef tonemap_fidelityfx_lpm
// AMD FidelityFX Luma Preserving Mapper
// Site: https://gpuopen.com/fidelityfx-lpm/
// Source: https://github.com/GPUOpen-Effects/FidelityFX-LPM
vec3 tonemap_fidelityfx_lpm(vec3 rgb) {
	//rgb = fidelityfx_lpm(rgb);
	return rgb;
}
#endif

vec3 tonemap_hatchling(vec3 rgb) {
	//rgb *= 1.6;
	/*vec3 lin_rgb = (rgb * working_to_display_color);
	lin_rgb = hatchling_tonemap(lin_rgb);
	rgb = Linear3(lin_rgb) * display_to_working_color;*/
	return rgb;
}


float vignette(vec2 uv) {
    const float vignette_size = 16.0;
    const float vignette_intensity = 0.08 * VIGNETTE_INTENSITY;

	float darkness_pulse = 1.0 - dampen(abs(cos(2.0 * frameTimeCounter)));

    float vignette = vignette_size * (uv.x * uv.y - uv.x) * (uv.x * uv.y - uv.y);
          vignette = pow(vignette, vignette_intensity + 0.1 * biome_cave + 0.3 * blindness + 0.2 * darkness_pulse * darknessFactor);

    return vignette;
}

vec3 pre_tonemap(vec3 rgb) {
#ifdef TONEMAP_HUE_FADE
	rgb = hdr_hue_fade(rgb);
#endif
	return rgb;
}

void main() {
	ivec2 texel = ivec2(gl_FragCoord.xy);

	scene_color = texelFetch(colortex5, texel, 0).rgb;

	float exposure = texelFetch(colortex5, ivec2(0), 0).a;

#ifdef BLOOM
	vec3 fog_bloom;
	vec3 bloom = get_bloom(fog_bloom);
	float bloom_intensity = 0.12 * BLOOM_INTENSITY;

	scene_color = mix(scene_color, bloom, bloom_intensity);

#ifdef BLOOMY_FOG
	float fog_transmittance = texture(colortex3, uv * taau_render_scale).x;
	scene_color = mix(fog_bloom, scene_color, pow(fog_transmittance, BLOOMY_FOG_INTENSITY));
#endif
#endif

	scene_color *= exposure;

#ifdef VIGNETTE
	scene_color *= vignette(uv);
#endif

	scene_color = grade_input(scene_color);

	scene_color = pre_tonemap(scene_color);

/* "/include/tonemapping/zcam_justjohn.glsl" */
#if TONEMAP_COLORTEST == COLORTEST_JJS_ZCAM
	vec2 position = vec2(uv.x - frameTimeCounter * 0.2, uv.y);
	vec3 ICh = vec3(exp(position.y * 3.0) - 1.0, 0.07, position.x * 5.0);
    vec3 sRGB = max(vec3(0.0), XYZ_to_sRGB * ICh_to_XYZ(ICh));
	scene_color = sRGB * display_to_working_color;
#elif TONEMAP_COLORTEST == COLORTEST_BOTTOSSON
#define COLORTEST_BOTTOSSON_ROT_SPEED 0.2
#define COLORTEST_BOTTOSSON_BRIGHT_SPEED 0.05
#define COLORTEST_BOTTOSSON_START_BRIGHT 1./sqrt(2.0)
#define COLORTEST_BOTTOSSON_AMOUNT 12
#define COLORTEST_BOTTOSSON_BRIGHTNESS 2.0

	vec2 pos = (uv * 2.0 - 1.0) * vec2(1.0 * aspectRatio, 1.0);
	pos = vec2(length(pos), atan(pos.y, pos.x) / tau);
	pos.y += mod(frameTimeCounter * COLORTEST_BOTTOSSON_ROT_SPEED, 1.0);
	pos.y *= tau;
	pos = vec2(cos(pos.y), sin(pos.y)) * pos.x;
	pos.y -= 0.22;

	float time = frameTimeCounter * COLORTEST_BOTTOSSON_BRIGHT_SPEED + COLORTEST_BOTTOSSON_START_BRIGHT * 0.5;

	scene_color = bottosson_color_test(vec3(0.0), pos, time, COLORTEST_BOTTOSSON_BRIGHTNESS, COLORTEST_BOTTOSSON_AMOUNT) * display_to_working_color;
	//scene_color = vec3(aspectRatio-length(pos), clamp01(abs(pos.x)), clamp01(pos.y));
#endif

#ifdef TONEMAP_COMPARISON
	scene_color = uv.x < 0.5 ? tonemap_left(scene_color) : tonemap_right(scene_color);
#else
	scene_color = tonemap(scene_color);
#endif

	scene_color = clamp01(scene_color * working_to_display_color);
	scene_color = grade_output(scene_color);

// Tonemap plot
#if TONEMAP_PLOT == TONEMAP_PLOT_WHITE
	const float scale = 2.0;
	vec2 uv_scaled = uv * scale * vec2(1.0, 1.0 / aspectRatio);
	float x = uv_scaled.x;
	float y = tonemap(pre_tonemap(vec3(x))).x;

	if (abs(uv_scaled.x - 1.0) < 0.001 * scale) scene_color = vec3(1.0, 0.0, 0.0);
	if (abs(uv_scaled.y - 1.0) < 0.001 * scale) scene_color = vec3(1.0, 0.0, 0.0);
	if (abs(uv_scaled.y - y) < 0.001 * scale) scene_color = vec3(1.0);

#elif TONEMAP_PLOT == TONEMAP_PLOT_COLOR
	const float scale = 2.0;
	const float line  = scale * 0.0015;
	const float line2 = scale * 0.002;
	vec2 uv_scaled = uv * scale * vec2(1.0, 1.0 / aspectRatio);
	float x = uv_scaled.x;
    mat4x3 c = mat4x3(vec3(x, 0.0, 0.0), vec3(0.0, x, 0.0), vec3(0.0, 0.0, x), vec3(x));
	mat4x3 y = mat4x3(
		tonemap(pre_tonemap(c[0] * display_to_working_color)),
		tonemap(pre_tonemap(c[1] * display_to_working_color)),
		tonemap(pre_tonemap(c[2] * display_to_working_color)),
		tonemap(pre_tonemap(c[3] * display_to_working_color))
	);

	if (abs(uv_scaled.x - 1.0) < line) scene_color = vec3(0.0);
	if (abs(uv_scaled.y - 1.0) < line) scene_color = vec3(0.0);

	/*if (abs(uv_scaled.y - y[0].r) < line2) scene_color += c[0];
	if (abs(uv_scaled.y - y[1].g) < line2) scene_color += c[1];
	if (abs(uv_scaled.y - y[2].b) < line2) scene_color += c[2];
	if (abs(uv_scaled.y - y[3].x) < line2) scene_color += c[3];*/

	if (abs(uv_scaled.y - y[3].x) < line) scene_color = y[3] * working_to_display_color;

	if (abs(uv_scaled.y - y[0].r) < line
	 || abs(uv_scaled.y - y[1].g) < line
	 || abs(uv_scaled.y - y[2].b) < line) scene_color = vec3(0.0);

	if (abs(uv_scaled.y - y[0].r) < line) scene_color += y[0] * working_to_display_color;
	if (abs(uv_scaled.y - y[1].g) < line) scene_color += y[1] * working_to_display_color;
	if (abs(uv_scaled.y - y[2].b) < line) scene_color += y[2] * working_to_display_color;

#endif
}
