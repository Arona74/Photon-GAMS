/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/c1_blend_layers
  Apply volumetric fog

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

out vec2 uv;

flat out vec3 ambient_color;
flat out vec3 light_color;

#if defined WORLD_OVERWORLD
#include "/include/fog/overworld/coeff_struct.glsl"
flat out AirFogCoefficients air_fog_coeff;
#endif

// ------------
//   Uniforms
// ------------

uniform sampler2D colortex4; // Sky map, lighting color palette

uniform float rainStrength;
uniform float sunAngle;

uniform int worldTime;
uniform int worldDay;

uniform vec3 sun_dir;

uniform float wetness;

uniform float eye_skylight;

uniform float biome_temperate;
uniform float biome_arid;
uniform float biome_snowy;
uniform float biome_taiga;
uniform float biome_jungle;
uniform float biome_swamp;
uniform float biome_may_rain;
uniform float biome_may_snow;
uniform float biome_temperature;
uniform float biome_humidity;

uniform float time_sunrise;
uniform float time_noon;
uniform float time_sunset;
uniform float time_midnight;

uniform float desert_sandstorm;

#if defined WORLD_OVERWORLD
#include "/include/fog/overworld/coeff.glsl"
#endif

uniform vec2 view_res;

void main() {
	uv = gl_MultiTexCoord0.xy;

	int lighting_color_x = SKY_MAP_LIGHT_X;
	light_color   = texelFetch(colortex4, ivec2(lighting_color_x, 0), 0).rgb;
	ambient_color = texelFetch(colortex4, ivec2(lighting_color_x, 1), 0).rgb;

#if defined WORLD_OVERWORLD
	air_fog_coeff = calculate_air_fog_coefficients();
#endif

	vec2 vertex_pos = gl_Vertex.xy * taau_render_scale;
	gl_Position = vec4(vertex_pos * 2.0 - 1.0, 0.0, 1.0);
}