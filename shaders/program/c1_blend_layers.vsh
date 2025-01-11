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

uniform sampler2D colortex4; // Sky map, lighting colors

uniform vec2 view_res;

void main() {
	uv = gl_MultiTexCoord0.xy;

	int lighting_color_x = SKY_MAP_LIGHT_X;
	light_color   = texelFetch(colortex4, ivec2(lighting_color_x, 0), 0).rgb;
	ambient_color = texelFetch(colortex4, ivec2(lighting_color_x, 1), 0).rgb;

	vec2 vertex_pos = gl_Vertex.xy * taau_render_scale;
	gl_Position = vec4(vertex_pos * 2.0 - 1.0, 0.0, 1.0);
}

