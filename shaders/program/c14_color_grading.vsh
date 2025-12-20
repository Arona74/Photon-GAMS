/*
--------------------------------------------------------------------------------

  Photon Shader by SixthSurge

  program/c14_color_grading:
  Apply bloom, color grading and tone mapping then convert to rec. 709

--------------------------------------------------------------------------------
*/

#include "/include/global.glsl"

#ifdef LENS_FLARE
uniform vec3 upPosition;
uniform vec3 sunPosition;
#endif

out vec2 uv;

#ifdef LENS_FLARE
flat out vec3 upVec, sunVec;
#endif

void main() {
	uv = gl_MultiTexCoord0.xy;

#ifdef LENS_FLARE
  upVec = normalize(upPosition);
  sunVec = normalize(sunPosition);
#endif

	gl_Position = vec4(gl_Vertex.xy * 2.0 - 1.0, 0.0, 1.0);
}
