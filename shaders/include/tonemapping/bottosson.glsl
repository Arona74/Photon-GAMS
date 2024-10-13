#ifndef INCLUDE_TONEMAPPING_BOTTOSSON
#define INCLUDE_TONEMAPPING_BOTTOSSON

// Source: https://www.shadertoy.com/view/fsXcz4

// Copyright(c) 2022 Björn Ottosson
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this softwareand associated documentation files(the "Software"), to deal in
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and /or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions :
// The above copyright noticeand this permission notice shall be included in all
// copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

const float softness_scale = BOTTOSSON_SOFTNESS;  // controls softness of RGB clipping
const float offset         = BOTTOSSON_LIGHT_SAT; // controls how colors desaturate as they brighten. 0 results in that colors never fluoresce, 1 in very saturated colors
const float chroma_scale   = BOTTOSSON_CHROMA;  // overall scale of chroma

// Origin: https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
// Using this since it was easy to differentiate, same technique would work for any curve
vec3 s_curve(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    x = max0(x);
    return clamp01((x * (a * x + b)) / (x * (c * x + d) + e));
}

// derivative of s-curve
vec3 d_s_curve(vec3 x) {
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;

    x = max0(x);
    vec3 r = (x * (c * x + d) + e);
    return (a * x * (d * x + 2.0 * e) + b * (e - c * x*x)) / (r*r);
}

vec3 btsn_tonemap_per_channel(vec3 c) {
    return s_curve(c);
}

vec2 findCenterAndPurity(vec3 x) {
    // Matrix derived for (c_smooth+s_smooth) to be an approximation of the macadam limit
    // this makes it some kind of g0-like estimate
    mat3 M = mat3(
        2.26775149, -1.43293879,  0.1651873,
        -0.98535505,  2.1260072, -0.14065215,
        -0.02501605, -0.26349465,  1.2885107
    );

    x = x * M;

    float x_min = min_of(x);
    float x_max = max_of(x);

    float c = (x_max + x_min) * 0.5;
    float s = (x_max - x_min);

    // math trickery to create values close to c and s, but without producing hard edges
    vec3 y = (x - c) / s;
    float c_smooth = c + dot(y*y*y, vec3(1.0/3.0)) * s;
    float s_smooth = sqrt(dot(x - c_smooth,x - c_smooth) / 2.0);
    return vec2(c_smooth, s_smooth);
}

vec3 toLms(vec3 c) {
    mat3 rgbToLms = mat3(
        0.4122214708, 0.5363325363, 0.0514459929,
        0.2119034982, 0.6806995451, 0.1073969566,
        0.0883024619, 0.2817188376, 0.6299787005
    );

    vec3 lms_ = c * rgbToLms;
    return sign(lms_) * pow(abs(lms_), vec3(1.0/3.0));
}

float calculateC(vec3 lms) {
    // Most of this could be precomputed
    // Creating a transform that maps R,G,B in the target gamut to have same distance from grey axis

    vec3 lmsR = toLms(vec3(1.0,0.0,0.0));
    vec3 lmsG = toLms(vec3(0.0,1.0,0.0));
    vec3 lmsB = toLms(vec3(0.0,0.0,1.0));

    vec3 uDir = (lmsR - lmsG) / sqrt(2.0);
    vec3 vDir = (lmsR + lmsG - 2.0 * lmsB) / sqrt(6.0);

    mat3 to_uv = inverse(mat3(
        1.0, uDir.x, vDir.x,
        1.0, uDir.y, vDir.y,
        1.0, uDir.z, vDir.z
    ));

    vec3 _uv = lms * to_uv;

    return sqrt(sqr(_uv.y) + sqr(_uv.z));

    float a = 1.9779984951f * lms.x - 2.4285922050f * lms.y + 0.4505937099f * lms.z;
    float b = 0.0259040371f * lms.x + 0.7827717662f * lms.y - 0.8086757660f * lms.z;

    return sqrt(a*a + b*b);
}

vec2 calculateMC(vec3 c) {
    vec3 lms = toLms(c);
    float M = findCenterAndPurity(lms).x;
    return vec2(M, calculateC(lms));
}

vec2 expandShape(vec3 rgb, vec2 ST) {
    vec2 MC = calculateMC(rgb);
    vec2 STnew = vec2(MC.x / MC.y, (1.0 - MC.x) / MC.y);
    STnew = (STnew + 3.0 * sqr(STnew) * MC.y);

    return vec2(min(ST.x, STnew.x), min(ST.y, STnew.y));
}

float expandScale(vec3 rgb, vec2 ST, float scale) {
    vec2 MC = calculateMC(rgb);
    float Cnew = (rcp((ST.x / MC.x) + (ST.y / (1.0 - MC.x))));

    return max(MC.y/Cnew, scale);
}

vec2 approximateShape() {
    float m = -softness_scale * 0.2;
    float s = 1.0; //+ (softness_scale * 0.2 + softness_scale * 0.8);

    vec2 ST = vec2(1000.0, 1000.0);
    ST = expandShape(m + s * vec3(1.0, 0.0, 0.0), ST);
    ST = expandShape(m + s * vec3(1.0, 1.0, 0.0), ST);
    ST = expandShape(m + s * vec3(0.0, 1.0, 0.0), ST);
    ST = expandShape(m + s * vec3(0.0, 1.0, 1.0), ST);
    ST = expandShape(m + s * vec3(0.0, 0.0, 1.0), ST);
    ST = expandShape(m + s * vec3(1.0, 0.0, 1.0), ST);

    float scale = 0.0;
    scale = expandScale(m + s * vec3(1.0, 0.0, 0.0), ST, scale);
    scale = expandScale(m + s * vec3(1.0, 1.0, 0.0), ST, scale);
    scale = expandScale(m + s * vec3(0.0, 1.0, 0.0), ST, scale);
    scale = expandScale(m + s * vec3(0.0, 1.0, 1.0), ST, scale);
    scale = expandScale(m + s * vec3(0.0, 0.0, 1.0), ST, scale);
    scale = expandScale(m + s * vec3(1.0, 0.0, 1.0), ST, scale);

    return ST/scale;
}

vec3 btsn_tonemap_hue_preserving(vec3 c) {
    mat3 toLms = mat3(
        0.4122214708, 0.5363325363, 0.0514459929,
        0.2119034982, 0.6806995451, 0.1073969566,
        0.0883024619, 0.2817188376, 0.6299787005
    );

    mat3 fromLms = mat3(
        +4.0767416621f , -3.3077115913, +0.2309699292,
        -1.2684380046f , +2.6097574011, -0.3413193965,
        -0.0041960863f , -0.7034186147, +1.7076147010
    );

    vec3 lms_ = c * toLms;
    vec3 lms = sign(lms_) * pow(abs(lms_), vec3(1.0/3.0));

    vec2 MP = findCenterAndPurity(lms);

    // apply tone curve

    // Approach 1: scale chroma based on derivative of chrome curve
#if BOTTOSSON_APPR == BOTTOSSON_APPR_CHROMA_SAT || BOTTOSSON_APPR == BOTTOSSON_APPR_CHROMA_LUM
    {
#if BOTTOSSON_APPR == BOTTOSSON_APPR_CHROMA_SAT
        float I = (MP.x+(1.0-offset)*MP.y);
#elif BOTTOSSON_APPR == BOTTOSSON_APPR_CHROMA_LUM
        // Remove comment to see what the results are with Oklab L
        float I = dot(lms, vec3(0.2104542553f, 0.7936177850f, - 0.0040720468f));
#endif

        lms = lms*I*I;
        I = I*I*I;
        vec3 dLms = lms - I;

        float Icurve = s_curve(vec3(I)).x;
        lms = 1.0f + chroma_scale * dLms * d_s_curve(vec3(I)) / Icurve;
        I = pow(Icurve, 1.0/3.0);

        lms = lms*I;
    }
#endif

    // Approach 2: Separate color into a whiteness/blackness part, apply scale to them independendtly
#if BOTTOSSON_APPR == BOTTOSSON_APPR_SEPARATE_BW
    {
        lms = chroma_scale*(lms - MP.x) + MP.x;

        float invBlackness = (MP.x + MP.y);
        float whiteness = (MP.x - MP.y);

        float invBlacknessC = pow(s_curve(vec3(cube(invBlackness))).x, 1.0/3.0);
        float whitenessC = pow(s_curve(vec3(cube(whiteness))).x, 1.0/3.0);

        lms = (invBlacknessC+whitenessC)/2.0 + (lms-(invBlackness+whiteness)/2.0)*(invBlacknessC-whitenessC)/(invBlackness-whiteness);
    }
#endif

    // compress to a smooth approximation of the target gamut
#ifdef BOTTOSSON_COMPRESS
    {
        float M = findCenterAndPurity(lms).x;
        vec2 ST = approximateShape(); // this can be precomputed, only depends on RGB gamut
        float C_smooth_gamut = (1.0)/((ST.x/(M)) + (ST.y/(1.0-M)));
        float C = calculateC(lms);

        lms = (lms-M)/sqrt(C*C/C_smooth_gamut/C_smooth_gamut+1.0) + M;
    }
#endif

    vec3 rgb = pow(lms, vec3(3.0)) * fromLms;

    return rgb;
}

vec3 btsn_soft_saturate(vec3 x, vec3 a) {
    a = clamp(a, 0.0, softness_scale);
    a = 1.0 + a;
    x = min(x, a);
    vec3 b = (a - 1.0) * sqrt(a / (2.0 - a));
    return 1.0 - (sqrt(sqr(x - a) + b*b) - b) / (sqrt(a*a + b*b) - b);
}

vec3 btsn_soft_clip_color(vec3 color) {
    // soft clip of rgb values to avoid artifacts of hard clipping
    // causes hues distortions, but is a smooth mapping
    // not quite sure this mapping is easy to invert, but should be possible to construct similar ones that do

    float grey = 0.2;

    vec3 x = color - grey;

    vec3 xsgn = sign(x);
    vec3 xscale = 0.5 + xsgn * (0.5 - grey);
    x /= xscale;

    float maxRGB = max_of(color);
    float minRGB = min_of(color);

    float softness_0 = maxRGB / (1.0 + softness_scale) * softness_scale;
    float softness_1 = (1.0 - minRGB) / (1.0 + softness_scale) * softness_scale;

    vec3 softness = vec3(0.5) * (softness_0 + softness_1 + xsgn * (softness_1 - softness_0));

    return grey + xscale * xsgn * btsn_soft_saturate(abs(x), softness);
}

vec3 bottosson_color_test(vec3 color, vec2 uv, float t, float b, uint n) {
    float n1 = float(n + 1u);
    float adjust = 1.0 / ((n1 * (n1 + 1.0) / 2.0) / 45.0);

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= i; j++) {
            float fi = float(i) / float(n);
            float fj = float(j) / float(n);

            float x = 1.5 * (fi - 0.5 * fj - 0.5);
            float y = 1.5 * sqrt(3.0) / 2.0 * (fj - 0.5);
            vec2 xyd = vec2(x, y) - uv;

            float d = min1(0.0006 / dot(xyd, xyd));

            vec3 c = vec3(fi-fj, fj, 1.0 - fi);

            color += pow(2.0,-3.0 * cos(tau * t)) * d * c * adjust * b;
        }
    }
	return color;
}

#endif  // INCLUDE_TONEMAPPING_BOTTOSSON
