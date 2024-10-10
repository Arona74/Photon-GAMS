// Procedural Nebula Generation for Minecraft Photon Shader

// Customizable settings
#ifndef NEBULA_ENABLED
	#define NEBULA_ENABLED // [on off]
#endif

#ifndef NEBULA_DENSITY
	#define NEBULA_DENSITY 0.50 // Default density, can be adjusted
#endif

#ifndef NEBULA_SCALE
	#define NEBULA_SCALE 0.10 // Default scale, can be adjusted
#endif

#ifndef NEBULA_BRIGHTNESS
	#define NEBULA_BRIGHTNESS 2.0 // Default brightness, can be adjusted
#endif

#ifndef NEBULA_COLOR1
	#define NEBULA_COLOR1 vec3(NEBULA_COLOR1_R, NEBULA_COLOR1_G, NEBULA_COLOR1_B) // Default color 1
#endif

#ifndef NEBULA_COLOR2
	#define NEBULA_COLOR2 vec3(NEBULA_COLOR2_R, NEBULA_COLOR2_G, NEBULA_COLOR2_B) // Default color 2
#endif

#ifndef NEBULA_COLOR3
	#define NEBULA_COLOR3 vec3(NEBULA_COLOR3_R, NEBULA_COLOR3_G, NEBULA_COLOR3_B) // Default color 3
#endif

#ifndef NEBULA_MOVEMENT_SPEED
	#define NEBULA_MOVEMENT_SPEED 0.05 // Default movement speed, can be adjusted
#endif

// Noise functions
float hash(float n) { return fract(sin(n) * 43758.5453123); }

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(
        mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
            mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
        mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
            mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z
    );
}

float fbm(vec3 x) {
    float v = 0.0;
    float a = 0.5;
    vec3 shift = vec3(100);
    for (int i = 0; i < 5; ++i) {
        v += a * noise(x);
        x = x * 2.0 + shift;
        a *= 0.5;
    }
    return v;
}

vec3 nebula_color(vec3 ray_dir, float time) {
    // Add time-based movement to the nebula
    vec3 p = ray_dir * 100.0 + vec3(time * NEBULA_MOVEMENT_SPEED);
    float density = fbm(p * NEBULA_SCALE) * NEBULA_DENSITY;
    
    // Create color variations using customizable colors
    vec3 nebula = mix(NEBULA_COLOR1, NEBULA_COLOR2, fbm(p * NEBULA_SCALE * 2.0));
    nebula = mix(nebula, NEBULA_COLOR3, fbm(p * NEBULA_SCALE * 3.0));
    
    // Add brightness variations
    float brightness = fbm(p * NEBULA_SCALE * 5.0) * NEBULA_BRIGHTNESS;
    nebula *= brightness;
    
    // Apply density
    nebula *= smoothstep(0.1, 0.6, density);
    
    return nebula;
}

vec3 draw_nebula(vec3 ray_dir, vec3 background) {
#ifdef NEBULA_ENABLED
#ifdef WORLD_OVERWORLD
    // Smooth transition from day to night
    float sunset_factor = smoothstep(0.1, -0.1, sun_dir.y);
    float night_factor = smoothstep(0.0, -0.1, sun_dir.y);
    float nebula_visibility = max(sunset_factor, night_factor);
    
    // Check if the ray direction is above the horizon
    float above_horizon = smoothstep(-0.01, 0.05, ray_dir.y);
    
    if (nebula_visibility > 0.0 && above_horizon > 0.0) {
        // Use frameTimeCounter for continuous movement
        float time = frameTimeCounter;
        vec3 nebula = nebula_color(ray_dir, time);
        
        // Blend nebula with background
        float blend_factor = smoothstep(0.0, 0.8, length(nebula)) * nebula_visibility * above_horizon;
        return mix(background, nebula, blend_factor);
    }
#endif
#endif

    return background;
}