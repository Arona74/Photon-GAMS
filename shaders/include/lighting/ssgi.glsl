uniform sampler2D colortex5;

const int   SSGI_SAMPLES = 8;
const int   SSGI_STEPS   = 16;
const float SSGI_RADIUS  = 8;
const float SSGI_THICKNESS = 8.0;
const float SSGI_INTENSITY = 1.0;

vec3 screen_to_view_space(vec2 uv, float depth) {
    return screen_to_view_space(combined_projection_matrix_inverse, vec3(uv, depth), false);
}

mat3 mat3_from_axis_angle(vec3 axis, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

vec3 ssgi(vec3 screen_pos, vec3 view_pos, vec3 view_normal, vec2 dither) {
    vec3 indirect_light = vec3(0.0);
    vec3 random_vec = normalize(vec3(dither * 2.0 - 1.0, 0.0));
    mat3 rotation_matrix = mat3_from_axis_angle(view_normal, random_vec.x * pi * 2.0);
    float step_size = SSGI_RADIUS / float(SSGI_STEPS);

    for (int i = 0; i < SSGI_SAMPLES; ++i) {
        float angle = (float(i) / float(SSGI_SAMPLES)) * pi * 2.0;
        vec3 ray_dir = rotation_matrix * normalize(vec3(cos(angle), sin(angle), 1.0));
        ray_dir = reflect(ray_dir, view_normal);

        for (int j = 1; j <= SSGI_STEPS; ++j) {
            vec3 current_step_pos = view_pos + ray_dir * step_size * float(j);
            vec4 projected_pos = gbufferProjection * vec4(current_step_pos, 1.0);
            projected_pos.xyz /= projected_pos.w;
            vec2 sample_uv = projected_pos.xy * 0.5 + 0.5;

            if (sample_uv.x < 0.0 || sample_uv.x > 1.0 || sample_uv.y < 0.0 || sample_uv.y > 1.0) {
                break;
            }

            float scene_depth = texture(depthtex1, sample_uv).x;
            vec3 scene_view_pos = screen_to_view_space(sample_uv, scene_depth);

            if (current_step_pos.z < scene_view_pos.z && scene_view_pos.z - current_step_pos.z < SSGI_THICKNESS) {
                vec3 hit_color = texture(colortex5, sample_uv).rgb;

                float falloff = 1.0 - smoothstep(0.0, SSGI_RADIUS, length(current_step_pos - view_pos));
                falloff *= falloff;

                indirect_light += hit_color * falloff;
                break;
            }
        }
    }

    return indirect_light / float(SSGI_SAMPLES) * SSGI_INTENSITY;
}