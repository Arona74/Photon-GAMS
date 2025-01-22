#if !defined INCLUDE_MISC_MATERIAL
#define INCLUDE_MISC_MATERIAL

#include "/include/tonemapping/aces/matrices.glsl"
#include "/include/utility/color.glsl"

const float air_n   = 1.000293; // for 0°C and 1 atm
const float water_n = 1.333;    // for 20°C

struct Material {
	vec3 albedo;
	vec3 emission;
	vec3 f0;
	vec3 f82; // hardcoded metals only
	float roughness;
	float sss_amount;
	float sheen_amount; // SSS "sheen" for tall grass
	float porosity;
	float ssr_multiplier;
	bool is_metal;
	bool is_hardcoded_metal;
};

const Material water_material = Material(vec3(0.0), vec3(0.0), vec3(0.02), vec3(0.0), 0.002, 1.0, 0.0, 0.0, 1.0, false, false);

#if TEXTURE_FORMAT == TEXTURE_FORMAT_LAB
void decode_specular_map(vec4 specular_map, inout Material material) {
	// f0 and f82 values for hardcoded metals from Jessie LC (https://github.com/Jessie-LC)
	const vec3[] metal_f0 = vec3[](
		vec3(0.78, 0.77, 0.74), // Iron
		vec3(1.00, 0.90, 0.61), // Gold
		vec3(1.00, 0.98, 1.00), // Aluminum
		vec3(0.77, 0.80, 0.79), // Chrome
		vec3(1.00, 0.89, 0.73), // Copper
		vec3(0.79, 0.87, 0.85), // Lead
		vec3(0.92, 0.90, 0.83), // Platinum
		vec3(1.00, 1.00, 0.91)  // Silver
	);
	const vec3[] metal_f82 = vec3[](
		vec3(0.74, 0.76, 0.76),
		vec3(1.00, 0.93, 0.73),
		vec3(0.96, 0.97, 0.98),
		vec3(0.74, 0.79, 0.78),
		vec3(1.00, 0.90, 0.80),
		vec3(0.83, 0.80, 0.83),
		vec3(0.89, 0.87, 0.81),
		vec3(1.00, 1.00, 0.95)
	);

	material.roughness = sqr(1.0 - specular_map.r);
	material.emission = max(material.emission, material.albedo * specular_map.a * float(specular_map.a != 1.0));

	if (specular_map.g < 229.5 / 255.0) {
		#ifdef F0_FALLBACK
		if (specular_map.g < 0.5 / 255.0) specular_map.g = 0.04;
		#endif

		// Dielectrics
		material.f0 = max(material.f0, specular_map.g);

		float has_sss = step(64.5 / 255.0, specular_map.b);
		material.sss_amount = max(material.sss_amount, linear_step(64.0 / 255.0, 1.0, specular_map.b * has_sss));
		material.porosity = linear_step(0.0, 64.0 / 255.0, max0(specular_map.b - specular_map.b * has_sss));
	} else if (specular_map.g < 237.5 / 255.0) {
		// Hardcoded metals
		uint metal_id = clamp(uint(255.0 * specular_map.g) - 230u, 0u, 7u);

		material.f0 = metal_f0[metal_id];
		material.f82 = metal_f82[metal_id];
		material.is_metal = true;
		material.is_hardcoded_metal = true;
	} else {
		// Albedo metal
		material.f0 = material.albedo;
		material.is_metal = true;
	}

	material.ssr_multiplier = step(0.01, (material.f0.x - material.f0.x * material.roughness * SSR_ROUGHNESS_THRESHOLD)); // based on Kneemund's method
}
#elif TEXTURE_FORMAT == TEXTURE_FORMAT_OLD
void decode_specular_map(vec4 specular_map, inout Material material) {
	material.roughness = sqr(1.0 - specular_map.r);
	material.is_metal  = specular_map.g > 0.5;
	material.f0        = material.is_metal ? material.albedo : material.f0;
	material.emission  = max(material.emission, material.albedo * specular_map.b);

	material.ssr_multiplier = step(0.01, (material.f0.x - material.f0.x * material.roughness * SSR_ROUGHNESS_THRESHOLD)); // based on Kneemund's method
}
#endif

void decode_specular_map(vec4 specular_map, inout Material material, out bool parallax_shadow) {
#if defined POM && defined POM_SHADOW
		// Specular map alpha >= 0.5 => parallax shadow
		parallax_shadow = specular_map.a >= 0.5;
		specular_map.a = fract(specular_map.a * 2.0);
#endif

		decode_specular_map(specular_map, material);
}

Material material_from(vec3 albedo_srgb, uint material_mask, vec3 world_pos, vec3 normal, inout vec2 light_levels) {
	vec3 block_pos = fract(world_pos);

	// Create material with default values

	Material material;
	material.albedo             = srgb_eotf_inv(albedo_srgb) * rec709_to_rec2020;
	material.emission           = vec3(0.0);
	material.f0                 = vec3(0.02);
	material.f82                = vec3(0.0);
	material.roughness          = 1.0;
	material.sss_amount         = 0.0;
	material.sheen_amount       = 0.0;
	material.porosity           = 0.0;
	material.ssr_multiplier     = 0.0;
	material.is_metal           = false;
	material.is_hardcoded_metal = false;

	// Hardcoded materials for specific blocks
	// Using binary split search to minimise branches per fragment (TODO: measure impact)

	vec3 hsl = rgb_to_hsl(albedo_srgb);
	vec3 albedo_sqrt = sqrt(material.albedo);

	if (material_mask < 32u) { // 0-32
		if (material_mask < 16u) { // 0-16
			if (material_mask < 8u) { // 0-8
				if (material_mask < 4u) { // 0-4
					if (material_mask < 2u) { // 0-2
						if (material_mask == 0u) { // 0
							#ifdef HARDCODED_SPECULAR
							// Default
							float smoothness = 0.33 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 1
							// Water
						}
					} else { // 2-4
						if (material_mask == 2u) { // 2
							#ifdef HARDCODED_SSS
							// Small plants
							material.sss_amount = 0.5;
							material.sheen_amount = 0.5;
							#endif
						} else { // 3
							#ifdef HARDCODED_SSS
							// Tall plants (lower half)
							material.sss_amount = 0.5;
							material.sheen_amount = 0.5;
							#endif
						}
					}
				} else { // 4-8
					if (material_mask < 6u) { // 4-6
						if (material_mask == 4u) { // 4
							#ifdef HARDCODED_SSS
							// Tall plants (upper half)
							material.sss_amount = 0.5;
							material.sheen_amount = 0.5;
							#endif
						} else { // 5
							// Leaves
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * smoothstep(0.16, 0.5, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.sheen_amount = 0.5;
							#endif

							#ifdef HARDCODED_SSS
							material.sss_amount = 1.0;
							#endif
						}
					} else { // 6-8
						if (material_mask == 6u) { // 6
							#ifdef HARDCODED_SPECULAR
							// Grass, stone, spruce and dark oak planks
							float smoothness = 0.4 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif

							#ifdef HARDCODED_POROSITY
							material.porosity = 0.2;
							#endif
						} else { // 7
							// Sand
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.8 * linear_step(0.81, 0.96, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif

							#ifdef HARDCODED_POROSITY
							material.porosity = 1.35;
							#endif
						}
					}
				}
			} else { // 8-16
				if (material_mask < 12u) { // 8-12
					if (material_mask < 10u) { // 8-10
						if (material_mask == 8u) { // 8
							// Ice
							#ifdef HARDCODED_SPECULAR
							float smoothness = pow4(linear_step(0.4, 0.8, hsl.z)) * 0.6;
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 1.0;
							#endif

							#ifdef HARDCODED_SSS
							// Strong SSS
							material.sss_amount = 0.75;
							#endif

							#ifdef HARDCODED_POROSITY
							material.porosity = 0.2;
							#endif
						} else { // 9
							// Red sand, birch planks
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * linear_step(0.61, 0.85, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif

							#ifdef HARDCODED_POROSITY
							material.porosity = 0.5;
							#endif
						}
					} else { // 10-12
						if (material_mask == 10u) { // 10
							// Oak, jungle and acacia planks, granite and diorite
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * linear_step(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif

							#ifdef HARDCODED_POROSITY
							material.porosity = 0.2;
							#endif
						} else { // 11
							//Nether bricks
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * linear_step(0.13, 0.4, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = material.albedo;
							material.ssr_multiplier = 0.04;
							#endif
						}
					}
				} else { // 12-16
					if (material_mask < 14u) { // 12-14
						if (material_mask == 12u) { // 12
							// Metals
							#ifdef HARDCODED_SPECULAR
							float smoothness = 1.13 * sin(linear_step(0.0, 1.0, hsl.z));
							material.roughness = min(sin(1.0 - smoothness), 0.65);
							material.f0 = sqr(material.albedo) * 1.2;
							material.is_metal = true;
							material.ssr_multiplier = 1.0;
							#endif
						} else { // 13
							// Gems
							#ifdef HARDCODED_SPECULAR
							float smoothness = 1.1 * sin(sqrt(linear_step(0.1, 0.9, hsl.z)));
							material.roughness = min(sin(sqrt(1.0 - smoothness)), 0.7);
							material.f0 = sqr(material.albedo) * 1.8 + vec3(0.01);
							material.is_metal = true;
							material.ssr_multiplier = 1.0;
							#endif
						}
					} else { // 14-16
						if (material_mask == 14u) { // 14
							#ifdef HARDCODED_SSS
							// Strong SSS
							material.sss_amount = 0.6;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
							#ifdef HARDCODED_POROSITY
							material.porosity = 1.5;
							#endif
						} else { // 15
							#ifdef HARDCODED_SSS
							// Weak SSS
							material.sss_amount = 0.22;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
							#ifdef HARDCODED_POROSITY
							material.porosity = 1.25;
							#endif
						}
					}
				}
			}
		} else { // 16-32
			if (material_mask < 24u) { // 16-24
				if (material_mask < 20u) { // 16-20
					if (material_mask < 18u) { // 16-18
						if (material_mask == 16u) { // 16
							#ifdef HARDCODED_EMISSION
							// Chorus plant
							material.emission  = 0.25 * albedo_sqrt * pow4(hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 17
							#ifdef HARDCODED_SPECULAR
							// End stone
							float smoothness = 0.4 * linear_step(0.61, 0.85, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 1.0;
							#endif
						}
					} else { // 18-20
						if (material_mask == 18u) { // 18
							// Metals
							#ifdef HARDCODED_SPECULAR
							float smoothness = 1.13 * sin(linear_step(0.0, 1.0, hsl.z));
							material.roughness = min(sin(1.0 - smoothness), 0.88);
							material.f0 = sqr(material.albedo);
							material.is_metal = true;
							material.ssr_multiplier = 1.0;
							#endif
						} else { // 19
							// Warped stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.yz, 1.0 - block_pos.yz),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.x))
							);
							float blue = isolate_hue(hsl, 200.0, 60.0);
							material.emission = albedo_sqrt * hsl.y * blue * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 20-24
					if (material_mask < 22u) { // 20-22
						if (material_mask == 20u) { // 20
							// Warped stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.xz, 1.0 - block_pos.xz),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.y))
							);
							float blue = isolate_hue(hsl, 200.0, 60.0);
							material.emission = albedo_sqrt * hsl.y * blue * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 21
							// Warped stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.xy, 1.0 - block_pos.xy),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.z))
							);
							float blue = isolate_hue(hsl, 200.0, 60.0);
							material.emission = albedo_sqrt * hsl.y * blue * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 22-24
						if (material_mask == 22u) { // 22
							// Warped hyphae
							#ifdef HARDCODED_EMISSION
							float blue = isolate_hue(hsl, 200.0, 60.0);
							material.emission = albedo_sqrt * hsl.y * blue;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 23
							// Crimson stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.yz, 1.0 - block_pos.yz),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.x))
							);
							material.emission = albedo_sqrt * linear_step(0.33, 0.5, hsl.z) * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				}
			} else { // 24-32
				if (material_mask < 28u) { // 24-28
					if (material_mask < 26u) { // 24-26
						if (material_mask == 24u) { // 24
							// Crimson stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.xz, 1.0 - block_pos.xz),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.y))
							);
							material.emission = albedo_sqrt * linear_step(0.33, 0.5, hsl.z) * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 25
							// Crimson stem
							#ifdef HARDCODED_EMISSION
							float emission_amount = mix(
								1.0,
								float(any(lessThan(
									vec4(block_pos.xy, 1.0 - block_pos.xy),
									vec4(rcp(16.0) - 1e-3)
								))),
								step(0.5, abs(normal.z))
							);
							material.emission = albedo_sqrt * linear_step(0.33, 0.5, hsl.z) * emission_amount;
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 26-28
						if (material_mask == 26u) { // 26
							// Crimson hyphae
							#ifdef HARDCODED_EMISSION
							material.emission = albedo_sqrt * linear_step(0.33, 0.5, hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
		     				float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 27
							// Particles
							#ifdef HARDCODED_PARTICLES
							vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
							float l = 0.5 * (min_of(ap1) + max_of(ap1));
							float redness = ap1.r * rcp(ap1.g + ap1.b);
							material.emission = PARTICLES_BRIGHTNESS * (0.47 * material.albedo * step(0.378, redness * l) + 2.5 * albedo_sqrt * (0.2 + 0.8 * isolate_hue(hsl, 15.0, 15.0)) * step(0.871, hsl.y) * step(0.35, hsl.z));

							#if defined WORLD_END && defined END_COLORED_LIGHTING
							material.emission -= 0.2 * albedo_sqrt * linear_step(0.78, 0.85, hsl.z);
							material.emission += 0.3 * sqr(sqr(vec3(0.70, 0.10, 1.0) * 2)) * linear_step(0.78, 0.85, hsl.z) * float(albedo_sqrt);
							#endif

							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * linear_step(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 28-32
					if (material_mask < 30) { // 28-30
						if (material_mask == 28u) { // 28  Ore
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
							#ifdef HARDCODED_ORE
							material.emission = 0.07 * ORE_BRIGHTNESS * albedo_sqrt * albedo_sqrt * (0.9 + 0.1 * isolate_hue(hsl, 15.0, 0.1)) * step(0.01, hsl.y) * step(0.00001, hsl.x) * step(0.35, hsl.z);
							#endif

						} else { // 29 Diamond ore
						#ifdef HARDCODED_SPECULAR
						float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
						material.roughness = sqr(1.0 - smoothness);
						material.f0 = vec3(0.02);
						#endif
						#ifdef HARDCODED_ORE
						vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
						float l = 0.5 * (min_of(ap1) + max_of(ap1));
						float redness = ap1.r * rcp(ap1.g + ap1.b);
						material.emission = 0.5 * ORE_BRIGHTNESS * material.albedo * step(0.40, redness * l) + 0.3 * albedo_sqrt * albedo_sqrt * (0.2 + 0.8 * isolate_hue(hsl, 15.0, 15.0)) * step(0.10, hsl.y) * step(0.76, hsl.z);
						#endif
						}
					} else { // 30-32
						if (material_mask == 30u) { // 30
						#ifdef HARDCODED_SSS
							// SSS (not transparent blocks)
							material.sss_amount = 0.22;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * linear_step(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif

						} else { // 31 Polished blocks
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 0.1;
							#endif
						}
					}
				}
			}
		}
	} else if (material_mask < 64u) { // 32-64
		if (material_mask < 48u) { // 32-48
			if (material_mask < 40u) { // 32-40
				if (material_mask < 36u) { // 32-36
					if (material_mask < 34u) { // 32-34
						if (material_mask == 32u) { // 32
							#ifdef HARDCODED_EMISSION
							// Strong white light
							material.emission = 1.00 * albedo_sqrt * (0.1 + 0.9 * cube(hsl.z));
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 33
							#ifdef HARDCODED_EMISSION
							// Medium white light
							material.emission = 0.66 * albedo_sqrt * linear_step(0.75, 0.9, hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 34-36
						if (material_mask == 34u) { // 34
							#ifdef HARDCODED_EMISSION
							// Weak white light
							material.emission = vec3(1.0) * albedo_sqrt * 0.45 - 0.449 * albedo_sqrt * step(0.20, hsl.x) * step(0.10, hsl.y) * step(0.15, hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 35
							#ifdef HARDCODED_EMISSION
							// Strong golden light
							material.emission  = 0.85 * albedo_sqrt * hsl.z * linear_step(0.4, 0.6, 0.2 * hsl.y + 0.55 * hsl.z);
							light_levels.x *= 0.85;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 36-40
					if (material_mask < 38u) { // 36-38
						if (material_mask == 36u) { // 36
							#ifdef HARDCODED_EMISSION
							// Medium golden light
							material.albedo = material.albedo * vec3(0.0) + float(material.albedo) * vec3(MEDIUM_GOLDEN_LIGHT_R, MEDIUM_GOLDEN_LIGHT_G, MEDIUM_GOLDEN_LIGHT_B);
							material.emission  = 2.0 * sqr(albedo_sqrt) * linear_step(0.78, 0.85, hsl.z);

							#if defined WORLD_END && defined END_COLORED_LIGHTING
							material.albedo = material.albedo * vec3(0.0) + float(material.albedo) * vec3(0.70, 0.10, 1.00);
							material.emission  = 2.0 * sqr(albedo_sqrt) * linear_step(0.78, 0.85, hsl.z);
							#endif
							light_levels.x *= 0.85;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 37
							#ifdef HARDCODED_EMISSION
							// Weak golden light
							float blue = isolate_hue(hsl, 200.0, 30.0);
							material.emission = 1.4 * albedo_sqrt * linear_step(0.47, 0.50, 0.2 * hsl.y + 0.49 * hsl.z + 0.1 * blue);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 38-40
						if (material_mask == 38u) { // 38
							#ifdef HARDCODED_EMISSION
							// Redstone components
							vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
							float l = 0.5 * (min_of(ap1) + max_of(ap1));
							float redness = ap1.r * rcp(ap1.g + ap1.b);
							material.emission = 0.33 * material.albedo * step(0.45, redness * l);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 39
							#ifdef HARDCODED_EMISSION
							// Lava
							material.albedo = float(material.albedo) * vec3(LAVA_LIGHT_R, LAVA_LIGHT_G, LAVA_LIGHT_B);
							material.emission = LAVA_EMISSION_I * (12 * albedo_sqrt * smoothstep(1.0, 0.6 ,hsl.y) * step(0.1, hsl.y) - 0.3 * albedo_sqrt * smoothstep(0.6, 1.0 ,hsl.y));

							#if defined WORLD_END && defined END_COLORED_LIGHTING
							material.albedo = float(material.albedo) * vec3(0.70, 0.10, 1.00);
							material.emission = LAVA_EMISSION_I * (12 * albedo_sqrt * smoothstep(1.0, 0.6 ,hsl.y) * step(0.1, hsl.y) - 0.3 * albedo_sqrt * smoothstep(0.6, 1.0 ,hsl.y));
							#endif

							#endif
						}
					}
				}
			} else { // 40-48
				if (material_mask < 44u) { // 40-44
					if (material_mask < 42u) { // 40-42
						if (material_mask == 40u) { // 40
							#ifdef HARDCODED_EMISSION
							// Medium orange emissives
							material.albedo = float(material.albedo) * vec3(MEDIUM_ORANGE_LIGHT_R, MEDIUM_ORANGE_LIGHT_G, MEDIUM_ORANGE_LIGHT_B);
							//material.emission = 1.40 * sqr(albedo_srgb) * (0.1 + 0.9 * cube(hsl.z));
							material.emission = 1.00 * sqr(albedo_srgb) * (0.9 * cube(hsl.z)) + sqr(albedo_sqrt) * smoothstep(0.80, 1.00, hsl.z);

							#if defined WORLD_END && defined END_COLORED_LIGHTING
							material.albedo = float(material.albedo) * vec3(0.70, 0.10, 1.00);
							material.emission = 1.00 * sqr(albedo_srgb) * (0.9 * cube(hsl.z)) + sqr(albedo_sqrt) * smoothstep(0.80, 1.00, hsl.z);
							#endif
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 41
							#ifdef HARDCODED_EMISSION
							// Brewing stand
							material.emission  = 0.85 * albedo_sqrt * linear_step(0.77, 0.85, hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 42-44
						if (material_mask == 42u) { // 42
							#ifdef HARDCODED_EMISSION
							// Jack o' Lantern
							material.emission = 0.80 * albedo_sqrt * step(0.73, 0.8 * hsl.z);
							light_levels.x *= 0.85;
							#endif
							#ifdef HARDCODED_SSS
							material.sss_amount = 0.2;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.4, 0.8, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 43
							#ifdef HARDCODED_EMISSION
							// Soul lights
							float blue = isolate_hue(hsl, 200.0, 28.0);
							material.emission = 0.66 * albedo_sqrt * linear_step(0.8, 1.0, blue + hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 44-48
					if (material_mask < 46u) { // 44-46
						if (material_mask == 44u) { // 44
							#ifdef HARDCODED_EMISSION
							// Beacon
							material.emission = step(0.2, hsl.z) * albedo_sqrt * step(max_of(abs(block_pos - 0.5)), 0.4);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 45
							#ifdef HARDCODED_EMISSION
							// End portal frame
							material.emission = 0.33 * material.albedo * isolate_hue(hsl, 120.0, 50.0);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 46-48
						if (material_mask == 46u) { // 46
							#ifdef HARDCODED_EMISSION
							// Sculk
							material.emission = 0.2 * material.albedo * isolate_hue(hsl, 200.0, 40.0) * smoothstep(0.5, 0.7, hsl.z) * (1.0 - linear_step(0.0, 20.0, distance(world_pos, cameraPosition)));
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 47
							//Pink glow
							#ifdef HARDCODED_EMISSION
							//material.emission = vec3(0.75) * 1.5 * isolate_hue(hsl, 325.0, 36.7) + 2.5 * albedo_sqrt * pow4(hsl.z) * step(0.53, hsl.z);
							material.emission = vec3(material.albedo);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * linear_step(0.03, 0.01, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 0.8;
							#endif
						}
					}
				}
			}
		} else { // 48-64
			if (material_mask < 56u) { // 48-56
				if (material_mask < 52u) { // 48-52
					if (material_mask < 50u) { // 48-50
						if (material_mask == 48u) { // 48
							material.emission = 0.5 * albedo_sqrt * linear_step(0.5, 0.6, hsl.z);
						} else { // 49
							#ifdef HARDCODED_EMISSION
							// Nether mushrooms
							material.emission = 0.80 * albedo_sqrt * step(0.73, 0.1 * hsl.y + 0.7 * hsl.z);
							#endif
							#ifdef HARDCODED_SSS
							material.sss_amount = 0.25;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 50-52
						if (material_mask == 50u) { // 50
							#ifdef HARDCODED_EMISSION
							// Redstone wire
							vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
							float l = 0.5 * (min_of(ap1) + max_of(ap1));
							float redness = ap1.r * rcp(ap1.g + ap1.b);
							material.emission = (0.2 * material.albedo * step(0.28, redness * l) + 0.3 * (0.38 + 0.02 * isolate_hue(hsl, 15.0, 15.0)) * step(0.771, hsl.y) * step(0.0, hsl.z));
							#endif
							#ifdef HARDCODED_SSS
							material.sss_amount = 0.15;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.75 * sin(smoothstep(0.2, 1.0, hsl.z));
							material.roughness = min(sin(1.1 - smoothness), 0.95);
							material.f0 = material.albedo * vec3(0.168, 0.043, 0.043);
							#endif
						} else { // 51
							#ifdef HARDCODED_EMISSION
							// Ochre froglight
							material.emission = 0.40 * albedo_sqrt * (0.1 + 0.9 * cube(hsl.z));
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 52-56
					if (material_mask < 54u) { // 52-54
						if (material_mask == 52u) { // 52
							#ifdef HARDCODED_EMISSION
							// Verdant froglight
							material.emission = 0.40 * albedo_sqrt * (0.1 + 0.9 * cube(hsl.z));
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 53
							#ifdef HARDCODED_EMISSION
							// Pearlescent froglight
							material.emission = 0.40 * albedo_sqrt * (0.1 + 0.9 * cube(hsl.z));
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 54-56
						if (material_mask == 54u) { // 54
							#ifdef HARDCODED_EMISSION
							// ench table
							float blue = isolate_hue(hsl, 200.0, 120.0);
							material.emission = 0.66 * albedo_sqrt * linear_step(0.97, 1.0, blue + hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.5 * linear_step(0.03, 0.01, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 0.8;
							#endif
						} else { // 55
							#ifdef HARDCODED_EMISSION
							// Amethyst cluster
							material.emission = vec3(0.35) * 2.0 * albedo_srgb * sin((0.1 + 0.9 * hsl.z)) * albedo_sqrt * pow4(hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * linear_step(0.1, 0.2, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							material.ssr_multiplier = 1.0;
							#endif
							#ifdef HARDCODED_SSS
							material.sss_amount = 0.3;
							#endif
						}
					}
				}
			} else { // 56-64
				if (material_mask < 60u) { // 56-60
					if (material_mask < 58u) { // 56-58
						if (material_mask == 56u) { // 56
							#ifdef HARDCODED_EMISSION
							// Calibrated sculk sensor
							material.emission  = 0.2 * material.albedo * isolate_hue(hsl, 200.0, 40.0) * smoothstep(0.5, 0.7, hsl.z) * (1.0 - linear_step(0.0, 20.0, distance(world_pos, cameraPosition)));
							material.emission += vec3(0.20) * (0.1 + 0.9 * hsl.z) * step(0.5, isolate_hue(hsl, 270.0, 50.0) + 0.55 * hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						} else { // 57
							#ifdef HARDCODED_EMISSION
							// Active sculk sensor
							material.emission = vec3(0.20) * (0.1 + 0.9 * hsl.z);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					} else { // 58-60
						if (material_mask == 58u) { // 58
							#ifdef HARDCODED_EMISSION
							// Redstone block
							vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
							float l = 0.5 * (min_of(ap1) + max_of(ap1));
							float redness = ap1.r * rcp(ap1.g + ap1.b);
							material.emission = (0.3 * material.albedo * step(0.3, redness * l) + 0.3 * (0.38 + 0.02 * isolate_hue(hsl, 15.0, 15.0)) * step(0.771, hsl.y) * step(0.0, hsl.z));
							#endif
							#ifdef HARDCODED_SSS
							material.sss_amount = 0.15;
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.75 * sin(smoothstep(0.2, 1.0, hsl.z));
							material.roughness = min(sin(1.1 - smoothness), 0.95);
							material.f0 = material.albedo * vec3(0.168, 0.043, 0.043);
							#endif
						} else { // 59
							//Purple weak light
							#ifdef HARDCODED_EMISSION
							material.emission = 0.80 * albedo_sqrt * albedo_srgb * step(0.8, 0.1 * hsl.y + 0.7 * hsl.z + 0.19 * hsl.x);
							#endif
							#ifdef HARDCODED_SPECULAR
							float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
							material.roughness = sqr(1.0 - smoothness);
							material.f0 = vec3(0.02);
							#endif
						}
					}
				} else { // 60-64
					if (material_mask < 62u) { // 60-62
						if (material_mask == 60u) { // 60

						} else { // 61

						}
					} else { // 62-64
						if (material_mask == 62u) { // 62
							// Nether portal
							material.emission = vec3(1.0);
						} else {  // 63
							// End portal
							material.emission = vec3(1.0);
						}
					}
				}
			}
		}
	}

		if(material_mask == 80u) { // 80
			// Powered lightning rod
			material.emission = vec3(1.0);
		}


	// 			if (material_mask > 191u && material_mask < 240u || material_mask > 63u && material_mask < 80u) { // 264 - 331
	// 	// 343
	// 	#ifdef HARDCODED_EMISSION
	// 	material.emission = material.albedo * isolate_hue(hsl, 0, 360);
	// 	#endif
	// 	#ifdef HARDCODED_SPECULAR
	// 	float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
	// 	material.roughness = sqr(1.0 - smoothness);
	// 	material.f0 = vec3(0.02);
	// 	#endif
	// }

	if (material_mask == 192u || material_mask == 208u || material_mask == 64u || material_mask == 224u) { // 264 - 331
		// Red
		#ifdef HARDCODED_EMISSION
		material.emission = albedo_sqrt * isolate_hue(hsl,-40.0, 65);
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

	if (material_mask == 193u || material_mask == 209u || material_mask == 65u || material_mask == 225u) { // 264 - 331
		// Orange
		#ifdef HARDCODED_EMISSION
		material.emission = material.albedo * isolate_hue(hsl, 15.0, 50) * linear_step(0.2, 1.0, hsl.z);
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask == 194u || material_mask == 210u || material_mask == 66u || material_mask == 226u) { // 264 - 331
		// Yellow
		#ifdef HARDCODED_EMISSION
		material.emission = max(material.albedo * isolate_hue(hsl, 45.0, 20), material.albedo * linear_step(0.7, 0.8, hsl.z));
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask == 195u || material_mask == 211u || material_mask == 67u || material_mask == 227u) { // 264 - 331
		// Brown
		#ifdef HARDCODED_EMISSION
		material.emission = max(material.albedo * isolate_hue(hsl, 45.0, 25), material.albedo * linear_step(0.7, 0.8, hsl.z));
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask == 196u || material_mask == 212u || material_mask == 68u || material_mask == 228u || material_mask == 197u || material_mask == 213u || material_mask == 69u || material_mask == 229u) { // 264 - 331
		// Green, Lime
		#ifdef HARDCODED_EMISSION
		material.emission = max(material.albedo * isolate_hue(hsl, 105.0, 40), material.albedo * linear_step(0.5, 0.8, hsl.z));
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask > 197u && material_mask < 201u || material_mask > 213u && material_mask < 217u || material_mask > 69u && material_mask < 73u || material_mask > 229u && material_mask < 233u) { // 264 - 331
		// Blue, Light blue, Cyan
		#ifdef HARDCODED_EMISSION
		material.emission = max(material.albedo * isolate_hue(hsl, 230.0, 70), material.albedo * linear_step(0.7, 0.8, hsl.z));
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask > 200u && material_mask < 203u || material_mask > 216u && material_mask < 219u || material_mask > 72u && material_mask < 75u || material_mask > 232u && material_mask < 235u) { // 264 - 331
		// Purple, Magenta
		#ifdef HARDCODED_EMISSION
		material.emission = material.albedo * albedo_srgb * isolate_hue(hsl, 300, 42);
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

			if (material_mask == 203u || material_mask == 219u || material_mask == 75u || material_mask == 235u) { // 264 - 331
		//Pink
		#ifdef HARDCODED_EMISSION
		material.emission = material.albedo * albedo_srgb * isolate_hue(hsl, 300, 50);
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

		if (material_mask > 203u && material_mask < 208u || material_mask > 219u && material_mask < 224u || material_mask > 75u && material_mask < 80u || material_mask > 235u && material_mask < 240u) { // 264 - 331
		// Black, White, Gray, Light gray
		#ifdef HARDCODED_EMISSION
		material.emission = max(material.albedo * isolate_hue(hsl, 50.0, 40), material.albedo * linear_step(0.7, 0.8, hsl.z));
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

	if (material_mask == 181u) {
		//Transparent blocks
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.95, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

	if (164u <= material_mask && material_mask < 180u) {
		// Stained glass, honey and slime
		#ifdef HARDCODED_SPECULAR
		material.f0 = vec3(0.04);
		material.roughness = 0.1;
		material.ssr_multiplier = 1.0;
		#endif

		#ifdef HARDCODED_SSS
		material.sss_amount = 0.5;
		#endif
	}

	if (material_mask == 180u) {
		// Obsidian
		#ifdef HARDCODED_SPECULAR
		float smoothness = linear_step(0.02, 0.4, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		material.ssr_multiplier = 0.8;
		#endif
	}

	if (material_mask > 184u && material_mask < 190u) { //Ore
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
		#ifdef HARDCODED_ORE
		vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
		float l = 0.5 * (min_of(ap1) + max_of(ap1));
		float redness = ap1.r * rcp(ap1.g + ap1.b);
		material.emission = 0.2 * ORE_BRIGHTNESS * albedo_sqrt * step(0.0, normalize(step(0.40, redness * l) + 0.3 * albedo_sqrt * albedo_sqrt * (0.2 + 0.8 * isolate_hue(hsl, 180.0, 180.0)) * step(0.10, hsl.y) * step(0.4, hsl.z)));
		#endif
	}

		if (material_mask == 191u) { //Redstone Ore
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
		#ifdef HARDCODED_ORE
		vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
		float l = 0.5 * (min_of(ap1) + max_of(ap1));
		float redness = ap1.r * rcp(ap1.g + ap1.b);
		material.emission = 0.2 * ORE_BRIGHTNESS * albedo_sqrt * step(0.0, normalize(step(0.40, redness * l) + 0.3 * albedo_sqrt * albedo_sqrt * (0.2 + 0.8 * isolate_hue(hsl, 180.0, 180.0)) * step(0.10, hsl.y) * step(0.4, hsl.z)));
		#endif
	}

	if (material_mask == 190u) { //Diamond ore
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
		#ifdef HARDCODED_ORE
		vec3 ap1 = material.albedo * rec2020_to_ap1_unlit;
		float l = 0.5 * (min_of(ap1) + max_of(ap1));
		float redness = ap1.r * rcp(ap1.g + ap1.b);
		material.emission = 0.2 * ORE_BRIGHTNESS * albedo_sqrt * step(0.0, normalize(step(0.40, redness * l) + 0.3 * albedo_sqrt * albedo_sqrt * (0.2 + 0.8 * isolate_hue(hsl, 180.0, 180.0)) * step(0.1, hsl.y) * step(0.76, hsl.z)));
		#endif
	}

	if (material_mask == 182u) {
		// End crystal
		#ifdef HARDCODED_EMISSION
		material.emission = isolate_hue(hsl, 310.0, 50.0) * 2.2 * albedo_sqrt * albedo_sqrt;
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 1.13 * sin(linear_step(0.0, 1.0, hsl.z));
		material.roughness = min(sin(1.0 - smoothness), 0.65);
		material.f0 = sqr(material.albedo) * 1.2;
		material.ssr_multiplier = 1.2;
		#endif
	}

	if (material_mask == 183u) {
		//
		#ifdef HARDCODED_EMISSION
		material.emission = 0.50 * albedo_sqrt * step(0.8, 0.1 * hsl.y + 0.7 * hsl.z + 0.3 * hsl.x);
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.4 * smoothstep(0.2, 0.6, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

	if (material_mask == 184u) {
		//
		#ifdef HARDCODED_EMISSION
		float blue = isolate_hue(hsl, 150.0, 60.0);
		material.emission = albedo_sqrt * hsl.y * blue;
		#endif
		#ifdef HARDCODED_SPECULAR
		float smoothness = 0.45 * smoothstep(0.01, 0.7, hsl.z);
		material.roughness = sqr(1.0 - smoothness);
		material.f0 = vec3(0.02);
		#endif
	}

	return material;
}

#endif // INCLUDE_MISC_MATERIAL
