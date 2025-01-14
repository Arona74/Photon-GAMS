#if !defined INCLUDE_VERTEX_DISPLACEMENT
#define INCLUDE_VERTEX_DISPLACEMENT

#if !defined PROGRAM_GBUFFERS_TERRAIN && !defined PROGRAM_SHADOW
	#undef WAVING_PLANTS
	#undef WAVING_LEAVES
#endif

#if !defined PROGRAM_GBUFFERS_WATER && !defined PROGRAM_SHADOW
	#undef WATER_DISPLACEMENT
#endif

#ifdef IS_IRIS
uniform vec3 eyePosition;
#else
#define eyePosition cameraPosition
#endif


#if defined WATER_DISPLACEMENT
float gerstner_wave(vec2 coord, vec2 wave_dir, float t, float noise, float wavelength) {
	// Gerstner wave function from Belmu in #snippets, modified
	const float g = 9.8;

	float k = tau / wavelength;
	float w = sqrt(g * k);

	float x = w * t - k * (dot(wave_dir, coord) + noise);

	return sqr(sin(x) * 0.5 + 0.5);
}

float get_water_displacement(vec3 world_pos, float skylight) {
	const float wave_frequency = 0.3 * WATER_WAVE_FREQUENCY;
	const float wave_speed     = 0.37 * WATER_WAVE_SPEED_STILL;
	const float wave_angle     = WATER_WAVE_ANGLE * degree;
	const float wavelength     = 1.0;
	const vec2  wave_dir       = vec2(cos(wave_angle), sin(wave_angle));

	float wave = gerstner_wave(world_pos.xy * wave_frequency, wave_dir, frameTimeCounter * wave_speed, 0.0, wavelength);
		  wave = (wave * 0.05 - 0.025) * (skylight * 0.9 + 0.1);

	return wave;
}
#endif

#if defined WAVING_PLANTS || defined WAVING_LEAVES
vec3 get_wind_displacement(vec3 world_pos, float wind_speed, float wind_strength, bool is_tall_plant_top_vertex) {
	const float wind_angle = 30.0 * degree;
	const vec2  wind_dir   = vec2(cos(wind_angle), sin(wind_angle));

	float t = wind_speed * frameTimeCounter;

	float gust_amount  = texture(noisetex, 0.05 * (world_pos.xz + wind_dir * t)).y;
		  gust_amount *= gust_amount;

	vec3 gust = vec3(wind_dir * gust_amount, 0.1 * gust_amount).xzy;

	world_pos = 32.0 * world_pos + 3.0 * t + vec3(0.0, golden_angle, 2.0 * golden_angle);
	vec3 wobble = sin(world_pos) + 0.5 * sin(2.0 * world_pos) + WAVING_CLEAR_STRENGTH * sin(WAVING_RAIN_STRENGTH * world_pos);

	if (is_tall_plant_top_vertex) { gust *= 2.0; wobble *= 0.5; }

	return wind_strength * (gust + 0.1 * wobble);
}
#endif

vec3 animate_vertex(vec3 world_pos, bool is_top_vertex, float skylight, uint material_mask) {
	float wind_speed = 0.3;
	float wind_strength = sqr(skylight) * (0.25 + 0.66 * rainStrength);

	// Displace plants close to the player
	vec3 to_player = eyePosition - world_pos;
	vec3 player_displacement = vec3(
		-6.0 * to_player.xz * exp2(-length(to_player * vec3(6.0, 2.0, 6.0))),
		0.0
	).xzy;

	switch (material_mask) {
#ifdef WATER_DISPLACEMENT
	case 1:
		world_pos.y += get_water_displacement(world_pos, skylight);
		return world_pos;
#endif

#ifdef WAVING_PLANTS
	case 2:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);

	case 3:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);

	case 4:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, is_top_vertex) + player_displacement);

	//Glowing Flowers
	case 192:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 65:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 66:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 69:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 70:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 71:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 72:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 73:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 74:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 75:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 76:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
	case 77:
		return world_pos + (get_wind_displacement(world_pos, wind_speed, wind_strength, false) + player_displacement) * float(is_top_vertex);
#endif

#ifdef WAVING_LEAVES
	case 5:
		return world_pos + get_wind_displacement(world_pos, wind_speed, wind_strength * 0.5, false);
#endif

	default:
		return world_pos;
	}
}

//Optifine #if fix, added boolean settings, removed 0.0 as SIZE sliders values to avoid comparing float
#if defined (WORLD_END) && defined (END_CURVATURE)
	#define CURVATURE_SIZE END_CURVATURE_SIZE
#elif defined(WORLD_NETHER) && defined (NETHER_CURVATURE)
	#define CURVATURE_SIZE NETHER_CURVATURE_SIZE
#elif defined(WORLD_OVERWORLD) && defined (OVERWORLD_CURVATURE)
	#define CURVATURE_SIZE OVERWORLD_CURVATURE_SIZE
#endif

vec3 world_curvature(vec3 scene_pos) {
	#if defined (WORLD_CURVATURE) && defined (CURVATURE_SIZE)
	//scene_pos += cameraPosition;
	#if CURVATURE_MODE == CURVATURE_MODE_SQUARED_DISTANCE

		scene_pos.y -= length_squared(scene_pos.xz) / CURVATURE_SIZE;

	#elif CURVATURE_MODE == CURVATURE_MODE_BASIC_SPHERICAL || CURVATURE_MODE == CURVATURE_MODE_LOG_SPHERICAL

		const float size = CURVATURE_SIZE;

		#if CURVATURE_MODE == CURVATURE_MODE_BASIC_SPHERICAL
		float h = scene_pos.y + size;
		#elif CURVATURE_MODE == CURVATURE_MODE_LOG_SPHERICAL
		float h = exp(scene_pos.y / size) * size;
		#endif

		vec2 azimuth = cartesian_to_polar(scene_pos.xz);
		float phi = azimuth.x / size;

		// Prevent looping with higher render distance
		if(phi > pi) return vec3(0, intBitsToFloat(0xff800000), 0); // negative infinity ¯\_(o_o)_/¯

		float h_sin_phi = h * sin(phi);
		scene_pos.y = h * cos(phi) /*+ 128.0*/ /* exp(128.0 / CURVATURE_SIZE)*/ - size;
		scene_pos.x = h_sin_phi * cos(azimuth.y);
		scene_pos.z = h_sin_phi * sin(azimuth.y);

	#endif
	//scene_pos -= cameraPosition;
#endif
	return scene_pos;
}

#endif // INCLUDE_VERTEX_DISPLACEMENT
