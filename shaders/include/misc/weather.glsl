#if !defined INCLUDE_MISC_WEATHER
#define INCLUDE_MISC_WEATHER

#include "/include/misc/weather_struct.glsl"
#include "/include/utility/color.glsl"
#include "/include/utility/random.glsl"

#define daily_weather_blend(weather_function) mix(weather_function(worldDay), weather_function(worldDay + 1), weather_mix_factor())

uint weather_day_index(int world_day) {
	// Start at noon
	world_day -= int(worldTime <= 6000);

	const uint day_count = 12;

#if WEATHER_DAY == -1
	uint day_index = uint(world_day);
	     day_index = (day_index) % day_count;
#else
	uint day_index = WEATHER_DAY;
#endif

	return day_index;
}

float weather_mix_factor() {
	return cubic_smooth(fract(float(worldTime) * rcp(24000.0) - 0.25));
}

float daily_weather_fogginess(int world_day) {
	const float[] fogginess = float[12](WEATHER_D0_FOGGINESS, WEATHER_D1_FOGGINESS, WEATHER_D2_FOGGINESS, WEATHER_D3_FOGGINESS, WEATHER_D4_FOGGINESS, WEATHER_D5_FOGGINESS, WEATHER_D6_FOGGINESS, WEATHER_D7_FOGGINESS, WEATHER_D8_FOGGINESS, WEATHER_D9_FOGGINESS, WEATHER_D10_FOGGINESS, WEATHER_D11_FOGGINESS);

	return fogginess[weather_day_index(world_day)];
}

// Clouds
void daily_weather_clouds(
	int world_day,
	out vec2 clouds_cumulus_coverage,
	out vec2 clouds_towering_cumulus_coverage,
	out vec2 clouds_thunderhead_coverage,
	out vec2 clouds_altocumulus_coverage,
	out vec2 clouds_cirrus_coverage,
	out float clouds_cumulus_congestus_amount,
	out float clouds_cumulonimbus_amount,
	out float clouds_stratus_amount
) {
	const uint day_count = 12;

	uint day_index = weather_day_index(world_day);

	switch (day_index) {
	case 0:
		clouds_cumulus_coverage            = vec2(WEATHER_D0_CLOUDS_CUMULUS_MIN, WEATHER_D0_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D0_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D0_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D0_CLOUDS_THUNDERHEAD_MIN, WEATHER_D0_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D0_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D0_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D0_CLOUDS_CIRRUS, WEATHER_D0_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D0_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D0_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D0_CLOUDS_STRATUS_AMOUNT;
		break;

	case 1:
		clouds_cumulus_coverage            = vec2(WEATHER_D1_CLOUDS_CUMULUS_MIN, WEATHER_D1_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D1_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D1_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D1_CLOUDS_THUNDERHEAD_MIN, WEATHER_D1_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D1_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D1_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D1_CLOUDS_CIRRUS, WEATHER_D1_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D1_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D1_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D1_CLOUDS_STRATUS_AMOUNT;
		break;

	case 2:
		clouds_cumulus_coverage            = vec2(WEATHER_D2_CLOUDS_CUMULUS_MIN, WEATHER_D2_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D2_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D2_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D2_CLOUDS_THUNDERHEAD_MIN, WEATHER_D2_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D2_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D2_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D2_CLOUDS_CIRRUS, WEATHER_D2_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D2_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D2_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D2_CLOUDS_STRATUS_AMOUNT;
		break;

	case 3:
		clouds_cumulus_coverage            = vec2(WEATHER_D3_CLOUDS_CUMULUS_MIN, WEATHER_D3_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D3_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D3_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D3_CLOUDS_THUNDERHEAD_MIN, WEATHER_D3_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D3_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D3_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D3_CLOUDS_CIRRUS, WEATHER_D3_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D3_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D3_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D3_CLOUDS_STRATUS_AMOUNT;
		break;

	case 4:
		clouds_cumulus_coverage            = vec2(WEATHER_D4_CLOUDS_CUMULUS_MIN, WEATHER_D4_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D4_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D4_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D4_CLOUDS_THUNDERHEAD_MIN, WEATHER_D4_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D4_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D4_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D4_CLOUDS_CIRRUS, WEATHER_D4_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D4_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D4_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D4_CLOUDS_STRATUS_AMOUNT;
		break;

	case 5:
		clouds_cumulus_coverage            = vec2(WEATHER_D5_CLOUDS_CUMULUS_MIN, WEATHER_D5_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D5_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D5_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D5_CLOUDS_THUNDERHEAD_MIN, WEATHER_D5_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D5_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D5_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D5_CLOUDS_CIRRUS, WEATHER_D5_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D5_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D5_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D5_CLOUDS_STRATUS_AMOUNT;
		break;

	case 6:
		clouds_cumulus_coverage            = vec2(WEATHER_D6_CLOUDS_CUMULUS_MIN, WEATHER_D6_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D6_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D6_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D6_CLOUDS_THUNDERHEAD_MIN, WEATHER_D6_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D6_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D6_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D6_CLOUDS_CIRRUS, WEATHER_D6_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D6_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D6_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D6_CLOUDS_STRATUS_AMOUNT;
		break;

	case 7:
		clouds_cumulus_coverage            = vec2(WEATHER_D7_CLOUDS_CUMULUS_MIN, WEATHER_D7_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D7_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D7_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D7_CLOUDS_THUNDERHEAD_MIN, WEATHER_D7_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D7_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D7_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D7_CLOUDS_CIRRUS, WEATHER_D7_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D7_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D7_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D7_CLOUDS_STRATUS_AMOUNT;
		break;

	case 8:
		clouds_cumulus_coverage            = vec2(WEATHER_D8_CLOUDS_CUMULUS_MIN, WEATHER_D8_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D8_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D8_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D8_CLOUDS_THUNDERHEAD_MIN, WEATHER_D8_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D8_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D8_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D8_CLOUDS_CIRRUS, WEATHER_D8_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D8_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D8_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D8_CLOUDS_STRATUS_AMOUNT;
		break;

	case 9:
		clouds_cumulus_coverage            = vec2(WEATHER_D9_CLOUDS_CUMULUS_MIN, WEATHER_D9_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D9_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D9_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D9_CLOUDS_THUNDERHEAD_MIN, WEATHER_D9_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D9_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D9_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D9_CLOUDS_CIRRUS, WEATHER_D9_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D9_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D9_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D9_CLOUDS_STRATUS_AMOUNT;
		break;

	case 10:
		clouds_cumulus_coverage            = vec2(WEATHER_D10_CLOUDS_CUMULUS_MIN, WEATHER_D10_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D10_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D10_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D10_CLOUDS_THUNDERHEAD_MIN, WEATHER_D10_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D10_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D10_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D10_CLOUDS_CIRRUS, WEATHER_D10_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D10_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D10_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D10_CLOUDS_STRATUS_AMOUNT;
		break;

	case 11:
		clouds_cumulus_coverage            = vec2(WEATHER_D11_CLOUDS_CUMULUS_MIN, WEATHER_D11_CLOUDS_CUMULUS_MAX);
		clouds_towering_cumulus_coverage   = vec2(WEATHER_D11_CLOUDS_TOWERING_CUMULUS_MIN, WEATHER_D11_CLOUDS_TOWERING_CUMULUS_MAX);
		clouds_thunderhead_coverage        = vec2(WEATHER_D11_CLOUDS_THUNDERHEAD_MIN, WEATHER_D11_CLOUDS_THUNDERHEAD_MAX);
		clouds_altocumulus_coverage     = vec2(WEATHER_D11_CLOUDS_ALTOCUMULUS_MIN, WEATHER_D11_CLOUDS_ALTOCUMULUS_MAX);
		clouds_cirrus_coverage          = vec2(WEATHER_D11_CLOUDS_CIRRUS, WEATHER_D11_CLOUDS_CIRROCUMULUS);
		clouds_cumulus_congestus_amount = WEATHER_D11_CLOUDS_CUMULUS_CONGESTUS_AMOUNT;
		clouds_cumulonimbus_amount      = WEATHER_D11_CLOUDS_CUMULONIMBUS_AMOUNT;
		clouds_stratus_amount           = WEATHER_D11_CLOUDS_STRATUS_AMOUNT;
		break;
	}
}

void clouds_weather_variation(
	out vec2 clouds_cumulus_coverage,
	out vec2 clouds_towering_cumulus_coverage,
	out vec2 clouds_thunderhead_coverage,
	out vec2 clouds_altocumulus_coverage,
	out vec2 clouds_cirrus_coverage,
	out float clouds_cumulus_congestus_amount,
	out float clouds_cumulonimbus_amount,
	out float clouds_stratus_amount
) {
	// Daily weather variation

#ifdef CLOUDS_DAILY_WEATHER
	vec2 coverage_cu_0, coverage_cu_1;
	vec2 coverage_tcu_0, coverage_tcu_1;
	vec2 coverage_th_0, coverage_th_1;
	vec2 coverage_ac_0, coverage_ac_1;
	vec2 coverage_ci_0, coverage_ci_1;
	float cu_con_0, cu_con_1;
	float cb_0, cb_1;
	float stratus_0, stratus_1;

	daily_weather_clouds(worldDay + 0, coverage_cu_0, coverage_tcu_0, coverage_th_0, coverage_ac_0, coverage_ci_0, cu_con_0, cb_0, stratus_0);
	daily_weather_clouds(worldDay + 1, coverage_cu_1, coverage_tcu_1, coverage_th_1, coverage_ac_1, coverage_ci_1, cu_con_1, cb_1, stratus_1);

	float mix_factor = weather_mix_factor();

	clouds_cumulus_coverage          = mix(coverage_cu_0, coverage_cu_1, mix_factor);
	clouds_towering_cumulus_coverage = mix(coverage_tcu_0, coverage_tcu_1, mix_factor);
	clouds_thunderhead_coverage      = mix(coverage_th_0, coverage_th_1, mix_factor);
	clouds_altocumulus_coverage     = mix(coverage_ac_0, coverage_ac_1, mix_factor);
	clouds_cirrus_coverage          = mix(coverage_ci_0, coverage_ci_1, mix_factor);
	clouds_cumulus_congestus_amount = mix(cu_con_0, cu_con_1, mix_factor);
	clouds_cumulonimbus_amount      = mix(cb_0, cb_1, mix_factor);
	clouds_stratus_amount           = mix(stratus_0, stratus_1, mix_factor);
#else
	clouds_cumulus_coverage          = vec2(0.4, 0.55);
	clouds_towering_cumulus_coverage = vec2(0.2, 0.35);
	clouds_thunderhead_coverage      = vec2(0.1, 0.25);
	clouds_altocumulus_coverage     = vec2(0.3, 0.5);
	clouds_cirrus_coverage          = vec2(0.4, 0.5);
	clouds_cumulus_congestus_amount = 0.0;
	clouds_cumulonimbus_amount      = 0.0;
	clouds_stratus_amount           = 0.0;
#endif

	// Weather influence

	clouds_cumulus_coverage = mix(clouds_cumulus_coverage, vec2(0.6, 0.8), wetness);
	clouds_towering_cumulus_coverage = mix(clouds_towering_cumulus_coverage, vec2(0.6, 0.8), wetness);
	clouds_thunderhead_coverage = mix(clouds_thunderhead_coverage, vec2(0.6, 0.8), wetness);
	clouds_altocumulus_coverage = mix(clouds_altocumulus_coverage, vec2(0.4, 0.9), wetness * 0.75);
	clouds_cirrus_coverage.x = mix(clouds_cirrus_coverage.x, 0.7, wetness * 0.50);
	clouds_cumulus_congestus_amount *= 1.0 - wetness;
	clouds_cumulonimbus_amount *= mix(1.0, 1.5 * wetness, wetness); // mix(clamp01(1.0 - 1.75 * wetness) ...
	clouds_stratus_amount = clamp01(clouds_stratus_amount + 0.7 * wetness);

	// User config values

	clouds_cumulus_coverage *= CLOUDS_CUMULUS_COVERAGE;
	clouds_towering_cumulus_coverage *= CLOUDS_TOWERING_CUMULUS_COVERAGE;
	clouds_thunderhead_coverage *= CLOUDS_THUNDERHEAD_COVERAGE;
	clouds_altocumulus_coverage *= CLOUDS_ALTOCUMULUS_COVERAGE;
	clouds_cirrus_coverage *= vec2(CLOUDS_CIRRUS_COVERAGE, CLOUDS_CIRROCUMULUS_COVERAGE);
	clouds_cumulus_congestus_amount *= CLOUDS_CUMULUS_CONGESTUS_COVERAGE;
	clouds_cumulonimbus_amount *= CLOUDS_CUMULONIMBUS_COVERAGE;

}

// [0] - bottom color
// [1] - top color
mat2x3 get_aurora_colors() {
	const mat2x3[] aurora_colors = mat2x3[](
		mat2x3( // 0
			vec3(AURORA_COLOR_0_BTM_R, AURORA_COLOR_0_BTM_G, AURORA_COLOR_0_BTM_B) * AURORA_COLOR_0_BTM_I,
			vec3(AURORA_COLOR_0_TOP_R, AURORA_COLOR_0_TOP_G, AURORA_COLOR_0_TOP_B) * AURORA_COLOR_0_TOP_I
		)
		, mat2x3( // 1
			vec3(AURORA_COLOR_1_BTM_R, AURORA_COLOR_1_BTM_G, AURORA_COLOR_1_BTM_B) * AURORA_COLOR_1_BTM_I,
			vec3(AURORA_COLOR_1_TOP_R, AURORA_COLOR_1_TOP_G, AURORA_COLOR_1_TOP_B) * AURORA_COLOR_1_TOP_I
		)
		, mat2x3( // 2
			vec3(AURORA_COLOR_2_BTM_R, AURORA_COLOR_2_BTM_G, AURORA_COLOR_2_BTM_B) * AURORA_COLOR_2_BTM_I,
			vec3(AURORA_COLOR_2_TOP_R, AURORA_COLOR_2_TOP_G, AURORA_COLOR_2_TOP_B) * AURORA_COLOR_2_TOP_I
		)
		, mat2x3( // 3
			vec3(AURORA_COLOR_3_BTM_R, AURORA_COLOR_3_BTM_G, AURORA_COLOR_3_BTM_B) * AURORA_COLOR_3_BTM_I,
			vec3(AURORA_COLOR_3_TOP_R, AURORA_COLOR_3_TOP_G, AURORA_COLOR_3_TOP_B) * AURORA_COLOR_3_TOP_I
		)
		, mat2x3( // 4
			vec3(AURORA_COLOR_4_BTM_R, AURORA_COLOR_4_BTM_G, AURORA_COLOR_4_BTM_B) * AURORA_COLOR_4_BTM_I,
			vec3(AURORA_COLOR_4_TOP_R, AURORA_COLOR_4_TOP_G, AURORA_COLOR_4_TOP_B) * AURORA_COLOR_4_TOP_I
		)
		, mat2x3( // 5
			vec3(AURORA_COLOR_5_BTM_R, AURORA_COLOR_5_BTM_G, AURORA_COLOR_5_BTM_B) * AURORA_COLOR_5_BTM_I,
			vec3(AURORA_COLOR_5_TOP_R, AURORA_COLOR_5_TOP_G, AURORA_COLOR_5_TOP_B) * AURORA_COLOR_5_TOP_I
		)
		, mat2x3( // 6
			vec3(AURORA_COLOR_6_BTM_R, AURORA_COLOR_6_BTM_G, AURORA_COLOR_6_BTM_B) * AURORA_COLOR_6_BTM_I,
			vec3(AURORA_COLOR_6_TOP_R, AURORA_COLOR_6_TOP_G, AURORA_COLOR_6_TOP_B) * AURORA_COLOR_6_TOP_I
		)
		, mat2x3( // 7
			vec3(AURORA_COLOR_7_BTM_R, AURORA_COLOR_7_BTM_G, AURORA_COLOR_7_BTM_B) * AURORA_COLOR_7_BTM_I,
			vec3(AURORA_COLOR_7_TOP_R, AURORA_COLOR_7_TOP_G, AURORA_COLOR_7_TOP_B) * AURORA_COLOR_7_TOP_I
		)
		, mat2x3( // 8
			vec3(AURORA_COLOR_8_BTM_R, AURORA_COLOR_8_BTM_G, AURORA_COLOR_8_BTM_B) * AURORA_COLOR_8_BTM_I,
			vec3(AURORA_COLOR_8_TOP_R, AURORA_COLOR_8_TOP_G, AURORA_COLOR_8_TOP_B) * AURORA_COLOR_8_TOP_I
		)
		, mat2x3( // 9
			vec3(AURORA_COLOR_9_BTM_R, AURORA_COLOR_9_BTM_G, AURORA_COLOR_9_BTM_B) * AURORA_COLOR_9_BTM_I,
			vec3(AURORA_COLOR_9_TOP_R, AURORA_COLOR_9_TOP_G, AURORA_COLOR_9_TOP_B) * AURORA_COLOR_9_TOP_I
		)
		, mat2x3( // 10
			vec3(AURORA_COLOR_10_BTM_R, AURORA_COLOR_10_BTM_G, AURORA_COLOR_10_BTM_B) * AURORA_COLOR_10_BTM_I,
			vec3(AURORA_COLOR_10_TOP_R, AURORA_COLOR_10_TOP_G, AURORA_COLOR_10_TOP_B) * AURORA_COLOR_10_TOP_I
		)
		, mat2x3( // 11
			vec3(AURORA_COLOR_11_BTM_R, AURORA_COLOR_11_BTM_G, AURORA_COLOR_11_BTM_B) * AURORA_COLOR_11_BTM_I,
			vec3(AURORA_COLOR_11_TOP_R, AURORA_COLOR_11_TOP_G, AURORA_COLOR_11_TOP_B) * AURORA_COLOR_11_TOP_I
		)
	);

#if AURORA_COLOR == -1
		const uint[] weights = uint[](
			AURORA_COLOR_0_WEIGHT, AURORA_COLOR_1_WEIGHT, AURORA_COLOR_2_WEIGHT, AURORA_COLOR_3_WEIGHT, AURORA_COLOR_4_WEIGHT,
			AURORA_COLOR_5_WEIGHT, AURORA_COLOR_6_WEIGHT, AURORA_COLOR_7_WEIGHT, AURORA_COLOR_8_WEIGHT, AURORA_COLOR_9_WEIGHT,
			AURORA_COLOR_10_WEIGHT, AURORA_COLOR_11_WEIGHT
		);
		const uint total_weight = AURORA_COLOR_0_WEIGHT + AURORA_COLOR_1_WEIGHT + AURORA_COLOR_2_WEIGHT + AURORA_COLOR_3_WEIGHT + AURORA_COLOR_4_WEIGHT +
			AURORA_COLOR_5_WEIGHT + AURORA_COLOR_6_WEIGHT + AURORA_COLOR_7_WEIGHT + AURORA_COLOR_8_WEIGHT + AURORA_COLOR_9_WEIGHT +
			AURORA_COLOR_10_WEIGHT + AURORA_COLOR_11_WEIGHT;
		mat2x3[total_weight] aurora_colors_weighted;
		for(uint i = 0u, index = 0u; i < weights.length(); i++) {
			for(uint j = 0u; j < weights[i]; j++, index++) {
				aurora_colors_weighted[index] = aurora_colors[i];
			}
		}

		uint day_index = uint(worldDay);
			 day_index = lowbias32(day_index) % aurora_colors_weighted.length();
		return aurora_colors_weighted[day_index];
#else
		return aurora_colors[uint(AURORA_COLOR)];
#endif

}

// 0.0 - no aurora
// 1.0 - full aurora
float get_aurora_amount() {
	float night = smoothstep(0.0, 0.2, -sun_dir.y);

#if   AURORA_NORMAL == AURORA_NEVER
	float aurora_normal = 0.0;
#elif AURORA_NORMAL == AURORA_RARELY
	float aurora_normal = float(lowbias32(uint(worldDay)) % 5 == 1);
#elif AURORA_NORMAL == AURORA_ALWAYS
	float aurora_normal = 1.0;
#endif

#if   AURORA_SNOW == AURORA_NEVER
	float aurora_snow = 0.0;
#elif AURORA_SNOW == AURORA_RARELY
	float aurora_snow = float(lowbias32(uint(worldDay)) % 5 == 1);
#elif AURORA_SNOW == AURORA_ALWAYS
	float aurora_snow = 1.0;
#endif

	return night * mix(aurora_normal, aurora_snow, biome_may_snow);
}

// 0.0 - no aurora
// 1.0 - full NLC
float get_nlc_amount() {
	float intensity = hash1(fract(float(worldDay) * golden_ratio));
	intensity = linear_step(CLOUDS_NOCTILUCENT_RARITY, 1.0, intensity);

	return dampen(intensity) * CLOUDS_NOCTILUCENT_INTENSITY;
}

DailyWeatherVariation get_daily_weather_variation() {
	DailyWeatherVariation daily_weather_variation;

	clouds_weather_variation(
		daily_weather_variation.clouds_cumulus_coverage,
		daily_weather_variation.clouds_towering_cumulus_coverage,
		daily_weather_variation.clouds_thunderhead_coverage,
		daily_weather_variation.clouds_altocumulus_coverage,
		daily_weather_variation.clouds_cirrus_coverage,
		daily_weather_variation.clouds_cumulus_congestus_amount,
		daily_weather_variation.clouds_cumulonimbus_amount,
		daily_weather_variation.clouds_stratus_amount
	);

	daily_weather_variation.fogginess = daily_weather_blend(daily_weather_fogginess);
	daily_weather_variation.nlc_amount = get_nlc_amount();
	daily_weather_variation.aurora_amount = get_aurora_amount();
	daily_weather_variation.aurora_colors = get_aurora_colors();

	return daily_weather_variation;
}

#endif // INCLUDE_MISC_WEATHER