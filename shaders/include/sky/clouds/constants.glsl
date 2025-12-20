#if !defined INCLUDE_SKY_CLOUDS_CONSTANTS
#define INCLUDE_SKY_CLOUDS_CONSTANTS

#include "/include/sky/atmosphere.glsl"

#ifdef CLOUDS_GAMS_CUMULUS
const float clouds_cumulus_radius          = planet_radius + CLOUDS_CUMULUS_ALTITUDE * 2.45;
#else
const float clouds_cumulus_radius          = planet_radius + CLOUDS_CUMULUS_ALTITUDE;
#endif
const float clouds_cumulus_thickness       = CLOUDS_CUMULUS_ALTITUDE * CLOUDS_CUMULUS_THICKNESS;
const float clouds_cumulus_top_radius      = clouds_cumulus_radius + clouds_cumulus_thickness;

const float clouds_altocumulus_radius      = planet_radius + CLOUDS_ALTOCUMULUS_ALTITUDE;
const float clouds_altocumulus_thickness   = CLOUDS_ALTOCUMULUS_ALTITUDE * CLOUDS_ALTOCUMULUS_THICKNESS;
const float clouds_altocumulus_top_radius  = clouds_altocumulus_radius + clouds_altocumulus_thickness;

const float clouds_cirrus_radius           = planet_radius + CLOUDS_CIRRUS_ALTITUDE;
const float clouds_cirrus_thickness        = CLOUDS_CIRRUS_ALTITUDE * CLOUDS_ALTOCUMULUS_THICKNESS;
const float clouds_cirrus_top_radius       = clouds_cirrus_radius + clouds_cirrus_thickness;
const float clouds_cirrus_extinction_coeff = 0.15;
const float clouds_cirrus_scattering_coeff = clouds_cirrus_extinction_coeff;

const float clouds_noctilucent_altitude    = 80000.0;
const float clouds_noctilucent_radius      = planet_radius + clouds_noctilucent_altitude;

// GAMS Clouds
const float clouds_cumulus_congestus_radius           = planet_radius + CLOUDS_CUMULUS_CONGESTUS_ALTITUDE;
const float clouds_cumulus_congestus_thickness        = CLOUDS_CUMULUS_CONGESTUS_ALTITUDE * CLOUDS_CUMULUS_CONGESTUS_THICKNESS;
const float clouds_cumulus_congestus_top_radius       = clouds_cumulus_congestus_radius + clouds_cumulus_congestus_thickness;
const float clouds_cumulus_congestus_distance         = CLOUDS_CUMULUS_CONGESTUS_DISTANCE;
const float clouds_cumulus_congestus_end_distance     = CLOUDS_CUMULUS_CONGESTUS_END_DISTANCE;
float clouds_cumulus_congestus_extinction_coeff       = 0.08;
float clouds_cumulus_congestus_scattering_coeff       = clouds_cumulus_congestus_extinction_coeff * (1.0 - 0.2 * rainStrength);

const float clouds_cumulonimbus_radius                = planet_radius + CLOUDS_CUMULONIMBUS_ALTITUDE;
const float clouds_cumulonimbus_thickness             = /*CLOUDS_CUMULONIMBUS_ALTITUDE*/ 1146 * CLOUDS_CUMULONIMBUS_THICKNESS;
const float clouds_cumulonimbus_top_radius            = clouds_cumulonimbus_radius + clouds_cumulonimbus_thickness * CLOUDS_CUMULONIMBUS_END_DISTANCE * 0.00004;
float clouds_cumulonimbus_distance                    = mix(CLOUDS_CUMULONIMBUS_DISTANCE, 0.0, rainStrength);
const float clouds_cumulonimbus_end_distance          = CLOUDS_CUMULONIMBUS_END_DISTANCE;
const float clouds_cumulonimbus_blend_distance        = 6.0;
const float clouds_cumulonimbus_extinction_coeff      = 0.05 * CLOUDS_CUMULONIMBUS_DENSITY;
float clouds_cumulonimbus_scattering_coeff            = clouds_cumulonimbus_extinction_coeff * (1.0 - 0.33 * rainStrength);

const float clouds_thunderhead_radius                 = planet_radius + CLOUDS_THUNDERHEAD_ALTITUDE;
const float clouds_thunderhead_thickness              = CLOUDS_THUNDERHEAD_ALTITUDE * CLOUDS_THUNDERHEAD_THICKNESS;
const float clouds_thunderhead_top_radius             = clouds_thunderhead_radius + clouds_thunderhead_thickness * 25.0;

const float clouds_towering_cumulus_radius            = planet_radius + CLOUDS_TOWERING_CUMULUS_ALTITUDE * 2.3;
const float clouds_towering_cumulus_thickness         = CLOUDS_TOWERING_CUMULUS_ALTITUDE * CLOUDS_TOWERING_CUMULUS_THICKNESS;
const float clouds_towering_cumulus_top_radius        = clouds_towering_cumulus_radius + clouds_towering_cumulus_thickness * 10.0;

// GAMS Clouds - Old daily weather variables set in Settings, default use D0 values
vec2 clouds_towering_cumulus_coverage                 = vec2(CLOUDS_TOWERING_CUMULUS_MIN, CLOUDS_TOWERING_CUMULUS_MAX);
vec2 clouds_thunderhead_coverage                      = vec2(CLOUDS_THUNDERHEAD_MIN, CLOUDS_THUNDERHEAD_MAX);

#endif // INCLUDE_SKY_CLOUDS_CONSTANTS
