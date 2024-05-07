#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

///////////////////////////////////////////////////////////////////////////////
// Material
///////////////////////////////////////////////////////////////////////////////
uniform vec3 material_color;
uniform float material_reflectivity;
uniform float material_metalness;
uniform float material_fresnel;
uniform float material_shininess;
uniform float material_emission;

uniform int has_emission_texture;
uniform int has_color_texture;
layout(binding = 0) uniform sampler2D colorMap;
layout(binding = 5) uniform sampler2D emissiveMap;

///////////////////////////////////////////////////////////////////////////////
// Environment
///////////////////////////////////////////////////////////////////////////////
layout(binding = 6) uniform sampler2D environmentMap;
layout(binding = 7) uniform sampler2D irradianceMap;
layout(binding = 8) uniform sampler2D reflectionMap;
uniform float environment_multiplier;

///////////////////////////////////////////////////////////////////////////////
// Light source
///////////////////////////////////////////////////////////////////////////////
uniform vec3 point_light_color = vec3(1.0, 1.0, 1.0);
uniform float point_light_intensity_multiplier = 50.0;

///////////////////////////////////////////////////////////////////////////////
// Constants
///////////////////////////////////////////////////////////////////////////////
#define PI 3.14159265359

///////////////////////////////////////////////////////////////////////////////
// Input varyings from vertex shader
///////////////////////////////////////////////////////////////////////////////
in vec2 texCoord;
in vec3 viewSpaceNormal;
in vec3 viewSpacePosition;

///////////////////////////////////////////////////////////////////////////////
// Input uniform variables
///////////////////////////////////////////////////////////////////////////////
uniform mat4 viewInverse;
uniform vec3 viewSpaceLightPosition;

///////////////////////////////////////////////////////////////////////////////
// Output color
///////////////////////////////////////////////////////////////////////////////
layout(location = 0) out vec4 fragmentColor;



vec3 calculateDirectIllumiunation(vec3 wo, vec3 n)
{
	vec3 wi = normalize(viewSpaceLightPosition - viewSpacePosition);
	float wi_dot_n = dot(n, wi);
	if(wi_dot_n < 0.0)
	{
		return vec3(0.0);
	}
	float d = length(viewSpaceLightPosition - viewSpacePosition);
	vec3 Li = point_light_color * point_light_intensity_multiplier / (d * d);
	vec3 diffuse = material_color * wi_dot_n * Li * 1.0/PI;	

    // brdf
	vec3 wh = normalize(wi + wo);
	float wh_dot_wi = dot(wh, wi);
	float F_wi = material_fresnel + (1.0 - material_fresnel) * pow(1.0 - wh_dot_wi, 5.0);
	float wh_dot_n = dot(wh, n);
	float D_wi = (material_shininess + 2.0) / (2.0 * PI) * pow(wh_dot_n, material_shininess);
	float wo_dot_n = dot(wo, n);
	float wh_dot_wo = dot(wh, wo);
	float G_wi = 2.0 * wh_dot_n * wi_dot_n / wh_dot_wo;
	float G_wo = 2.0 * wh_dot_n * wo_dot_n / wh_dot_wo;
	float G = min(1.0, min(G_wi, G_wo));
	float denominator = 4.0 * clamp(wi_dot_n * wo_dot_n, 0.0001, 1.0);
	float brdf = F_wi * D_wi * G / denominator;

	vec3 dielectric_term = brdf * wi_dot_n * Li + (1.0 - F_wi) * diffuse;
	vec3 metal_term = brdf * material_color * wi_dot_n * Li;
	vec3 microfacet_term = material_metalness * metal_term + (1.0 - material_metalness) * dielectric_term;
	return microfacet_term;
}

vec3 calculateIndirectIllumination(vec3 wo, vec3 n)
{
	vec3 world_normal = vec3(viewInverse * vec4(n, 0.0));

	float theta = acos(max(-1.0f, min(1.0f, world_normal.y)));
	float phi = atan(world_normal.z, world_normal.x);
	if(phi < 0.0f)
		phi = phi + 2.0f * PI;
	vec2 lookup = vec2(phi / (2.0 * PI), theta / PI);
	vec3 Li = environment_multiplier * texture(irradianceMap, lookup).rgb;
	vec3 diffuse_term = material_color * (1.0 / PI) * Li;

	vec3 wi = normalize(reflect(-wo, n));
	vec3 wr = normalize(vec3(viewInverse * vec4(wi, 0.0)));
	theta = acos(max(-1.0f, min(1.0f, wr.y)));
	phi = atan(wr.z, wr.x);
	if(phi < 0.0f)
		phi = phi + 2.0f * PI;
	lookup = vec2(phi / (2.0 * PI), theta / PI);
	float roughness = sqrt(sqrt(2.0 / (material_shininess + 2.0)));
	Li = environment_multiplier * textureLod(reflectionMap, lookup, roughness * 7.0).rgb;
	vec3 wh = normalize(wi + wo);
	float wodotwh = max(0.0, dot(wo, wh));
	float F = material_fresnel + (1.0 - material_fresnel) * pow(1.0 - wodotwh, 5.0);
	vec3 dielectric_term = F * Li + (1.0 - F) * diffuse_term;
	vec3 metal_term = F * material_color * Li;

	vec3 microfacet_term = material_metalness * metal_term + (1.0 - material_metalness) * dielectric_term;

	return microfacet_term;

}

void main()
{
	float visibility = 1.0;
	float attenuation = 1.0;

	vec3 wo = -normalize(viewSpacePosition);
	vec3 n = normalize(viewSpaceNormal);

	// Direct illumination
	vec3 direct_illumination_term = visibility * calculateDirectIllumiunation(wo, n);

	// Indirect illumination
	vec3 indirect_illumination_term = calculateIndirectIllumination(wo, n);

	///////////////////////////////////////////////////////////////////////////
	// Add emissive term. If emissive texture exists, sample this term.
	///////////////////////////////////////////////////////////////////////////
	vec3 emission_term = material_emission * material_color;
	if(has_emission_texture == 1)
	{
		emission_term = texture(emissiveMap, texCoord).xyz;
	}

	vec3 shading = direct_illumination_term + indirect_illumination_term + emission_term;

	fragmentColor = vec4(shading, 1.0);
	return;
}
