#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

in vec3 vertexWorldSpacePos;    // world space position
in vec3 vertexTextureSpacePos;  // texture space position

uniform vec3 worldSpaceLightPosition;
uniform float pointLightIntensity;
uniform vec3 worldSpaceCameraPosition;
uniform float factor;
uniform int enable_light_tracing;

layout(binding = 0) uniform sampler3D densityTex;
layout(binding = 1) uniform sampler3D occupiedTex;

layout(location = 0) out vec4 fragmentColor;

int numSamples = 50;
float step = 0.02;
float modelScale = 20;

// parallel light, from the right
vec3 light_dir;

float getTray(vec3 pos, vec3 dir) {
    float Tl = 1.0;
    for(int s = 0; s < numSamples; ++s) {
        pos += dir * step;
        float density = texture(densityTex, pos).x;
        if (density < 0.01)
            continue;

        Tl *= exp(-factor * density * step);
        if(Tl <= 0.01)
            break;
    }
    return Tl;
}

void main() {

    vec3 pos = vertexTextureSpacePos.xyz;

    vec3 wpos = vertexWorldSpacePos.xyz;

    vec3 eyeDir = normalize(wpos - worldSpaceCameraPosition);
    // vec3 eyeDir = vec3(0.0, 0.0, -1.0);

    float T = 1.0;
    vec3 Lo = vec3(0.0);
    for(int i = 0; i < numSamples; ++i) {
        float density = texture(densityTex, pos).x;
        // skip empty space
        if(density > 0.01) {
            T *= exp(-factor * density * step);
            if(T <= 0.01)
                break;

            light_dir = normalize(worldSpaceLightPosition - wpos);
            float dist = length(worldSpaceLightPosition - wpos);

            float Tray = getTray(pos, light_dir);

            float Tvox = exp(-factor * density * step);

            float Lvox = pointLightIntensity * (1 - Tvox) * Tray/ dist / dist;

            Lo += Lvox * T;
        }

        pos += eyeDir * step;
        wpos += modelScale * eyeDir * step;
    }

    fragmentColor = vec4(Lo, 1-T);
}
