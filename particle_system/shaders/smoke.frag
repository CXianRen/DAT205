#version 420

// required by GLSL spec Sect 4.5.3 (though nvidia does not, amd does)
precision highp float;

in vec3 vertexWorldSpacePos;    // world space position
in vec3 vertexTextureSpacePos;  // texture space position

uniform vec3  worldSpaceLightPosition;
uniform float pointLightIntensity;
uniform vec3  worldSpaceCameraPosition;

layout(binding = 0) uniform sampler3D densityTex;


layout(location = 0) out vec4 fragmentColor;


void main()
{
 
    float numSamples = 64;
    float numLightSamples = numSamples;

    float scale = 2.0 / 64;
    float lscale = scale;
 
    vec3 pos = vertexTextureSpacePos; 
    vec3 wpos = vertexWorldSpacePos.xyz; 

    vec3 eyeDir = normalize(wpos - worldSpaceCameraPosition)*scale;
 
    // transmittance
    float T = 1.0;
    // in-scattered radiance
    vec3 Lo = vec3(0.0);
    for (int i=0; i < numSamples; ++i)
    {
        // sample density
        float density = texture(densityTex, pos).x;
        // skip empty space
        if (density > 0.0)
        {
			// T = 0.0;
			// break;
            
			// attenuate ray-throughput
            T *= 1.0 - density*scale*5;
            if (T <= 0.01)
            {
                break;
            }
            // // point light dir in texture space
            vec3 lightDir = normalize(worldSpaceLightPosition-wpos)*lscale;
 
            // sample light
            float Tl = 1.0; // transmittance along light ray
            vec3 lpos = pos + lightDir;
 
            for (int s=0; s < numLightSamples; ++s)
            {
                float ld = texture(densityTex, lpos).x;
                Tl *= 1.0-5.0f*lscale*ld;
 
                if (Tl <= 0.01)
                    break;
 
                lpos += lightDir;
            }
 
            vec3 Li = pointLightIntensity * vec3(1.0f) * Tl;
 
            Lo += Li*T*density*scale;
        }
 
        pos += eyeDir;
        wpos += eyeDir;
    }

	fragmentColor = vec4(Lo,1.0-T);
}
