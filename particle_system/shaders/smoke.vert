#version 420

layout(location = 0) in vec3 position;

uniform mat4 modelMatrix;
uniform mat4 modelViewProjectionMatrix;

// output variables
out vec3 vertexWorldSpacePos;     
out vec3 vertexTextureSpacePos; 

void main()
{
	// first transform to the center of the cube

	//  stretch the cube in y direction 
  	// pos = vec3(1.0,1.0,1.0) * pos.xyz;

	gl_Position = modelViewProjectionMatrix * vec4(position, 1.0);
	vertexWorldSpacePos = (modelMatrix * vec4(position, 1.0)).xyz;
	vertexTextureSpacePos = position;
}
