#version 330

layout (location = 0) in vec4 Position; //VERTICES
layout (location = 5) in vec3 Size; //ATTRIBUTE_0
layout (location = 4) in vec3 Color; //Color

out vec3 pass_voxelSize; 
out vec3 pass_color; 

void main() 
{
    gl_Position = Position;
	pass_voxelSize = Size; 
	pass_color = Color;
}