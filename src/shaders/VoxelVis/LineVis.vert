#version 330

layout (location = 0) in vec4 Position; //VERTICES
layout (location = 4) in vec3 Color; //Color

out vec3 pass_color; 

uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

void main() 
{
    gl_Position = projectionMatrix * viewMatrix * Position;
    pass_color = Color;
}