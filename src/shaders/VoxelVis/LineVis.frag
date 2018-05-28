#version 330 core

in vec3 pass_color;
out vec4 fragmentColor;

void main()
{
    fragmentColor = vec4( pass_color, 1.f);
}