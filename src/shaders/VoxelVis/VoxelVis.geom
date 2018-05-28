#version 330 

layout (points) in;
layout (line_strip, max_vertices = 24) out;

in vec3 pass_voxelSize[];
in vec3 pass_color[];

out vec4 color;

uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;

void main(void)
  {
    // Normals for the triangle vertices    
    mat4 VP = projectionMatrix * viewMatrix;
	vec4 tmp_color = vec4( pass_color[0], 1.0f);

	vec4 p = gl_in[0].gl_Position;
	vec3 size = pass_voxelSize[0];
	vec4 A = VP * ( vec4( -size.x,  size.y,  size.z, 0.f) + p);
    vec4 B = VP * ( vec4( -size.x, -size.y,  size.z, 0.f) + p);
    vec4 C = VP * ( vec4(  size.x, -size.y,  size.z, 0.f) + p);
    vec4 D = VP * ( vec4(  size.x,  size.y,  size.z, 0.f) + p);
    vec4 E = VP * ( vec4(  size.x,  size.y, -size.z, 0.f) + p);
    vec4 F = VP * ( vec4(  size.x, -size.y, -size.z, 0.f) + p);
    vec4 G = VP * ( vec4( -size.x, -size.y, -size.z, 0.f) + p);
    vec4 H = VP * ( vec4( -size.x,  size.y, -size.z, 0.f) + p);

	gl_Position = A; color = tmp_color; EmitVertex ();
	gl_Position = B; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = B; color = tmp_color; EmitVertex ();
	gl_Position = C; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = C; color = tmp_color; EmitVertex ();
	gl_Position = D; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = D; color = tmp_color; EmitVertex ();
	gl_Position = A; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = A; color = tmp_color; EmitVertex ();
	gl_Position = H; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = D; color = tmp_color; EmitVertex ();
	gl_Position = E; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = C; color = tmp_color; EmitVertex ();
	gl_Position = F; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = B; color = tmp_color; EmitVertex ();
	gl_Position = G; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = E; color = tmp_color; EmitVertex ();
	gl_Position = F; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = F; color = tmp_color; EmitVertex ();
	gl_Position = G; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = G; color = tmp_color; EmitVertex ();
	gl_Position = H; color = tmp_color; EmitVertex ();
	EndPrimitive ();

	gl_Position = H; color = tmp_color; EmitVertex ();
	gl_Position = E; color = tmp_color; EmitVertex ();
	EndPrimitive ();
}
