#version 330 core

in vec3 passPosition;
in vec3 passNormal;
in vec2 passTcoord;

uniform mat4 viewMatrix;

uniform struct 
{
	vec3 kd, ks;
	float kt;
    vec3 diffColor;
	vec3 specColor;
    float shininess;
} mat;

uniform struct LIGHT
{
  vec4 pos;    //pos.w = 0 dir. light, pos.w = 1 point light
  vec3 col;
  vec3 spot_direction; 
  float spot_exponent; 
  float spot_cutoff;  // no spotlight if cutoff = 0
} light[8]; //MAX_LIGHTS

uniform struct FogParams
{ 
   vec3 col; // Fog color
   float start; // This is only for linear fog
   float end; // This is only for linear fog
   float density; // For exp and exp2 equation  
   int mode; // 1 = linear, 2 = exp, 3 = exp2
} fog;

uniform int numLights;

uniform vec3 lightAmbient;

uniform int useColorTexture;
uniform sampler2D colorTexture; 

out vec4 fragmentColor;

float getFogFactor(FogParams fog, float z) 
{ 
   float f = 0.0; 
   if (fog.mode == 1) 
      f = (fog.end-z)/(fog.end-fog.start); 
   else if (fog.mode == 2) 
      f = exp(-fog.density*z); 
   else if( fog.mode == 3) 
      f = exp(-pow(fog.density*z, 2.0)); 
       
   f = clamp(f, 0.0, 1.0); 
     
   return f; 
}

void main() 
{ 
	vec3 lightVector;
	float spot;
	vec3 diffuse_color;
	float diffuse_alpha;

	// diffuse color 
	if (useColorTexture != 0)
	{
        diffuse_color = texture( colorTexture, passTcoord ).rgb;
        //diffuse_alpha = texture( colorTexture, passTcoord ).a;
		diffuse_alpha = 1.f - mat.kt;
    }

	else
	{
        diffuse_color = mat.diffColor;
        diffuse_alpha = 1.f - mat.kt;
    }

	//ambient light
	fragmentColor.rgb = mat.kd * diffuse_color * lightAmbient;
	
	for ( int i = 0; i < numLights; i++)
	{
        // Diffuse
		vec3 light_camcoord = (viewMatrix * light[i].pos).xyz;
		if (light[i].pos.w > 0.001f)
			lightVector = normalize( light_camcoord - passPosition);
        else
			lightVector = normalize(light_camcoord);
        float cos_phi = max( dot( passNormal, lightVector), 0.0f);
        
		// Specular
		vec3 eye = normalize( -passPosition);
        vec3 reflection = normalize( reflect( -lightVector, passNormal));
        float cos_psi_n = pow( max( dot( reflection, eye), 0.0f), mat.shininess);
        
		// Spotlight	
		if (light[i].spot_cutoff < 0.001f)
			spot = 1.0;
        else
		{
            float cos_phi_spot = max( dot( -lightVector, mat3(viewMatrix) * light[i].spot_direction), 0.0f);
            if( cos_phi_spot >= cos( light[i].spot_cutoff))
				spot = pow( cos_phi_spot, light[i].spot_exponent);
            else
				spot = 0.0f;
        }

		// All together 
		fragmentColor.rgb += mat.kd * spot * diffuse_color * cos_phi * light[i].col;
		fragmentColor.rgb += mat.ks * spot * mat.specColor * cos_psi_n * light[i].col;
	}

	//Add fog
	if (fog.mode != 0)
	{
		float f = getFogFactor(fog, length( passPosition));
        fragmentColor.rgb = f * fragmentColor.rgb + (1-f) * fog.col;
	}

    fragmentColor.a = diffuse_alpha;
}