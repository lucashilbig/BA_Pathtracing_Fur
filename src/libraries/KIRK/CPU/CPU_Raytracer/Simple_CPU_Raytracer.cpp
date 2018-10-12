#include "Simple_CPU_Raytracer.h"

KIRK::CPU::SimpleCPURaytracer::SimpleCPURaytracer() : KIRK::GuiElement()
{
	m_num_supersamples = 3;
	m_num_blursamples = 5;
	m_num_lightsamples = 4;
	m_num_poissondisksamples = 2;
	m_adaptive_depth = 2;
	m_max_adaptive_difference = 0.5f;

	m_num_threads = KIRK::ThreadManager::maxThreads(); //Maximum thread count
};

KIRK::CPU::SimpleCPURaytracer::SimpleCPURaytracer(std::shared_ptr<KIRK::CPU::Scene> cpuscene, int depth) : SimpleCPURaytracer()
{
	SimpleCPURaytracer::init(cpuscene);
}

KIRK::CPU::SimpleCPURaytracer::~SimpleCPURaytracer()
{
}

void KIRK::CPU::SimpleCPURaytracer::render()
{

	KIRK::ThreadManager::for_loop_double(0, 0, m_render_texture->getSize().x, m_render_texture->getSize().y, [=](int x, int y)
	{

		KIRK::Color::RGBA color;
		if (isEnabled(RTFLAG_USE_SUPERSAMPLING))
		{
			color = superSampling(x, y);
		}
		else if (isEnabled(RTFLAG_USE_ADAPTIVE_SAMPLING))
		{
			color = adaptiveSampling(x, y);
		}
		else if (isEnabled(RTFLAG_USE_POISSONDISK_SAMPLING))
		{
			color = poissonDiskSampling(x, y);
		}
		else
		{
			KIRK::Ray ray = m_cpuscene->getActiveCamera().getRayFromPixel(x, y);
			color = trace(ray, 0, 1.0);
		}

		m_render_texture->setPixel(x, y, color);

	}, m_num_threads, isEnabled(RTFLAG_USE_PROGRESSBAR));
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::trace(Ray &ray, int level, float weight)
{
	//depth of field
	if (level < 1 && isEnabled(RTFLAG_USE_DOF))
	{
		return depthOfField(level, weight, ray.m_direction);
	}

	//Check for an intersection.
	Intersection hit(ray);
	if (m_cpuscene->getDataStructure().closestIntersection(&hit))
	{
		//if we have marschnerHairBSDF we use the specific shade method
		if (hit.m_object->getMaterial()->m_current_bsdf == 6)
			return shadeMarschnerHair(hit, level, weight);
		else//we use standart shade method
			return shade(hit, level, weight);
	}

	// Background color
	return (/*Show in background*/true || level != 0) ? m_cpuscene->getEnvironment().getColor(ray.m_direction) : KIRK::Color::BLACK;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::lightShading(const glm::vec3 &pos, const glm::vec3 &norm, const glm::vec2 &texcoord, const glm::vec3 &view, KIRK::Material *material,
	const KIRK::Color::RGBA &diff_color)
{
	KIRK::Color::RGBA color = KIRK::Color::RGBA(0.0f);

	// 2# ambient part of the illumination -> "global illumination"
	KIRK::Color::RGBA ambient_color(0.0f);
	ambient_color += m_cpuscene->getEnvironment().getAmbientLight()  * diff_color;

	// The following value is kind of a cheaty way to remove one if-clause.
	// if useSoftShadows is true, m_lightsamples will stay unchanged
	// otherwise, it will be clamped to 1 sample, because we don't need more.
	float clamped_samples = (m_num_lightsamples - 1)*isEnabled(RTFLAG_USE_SOFT_SHADOWS) + 1;
	float light_count = m_cpuscene->getLights().size();

	// recompute normal to consider both sides of polygons correctly in lighting
	float m_dot = (glm::dot(norm, view));
	glm::vec3 norm_view = norm;
	if (m_dot >= 1e-5f)
	{
		//Flip towards view
		norm_view = -norm;
	}

	//direct illumination for each lightsource
	std::vector<KIRK::Color::RGBA> directIllum(light_count, KIRK::Color::RGBA(0.0f));
	std::vector<float> attenuations(light_count, 0.0f);

	for (unsigned int i = 0; i < light_count; i++)
	{

		Ray lightToHit = m_cpuscene->getLights()[i]->calcLightdir(pos, attenuations[i], false);

		if (attenuations[i] > 0.0f && glm::dot(norm_view, glm::normalize(lightToHit.m_direction)) >= 0.0f)
		{
			glm::vec3 n_lightdir = glm::normalize(lightToHit.m_direction);
			// angle between normalvector and light direction
			float cos_phi = glm::max(glm::dot(norm, n_lightdir), 0.0f);

			// Kd value * angle * material color * illumination from source
			directIllum[i] += cos_phi * diff_color * m_cpuscene->getLights()[i]->m_color * attenuations[i];

			// specular shading (phong)
			glm::vec3 reflection = glm::reflect(n_lightdir, norm);
			float cos_psi_n = pow(glm::max(glm::dot(reflection, view), 0.0f), material->shininess());
			directIllum[i] += material->fetchParameterFloat<MatParamType::REFLECTIVITY>(texcoord) * cos_psi_n * material->fetchParameterColor<MatParamType::SPECULAR>(texcoord) * m_cpuscene->getLights()[i]->m_color * attenuations[i];
		}
	}

	// shadow weight for each lightsource and ambient light (AO)
	std::vector<int> shadowWeight(light_count, 0);
	int ambientShadowWeight = 0;
	for (unsigned int sample = 0; sample < clamped_samples; sample++)
	{
		// ambient occlusion
		if (isEnabled(RTFLAG_USE_SOFT_SHADOWS))
		{

			Ray ambient_occlusion_test(pos + norm_view * 1e-3f, norm_view);
			ambient_occlusion_test.jitterBy(3);

			// shadow test
			if (m_cpuscene->getDataStructure().isIntersection(&ambient_occlusion_test, 0.4f))
			{
				ambientShadowWeight++;
			}
		}

		for (unsigned int i = 0; i < light_count; i++)
		{
			// If there is no light here we can skip shadow calculation
			if (attenuations[i] <= 0.0f)
			{
				continue;
			}

			float bias = 1e-2f;
			float attenuation; // unused here

			glm::vec3 origin = pos + bias * norm_view;
			Ray lightToHit = m_cpuscene->getLights()[i]->calcLightdir(origin, attenuation, isEnabled(RTFLAG_USE_SOFT_SHADOWS));

			// shadow test
			if (glm::dot(norm_view, glm::normalize(lightToHit.m_direction)) < 0 || m_cpuscene->getDataStructure().isIntersection(&lightToHit, 1.0f))
			{
				shadowWeight[i]++;
			}
		}
	}
	// calculate final ambient lighting after AO
	float ambient_fac = static_cast<float>(ambientShadowWeight) / clamped_samples;
	color += (1 - ambient_fac) * ambient_color;

	// calculate final direct lighting for each lightsource considering (soft) shadows
	for (int i = 0; i < light_count; i++)
	{
		float fac = static_cast<float>(shadowWeight[i]) / clamped_samples;
		color += (1 - fac) * directIllum[i];
	}

	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::reflection(const glm::vec3 &pos, const glm::vec3 &view, const glm::vec3 &norm, float falloff, float roughness, int level, float weight)
{
	glm::vec3 specdir, origin;
	KIRK::Color::RGBA color;

	float m_dot = (glm::dot(norm, view));
	glm::vec3 norm_view = norm;
	if (std::abs(m_dot) >= 1e-5f)
	{
		// flip normal towards view but still parallel to the old normal.
		norm_view = -glm::normalize(m_dot*norm);
	}
	// get the diretion of the reflecting ray
	specdir = glm::normalize(glm::reflect(view, norm_view));
	// place it directly in front of our hit
	// We'll take a fairly large epsilon to compensate for inaccurate smoothed normals on flat triangles
	origin = pos + 1e-2 * norm_view;

	// create a ray with the calculated origin and direction
	Ray reflect(origin, specdir);
	reflect.jitterBy(roughness);
	// trace it and add the color if we hit something
	color = weight * trace(reflect, level + 1, falloff * weight);

	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::refraction(Intersection &hit, const glm::vec3 &pos, const glm::vec3 &view, const glm::vec3 &norm, float falloff, float roughness,
	float ior, int level, float weight)
{
	glm::vec3 transdir, origin;
	KIRK::Color::RGBA color(0.0f);

	if (hit.m_enter)
		transdir = glm::refract(view, norm, 1.0f / ior);
	else
		transdir = glm::refract(view, -norm, ior);

	if (transdir != glm::vec3(0) && !(std::isnan(transdir.x))) //we got no total reflection
	{
		transdir = glm::normalize(transdir);
		origin = pos + cRayEpsilon * transdir;
		Ray refract(origin, transdir);
		refract.jitterBy(roughness);
		color = weight * trace(refract, level + 1, falloff * weight);
	}
	else //we got total reflection
	{
		color = reflection(pos, view, norm, falloff, roughness, level, weight);
	}

	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::depthOfField(int level, float weight, const glm::vec3 &view)
{
	KIRK::Color::RGBA color(0.f);

	// We can add separate DOF-samples to accelerate simple DOF applications
	for (int di = 0; di < m_num_blursamples; di++)
	{
		//KIRK::Camera::PixelRay dof_ray = getCamera()->transformToDof(view);
		Ray nray = m_cpuscene->getActiveCamera().transformToDof(view);
		color += trace(nray, level + 1, weight) / (float)(m_num_blursamples);
	}
	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::superSampling(int x, int y)
{
	KIRK::Color::RGBA color(0);
	float sampleWeight = m_num_supersamples * m_num_supersamples;
	float steps_per_pixel = 10.f; 	float step = 1.0f / steps_per_pixel;
	for (int j = 0; j < m_num_supersamples; j++)
	{
		for (int i = 0; i < m_num_supersamples; i++)
		{
			// This approach just uniformly adds subpixel-samples at discrete places
			float delta_x = i / (double)m_num_supersamples;
			float delta_y = j / (double)m_num_supersamples;
			//"Squeeze" deltas into pixel a bit.
			delta_x *= (steps_per_pixel - 2) / steps_per_pixel;
			delta_x += step;
			delta_y *= (steps_per_pixel - 2) / steps_per_pixel;
			delta_y += step;
			KIRK::Ray ray(m_cpuscene->getActiveCamera().getRayFromPixel(x, y, delta_x, delta_y));
			color += trace(ray, 0, 1.0) / sampleWeight;
		}
	}
	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::adaptiveSampling(int x, int y)
{
	// Avoid information loss due to redundant calculation of border samples
	float step = 1.0f / 20.0f;
	float x1 = x + step;
	float y1 = y + step;
	float x2 = x + 1 - step;
	float y2 = y + 1 - step;
	// 4 corners of pixel, not the actual corners because of step (on purpose)
	Ray ray1 = m_cpuscene->getActiveCamera().getRayFromPixel(x1, y2);
	Ray ray2 = m_cpuscene->getActiveCamera().getRayFromPixel(x2, y2);
	Ray ray3 = m_cpuscene->getActiveCamera().getRayFromPixel(x1, y1);
	Ray ray4 = m_cpuscene->getActiveCamera().getRayFromPixel(x2, y1);
	Color::RGBA color1 = trace(ray1, 0, 1.0);
	Color::RGBA color2 = trace(ray2, 0, 1.0);
	Color::RGBA color3 = trace(ray3, 0, 1.0);
	Color::RGBA color4 = trace(ray4, 0, 1.0);
	return adaptiveSamplingRecursive(ray1.m_direction, ray2.m_direction, ray3.m_direction, ray4.m_direction,
		color1, color2, color3, color4, ray1.m_origin, 0);
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::adaptiveSamplingRecursive(glm::vec3 &dir1, glm::vec3 &dir2, glm::vec3 &dir3, glm::vec3 &dir4,
	Color::RGBA color1, Color::RGBA color2, Color::RGBA color3, Color::RGBA color4, glm::vec3 &position, int depth)
{

	if (depth < m_adaptive_depth)
	{
		float diff1 = glm::length(color1 - color2);
		float diff2 = glm::length(color1 - color3);
		float diff3 = glm::length(color1 - color4);
		float diff4 = glm::length(color2 - color3);
		float diff5 = glm::length(color2 - color4);
		float diff6 = glm::length(color3 - color4);
		// check if there are rapid changes between any two of the four samples
		if (diff1 > m_max_adaptive_difference
			|| diff2 > m_max_adaptive_difference
			|| diff3 > m_max_adaptive_difference
			|| diff4 > m_max_adaptive_difference
			|| diff5 > m_max_adaptive_difference
			|| diff6 > m_max_adaptive_difference)
		{
			// generate new samples and get their color
			glm::vec3 n1 = (dir1 + dir2) / 2.0f;
			glm::vec3 n2 = (dir1 + dir3) / 2.0f;
			glm::vec3 n3 = (dir1 + dir4) / 2.0f;
			glm::vec3 n4 = (dir2 + dir4) / 2.0f;
			glm::vec3 n5 = (dir3 + dir4) / 2.0f;
			Ray rn1(position, n1);
			Ray rn2(position, n2);
			Ray rn3(position, n3);
			Ray rn4(position, n4);
			Ray rn5(position, n5);
			Color::RGBA cn1 = trace(rn1, 0, 1.0);
			Color::RGBA cn2 = trace(rn2, 0, 1.0);
			Color::RGBA cn3 = trace(rn3, 0, 1.0);
			Color::RGBA cn4 = trace(rn4, 0, 1.0);
			Color::RGBA cn5 = trace(rn5, 0, 1.0);
			// continue recursively for the subdivided square
			color1 = adaptiveSamplingRecursive(dir1, n1, n2, n3, color1, cn1, cn2, cn3, position, depth + 1);
			color2 = adaptiveSamplingRecursive(n1, dir2, n3, n4, cn1, color2, cn3, cn4, position, depth + 1);
			color3 = adaptiveSamplingRecursive(n2, n3, dir3, n5, cn2, cn3, color3, cn5, position, depth + 1);
			color4 = adaptiveSamplingRecursive(n3, n4, n5, dir4, cn3, cn4, cn5, color4, position, depth + 1);
		}
	}
	return (color1 + color2 + color3 + color4) / 4.0f;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::poissonDiskSampling(int x, int y)
{
	KIRK::Color::RGBA color(0);

	for (const auto& p : m_poissonDisks[m_num_poissondisksamples - 2])
	{
		Ray ray = m_cpuscene->getActiveCamera().getRayFromPixel(x, y, p.x, p.y);

		color = color + (trace(ray, 0, 1.0) / static_cast<float>(m_num_poissondisksamples));
	}
	return color;
}


KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::shade(Intersection &hit, int level, float weight)
{
	// get the objected we hit. static_cast cause we checked correct object type in render method already
	Triangle *obj = static_cast<KIRK::Triangle*>(hit.m_object);

	// calculation of the normal
	obj->calcNormal(&hit);
	obj->calcTcoord(&hit);

	// get the position, normal, and direction of the intersection
	glm::vec3 pos = hit.m_location;
	glm::vec3 norm = hit.m_normal;
	glm::vec3 view = glm::normalize(hit.m_ray.m_direction);
	glm::vec2 texcoord = hit.m_texcoord;

	// Material
	KIRK::Material current_material = *obj->getMaterial();//*m_scene->m_materials[0];// = KIRK::CVKConverterUtils::convertMaterial(obj->getMaterial());

	KIRK::Color::RGBA color(0.0f);

	////////////////////////////////////////////////////////////////////
	//  HOW THE FOLLOWING STEPS WORK:
	//  1) We use the diffuse and phong shaded color as a base.
	//  2) If there is transparency, then we apply a refractive transparent
	//     shading onto the color.
	//  3) Finally, we add a fresnel oriented reflection component.
	////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////
	//  DIFFUSION:
	//  We just use a simple Phong Shading model to add onto the default
	//  diffuse color.
	////////////////////////////////////////////////////////////////////
	// direct illumination for each lightsource
	// Simple diffuse and specular shading.
	// We don't distinguish between soft and direct shading anymore, that can be done in the method itself, so we can save quite some tens of lines of code.
	// In combination with diff_color: LambertianReflectionBSDF

	KIRK::Color::RGBA base_color = current_material.fetchParameterColor<KIRK::MatParamType::DIFFUSE>(hit.m_texcoord);//diffuseShading(hit);

	color = lightShading(pos, norm, texcoord, view, &current_material, base_color);

	float m_dot = (glm::dot(norm, view));
	glm::vec3 norm_view = norm;
	if (std::abs(m_dot) >= 1e-5f)
	{
		//Flip towards view
		norm_view = -glm::normalize(m_dot*norm);
	}

	float angle = glm::angle(-view, norm_view);
	if (angle > 180)
		angle = angle - glm::radians(180.f);
	if (angle > 90)
		angle = angle - (glm::radians(90.f) - angle);
	float r_0 = glm::pow((1 - 1.56f) / (1 + 1.56f), 2);
	float r_theta = r_0 + (1 - r_0)*glm::pow((1 - glm::cos(angle)), 5);

	float reflectivity = current_material.fetchParameterFloat<KIRK::MatParamType::REFLECTIVITY>(texcoord);
	float transparency = current_material.fetchParameterFloat<KIRK::MatParamType::TRANSPARENCY>(texcoord);
	float roughness = current_material.fetchParameterFloat<KIRK::MatParamType::ROUGHNESS>(texcoord);
	Color::RGBA volume = current_material.fetchParameterColor<KIRK::MatParamType::VOLUME>(texcoord);
	Color::RGBA specular = current_material.fetchParameterColor<KIRK::MatParamType::SPECULAR>(texcoord);

	float fresnel = glm::clamp(reflectivity*reflectivity - (transparency*transparency) + r_theta * (reflectivity), 0.f, 1.f);

	//return KIRK::Color::color_rgba(fresnel);

	//return KIRK::Color::color_rgba(fresnel);

	////////////////////////////////////////////////////////////////////
	//  REFRACTION:
	//  The color will be mixed with another adapted fresnel
	//  term, now dependant on the base transparency.
	////////////////////////////////////////////////////////////////////
	float fresnel_transparency = transparency * (1 - fresnel);
	if ((fresnel_transparency * weight > Minweight) && (level < m_depth) && isEnabled(RTFLAG_USE_REFRACTIONS))
	{
		color = KIRK::Color::mix(color, 1.f*volume*refraction(hit, pos, view, norm, fresnel_transparency, roughness, current_material.m_ior, level, 1), transparency);
	}

	////////////////////////////////////////////////////////////////////
	//  REFLECTION:
	//  The diffuse shaded part will be mixed with an adapted
	//  fresnel value. This one will be offset by a base reflectivity.
	////////////////////////////////////////////////////////////////////
	if ((fresnel * weight > Minweight) && (level < m_depth) && isEnabled(RTFLAG_USE_REFLECTIONS))
	{
		color = KIRK::Color::mix(color, specular*reflection(pos, view, norm, fresnel, roughness, level, weight), fresnel);
	}

	return color;
}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::shadeMarschnerHair(Intersection &hit, int level, float weight)
{
	// get the objected we hit. static_cast cause we checked correct object type in render method already
	KIRK::Object *cylinder_obj = hit.m_object;

	// calculation of the normal
	cylinder_obj->calcNormal(&hit);
	cylinder_obj->calcTcoord(&hit);

	// get the position, normal, and direction of the intersection
	glm::vec3 pos = hit.m_location;
	glm::vec3 normal = hit.m_normal;
	glm::vec2 texcoord = hit.m_texcoord;
	Material* material = hit.m_object->getMaterial();
	glm::vec3 norm_in_ray = glm::normalize(hit.m_ray.m_direction);//normalized input ray

	/*******************************************************/
	/******************* LIGHT SAMPLING ********************/
	/*******************************************************/
	//Random Value between 0 and 1, to select light source randomly
	std::random_device rd;
	std::mt19937 gen = std::mt19937(rd());
	std::uniform_real_distribution<> dist_01{ 0, 1 };

	int lightIndex = dist_01(gen) * m_cpuscene->getLights().size();
	auto light = m_cpuscene->getLights()[lightIndex].get();

	float attenuation;
	KIRK::Ray hitToLightRay = light->calcLightdir(hit.m_location, attenuation, true);//calculates ray from intersection point towards light source
	/*******************************************************/  

	//calculation of tangent vector
	glm::vec3 c1 = glm::cross(normal, glm::vec3(0.0, 0.0, 1.0));
	glm::vec3 c2 = glm::cross(normal, glm::vec3(0.0, 1.0, 0.0));
	glm::vec3 tangent = (glm::length(c1) > glm::length(c2)) ? glm::normalize(c1) : glm::normalize(c2);

	//move light ray to cylinders (fibers) local space.
	/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model
	and we used the local_output_ray to store our ray towards the light source we want to sample*/
	glm::vec3 hit_to_light = Math::worldToLocal(glm::normalize(hitToLightRay.m_direction), cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

	//calculate spherical coordinates of light ray. NOTE: std::atan2 is atan(y/x) using the signs of arguments to determine the correct quadrant
	float theta_i = std::atan2(hit_to_light.x, hit_to_light.z);//Angle between light ray and fibers v-w(normal-)plane(Inclinations with respect of the fibers normalplane)
	float phi_i = std::atan2(std::hypot(hit_to_light.x, hit_to_light.z), hit_to_light.y);//azimuth angle. Around x(fibers u-)-axis
	
	//values of scattering function parts
	glm::vec3 scat_r(0.f);
	glm::vec3 scat_tt(0.f);
	glm::vec3 scat_trt(0.f);

	//random generator for lobe alpha shift and beta width
	std::uniform_real_distribution<float> dist(5.0f, std::nextafter(10.0f, DBL_MAX));//std::nextafter so we get the [5,10] interval instead of [5,10)
	float alpha_r = -1.0f * dist(gen);//longitudinal shift: R lobe. Suggested value from marschner hair paper between -10 and -5 degrees
	float beta_r = dist(gen);//longitudinal width (stdev.): R lobe. Suggested value from marschner hair paper between 5 and 10 degrees		

	////////////////////////////////////////////////////////////////////
	//  HOW THE FOLLOWING STEPS WORK:
	//  1) We trace a Ray for every Cylinder-Path. R, TT and TRT
	//  2) We sum the Scattering values from each Path to the final Scattering value.
	//  3) Finally, we calculate the color with the scattering value and the direct Light at the intersection point
	////////////////////////////////////////////////////////////////////	

	//////////////////////////////////////////
	////
	//// Cylinder Path: R
	////
	//////////////////////////////////////////
	{
		//reflect input ray on the surface
		glm::vec3 out_ray = glm::reflect(norm_in_ray, glm::faceforward(normal, norm_in_ray, normal));
		//rotate towards normal to account for the tilted fiber surface
		out_ray = glm::vec3(glm::vec4(out_ray, 0.f) * glm::rotate(glm::radians(2 * alpha_r), tangent));

		//move output ray to cylinders (fibers) local space.
		/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model*/
		glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(out_ray), cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

		//--------------------------------------------------------------------
		// M_r(theta_h) : Marschner marginal, longitudinal	scattering function(M)  			
		//--------------------------------------------------------------------

		//calculate parameters for M
		float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane
		float theta_h = (theta_r + theta_i) / 2;//theta half angle
		float theta_d = (theta_r - theta_i) / 2;//theta difference angle

		//M_r(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation. TODO: Try Logistic function as suggested in pixar paper
		float m_r = BSDFHelper::normal_gauss_pdf(theta_h - glm::radians(alpha_r), 0.0f, beta_r);

		//--------------------------------------------------------------------
		// N_r(phi) : Marschner conditional, azimuthal scattering function(N)
		//--------------------------------------------------------------------

		//calculate parameters for N
		float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
		float phi = phi_r - phi_i;
		float h = glm::sin(phi) * -0.5f;//root of the approximation for phi^ in marschner paper (Equation 10) with p = 0

		//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
		float dphi_dh = -2.f / glm::sqrt(1 - h * h);

		//help variables for bravais calculation
		float cos_theta_d = glm::cos(theta_d);
		float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
		//Bravais (virtual index of reflection eta_one, eta_two) calculation
		float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
		float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)

		//calculate attenuation factor with fresnel
		float fresnel = BSDFHelper::dialectricFresnel(glm::cos(glm::asin(h)), bravais_first, bravais_sec);

		//final term for N_r(phi). Marschner Equation 8
		float n_r = 0.5f * fresnel / glm::abs(2 * dphi_dh);
		
		//--------------------------------------------------------------------

		//Final value for the combined scattering function
		scat_r = glm::vec3(m_r * n_r / glm::pow2(cos_theta_d));
	}

	//////////////////////////////////////////
	////
	//// Cylinder Path: TT
	////
	//////////////////////////////////////////
	{
		//refract input ray on surface and put it into new intersection
		glm::vec3 t_dir = glm::refract(norm_in_ray, glm::faceforward(normal, norm_in_ray, normal), 1.0f / material->m_ior);
		glm::vec3 t_normal;
		Intersection t_hit(Ray(pos + 1e-4f * t_dir, t_dir));
		//get the intersection point on the second wall
		if (m_cpuscene->getDataStructure().closestIntersection(&t_hit))
		{
			cylinder_obj->calcNormal(&t_hit);// calculation of the normal
			t_normal = t_hit.m_normal;
		}
		else { t_normal = -normal; }

		//refract on the second wall of the cylinder to get output ray
		glm::vec3 out_ray = glm::refract(glm::normalize(t_dir), glm::faceforward(t_normal, glm::normalize(t_dir), t_normal), 1.f / material->m_ior);

		//calculation of tangent vector of second intersection
		glm::vec3 c3 = glm::cross(t_normal, glm::vec3(0.0, 0.0, 1.0));
		glm::vec3 c4 = glm::cross(t_normal, glm::vec3(0.0, 1.0, 0.0));
		glm::vec3 t_tangent = (glm::length(c3) > glm::length(c4)) ? glm::normalize(c3) : glm::normalize(c4);
		//rotate towards normal to account for the tilted fiber surface
		out_ray = glm::vec3(glm::vec4(out_ray, 0.f) * glm::rotate(glm::radians(-alpha_r / 2), t_tangent));

		//move output ray to cylinders (fibers) local space.
		/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model*/
		glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(out_ray), cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

		//--------------------------------------------------------------------
		// M_tt(theta_h) : Marschner marginal, longitudinal	scattering function(M)  
		//--------------------------------------------------------------------
		//calculate parameters for M
		float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane
		float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
		float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle
		float gaussian_x = theta_h - glm::radians((-alpha_r / 2.f));//longitudinal shift for TT-Path is smaller than for R-Path

		//M_tt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
		float m_tt = BSDFHelper::normal_gauss_pdf(gaussian_x, 0.0f, beta_r / 2.f);

		//--------------------------------------------------------------------
		// N_tt(phi) : Marschner conditional, azimuthal scattering function(N)
		//--------------------------------------------------------------------

		//calculate parameters for N
		float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
		float phi = phi_r - phi_i;


		//calculation of h from marschner but the concrete equation is from d'Eon paper (Above Equation 9)
		float a = 1 / material->m_ior;
		float nenner = glm::sqrt(1 + glm::pow2(a) - 2 * a * glm::sign(phi) * glm::sin(phi / 2.f));
		float h = (glm::sign(phi) * glm::cos(phi / 2.f)) / nenner;
		float gamma_i = glm::asin(h);

		//Bravais (virtual index of reflection eta_one, eta_two) calculation
		float cos_theta_d = glm::cos(theta_d);
		float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
		float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta). From Marschner Paper appendix
		float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)
		float c = glm::asin(1.f / bravais_first);

		//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
		float dphi_dh = 1.f / (1 / glm::sqrt(1 - h * h) * ((-(24 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i)) + (6 * c / glm::pi<float>() - 2)));

		//calculate fresnel for attenuation factor
		float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);

		//new color absorption coefficient
		glm::vec3 sigma_a = KIRK::MarschnerHairBSDF::calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
		//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));

		//attenuation factor, which contains color absorbtion. Marschner above Equation 8
		glm::vec3 att_color = glm::pow2(1 - fresnel) * KIRK::MarschnerHairBSDF::T(sigma_a, glm::asin(h / bravais_first));

		//final term for N_tt(phi). Marschner Equation 8.
		glm::vec3 n_tt = 0.5f * att_color * glm::abs(2 * dphi_dh);
		
		//--------------------------------------------------------------------

		//add tt-path value to final scattering function value
		scat_tt = (m_tt * n_tt / glm::pow2(cos_theta_d));
	}

	//////////////////////////////////////////
	////
	//// Cylinder Path: TRT
	////
	//////////////////////////////////////////
	{
		//refract input ray on surface
		glm::vec3 t_dir = glm::refract(norm_in_ray, glm::faceforward(normal, norm_in_ray, normal), 1.0f / material->m_ior);
		glm::vec3 t_normal;
		//New Intersection for T path
		Intersection t_hit(Ray(pos + 1e-4f * t_dir, t_dir));
		//get the intersection point on the second wall
		if (m_cpuscene->getDataStructure().closestIntersection(&t_hit))
		{
			cylinder_obj->calcNormal(&t_hit);// calculation of the normal
			t_normal = t_hit.m_normal;
		}
		else { t_normal = -normal; }

		//Reflect on the second cylinder wall 
		glm::vec3 tr_dir = glm::reflect(glm::normalize(t_dir), faceforward(t_normal, glm::normalize(t_dir), t_normal));
		glm::vec3 tr_normal;
		//New Intersection for TR path
		Intersection tr_hit(Ray(t_hit.m_location + 1e-4f * tr_dir, tr_dir));
		//get the intersection for the second point on the first wall
		if (m_cpuscene->getDataStructure().closestIntersection(&tr_hit))
		{
			cylinder_obj->calcNormal(&tr_hit);// calculation of the normal
			tr_normal = tr_hit.m_normal;
		}
		else { tr_normal = normal; }

		//finaly refract on the first wall of the cylinder to get output ray
		glm::vec3 out_ray = glm::refract(glm::normalize(tr_dir), glm::faceforward(tr_normal, glm::normalize(tr_dir), tr_normal), 1.0f / material->m_ior);

		//calculation of tangent vector of second intersection
		glm::vec3 c3 = glm::cross(tr_normal, glm::vec3(0.0, 0.0, 1.0));
		glm::vec3 c4 = glm::cross(tr_normal, glm::vec3(0.0, 1.0, 0.0));
		glm::vec3 tr_tangent = (glm::length(c3) > glm::length(c4)) ? glm::normalize(c3) : glm::normalize(c4);
		//rotate towards normal to account for the tilted fiber surface
		out_ray = glm::vec3(glm::vec4(out_ray, 0.f) * glm::rotate(glm::radians(-3.f * alpha_r / 2.f), tr_tangent));

		//move output ray to cylinders (fibers) local space.
		/*  Our v-axis(vector through the center of the cylinder(apexpoint - basepoint)) is u-axis in marschner model*/
		glm::vec3 out_ray_cyl = Math::worldToLocal(glm::normalize(out_ray), cylinder_obj->getV(), cylinder_obj->getU(), cylinder_obj->getW());

		//--------------------------------------------------------------------
		// M_trt(theta_h) : Marschner marginal, longitudinal scattering function(M)
		//--------------------------------------------------------------------

		//calculate parameters for M
		float theta_r = std::atan2(out_ray_cyl.x, out_ray_cyl.z);//Angle between reflected output ray and fibers u-w(normal-) plane
		float theta_h = (theta_r + theta_i) / 2.f;//theta half angle
		float theta_d = (theta_r - theta_i) / 2.f;//theta difference angle
		float gaussian_x = theta_h - glm::radians((-3 * alpha_r / 2.f));//longitudinal shift for TT-Path is smaller than for R-Path

		//M_trt(theta_h) -> gaussian function with zero-mean and our fibers material standart derivation.
		float m_trt = BSDFHelper::normal_gauss_pdf(gaussian_x, 0.0f, 2.f * beta_r);

		//--------------------------------------------------------------------
		// N_trt(phi) : Marschner conditional, azimuthal scattering function(N) [Chapter 5.2.2]
		//--------------------------------------------------------------------

		//calculate parameters for N
		float phi_r = std::atan2(std::hypot(out_ray_cyl.x, out_ray_cyl.z), out_ray_cyl.y);
		float phi = phi_r - phi_i;
		float w_c = 15.f;//azimuthal width of caustic. Between 10 and 25 degrees
		float k_g = 2.f;//glint scale factor. Between 0.5 and 5
		float t, phi_c, h;

		//help variables for bravais calculation
		float cos_theta_d = glm::cos(theta_d);
		float x1 = glm::sqrt(glm::pow2(material->m_ior) - glm::pow2(glm::sin(theta_d)));
		float bravais_first = x1 / cos_theta_d;//first bravais index (virtual eta)
		float bravais_sec = glm::pow2(material->m_ior) * cos_theta_d / x1;//second bravais index (virtual eta)
		float c = glm::asin(1.f / bravais_first);

		if (bravais_first < 2)//procedure from marschner paper chapter 5.2.2
		{
			float h_c = glm::sqrt((4 - glm::pow2(bravais_first)) / 3.f);//Marschner paper Equation 4
			phi_c = glm::radians(4 * glm::asin(h_c / bravais_first) - 2 * glm::asin(h_c) + 2 * glm::pi<float>());
			//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
			float dphi_dh_h_c = 1.f / (1 / glm::sqrt(1 - h_c * h_c) * ((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(glm::asin(h_c)) + (12 * c / glm::pi<float>() - 2))));
			h = glm::min(0.5f, 2 * glm::sqrt(2 * w_c / glm::abs(dphi_dh_h_c)));
			t = 1;
		}
		else
		{
			phi_c = 0;
			h = 0.5f;
			t = glm::smoothstep(2.f, 2.3f, bravais_first);
		}
		
		float gamma_i = glm::asin(h);
		//From ruenz12 bachelorthesis Equation 44, which derived it from marschner Equation 10.
		float dphi_dh = 1.f / (1 / glm::sqrt(1 - h * h) * ((-(48 * c / glm::pow3(glm::pi<float>())) * glm::pow2(gamma_i) + (12 * c / glm::pi<float>() - 2))));

		//calculate fresnel part of attenuation
		float fresnel = BSDFHelper::dialectricFresnel(glm::cos(gamma_i), bravais_first, bravais_sec);
		float gamma_t = glm::asin(h / bravais_first);
		float fresnel_exit = BSDFHelper::dialectricFresnel(glm::cos(gamma_t), 1 / bravais_first, 1 / bravais_sec);

		//new color absorption coefficient
		glm::vec3 sigma_a = KIRK::MarschnerHairBSDF::calcSigmaFromColor(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));
		//glm::vec3 sigma_a = glm::vec3(material->fetchParameterColor<MatParamType::SIGMA_A>(texcoord));

		//full attenuation
		glm::vec3 att_color = glm::pow2(1 - fresnel) * fresnel_exit * glm::pow2(KIRK::MarschnerHairBSDF::T(sigma_a, gamma_t));

		//calculate gauss values for power apporximation from section 5.2.2
		float gauss_0 = BSDFHelper::normal_gauss_pdf(0.f, 0.f, w_c);
		float gauss_phi_diff = BSDFHelper::normal_gauss_pdf(phi - phi_c, 0.f, w_c);
		float gauss_phi_sum = BSDFHelper::normal_gauss_pdf(phi + phi_c, 0.f, w_c);

		//final term for N_trt(phi). Marschner Equation 8
		glm::vec3 n_trt = 0.5f * att_color / glm::abs(2 * dphi_dh);
		n_trt *= (1 - t * gauss_phi_diff / gauss_0);
		n_trt *= (1 - t * gauss_phi_sum / gauss_0);
		n_trt += t * k_g * att_color * h * (gauss_phi_diff + gauss_phi_sum);

		//--------------------------------------------------------------------
		
		//return final scattering function
		scat_trt = m_trt * n_trt / glm::pow2(cos_theta_d);
	}

	/////////////////////////////////
	// Final Color calculations. We add a diffuse part to the marschner specular part as suggested by marschner.
	/////////////////////////////////

	//specular part of the color
	KIRK::Color::RGBA directLight = calcDirectLight(hit);
	KIRK::Color::RGBA specular(1000*scat_r + 10*scat_tt + 10*scat_trt, 1.0f);
	//KIRK::Color::RGBA specular(scat_r + scat_tt + scat_trt, 1.0f);

	//diffuse part of the color. Similar to Kajiya and Kay lighting
	//glm::vec4 diffuse = 0.4f * glm::sqrt(1 - sin_theta_i * sin_theta_i) * glm::vec4(material->fetchParameterColor<MatParamType::DIFFUSE>(texcoord));

	return directLight * KIRK::Color::RGBA(scat_r, 1.f) * glm::cos(theta_i);//Scattering Integral. MarschnerPaper Equation 1

}

KIRK::Color::RGBA KIRK::CPU::SimpleCPURaytracer::calcDirectLight(const KIRK::Intersection & hit)
{
	Color::RGBA directLight(0.0f);

	if (m_cpuscene->getLights().size() == 0)
		return directLight;

	/*******************************************************/
	/******************* LIGHT SAMPLING ********************/
	/*******************************************************/
	//Random Value between 0 and 1, to select light source randomly
	std::random_device rd;
	std::mt19937 gen = std::mt19937(rd());
	std::uniform_real_distribution<> dist{ 0, 1 };

	int lightIndex = dist(gen) * m_cpuscene->getLights().size();
	auto light = m_cpuscene->getLights()[lightIndex].get();

	float attenuation;
	KIRK::Ray hitToLight = light->calcLightdir(hit.m_location, attenuation, true);//calculates ray from intersection point towards light source
	glm::vec3 lightpos = hitToLight.m_origin + hitToLight.m_direction;
	hitToLight.m_origin += 1e-4f * glm::faceforward(hit.m_normal, hitToLight.m_origin - lightpos, hit.m_normal);
	hitToLight.m_direction = glm::normalize(hitToLight.m_direction);

	attenuation = 1;//set attenuation to 1 because we already calculated attenuation factor in the bsdf sample

	auto lightColor = light->m_color;

	if (light->m_color.x > 0 || light->m_color.y > 0 || light->m_color.z > 0)
	{
		lightColor *= attenuation//attenuation factor of the light. Calculated in light->calcLightdir()
			* std::abs(glm::dot(hitToLight.m_direction, hit.m_normal));//angle between ray towards light and normal

		if (lightColor.x > 0 || lightColor.y > 0 || lightColor.z > 0)
		{
			float t_max = glm::length(lightpos - hitToLight.m_origin);
			bool hitToLight_hasIntersection = m_cpuscene->getDataStructure().isIntersection(&hitToLight, t_max);
			if (!hitToLight_hasIntersection)
			{
				// Test for light Intersections
				for (int i = 0; i < m_cpuscene->getLights().size(); i++)
				{
					auto lightIntersectionTest = m_cpuscene->getLights()[i].get();
					float t;
					if (lightIntersectionTest->isIntersection(hitToLight, t) && (t < t_max))
					{
						hitToLight_hasIntersection = true;
						break;
					}
				}
			}
			lightColor *= !hitToLight_hasIntersection;
		}
		directLight += lightColor;
	}

	return directLight;
}



void KIRK::CPU::SimpleCPURaytracer::onGui()
{
	ImGui::Text("Tools:");
	ImGui::CheckboxFlags("Progress bar", &m_flags, RTFLAG_USE_PROGRESSBAR);
	ImGui::Spacing();
	ImGui::Text("Performance:");
	ImGui::SliderInt("Threads", &m_num_threads, 1, KIRK::ThreadManager::maxThreads());
	ImGui::Separator();
	ImGui::SliderInt("MaxDepth", &m_depth, 0, 20);
	ImGui::Spacing();
	ImGui::Text("Effects:");
	ImGui::CheckboxFlags("SoftShadows", &m_flags, RTFLAG_USE_SOFT_SHADOWS);
	ImGui::SliderInt("Light Samples", &m_num_lightsamples, 1, 256);
	ImGui::Spacing();
	ImGui::CheckboxFlags("Reflection", &m_flags, RTFLAG_USE_REFLECTIONS);
	ImGui::CheckboxFlags("Refraction", &m_flags, RTFLAG_USE_REFRACTIONS);
	ImGui::Separator();
	ImGui::CheckboxFlags("Super-Sampling", &m_flags, RTFLAG_USE_SUPERSAMPLING);
	ImGui::SliderInt("Samples", &m_num_supersamples, 2, 16);
	ImGui::CheckboxFlags("Adaptive-Sampling", &m_flags, RTFLAG_USE_ADAPTIVE_SAMPLING);
	ImGui::SliderInt("Adapt. Depth", &m_adaptive_depth, 1, 16);
	ImGui::SliderFloat("Adapt. Thresh.", &m_max_adaptive_difference, 0.01f, 1.f, "%.2f");
	ImGui::CheckboxFlags("Poisson Disk-Sampling", &m_flags, RTFLAG_USE_POISSONDISK_SAMPLING);
	ImGui::SliderInt("Samples", &m_num_poissondisksamples, 2, 16);
	ImGui::Separator();
	ImGui::CheckboxFlags("Depth of Field", &m_flags, RTFLAG_USE_DOF);
	ImGui::DragInt("DOF Samples", &m_num_blursamples, 1, 2, 500);

	//    ImGui::Separator();
	//    ImGui::Text("Experimental");
	//    ImGui::DragFloat("Focus", &getCamera()->m_focus_distance, 0.1, 0.1, 500);
}

