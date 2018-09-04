#ifndef LIBRARIES_CVK_RT_CVK_RT_SIMPLECPURAYTRACER_H_
#define LIBRARIES_CVK_RT_CVK_RT_SIMPLECPURAYTRACER_H_

//includes
#include <vector>
#include <thread>
#include <random>
#include <cmath>
#include "CVK/CVK_2/CVK_Perspective.h"
#include "KIRK/CPU/CPU_Datastructures/CPU_DataStructure.h"
#include "CPU_Scene.h"
#include "CPU_Raytracer.h"

#include "KIRK/CPU/CPU_Datastructures/CPU_NoDataStructure.h"
#include "KIRK/CPU/CPU_Datastructures/UniformGrid.h"
#include "KIRK/CPU/CPU_Datastructures/Octree.h"
#include "KIRK/CPU/CPU_Datastructures/CPU_BVH.h"

#include "KIRK/Utils/Threading.h"
#include "KIRK/Common/Camera.h"
#include "CVK/CVK_Utils/CVK_ConverterUtils.h"
#include "KIRK/Common/Material.h"
#include "KIRK/Common/Ray.h"
#include "KIRK/Common/Intersection.h"
#include "KIRK/Common/Environment.h"
#include "externals/PoissonDiskGenerator/PoissonGenerator.h"
#include "KIRK/Utils/Gui/Gui.h"

namespace KIRK
{
	namespace CPU
	{

		class SimpleCPURaytracer : public KIRK::GuiElement, public CPU_Raytracer
		{
		public:
			/** default constructor */
			SimpleCPURaytracer();

			SimpleCPURaytracer(std::shared_ptr<KIRK::CPU::Scene> cpuscene, int depth);

			/** destructor */
			~SimpleCPURaytracer();


			/** fill gui
			* fills gui with all controllable variables
			*/
			void onGui() override;

			void setFlags(unsigned int flags)
			{
				m_flags = flags;
			}

			void enable(unsigned int flag)
			{
				m_flags |= flag;
			}

			void disable(unsigned int flag)
			{
				m_flags &= ~flag;
			}

			bool isEnabled(unsigned int flag)
			{
				return m_flags & flag;
			}

			const static unsigned int RTFLAG_USE_PROGRESSBAR = 0x0001;
			const static unsigned int RTFLAG_USE_REFLECTIONS = 0x0002;
			const static unsigned int RTFLAG_USE_REFRACTIONS = 0x0004;
			const static unsigned int RTFLAG_USE_SOFT_SHADOWS = 0x0008;
			const static unsigned int RTFLAG_USE_DOF = 0x0010;
			const static unsigned int RTFLAG_USE_SUPERSAMPLING = 0x0020;
			const static unsigned int RTFLAG_USE_ADAPTIVE_SAMPLING = 0x0040;
			const static unsigned int RTFLAG_USE_POISSONDISK_SAMPLING = 0x0080;

			int m_num_supersamples;			//!< m_superSamples*m_superSamples rays per pixel
			int m_num_blursamples;			//!< rays per pixel
			int m_num_lightsamples;         //!< samples per light used for soft shadows
			int m_num_threads;         //!< CPU thread count to calculate on
			int m_num_poissondisksamples; //!< max samples for poisson disk sampling
			int m_adaptive_depth;      //!< max depth of the recursion from adaptiveSampling

			float m_max_adaptive_difference;	//!< Difference below which two colors are considered equal enough for adaptive sampling



		protected:
		private:
			unsigned int m_flags = 0;

			/** render image
			* ray traces the whole image, so loops over all pixels and traces each with a ray
			*/
			void render() override;

			//inline void showIterationInfo();
			/** trace ray
			 *
			 * traces the given ray and calls the intersection function of the datastructur
			 *
			 * @param ray ray for each pixel with direction and startposition
			 * @param level level of recursion depth, must be smaller than the maximum depth, else it terminates
			 * @param weight weight of each bounced ray, if it is too small, it terminates
			 * @return color color of background if we dont hit anything in the scene
			 */
			KIRK::Color::RGBA trace(Ray &ray, int level, float weight);

			/** direct illumination
			 *
			 * direct illumination for each lightsource 3# diffuse shading 4# specular shading
			 *
			 * @param pos position of intersection
			 * @param norm normal of intersection
			 * @param view direction of intersection
			 * @param material material of the intersected object
			 * @param diff_color diffuse color of the intersected object, calculated in JustColor
			 * @return color with direct illumination
			 */
			KIRK::Color::RGBA lightShading(const glm::vec3 &pos, const glm::vec3 &norm, const glm::vec2 &texcoord, const glm::vec3 &view, KIRK::Material *material,
				const KIRK::Color::RGBA &diff_color);

			/** reflections
			 *
			 * calculates for a given intersection the reflection color recursive
			 *
			 * @param pos position of intersection
			 * @param view direction of intersecting ray
			 * @param norm normal of intersection
			 * @param falloff (here: specular reflectivity of the material) color weight loss factor after each bounce
			 * @param level level of recursion depth, must be smaller than the maximum depth, else it terminates
			 * @param weight weight of each bounced ray, if it is too small, it terminates
			 * @return color of the reflection
			 */
			KIRK::Color::RGBA reflection(const glm::vec3 &pos, const glm::vec3 &view, const glm::vec3 &norm, float falloff, float roughness, int level, float weight);

			/** refractions
			 *
			 * calculates for a given intersection the refraction color recursive
			 *
			 * @param hit intersection of ray with the object
			 * @param pos position of intersection
			 * @param view direction of intersection
			 * @param norm normal of intersection
			 * @param falloff (here: transparency of the material) color weight loss factor after each bounce
			 * @param ior ior is the index of refraction
			 * @param level level of recursion depth, must be smaller than the maximum depth, else it terminates
			 * @param weight weight of each bounced ray, if it is too small, it terminates
			 * @return color of the refraction
			 */
			KIRK::Color::RGBA refraction(Intersection &hit, const glm::vec3 &pos, const glm::vec3 &view, const glm::vec3 &norm, float falloff, float roughness, float ior, int level,
				float weight);

			/** depth of field
			 *
			 * calculate distance to focus plane, by that generate random new rays for each pixel
			 *
			 * @param level
			 * @param weight
			 * @param view
			 * @return color
			 */
			KIRK::Color::RGBA depthOfField(int level, float weight, const glm::vec3 &view);

			/** super-sampling
			 *
			 * generate A*A rays per pixel
			 *
			 * @param x
			 * @param y
			 */
			KIRK::Color::RGBA superSampling(int x, int y);

			/** adaptive sampling
			 *
			 * generate A*A rays per pixel if there was a change in material or background
			 *
			 * @param x
			 * @param y
			 */
			KIRK::Color::RGBA adaptiveSampling(int x, int y);
			KIRK::Color::RGBA adaptiveSamplingRecursive(glm::vec3 &dir1, glm::vec3 &dir2, glm::vec3 &dir3, glm::vec3 &dir4,
				Color::RGBA color1, Color::RGBA color2, Color::RGBA color3, Color::RGBA color4, glm::vec3 &position, int depth);

			/** poisson disk sampling
			*
			*
			*
			* @param x
			* @param y
			*/
			KIRK::Color::RGBA poissonDiskSampling(int x, int y);


			/** shade
			 *
			 * Applies several shading methods depending on user settings.
			 *
			 * @param hit
			 * @param level
			 * @param weight
			 * @return color The final color.
			 */
			KIRK::Color::RGBA shade(Intersection &hit, int level, float weight);

			/** shade
			*
			* Applies the marschner Hair shading model
			*
			* @param hit
			* @param level
			* @param weight
			* @return color The final color.
			*/
			KIRK::Color::RGBA shadeMarschnerHair(Intersection &hit, int level, float weight);

			// Predefined PoissonDisks for poissonDisk Sampling
			std::vector<std::vector<PoissonGenerator::sPoint>> m_poissonDisks{
				{ { 0.876225f,0.73452f },{ 0.270721f,0.308664f } },
				{ { 0.169044f,0.162109f },{ 0.937509f,0.0716912f },{ 0.53597f,0.728993f } },
				{ { 0.992372f,0.949411f },{ 0.0445899f,0.692603f },{ 0.380319f,0.210391f },{ 0.87717f,0.270143f } },
				{ { 0.620592f,0.470213f },{ 0.0960264f,0.985463f },{ 0.060037f,0.410037f },{ 0.390997f,0.0112952f },{ 0.647264f,0.996541f } },
				{ { 0.551783f,0.451631f },{ 0.25204f,0.0159029f },{ 0.573596f,0.926213f },{ 0.0997206f,0.670109f },{ 0.980341f,0.145409f },{ 0.97494f,0.576453f } },
				{ { 0.214503f,0.0439718f },{ 0.545092f,0.241408f },{ 0.428024f,0.670745f },{ 0.0782009f,0.499675f },{ 0.940662f,0.519019f },{ 0.874426f,0.894647f },{ 0.886778f,0.00616533f } },
				{ { 0.104548f,0.530203f },{ 0.578633f,0.682645f },{ 0.314238f,0.962926f },{ 0.600519f,0.200092f },{ 0.254885f,0.043577f },{ 0.974884f,0.141768f },{ 0.994096f,0.534767f },{ 0.92046f,0.883586f } },
				{ { 0.901776f,0.792632f },{ 0.589908f,0.386356f },{ 0.330919f,0.980232f },{ 0.964372f,0.357401f },{ 0.427321f,0.0535911f },{ 0.0559535f,0.0748883f },{ 0.047981f,0.463666f },{ 0.362323f,0.64631f },{ 0.0106121f,0.795087f } },
				{ { 0.933456f,0.396886f },{ 0.617932f,0.370689f },{ 0.700235f,0.760013f },{ 0.798147f,0.0315559f },{ 0.387831f,0.0311461f },{ 0.283536f,0.408451f },{ 0.171462f,0.998523f },{ 0.960136f,0.994702f },{ 0.369335f,0.7431f },{ 0.00755554f,0.65519f } },
				{ { 0.581315f,0.292595f },{ 0.228379f,0.621364f },{ 0.588956f,0.702809f },{ 0.141755f,0.000873208f },{ 0.739932f,0.0314237f },{ 0.965892f,0.722718f },{ 0.999566f,0.198001f },{ 0.273087f,0.29625f },{ 0.393969f,0.973714f },{ 0.774489f,0.956949f },{ 0.0772681f,0.960969f } },
				{ { 0.912188f,0.802024f },{ 0.961714f,0.51178f },{ 0.734641f,0.31379f },{ 0.366667f,0.812687f },{ 0.493309f,0.514161f },{ 0.658913f,0.952572f },{ 0.209367f,0.294254f },{ 0.996308f,0.156222f },{ 0.733691f,0.0240584f },{ 0.127163f,0.596387f },{ 0.280606f,0.00790405f },{ 0.0859847f,0.934488f } },
				{ { 0.314857f,0.34218f },{ 0.0389938f,0.379982f },{ 0.144414f,0.838981f },{ 0.580802f,0.252321f },{ 0.0522627f,0.0944282f },{ 0.591252f,0.805347f },{ 0.937681f,0.37916f },{ 0.819281f,0.979265f },{ 0.98659f,0.706316f },{ 0.368842f,0.636209f },{ 0.331632f,0.0554404f },{ 0.724688f,0.559545f },{ 0.909652f,0.0450434f } },
				{ { 0.287403f,0.450963f },{ 0.36062f,0.839246f },{ 0.698174f,0.535026f },{ 0.522492f,0.32128f },{ 0.068157f,0.679375f },{ 0.278968f,0.138681f },{ 0.89391f,0.774232f },{ 0.901022f,0.0746571f },{ 0.964459f,0.361234f },{ 0.615598f,0.0303879f },{ 0.00361177f,0.00946213f },{ 0.616477f,0.982358f },{ 0.00689887f,0.287614f },{ 0.0819751f,0.953541f } },
				{ { 0.364089f,0.902787f },{ 0.519908f,0.581197f },{ 0.107378f,0.719646f },{ 0.800155f,0.674384f },{ 0.855714f,0.980433f },{ 0.922872f,0.183384f },{ 0.620388f,0.866241f },{ 0.674646f,0.338558f },{ 0.305266f,0.427158f },{ 0.442797f,0.195435f },{ 0.0339634f,0.447619f },{ 0.13459f,0.197146f },{ 0.0296116f,0.983108f },{ 0.678076f,0.0396471f },{ 0.953666f,0.463888f } },
				{ { 0.386733f,0.247876f },{ 0.745081f,0.504502f },{ 0.0385348f,0.238517f },{ 0.522442f,0.673943f },{ 0.496945f,0.00602357f },{ 0.141586f,0.502205f },{ 0.187005f,0.0286478f },{ 0.708468f,0.198206f },{ 0.972764f,0.337947f },{ 0.948924f,0.0754391f },{ 0.700029f,0.976794f },{ 0.992974f,0.68738f },{ 0.0380665f,0.900975f },{ 0.346203f,0.999184f },{ 0.983111f,0.979181f },{ 0.278199f,0.749605f } },
				{ { 0.992993f,0.40671f },{ 0.769727f,0.241442f },{ 0.915908f,0.874689f },{ 0.619215f,0.489455f },{ 0.922161f,0.0514146f },{ 0.627607f,0.0139559f },{ 0.50087f,0.265794f },{ 0.250488f,0.00668951f },{ 0.252711f,0.272758f },{ 0.174894f,0.656821f },{ 0.00481866f,0.414715f },{ 0.0273106f,0.155253f },{ 0.448527f,0.715451f },{ 0.061912f,0.876643f },{ 0.322863f,0.975321f },{ 0.582286f,0.954589f },{ 0.822804f,0.649423f } },
				{ { 0.721329f,0.749725f },{ 0.786439f,0.990406f },{ 0.501183f,0.454394f },{ 0.376124f,0.751396f },{ 0.986374f,0.371044f },{ 0.726934f,0.359143f },{ 0.952671f,0.612057f },{ 0.194439f,0.946106f },{ 0.260489f,0.505897f },{ 0.0110943f,0.468902f },{ 0.105843f,0.713966f },{ 0.554982f,0.944148f },{ 0.138049f,0.248563f },{ 0.00831037f,0.00314012f },{ 0.387097f,0.220768f },{ 0.276649f,0.00226183f },{ 0.700076f,0.0875657f },{ 0.959238f,0.0826479f } },
				{ { 0.758271f,0.513327f },{ 0.918706f,0.302913f },{ 0.975563f,0.638283f },{ 0.647543f,0.105521f },{ 0.320428f,0.512952f },{ 0.524981f,0.876343f },{ 0.506304f,0.350469f },{ 0.753407f,0.756016f },{ 0.99111f,0.890908f },{ 0.974186f,0.0295781f },{ 0.320497f,0.764779f },{ 0.000435591f,0.348967f },{ 0.0734909f,0.622039f },{ 0.262959f,0.243663f },{ 0.0423637f,0.0415649f },{ 0.0942528f,0.942575f },{ 0.353191f,0.00810907f },{ 0.32676f,0.99417f },{ 0.548587f,0.613926f } },
				{ { 0.991927f,0.295881f },{ 0.727936f,0.236153f },{ 0.918267f,0.537868f },{ 0.677784f,0.491502f },{ 0.665087f,0.0162864f },{ 0.888411f,0.973133f },{ 0.665411f,0.772168f },{ 0.424267f,0.749233f },{ 0.452013f,0.52692f },{ 0.476356f,0.186048f },{ 0.982725f,0.75268f },{ 0.270039f,0.931888f },{ 0.534102f,0.980454f },{ 0.210417f,0.366478f },{ 0.241046f,0.0814072f },{ 0.0368737f,0.216178f },{ 0.00158063f,0.650603f },{ 0.00104989f,0.996817f },{ 0.228406f,0.618655f },{ 0.886539f,0.0610406f } },
				{ { 0.961898f,0.779454f },{ 0.637933f,0.854858f },{ 0.739266f,0.655549f },{ 0.993227f,0.53135f },{ 0.822513f,0.389178f },{ 0.803589f,0.996991f },{ 0.614648f,0.459445f },{ 0.293008f,0.731644f },{ 0.448355f,0.9856f },{ 0.766285f,0.131903f },{ 0.525378f,0.130455f },{ 0.986668f,0.169667f },{ 0.510398f,0.662312f },{ 0.33137f,0.263844f },{ 0.288677f,0.481056f },{ 0.0288707f,0.723423f },{ 0.100514f,0.149707f },{ 0.0153515f,0.458671f },{ 0.329768f,0.00551814f },{ 0.227408f,0.997598f },{ 0.014315f,0.948195f } },
				{ { 0.861448f,0.39533f },{ 0.872867f,0.135383f },{ 0.645922f,0.377428f },{ 0.923856f,0.651733f },{ 0.500197f,0.61749f },{ 0.71079f,0.771096f },{ 0.493734f,0.198869f },{ 0.323529f,0.776439f },{ 0.993281f,0.963391f },{ 0.545224f,0.921555f },{ 0.204907f,0.312797f },{ 0.18371f,0.0042868f },{ 0.414627f,0.000812635f },{ 0.319568f,0.992993f },{ 0.648935f,0.00990777f },{ 0.421049f,0.416701f },{ 0.0328644f,0.532359f },{ 0.259562f,0.569055f },{ 0.106114f,0.779046f },{ 0.0369884f,0.163719f },{ 0.050361f,0.99389f },{ 0.776355f,0.989412f } },
				{ { 0.539683f,0.339726f },{ 0.760168f,0.37798f },{ 0.28764f,0.333862f },{ 0.529721f,0.578776f },{ 0.29327f,0.0998164f },{ 0.553404f,0.0160797f },{ 0.740724f,0.160993f },{ 0.790615f,0.65005f },{ 0.201719f,0.545331f },{ 0.934776f,0.250462f },{ 0.966418f,0.484367f },{ 0.975057f,0.0446468f },{ 0.983336f,0.732129f },{ 0.739861f,0.999939f },{ 0.334618f,0.721845f },{ 0.257085f,0.924276f },{ 0.00376096f,0.378338f },{ 0.000503615f,0.643494f },{ 0.586293f,0.787791f },{ 0.948909f,0.952471f },{ 0.485737f,0.988914f },{ 0.0231779f,0.124771f },{ 0.0264208f,0.969608f } },
				{ { 0.760665f,0.0787493f },{ 0.790236f,0.422466f },{ 0.544694f,0.0472092f },{ 0.927503f,0.197833f },{ 0.637532f,0.275914f },{ 0.415413f,0.246208f },{ 0.991211f,0.0026554f },{ 0.237417f,0.566387f },{ 0.128003f,0.233702f },{ 0.338045f,0.0210838f },{ 0.616506f,0.548254f },{ 0.847296f,0.765835f },{ 0.988881f,0.501192f },{ 0.627354f,0.765245f },{ 0.74023f,0.96885f },{ 0.433165f,0.973701f },{ 0.345819f,0.740165f },{ 0.993057f,0.916932f },{ 0.0244333f,0.742368f },{ 0.0315515f,0.510831f },{ 0.00668293f,0.0304762f },{ 0.189864f,0.882306f },{ 0.0100299f,0.992368f },{ 0.411639f,0.451182f } },
				{ { 0.134175f,0.980575f },{ 0.15416f,0.639774f },{ 0.348765f,0.713901f },{ 0.349891f,0.994946f },{ 0.00372462f,0.802647f },{ 0.552245f,0.934461f },{ 0.556242f,0.701996f },{ 0.290202f,0.273464f },{ 0.0259436f,0.278417f },{ 0.476911f,0.438297f },{ 0.0111183f,0.00853717f },{ 0.349266f,0.0439037f },{ 0.0132533f,0.497495f },{ 0.275364f,0.477766f },{ 0.680351f,0.529737f },{ 0.834442f,0.668351f },{ 0.786403f,0.927282f },{ 0.525136f,0.242398f },{ 0.997553f,0.327382f },{ 0.796054f,0.196073f },{ 0.550889f,0.0157787f },{ 0.996704f,0.998712f },{ 0.991598f,0.797753f },{ 0.996261f,0.540838f },{ 0.916177f,0.0060515f } }
			};
		};
	};
}

#endif /* LIBRARIES_CVK_RT_CVK_RT_SIMPLECPURAYTRACER_H_ */
