#ifndef __CVK_DEFS_H
#define __CVK_DEFS_H

#include <vector>
#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_precision.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtx/component_wise.hpp>

//OpenGL
#define VERTICES 0
#define NORMALS 1
#define TEXTURECOORDS 2
#define TANGENTS 3
#define COLOUR 4
#define ATTRIBUTE_0 5
#define ATTRIBUTE_1 6

#define COLOR_TEXTURE_UNIT GL_TEXTURE0
#define NORMAL_TEXTURE_UNIT GL_TEXTURE1
#define SHADOW_TEXTURE_UNIT GL_TEXTURE2
#define CUBE_MAP_TEXTURE_UNIT GL_TEXTURE3

#define VERTEX_SHADER_BIT 1
#define TESS_CONTROL_BIT 2
#define TESS_EVAL_BIT 4
#define GEOMETRY_SHADER_BIT 8
#define FRAGMENT_SHADER_BIT 16
#define COMPUTE_SHADER_BIT 32

#define MAT_AMBIENT_BIT 1
#define MAT_DIFFUSE_BIT 2
#define MAT_SPECULAR_BIT 4
#define MAT_TRANSPARENCY_BIT 8

#define NO_FOG 0
#define FOG_LINEAR 1
#define FOG_EXP 2
#define FOG_EXP2 3

#define MAX_LIGHTS 8

#define INVALID_OGL_VALUE 0xFFFFFFFF

//Geometry Objects
#define CVK_CONE 0
#define CVK_CUBE 1
#define CVK_GEOMETRY 2
#define CVK_LINES 3
#define CVK_PLANE 4
#define CVK_POINTS 5
#define CVK_SPHERE 6
#define CVK_TEAPOT 7
#define CVK_TRIANGLE 8
#define CVK_WIRECUBE 9

//Ray Tracing
#define EPS 1E-5
//#define Rayeps 1e-4f
#define Minweight 0.01f

//Colours
#define darkgrey glm::vec3( 0.2, 0.2, 0.2)
#define bluegrey glm::vec3( 0.8, 0.8, 1.0)
#define bronze glm::vec3( 0.59, 0.35, 0.22)

#define Aquamarine glm::vec3( .439216, .858824, .576471)
#define Black glm::vec3( 0, 0, 0)
#define Blue glm::vec3( 0, 0, 1)
#define BlueViolet glm::vec3( .623529, .372549, .623529)
#define Brown glm::vec3( .647059, .164706, .164706)
#define CadetBlue glm::vec3( .372549, .623529, .623529)
#define Coral glm::vec3( 1, .498039, 0)
#define CornflowerBlue glm::vec3( .258824, .258824, .435294)
#define Cyan glm::vec3( 0, 1, 1)
#define DarkGreen glm::vec3( .184314, .309804, .184314)
#define DarkOliveGreen glm::vec3( .309804, .309804, .184314)
#define DarkOrchid glm::vec3( .6, .196078, .8)
#define DarkSlateBlue glm::vec3( .419608, .137255, .556863)
#define DarkSlateGray glm::vec3( .184314, .309804, .309804)
#define DarkSlateGrey glm::vec3( .184314, .309804, .309804)
#define DarkTurquoise glm::vec3( .439216, .576471, .858824)
#define DimGray glm::vec3( .329412, .329412, .329412)
#define DimGrey glm::vec3( .329412, .329412, .329412)
#define Firebrick glm::vec3( .556863, .137255, .137255)
#define ForestGreen glm::vec3( .137255, .556863, .137255)
#define Gold glm::vec3( .8, .498039, .196078)
#define Goldenrod glm::vec3( .858824, .858824, .439216)
#define Gray glm::vec3( .752941, .752941, .752941)
#define Green glm::vec3( 0, 1, 0)
#define GreenYellow glm::vec3( .576471, .858824, .439216)
#define Grey glm::vec3( .752941, .752941, .752941)
#define IndianRed glm::vec3( .309804, .184314, .184314)
#define Khaki glm::vec3( .623529, .623529, .372549)
#define LightBlue glm::vec3( .74902, .847059, .847059)
#define LightGray glm::vec3( .658824, .658824, .658824)
#define LightGrey glm::vec3( .658824, .658824, .658824)
#define LightSteelBlue glm::vec3( .560784, .560784, .737255)
#define LimeGreen glm::vec3( .196078, .8, .196078)
#define Magenta glm::vec3( 1, 0, 1)
#define Maroon glm::vec3( .556863, .137255, .419608)
#define MediumAquamarine glm::vec3( .196078, .8, .6)
#define MediumBlue glm::vec3( .196078, .196078, .8)
#define MediumForestGreen glm::vec3( .419608, .556863, .137255)
#define MediumGoldenrod glm::vec3( .917647, .917647, .678431)
#define MediumOrchid glm::vec3( .576471, .439216, .858824)
#define MediumSeaGreen glm::vec3( .258824, .435294, .258824)
#define MediumSlateBlue glm::vec3( .498039, 0, 1)
#define MediumSpringGreen glm::vec3( .498039, 1, 0)
#define MediumTurquoise glm::vec3( .439216, .858824, .858824)
#define MediumVioletRed glm::vec3( .858824, .439216, .576471)
#define MidnightBlue glm::vec3( .184314, .184314, .309804)
#define Navy glm::vec3( .137255, .137255, .556863)
#define NavyBlue glm::vec3( .137255, .137255, .556863)
#define Orange glm::vec3( .8, .196078, .196078)
#define OrangeRed glm::vec3( 1, 0, .498039)
#define Orchid glm::vec3( .858824, .439216, .858824)
#define PaleGreen glm::vec3( .560784, .737255, .560784)
#define Pink glm::vec3( .737255, .560784, .560784)
#define Plum glm::vec3( .917647, .678431, .917647)
#define Red glm::vec3( 1, 0, 0)
#define Salmon glm::vec3( .435294, .258824, .258824)
#define SeaGreen glm::vec3( .137255, .556863, .419608)
#define Sienna glm::vec3( .556863, .419608, .137255)
#define SkyBlue glm::vec3( .196078, .6, .8)
#define SlateBlue glm::vec3( 0, .498039, 1)
#define SpringGreen glm::vec3( 0, 1, .498039)
#define SteelBlue glm::vec3( .137255, .419608, .556863)
#define Thistle glm::vec3( .847059, .74902, .847059)
#define Turquoise glm::vec3( .678431, .917647, .917647)
#define Violet glm::vec3( .309804, .184314, .309804)
#define VioletRed glm::vec3( .8, .196078, .6)
#define Wheat glm::vec3( .847059, .847059, .74902)
#define White glm::vec3( .988235, .988235, .988235)
#define Yellow glm::vec3( 1, 1, 0)
#define YellowGreen glm::vec3( .6, .8, .196078)
#define aquamarine glm::vec3( .439216, .858824, .576471)
#define black glm::vec3( 0, 0, 0)
#define blue glm::vec3( 0, 0, 1)
#define blue_violet glm::vec3( .623529, .372549, .623529)
#define brown glm::vec3( .647059, .164706, .164706)
#define cadet_blue glm::vec3( .372549, .623529, .623529)
#define coral glm::vec3( 1, .498039, 0)
#define cornflower_blue glm::vec3( .258824, .258824, .435294)
#define cyan glm::vec3( 0, 1, 1)
#define dark_green glm::vec3( .184314, .309804, .184314)
#define dark_olive_green glm::vec3( .309804, .309804, .184314)
#define dark_orchid glm::vec3( .6, .196078, .8)
#define dark_slate_blue glm::vec3( .419608, .137255, .556863)
#define dark_slate_gray glm::vec3( .184314, .309804, .309804)
#define dark_slate_grey glm::vec3( .184314, .309804, .309804)
#define dark_turquoise glm::vec3( .439216, .576471, .858824)
#define dim_gray glm::vec3( .329412, .329412, .329412)
#define dim_grey glm::vec3( .329412, .329412, .329412)
#define firebrick glm::vec3( .556863, .137255, .137255)
#define forest_green glm::vec3( .137255, .556863, .137255)
#define gold glm::vec3( .8, .498039, .196078)
#define goldenrod glm::vec3( .858824, .858824, .439216)
#define gray glm::vec3( .752941, .752941, .752941)
#define green glm::vec3( 0, 1, 0)
#define green_yellow glm::vec3( .576471, .858824, .439216)
#define grey glm::vec3( .752941, .752941, .752941)
#define indian_red glm::vec3( .309804, .184314, .184314)
#define khaki glm::vec3( .623529, .623529, .372549)
#define light_blue glm::vec3( .74902, .847059, .847059)
#define light_gray glm::vec3( .658824, .658824, .658824)
#define light_grey glm::vec3( .658824, .658824, .658824)
#define light_steel_blue glm::vec3( .560784, .560784, .737255)
#define lime_green glm::vec3( .196078, .8, .196078)
#define magenta glm::vec3( 1, 0, 1)
#define maroon glm::vec3( .556863, .137255, .419608)
#define medium_aquamarine glm::vec3( .196078, .8, .6)
#define medium_blue glm::vec3( .196078, .196078, .8)
#define medium_forest_green glm::vec3( .419608, .556863, .137255)
#define medium_goldenrod glm::vec3( .917647, .917647, .678431)
#define medium_orchid glm::vec3( .576471, .439216, .858824)
#define medium_sea_green glm::vec3( .258824, .435294, .258824)
#define medium_slate_blue glm::vec3( .498039, 0, 1)
#define medium_spring_green glm::vec3( .498039, 1, 0)
#define medium_turquoise glm::vec3( .439216, .858824, .858824)
#define medium_violet_red glm::vec3( .858824, .439216, .576471)
#define midnight_blue glm::vec3( .184314, .184314, .309804)
#define navy glm::vec3( .137255, .137255, .556863)
#define navy_blue glm::vec3( .137255, .137255, .556863)
#define orange glm::vec3( .8, .196078, .196078)
#define orange_red glm::vec3( 1, 0, .498039)
#define orchid glm::vec3( .858824, .439216, .858824)
#define pale_green glm::vec3( .560784, .737255, .560784)
#define pink glm::vec3( .737255, .560784, .560784)
#define plum glm::vec3( .917647, .678431, .917647)
#define red glm::vec3( 1, 0, 0)
#define salmon glm::vec3( .435294, .258824, .258824)
#define sea_green glm::vec3( .137255, .556863, .419608)
#define sienna glm::vec3( .556863, .419608, .137255)
#define sky_blue glm::vec3( .196078, .6, .8)
#define slate_blue glm::vec3( 0, .498039, 1)
#define spring_green glm::vec3( 0, 1, .498039)
#define steel_blue glm::vec3( .137255, .419608, .556863)
#define thistle glm::vec3( .847059, .74902, .847059)
#define turquoise glm::vec3( .678431, .917647, .917647)
#define violet glm::vec3( .309804, .184314, .309804)
#define violet_red glm::vec3( .8, .196078, .6)
#define wheat glm::vec3( .847059, .847059, .74902)
#define white glm::vec3( .988235, .988235, .988235)
#define yellow glm::vec3( 1, 1, 0)
#define yellow_green glm::vec3( .6, .8, .196078)

#endif /* __CVK_DEFS_H */

