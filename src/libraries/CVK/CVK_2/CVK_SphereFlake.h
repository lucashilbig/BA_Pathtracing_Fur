#ifndef __CVK_SPHEREFLAKE_H
#define __CVK_SPHEREFLAKE_H

#include "CVK_Node.h"
#include "CVK_Sphere.h"

namespace CVK {
class SphereFlake : public CVK::Sphere
{
public:
    SphereFlake(int recursionDepth, std::shared_ptr<Material> mat, std::shared_ptr<Node> scene);
    SphereFlake(int recursionDepth, std::shared_ptr<Material> mat, std::shared_ptr<Node> scene, glm::mat4 ModelMatrix, float radius);

    ~SphereFlake();

private:
    void generateSphereFlake(int ctr, int recursionDepth, glm::mat4 ModelMatrix, float radius);
    std::shared_ptr<Material> m_mat;
	std::shared_ptr<Node> m_scene;
    static std::shared_ptr<CVK::Geometry> s_sphere;
};

};

#endif /* __CVK_SPHEREFLAKE_H */

