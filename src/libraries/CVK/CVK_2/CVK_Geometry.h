#ifndef __CVK_GEOMETRY_H
#define __CVK_GEOMETRY_H

#include "CVK_Defs.h"
#include <fstream>

namespace CVK {

class Geometry
{
public:
    /**
     * Standard Constructor for Geometry
     */
    Geometry();
    ~Geometry();
    /**
     * @brief Getter for geometry type
     * @return the geometry type of this object
     */
    int getGeoType();

    void createBuffers();
    /**
     * Renders this Geometry by using its VAO as Triangles in OpenGL. Can be overridden in
     * subclasses, for example to render lines or points.
     * @brief Render this object
     */
    virtual void render();

    /**
     * Computes the tangents for all vertices of this Geometry. Uses the values of the vertices
     * and the uv coordinates attribute. Stores data in attribute.
     * @brief Computes and stores tangents
     */
    void computeTangents();

    /**
     * @brief Standard Getter for vertices
     * @return the vertices of this object
     */
    std::vector<glm::vec4> *getVertices();
    /**
     * @brief Standard Getter for normals
     * @return the normals of this object as pointer
     */
    std::vector<glm::vec3> *getNormals();
    /**
     * @brief Standard Getter for uv coordinates
     * @return the uv coordinates of this object as pointer
     */
    std::vector<glm::vec2> *getUVs();
    /**
     * @brief Standard Getter for index list
     * @return the index list of this object as pointer
     */
    std::vector<unsigned int> *getIndex();
    /**
     * @brief Standard Getter for tangents
     * @return the tangents of this object as pointer
     */
    std::vector<glm::vec3> *getTangents();

    int getSizeInBytes();

    int getMaterialIndex() { return m_materialIndex; };

    void setMaterialIndex(int index) { m_materialIndex = index; };

protected:
    int m_sizeInBytes;
    int calcSize();

    int m_geotype; //!< Information for the geometry type (f.e. sphere)
    GLuint m_vao; //!< The OpenGL Vertex Array Object
    GLuint m_vertexbuffer; //!< A Vertex Buffer Object for storing vertex positions
    GLuint m_normalbuffer; //!< A Vertex Buffer Object for storing vertex normals
    GLuint m_uvbuffer; //!< A Vertex Buffer Object for storing vertex uv coordinates
    GLuint m_indexlist; //!< A Vertex Buffer Object for storing vertex indices
    GLuint m_tangentbuffer; //!< A Vertex Buffer Object for storing vertex tangents

    int m_points; //!< Number of all vertices
    int m_indices; //!< Number of all indices

    std::vector<glm::vec4> m_vertices; //!< A list of all vertex positions
    std::vector<glm::vec3> m_normals; //!< A list of all vertex normals
    std::vector<glm::vec2> m_uvs; //!< A list of all vertex uv coordinates
    std::vector<unsigned int> m_index; //!< A list of all vertex indices
    std::vector<glm::vec3> m_tangents; //!< A list of all vertex tangents

    int m_materialIndex;
};

};

#endif /* __CVK_GEOMETRY_H */

