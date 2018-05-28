#ifndef __SHADERMINIMAL_H
#define __SHADERMINIMAL_H

#include "CVK_Defs.h"
#include "CVK_ShaderSet.h"
#include "CVK_State.h"
#include "CVK_Node.h"

namespace CVK {

/**
 * Minimal shader class impementation using the ShaderSet. The model, view and projection
 * matrices are set. Uses camera bound to State. Has to collect uniform locations for variables from shader first.
 * @brief Minimal shader that sets model, view and projection matrix
 * @see State
 */
class ShaderMinimal : public CVK::ShaderSet
{
public:

    /**
    * Standard constructor for ShaderMinimal
    */
    ShaderMinimal() {}

    /**
    * Constructor for ShaderMinimal with given parameters. Collects uniform locations for
    * all used variables from Shader Program.
    * @param shader_mask Describes which shader files are used
    * @param shader_paths Array of paths to shader files
    */
    ShaderMinimal(GLuint shader_mask, const char **shader_paths);
    /**
    * Initialize this shader with given parameters.
    * @param shader_mask Describes which shader files are used
    * @param shader_paths Array of paths to shader files
    */
    void setShaderPaths(GLuint shader_mask, const char **shader_paths);
    /**
     * @brief Standard Setter for model matrix in shader
     * @param modelmatrix The new model matrix of this object
     */
    void updateModelMatrix(glm::mat4 modelmatrix);
    /**
     * Sets scene dependent variables in Shader. Namely view and projection matrix by using Camera in State.
     * @brief Sets scene variables
     * @see State
     * @see Camera
     */
    virtual void update();
    /**
    * Sets node dependent variables in Shader. None used in ShaderMinimal, but needed for subclasses.
    * @brief Sets node variables
    * @see Node
    */
    virtual void update(CVK::Node *node);

private:
    GLuint m_modelMatrixID;
    GLuint m_viewMatrixID, m_projectionMatrixID;

};

};

#endif /*__SHADERMINIMAL_H*/