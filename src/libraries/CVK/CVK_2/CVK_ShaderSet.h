#ifndef __CVK_SHADERSET_H
#define __CVK_SHADERSET_H

#include <map>
#include <memory>

#include "CVK_Defs.h"
#include "CVK_State.h"

namespace CVK {

/**
 * Base class for Shader. Loads the shader files and compiles them with OpenGL to shader programs.
 * @brief Class for loading and compiling shader files
 */
class ShaderSet
{
public:
    /**
     * Standard Constructor for ShaderSet. Normally not used.
     */
    ShaderSet();
    /**
     * Constructor for ShaderSet with given parameters
     * @param shader_mask Describes which shader files are used
     * @param ShaderNames Array of paths to shader files
     */
    ShaderSet(GLuint shader_mask, const char **ShaderNames);
    /**
     * Standard Destructor for ShaderSet
     */
    ~ShaderSet();

    /**
     * Loads the shader files to strings and compiles them with OpenGL
     * @brief Loading and compiling shader files
     * @param shader_mask Describes which shader files are used
     * @param ShaderNames Array of paths to shader files
     */
    void GenerateShaderProgramm(GLuint shader_mask, const char **ShaderNames);
    /**
     * @brief Standard Getter for OpenGL program id
     * @return The OpenGL program id of this object
     */
    GLuint getProgramID();
    /**
     * Binds the OpenGL shader program with OpenGL so that it can be used for rendering
     * @brief Binds shader program
     */
    void useProgram();
    /**
     * Convenience method to bind an OpenGL texture object and use it in the shader at given location.
     * Has to be used in shader implementation
     * @brief Binds the given texture to use in shader
     * @param num The location, where to bind the texture
     * @param texture The OpenGL texture object to bind
     */
    void setTextureInput(int num, GLuint texture);
    /**
    * Convenience method to set a float variable to the given variable name within a shader
    * @brief Sets the value of the float variable in the shader
    * @param variableName The name of the variable in the shader
    * @param value The float value of the variable
    */
    void setValue(const char *variableName, float value);

    /**
    * Method to set a variable to the given variable name within a shader
    * @brief Sets the value of the variable in the shader
    * @param variableName The name of the variable in the shader
    * @param value The float value of the variable
    */
    void setUniform(const char *variableName, float value) { glUniform1f(getLocation(variableName), value); }

    void setUniform(const char *variableName, int value) { glUniform1i(getLocation(variableName), value); }

    void setUniform(const char *variableName, glm::vec2 &value)
    {
        glUniform2f(getLocation(variableName), value.x, value.y);
    }

    void setUniform(const char *variableName, glm::vec3 &value)
    {
        glUniform3f(getLocation(variableName), value.x, value.y, value.z);
    }

    void setUniform(const char *variableName, glm::vec4 &value)
    {
        glUniform4f(getLocation(variableName), value.x, value.y, value.z, value.w);
    }

    void setUniform(const char *variableName, glm::ivec2 &value)
    {
        glUniform2i(getLocation(variableName), value.x, value.y);
    }

    void setUniform(const char *variableName, glm::ivec3 &value)
    {
        glUniform3i(getLocation(variableName), value.x, value.y, value.z);
    }

    void setUniform(const char *variableName, glm::ivec4 &value)
    {
        glUniform4i(getLocation(variableName), value.x, value.y, value.z, value.w);
    }

    void setUniform(const char *variableName, glm::mat2 &value)
    {
        glUniformMatrix2fv(getLocation(variableName), 1, GL_FALSE, glm::value_ptr(value));
    }

    void setUniform(const char *variableName, glm::mat3 &value)
    {
        glUniformMatrix3fv(getLocation(variableName), 1, GL_FALSE, glm::value_ptr(value));
    }

    void setUniform(const char *variableName, glm::mat4 &value)
    {
        glUniformMatrix4fv(getLocation(variableName), 1, GL_FALSE, glm::value_ptr(value));
    }

    void setUniform(const char *variableName, int count, float value[])
    {
        glUniform1fv(getLocation(variableName), count, value);
    }

    void setUniform(const char *variableName, int count, glm::vec2 value[])
    {
        glUniform2fv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

    void setUniform(const char *variableName, int count, glm::vec3 value[])
    {
        glUniform3fv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

    void setUniform(const char *variableName, int count, glm::vec4 value[])
    {
        glUniform4fv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

    void setUniform(const char *variableName, int count, int value[])
    {
        glUniform1iv(getLocation(variableName), count, value);
    }

    void setUniform(const char *variableName, int count, glm::ivec2 value[])
    {
        glUniform2iv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

    void setUniform(const char *variableName, int count, glm::ivec3 value[])
    {
        glUniform3iv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

    void setUniform(const char *variableName, int count, glm::ivec4 value[])
    {
        glUniform4iv(getLocation(variableName), count, glm::value_ptr(value[0]));
    }

private:
    /**
    * Method to get a variables location within a shader
    * @brief Gets the location of the variable in the shader
    * @param variableName The name of the variable in the shader
    */
    GLuint getLocation(const char *variableName);
    /**
     * Checks if the shader can be compiled without errors. Prints result via standard output.
     * @brief Checks shader for errors
     * @param shaderID The shader to check
     */
    void checkShader(GLuint shaderID);
    /**
    * Checks if the shader program can be compiled without errors. Prints result via standard output.
    * @brief Checks shader program for errors
    * @param programID The shader program to check
    */
    void checkProgram(GLuint programID);
    /**
    * Loads the shader file and compiles it with OpenGL
    * @brief Loading and compiling of one shader file
    * @param shaderID OpenGL shader object id
    * @param fileName Path to shader source code file
    */
    void loadShaderSource(GLint shaderID, const char *fileName);

    GLuint m_shader_mask; //!< Which shaders does the program contain

protected:
    GLuint m_ProgramID; //!< The OpenGL shader program id
    std::vector<GLuint> m_textures; //!< Convenience list of OpenGL texture objects to set in shader program

    std::shared_ptr<std::map<std::string, GLuint>> m_uniform_locations; //!< Associative array of uniform location and variable names
};

};

#endif /*__CVK_SHADERSET_H*/
