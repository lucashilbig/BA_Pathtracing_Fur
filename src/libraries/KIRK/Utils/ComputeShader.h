
#ifndef __GLWRAP_H
#define __GLWRAP_H

#include <iostream>
#include <vector>
#include <functional>
#include <map>
#include <memory>

#include <glm/ext.hpp>
#include <GL/glew.h>

#include "KIRK/Utils/Log.h"
#include "KIRK/Common/Texture.h"

namespace KIRK
{
	/**
	 * \brief Wraps OpenGL's binding targets for buffers.
	 */
	enum class Target
	{
		ARRAY_BUFFER = GL_ARRAY_BUFFER,
		ATOMIC_COUNTER_BUFFER = GL_ATOMIC_COUNTER_BUFFER,
		COPY_READ_BUFFER = GL_COPY_READ_BUFFER,
		COPY_WRITE_BUFFER = GL_COPY_WRITE_BUFFER,
		DISPATCH_INDIRECT_BUFFER = GL_DISPATCH_INDIRECT_BUFFER,
		DRAW_INDIRECT_BUFFER = GL_DRAW_INDIRECT_BUFFER,
		ELEMENT_ARRAY_BUFFER = GL_ELEMENT_ARRAY_BUFFER,
		PIXEL_PACK_BUFFER = GL_PIXEL_PACK_BUFFER,
		PIXEL_UNPACK_BUFFER = GL_PIXEL_UNPACK_BUFFER,
		QUERY_BUFFER = GL_QUERY_BUFFER,
		SHADER_STORAGE_BUFFER = GL_SHADER_STORAGE_BUFFER,
		TEXTURE_BUFFER = GL_TEXTURE_BUFFER,
		TRANSFORM_FEEDBACK_BUFFER = GL_TRANSFORM_FEEDBACK_BUFFER,
		UNIFORM_BUFFER = GL_UNIFORM_BUFFER
	};

	/**
	 * \brief A class for wrapping OpenGL Buffers.
	 * \tparam T The Buffertarget, e.g. KIRK::Target::SHADER_STORAGE_BUFFER (default).
	 */
	template <Target T = Target::SHADER_STORAGE_BUFFER>
	class Buffer
	{
	public:
		Buffer()
		{
			glGenBuffers(1, &m_id);
		};

		~Buffer()
		{
			glDeleteBuffers(1, &m_id);
		};

		Buffer(Buffer&&) = default;
		Buffer& operator=(Buffer&& other) = default;
		// copy operations are automatically deleted now

		/**
		 * \brief Bind this Buffer.
		 */
		void bind() const
		{
			glBindBuffer(static_cast<GLenum>(T), m_id);
		}

		/**
		 * \brief executes glBindBufferBase on this Buffer.
		 * \param buffer_index forwarded to the index paramter of glBindBufferBase.
		 */
		void bindBufferBase(GLuint buffer_index) const
		{
			glBindBufferBase(static_cast<GLenum>(T), buffer_index, m_id);
		}

		/**
		 * \brief Maps the currently bound Buffer to a void-pointer for reading the data stored on the GPU.
		 * \param offset The offset where the mapping will begin, in bytes.
		 * \param length The length of the data that is to be mapped, in bytes.
		 * \param access OpenGL specific bitfield controlling access modes for glMapBufferRange.
		 * \return A void-pointer to the mapped data.
		 */
		static void* map(GLintptr offset, GLsizeiptr length, GLbitfield access) { return glMapBufferRange(static_cast<GLenum>(T), offset, length, access); };


		/**
		 * \brief Unmaps the currently bound Buffer. Should be called after finishing your work on the pointer returned by Buffer::map.
		 */
		static void unmap()
		{
			glUnmapBuffer(static_cast<GLenum>(T));
		};

		/**
		 * \brief Initializes a Buffer and it's data.
		 * \param size the size that your buffer will be, in bytes.
		 * \param data void-pointer to the data. Can be null. If not null, the pointer should point to memory of at least size bytes.
		 * \param usage What the buffer will be used for, see https://www.khronos.org/opengl/wiki/GLAPI/glBufferData#Description
		 */
		static void init(GLuint size, const void* data, GLenum usage)
		{
			glBufferData(static_cast<GLenum>(T), size, data, usage);
		};

		void clearToZero()
		{
			bind();
			uint32_t zero = 0;
			//For clearing to all-zero, it should not matter how the format is.
			glClearBufferData(static_cast<GLenum>(T), GL_R32F, GL_RED, GL_UNSIGNED_BYTE, &zero);
		}

	private:
		GLuint m_id = 0; //!< OpenGL-internal buffer-ID.
	};

	/**
	 * \brief A class that should simplify the usage of OpenGL compute shaders.
	 */
	class ComputeShader
	{
	public:
		using BufferMap = std::map<std::string, std::shared_ptr<Buffer<>>>;

		/**
		 * \brief Create a ComputeShader by loading as many files as you want to append them together to a big shader source.
		 * \param shaderFiles The Paths to your shader source files.
		 */
		explicit ComputeShader(std::vector<const char*> shaderFiles);

		/**
		 * \brief Create a ComputeShader by compiling an already existing shader source code.
		 * \param shaderSource The shader sourcecode.
		 */
		explicit ComputeShader(std::string shaderSource);

		~ComputeShader();

		/**
		 * \brief Short form for glUseProgram(GLuint);
		 */
		void use() const;

		/**
		 * \brief Runs the program with 3 dimensions and locks the GL_TEXTURE_UPDATE_BARRIER_BIT memory barrier.
		 * \param global_size_x x-Size of the set to be computed on.
		 * \param global_size_y y-Size of the set to be computed on.
		 * \param global_size_z z-Size of the set to be computed on.
		 * \param barriers What Barriers to set for glMemoryBarrier. Default: GL_ALL_BARRIER_BITS
		 */
		void dispatch3D(unsigned int global_size_x, unsigned int global_size_y, unsigned int global_size_z, GLbitfield barriers = GL_ALL_BARRIER_BITS) const;

		/**
		 * \brief Runs the program with 2 dimensions and locks the GL_TEXTURE_UPDATE_BARRIER_BIT memory barrier.
		 * \param global_size_x x-Size of the set to be computed on.
		 * \param global_size_y y-Size of the set to be computed on.
		 * \param barriers What Barriers to set for glMemoryBarrier. Default: GL_ALL_BARRIER_BITS
		 */
		void dispatch2D(unsigned int global_size_x, unsigned int global_size_y, GLbitfield barriers = GL_ALL_BARRIER_BITS) const;

		/**
		 * \brief Runs the program with one dimension and locks the GL_TEXTURE_UPDATE_BARRIER_BIT memory barrier.
		 * \param global_size Size of the whole set to be computed on.
		 * \param barriers What Barriers to set for glMemoryBarrier. Default: GL_ALL_BARRIER_BITS
		 */
		void dispatch1D(unsigned int global_size, GLbitfield barriers = GL_ALL_BARRIER_BITS) const;

		/**
		 * \brief Short form for glGetUniformLocation(GLuint, const GLchar*)
		 * \param var_name The name of the variable in the shader code.
		 * \return The location of the variable.
		 */
		GLuint getUniformLocation(const char* var_name) const;

		/**
		 * \brief Finds the uniform location of the texture, assigns and binds it to that one and uploads the texture data if desired.
		 * \param texture The texture to be bound.
		 * \param name The source uniform variable name.
		 * \param access Access parameter, e.g. GL_WRITE_ONLY or GL_READ_ONLY, matching to your declaration in shader-code.
		 * \param upload Whether you want to upload the texture data.
		 * \return The received uniform location and the texture unit (which will be the same value then.)
		 */
		GLuint assignUniformTexture(std::weak_ptr<Texture> texture, const GLchar* name, GLenum access, bool upload = true) const;

		/**
		 
		 @param location The texture to be bound.
		 @param name The source uniform variable name.
		 */;

		/**
		 * \brief Binds a Texture to a uniform location and uploads the texture data if desired.
		 * \param texture The texture to be bound.
		 * \param location The uniform location.
		 * \param access Access parameter, e.g. GL_WRITE_ONLY or GL_READ_ONLY
		 * \param upload Whether you want to upload the texture data.
		 */
		void assignUniformTexture(std::weak_ptr<Texture> texture, const GLuint location, GLenum access, bool upload = true) const;

		/**
		 * \brief Creates a new buffer and binds it to the shader.
		 * \param name A name for the new buffer. This has to correspond with the name it is referred by in the shader.
		 * \param buffer_size The size of the buffer.
		 * \param usage what the expected usage pattern for the buffer is. Only Performance related, see https://www.khronos.org/opengl/wiki/GLAPI/glBufferData#Description
		 * \param data The data pushed to the buffer, defaults to nullptr, if you do not need to upload any data, e.g. because the data is generated on GPU.
		 */
		void createBuffer(std::string name, GLuint buffer_size, GLenum usage, const void* data = nullptr);

		/**
		 * \brief Reuses an already existing buffer, binding it to name in this shader. This way you can reuse buffers, especially buffers from another ComputeShader.
		 * \param buffer_ptr shared_ptr to the Bufferobject, as retrieved by getBuffer.
		 * \param name The name that this buffer is recognized by in this shader.
		 * \return whether the operation was successfull
		 */
		bool assignBuffer(std::shared_ptr<Buffer<>> buffer_ptr, std::string name);

		/**
		 * \brief Assigns all buffers used in this Shader to target, using the same name
		 * \param target the ComputeShader to assign the buffers to
		 */
		void assignBuffersTo(ComputeShader& target);

		/**
		 * \brief Returns a shared_ptr to the Buffer recognized by name.
		 * \param name the name of the Buffer.
		 * \return shared_ptr to the Buffer.
		*/
		std::shared_ptr<Buffer<>> getBuffer(std::string name)
		{
			if(m_buffers.count(name) == 0)
			{
				m_buffers[name] = std::make_shared<Buffer<>>();
			}
			return m_buffers[name];
		};

		/**
		 * \brief Reads out a given output buffer.
		 * \param name The name of the buffer to read.
		 * \param offset The offset from where to begin reading, in bytes.
		 * \param size The number of bytes to read
		 * \param data Pointer to the data where the read will be stored.
		*/
		void readBuffer(std::string name, size_t offset, size_t size, void* data);

		/**
		 * \brief Writes into a given output buffer.
		 * \param name The name of the buffer to write to.
		 * \param name The offset where to start writing to, in bytes.
		 * \param size The number of bytes to write.
		 * \param data Pointer to the data that will be written to the buffer.
		*/
		void writeBuffer(std::string name, size_t offset, size_t size, void* data);

		void updateUniform(std::string name, bool value) const
		{
			use();
			glUniform1i(getUniformLocation(name.c_str()), value);
		}

		void updateUniform(std::string name, int value) const
		{
			use();
			glUniform1i(getUniformLocation(name.c_str()), value);
		}

		void updateUniform(std::string name, unsigned int value) const
		{
			use();
			glUniform1ui(getUniformLocation(name.c_str()), value);
		}

		void updateUniform(std::string name, float value) const
		{
			use();
			glUniform1f(getUniformLocation(name.c_str()), value);
		}

		void updateUniform(std::string name, double value) const
		{
			use();
			glUniform1f(getUniformLocation(name.c_str()), value);
		}

		void updateUniform(std::string name, GLuint64 value) const
		{
			use();
			glUniformui64vNV(getUniformLocation(name.c_str()), 1, &value);
		}

		void updateUniform(std::string name, glm::ivec2 vector) const
		{
			use();
			glUniform2iv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::ivec3 vector) const
		{
			use();
			glUniform3iv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::ivec4 vector) const
		{
			use();
			glUniform4iv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::vec2 vector) const
		{
			use();
			glUniform2fv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::vec3 vector) const
		{
			use();
			glUniform3fv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::vec4 vector) const
		{
			use();
			glUniform4fv(getUniformLocation(name.c_str()), 1, glm::value_ptr(vector));
		}

		void updateUniform(std::string name, glm::mat2 matrix) const
		{
			use();
			glUniformMatrix2fv(getUniformLocation(name.c_str()), 1, GL_FALSE, glm::value_ptr(matrix));
		}

		void updateUniform(std::string name, glm::mat3 matrix) const
		{
			use();
			glUniformMatrix3fv(getUniformLocation(name.c_str()), 1, GL_FALSE, glm::value_ptr(matrix));
		}

		void updateUniform(std::string name, glm::mat4 matrix) const
		{
			use();
			glUniformMatrix4fv(getUniformLocation(name.c_str()), 1, GL_FALSE, glm::value_ptr(matrix));
		}

		void updateUniform(std::string name, std::vector<glm::vec2> vector) const
		{
			use();
			glUniform2fv(getUniformLocation(name.c_str()), sizeof(vector), glm::value_ptr((&vector[0])[0]));
		}

		void updateUniform(std::string name, std::vector<glm::vec3> vector) const
		{
			use();
			glUniform3fv(getUniformLocation(name.c_str()), sizeof(vector), glm::value_ptr((&vector[0])[0]));
		}

		void updateUniform(std::string name, std::vector<glm::vec4> vector) const
		{
			use();
			glUniform4fv(getUniformLocation(name.c_str()), sizeof(vector), glm::value_ptr((&vector[0])[0]));
		}

		GLuint getProgramHandle() const { return m_program_handle; }
		GLuint getShaderHandle() const { return m_compute_shader; }
		glm::ivec3 getLocalSizes() const { return m_local_sizes; };

	private:

		/**
		 * \brief Upload the given shader source code for compilation purposes. Also sets local sizes if available.
		 * \param source The shader source code.
		 */
		void setShaderSource(std::string source);

		/**
		 * \brief Compiles the loaded sources.
		 */
		void compileSource() const;

		/**
		 * \brief Links the compiled shader program.
		 */
		void linkProgram() const;

		/**
		 * Checks whether the current Shader was compiled successfully and logs the result.
		 * Exits the program with error-code 40 if something went wrong.
		 */
		void checkShader() const;

		/**
		 * Checks whether the current Shader Program was linked successfully and logs the result.
		 * Exits the program with error-code 41 if something went wrong.
		 */
		void checkProgram() const;

		GLuint m_program_handle; //!< The shader program id.
		GLuint m_compute_shader; //!< The shader id. 
		BufferMap m_buffers;

		glm::ivec3 m_local_sizes = glm::ivec3(1); //!< Local sizes loaded from shader source layout tag (if available, 1 otherwise for each component).
	};

	constexpr const char* cDefaultVersion = "430"; //!< Default GLSL version code.

	/**
	 * \brief This class contains a couple of helper methods for compute shader source alterations and on-the-fly compilations.
	 */
	class ComputeShaderUtils
	{
	public:

		/**
		 * \brief Generates a GLSL version header.
		 * \param version the GLSL-Version-Code (default: 430)
		 * \return <code>#version (version)</code> (example: <code>#version 430</code>)
		 */
		static std::string generateHeader(const char* version = cDefaultVersion);

		/**
		 * \brief Creates a compute shader input("in") layout line.
		 * \param local_size_x The work group x dimension
		 * \param local_size_y The work group y dimension
		 * \param local_size_z The work group z dimension
		 * \return A string (with a new line at it's end) like <code>layout(local_size_x=1024, local_size_y=1, local_size_z=1) in;</code>
		 */
		static std::string generateLocalSizes(unsigned int local_size_x = 1, unsigned int local_size_y = 1, unsigned int local_size_z = 1);

		/**
		 * Loads a compute shader file from disk. All <code>#include</code> lines in the shadercode will be replaced by their corresponding file contents here.
		 * That way, you can load multiple files, while giving only one file name.
		 * \param file The shader file path.
		 * \param included Map of already included shader files.
		 * \return The final file content of the given file, as well as that of all referenced includes.
		 */
		static std::string loadFile(const char* file, std::map<std::string, int>& included);

		/**
		 * Loads a compute shader file from disk. All <code>#include</code> lines in the shadercode will be replaced by their corresponding file contents here.
		 * That way, you can load multiple files, while giving only one file name.
		 * \param file The shader file path.
		 * \return The final file content of the given file, as well as that of all referenced includes.
		 */
		static std::string loadFile(const char* file);

		/**
		 * \brief Loads some compute shader files from disk and merges them to a whole shader source. All <code>#include</code> lines will be replaced by their corresponding file contents here.
		 * \param files A vector of files to be merged into one source.
		 * \param included Map of already included shader files. It's best to leave this parameter out to pass an empty map.
		 * \return The final file content of the given files, as well as that of all referenced includes.
		 */
		static std::string loadFiles(std::vector<const char*> files, std::map<std::string, int>& included);

		/**
		 * \brief Loads some compute shader files from disk and merges them to a whole shader source. All <code>#include</code> lines will be replaced by their corresponding file contents here.
		 * \param files A vector of files to be merged into one source.
		 * \return The final file content of the given files, as well as that of all referenced includes.
		 */
		static std::string loadFiles(std::vector<const char*> files);

		/**
		 * \brief Takes the source string, finds the given function name and creates a <code>void main()</code> which only contains a call of the function.
		 * \param source The full shader source code containing at least the needed function and all buffers/uniforms it needs.
		 * \param function_name The name of the function to be called in main.
		 * \return The source with the appended main function.
		 */
		static std::string createMain(std::string source, const char* function_name);

		/**
		 * Takes the source string, finds all of the given function names and creates a <code>void main()</code> which contains one call of each of the function names
		 * in the order they are stored.
		 * \param source The full shader source code containing at least the needed functions and all buffers/uniforms they need.
		 * \param function_names A vector of all functions which will be called in the new main in the needed order. (0 will be called first, etc.)
		 * \return The source with the appended main function.
		 */
		static std::string createMain(std::string source, std::vector<const char*> function_names);

	private:
		/**
		 * \brief Finds and replaces all <code>#include</code>-lines in the given source by the files referenced there, if it's not been included yet.
		 * \param source The shader source. It will have replaced the included file contents.
		 * \param included A map containing all already included files.
		 */
		static void replaceIncludes(std::string& source, std::map<std::string, int>& included);
	};
}

#endif // __GLWRAP_H
