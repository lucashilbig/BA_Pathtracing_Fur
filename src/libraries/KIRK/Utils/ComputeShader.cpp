#include "ComputeShader.h"

#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <regex>


namespace KIRK
{
	ComputeShader::ComputeShader(std::vector<const char*> shaderFiles) : ComputeShader(ComputeShaderUtils::loadFiles(shaderFiles))
	{

	}

	ComputeShader::ComputeShader(std::string shaderSource)
	{
		//Create the handles
		m_program_handle = glCreateProgram();
		m_compute_shader = glCreateShader(GL_COMPUTE_SHADER);

		//Build program.
		setShaderSource(shaderSource);
		compileSource();
		linkProgram();
	}

	ComputeShader::~ComputeShader()
	{
		glDeleteShader(m_compute_shader);
		glDeleteProgram(m_program_handle);
	}

	void ComputeShader::setShaderSource(std::string source)
	{
		//Now loading and parsing the optimal local sizes.
		std::regex reg("(local_size_)([\\w])+\\s*=\\s*([0-9]+)");
		auto words_begin = std::sregex_iterator(source.begin(), source.end(), reg);
		auto words_end = std::sregex_iterator();

		for (std::sregex_iterator i = words_begin; i != words_end; ++i)
		{
			std::smatch match = *i;
			if (std::string(match[2]) == "x")
			{
				m_local_sizes.x = std::stoi(match[3]);
			}
			else if (std::string(match[2]) == "y")
			{
				m_local_sizes.y = std::stoi(match[3]);
			}
			else if (std::string(match[2]) == "z")
			{
				m_local_sizes.z = std::stoi(match[3]);
			}
		}

		LOG_DEBUG("Your value in x is %", m_local_sizes.x);
		LOG_DEBUG("Your value in y is %", m_local_sizes.y);
		LOG_DEBUG("Your value in z is %", m_local_sizes.z);


		const char* shader_source = source.c_str();
		const GLint source_size = strlen(shader_source);

		glShaderSource(m_compute_shader, 1, &shader_source, &source_size);
	}

	GLuint ComputeShader::assignUniformTexture(std::weak_ptr<Texture> texture, const GLchar* name, GLenum access, bool upload) const
	{
		use();
		GLuint texture_handle = getUniformLocation(name);
		assignUniformTexture(texture, texture_handle, access);
		return texture_handle;
	}

	void ComputeShader::assignUniformTexture(std::weak_ptr<Texture> texture, const GLuint location, GLenum access, bool upload) const
	{
		use();
		glUniform1i(location, location);

		//Both textures here have to have the same format (both RGBA, not one RGB and one RGBA).
		//Otherwise the output will be somewhat corrupt.
		//Source texture here is a simple image loaded from memory.
		if (upload)
			texture.lock()->glUpload();
		//Bind the texture to the according texture unit (see above) as read_only.
		texture.lock()->glBind(true, location, access);
	}

	GLuint ComputeShader::getUniformLocation(const char* var_name) const
	{
		return glGetUniformLocation(m_program_handle, var_name);
	}

	void ComputeShader::createBuffer(std::string name, GLuint buffer_size, GLenum usage, const void* data)
	{
		//Re-use existing buffer or create new one if not existing.
		auto buffer = getBuffer(name);
		if (assignBuffer(buffer, name))
		{
			buffer->bind();
			Buffer<>::init(buffer_size, data, usage);
		}
	}

	bool ComputeShader::assignBuffer(std::shared_ptr<Buffer<>> buffer_ptr, std::string name)
	{
		//Find buffer index and map it to the buffer_id.
		GLuint buffer_index = glGetProgramResourceIndex(m_program_handle, GL_SHADER_STORAGE_BLOCK, name.c_str());

		if (buffer_index != GL_INVALID_INDEX)
		{
			m_buffers[name] = buffer_ptr;
			//simply map the index to itself.
			glShaderStorageBlockBinding(m_program_handle, buffer_index, buffer_index);
			buffer_ptr->bindBufferBase(buffer_index);
			return true;
		}
		return false;
	}

	void ComputeShader::assignBuffersTo(ComputeShader &target)
	{
		for (auto &buffer : m_buffers)
		{
			target.assignBuffer(buffer.second, buffer.first);
		}
	}

	void ComputeShader::readBuffer(std::string name, size_t offset, size_t size, void* data)
	{
		//Bind and read.
		m_buffers[name]->bind();
		void* output = Buffer<>::map(offset, size, GL_MAP_READ_BIT);

		memcpy(data, output, size);

		Buffer<>::unmap();
	}

	void ComputeShader::writeBuffer(std::string name, size_t offset, size_t size, void* data)
	{
		//Bind and read.
		m_buffers[name]->bind();
		void* output = Buffer<>::map(offset, size, GL_MAP_WRITE_BIT);

		memcpy(output, data, size);

		Buffer<>::unmap();
	}

	void ComputeShader::use() const
	{
		glUseProgram(m_program_handle);
	}

	void ComputeShader::dispatch3D(unsigned int global_size_x, unsigned int global_size_y, unsigned int global_size_z, GLbitfield barriers) const
	{
		use();
		int x_count = global_size_x%m_local_sizes.x == 0 ? global_size_x / m_local_sizes.x : global_size_x / m_local_sizes.x + 1;
		int y_count = global_size_y%m_local_sizes.y == 0 ? global_size_y / m_local_sizes.y : global_size_y / m_local_sizes.y + 1;
		int z_count = global_size_z%m_local_sizes.z == 0 ? global_size_z / m_local_sizes.z : global_size_z / m_local_sizes.z + 1;
		LOG_DEBUG("Dispatching program on x=%, y=% and z=% work groups.", x_count, y_count, z_count);
		//Create as many 16x16 blocks as you can fit into the image and add one on
		//each dimension to include the borders under all circumstances.
		glDispatchCompute(x_count, y_count, z_count);
		//Lock OpenGL rendering. And Texture access.
		glMemoryBarrier(barriers);
	}

	void ComputeShader::dispatch2D(unsigned int global_size_x, unsigned int global_size_y, GLbitfield barriers) const
	{
		dispatch3D(global_size_x, global_size_y, 1, barriers);
	}

	void ComputeShader::dispatch1D(unsigned int global_size, GLbitfield barriers) const
	{
		dispatch2D(global_size, 1, barriers);
	}

	void ComputeShader::compileSource() const
	{
		glCompileShader(m_compute_shader);
		checkShader();
		glAttachShader(m_program_handle, m_compute_shader);
	}

	void ComputeShader::linkProgram() const
	{
		glLinkProgram(m_program_handle);
		checkProgram();
	}

	//checks a shader for compiler errors
	void ComputeShader::checkShader() const
	{
		GLint status;
		glGetShaderiv(m_compute_shader, GL_COMPILE_STATUS, &status);

		if (status == GL_FALSE)
		{
			GLint infoLogLength;
			glGetShaderiv(m_compute_shader, GL_INFO_LOG_LENGTH, &infoLogLength);

			GLchar *infoLog = new GLchar[infoLogLength + 1];
			glGetShaderInfoLog(m_compute_shader, infoLogLength, nullptr, infoLog);

			LOG_ERROR("Unable to compile ComputeShader: %", infoLog);
			delete[] infoLog;
			std::cout << "Press [ENTER] to exit...";
			std::cin.ignore();
			exit(40);
		}
		else
		{
			LOG_INFO("ComputeShader compiled!");
		}
	}

	//checks a program
	void ComputeShader::checkProgram() const
	{
		GLint status;
		glGetProgramiv(m_program_handle, GL_LINK_STATUS, &status);

		if (status == GL_FALSE)
		{
			GLint infoLogLength;
			glGetProgramiv(m_program_handle, GL_INFO_LOG_LENGTH, &infoLogLength);

			GLchar *infoLog = new GLchar[infoLogLength + 1];
			glGetProgramInfoLog(m_program_handle, infoLogLength, nullptr, infoLog);

			LOG_ERROR("Unable to link ComputeShader Program: %", infoLog);
			delete[] infoLog;
			std::cout << "Press [ENTER] to exit...";
			std::cin.ignore();
			exit(41);
		}
		else
		{
			LOG_INFO("ComputeShader Program linked!");
		}
	}

	std::string ComputeShaderUtils::generateHeader(const char* version)
	{
		return "#version " + std::string(version) + "\n";
	}

	std::string ComputeShaderUtils::generateLocalSizes(unsigned int local_size_x, unsigned int local_size_y, unsigned int local_size_z)
	{

		std::string lsx = "local_size_x = " + std::to_string(local_size_x);
		std::string lsy = "local_size_y = " + std::to_string(local_size_y);
		std::string lsz = "local_size_z = " + std::to_string(local_size_z);

		return "layout (" + lsx + "," + lsy + "," + lsz + ") in;\n";
	}

	std::string ComputeShaderUtils::loadFiles(std::vector<const char*> files, std::map<std::string, int> &included)
	{
		std::string fileContent = "";

		for (const char* file : files)
		{
			fileContent += loadFile(file, included);
		}

		return fileContent;
	}

	std::string ComputeShaderUtils::loadFiles(std::vector<const char*> files)
	{
		std::map<std::string, int> included = std::map<std::string, int>();
		return loadFiles(files, included);
	}

	std::string ComputeShaderUtils::loadFile(const char *file_name, std::map<std::string, int> &included)
	{
		if (included[std::string(file_name)] > 0) {
			//File already included in this run.
			return "";
		}

		/*if (std::string(file_name) == std::string(SHADERS_PATH "/compute/Pathtracer/bsdf_base.compute"))
		{
			return KIRK::GPU::BsdfBuilder::build();
		}*/

		std::string fileContent = "";
		std::string line;
		std::ifstream file(file_name);
		if (file.is_open())
		{
			while (!file.eof())
			{
				getline(file, line);
				fileContent += line + "\n";
			}
			file.close();
		}
		else
		{
			LOG_ERROR("Could not open Shader file: %", file_name);
		}

		replaceIncludes(fileContent, included);
		return fileContent;
	}

	std::string ComputeShaderUtils::loadFile(const char* file)
	{
		std::map<std::string, int> included = std::map<std::string, int>();
		return loadFile(file, included);
	}

	void ComputeShaderUtils::replaceIncludes(std::string &source, std::map<std::string, int> &included) {
		std::regex reg("#include \"([\\w|\\|/|.]*)\"");

		std::string s = source;

		std::sregex_iterator i = std::sregex_iterator(s.begin(), s.end(), reg);

		while (i != std::sregex_iterator()) {
			s = source;

			std::smatch match = *i;

			//Make path from relative to absolute.
			std::string include_line = std::string(SHADERS_PATH "/compute/") + std::string(match[1]);

			source = std::regex_replace(s, std::regex(std::string(match[0])), loadFile(include_line.c_str(), included));

			s = source;

			included[include_line]++;
			i = std::sregex_iterator(s.begin(), s.end(), reg);

		} ;
	}

	std::string ComputeShaderUtils::createMain(std::string function_file_content, std::vector<const char *> function_names)
	{
		std::string functions = "";
		for (const char* function_name : function_names)
		{
			std::regex reg("void ([\\w]*)");
			auto words_begin =
				std::sregex_iterator(function_file_content.begin(), function_file_content.end(), reg);
			auto words_end = std::sregex_iterator();

			bool found = false;

			for (std::sregex_iterator i = words_begin; i != words_end; ++i)
			{
				std::smatch match = *i;
				if (std::string(function_name) == std::string(match[1]))
				{
					found = true;
					break;
				}
			}

			if (!found)
			{
				LOG_ERROR("Could not find function %!", function_name);
				continue;
			}

			functions += "\n" + std::string(function_name) + "();";
		}

		function_file_content += "\nvoid main(){\n" + functions + "\n}\n";
		
		return function_file_content;
	}

	std::string ComputeShaderUtils::createMain(std::string function_file_content, const char* function_name)
	{
		//open file and "parse" input

		std::regex reg("void ([\\w]*)");
		auto words_begin =
			std::sregex_iterator(function_file_content.begin(), function_file_content.end(), reg);
		auto words_end = std::sregex_iterator();

		bool found = false;

		for (std::sregex_iterator i = words_begin; i != words_end; ++i)
		{
			std::smatch match = *i;
			if (std::string(function_name) == std::string(match[1]))
			{
				found = true;
				break;
			}
		}

		if (!found)
		{
			LOG_ERROR("Could not find function %!", function_name);
		}

		function_file_content += "\nvoid main(){\n" + std::string(function_name) + "();\n}";

		return function_file_content;
	}
}
