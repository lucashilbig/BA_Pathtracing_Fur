
#ifndef __BINARY_MODEL_UTILS_H__
#define __BINARY_MODEL_UTILS_H__

#include <string>
#include <map>
#include "../Common/Mesh.h"

namespace KIRK
{
	/**
	Loads a .b3df binary file from the given path and saves it as a std::shared_ptr<Mesh>.
	@param path The path containing the binary file.
	@return A std::shared_ptr<Mesh> containing all mesh parameters read from the binary file.
	*/
	std::shared_ptr<Mesh> loadBinaryGeometry(std::string path);

	/**
	Writes a .b3df binary file by accessing the parameters of a given Mesh.
	@param path The path the binary file is written to (appends .b3df at the end).
	@param geom A pointer to the Mesh with all the parameters relevant for the binary file.
	@param invertNormals An option for inverting the normals, set false by default.
	@param compress An option to compress the written file via lodepng. Set true by default.
	*/
	int writeBinaryGeometry(std::string path, const KIRK::Mesh* geom, bool invertNormals = false, bool compress = true);
}


#endif //__BINARY_MODEL_UTILS_H__
