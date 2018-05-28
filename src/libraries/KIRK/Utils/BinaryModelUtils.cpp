#include "BinaryModelUtils.h"
#include "externals/stb_image/stb_image.h"

std::shared_ptr<KIRK::Mesh> KIRK::loadBinaryGeometry(std::string path)
{
	LOG_INFO("Starting to read binary 3d file...");

	//Give ownership to the caller afterwards
	std::shared_ptr<KIRK::Mesh> geometry(new KIRK::Mesh());

	//do read
	int type, numPoints, numIndices, numFaces;
	std::ifstream file(path, std::ios::in | std::ios::binary);
	if (!file.is_open())
	{
		LOG_ERROR("Can't open file: '%' ", path.c_str());
		return 0;
	}
	file.seekg(0, file.end);
	long size = file.tellg();
	file.seekg(0, file.beg);

	std::vector<unsigned char> buffer(size);
	file.read((char*)buffer.data(), size);

	std::vector<unsigned char> uncompressed_data;

	if (strlen((char*)&buffer[0]) <= 1) //=> old Binary File
	{
		LOG_INFO("Uncompressed binary 3d file...");
		file.seekg(0, file.beg);
		file.read((char *)&type, sizeof(int));
		file.read((char *)&numPoints, sizeof(int));
		file.read((char *)&numIndices, sizeof(int));
	}
	else //=> new Binary File
	{
		lodepng::decompress(uncompressed_data, buffer);
		////////////////////////////////////////////////////////////////////////////////////////////////////
		//////
		//////	From the decompressed sequence of chars, we start at an index (e.g. 0). Then we just assume
		//////	that the pointer we have there is not a pointer to a char, but one to an int.
		//////	From that, we just access it via pointer de-referencing.
		//////
		////////////////////////////////////////////////////////////////////////////////////////////////////
		LOG_INFO("Compressed binary 3d file...");
		type = *((int*)&uncompressed_data.data()[0]);
		numPoints = *((int*)&uncompressed_data.data()[sizeof(int)]);
		numIndices = *((int*)&uncompressed_data.data()[2 * sizeof(int)]);
	}
	numFaces = numIndices / 3;

	LOG_INFO("Found type: %", type);
	LOG_INFO("Found % points.", numPoints);
	LOG_INFO("Found % indices.", numIndices);
	LOG_INFO("Found % faces.", numFaces); //needed for Iteration to load material indices

	int datNum = numPoints * 3;
	if (type & 1)
		datNum += numPoints * 3;
	if (type & 2)
		datNum += numPoints * 2;

	KIRK::Mesh::vertex vertex; //construct vertex object to initialize struct members
	if (strlen((char*)&buffer[0]) <= 1) //=> old Binary File
	{
		float *vData = new float[datNum];
		file.read((char *)vData, sizeof(float) * datNum);
		unsigned int *iData = new unsigned int[numIndices];
		file.read((char *)iData, sizeof(unsigned int) * numIndices);
		unsigned int *mData = new unsigned int[numFaces];
		file.read((char *)mData, sizeof(unsigned int) * numFaces);

		file.close();

		for (int t = 0; t < numPoints; t++)
		{
			vertex.position = glm::vec3(vData[t * 3], vData[t * 3 + 1], vData[t * 3 + 2]);
			geometry->m_vertices.push_back(vertex);
		}
		int o = numPoints * 3;
		if (type & 1)
		{
			for (int t = 0; t < numPoints; t++)
			{
				geometry->m_vertices[t].normal = glm::vec3(vData[t * 3 + o], vData[t * 3 + 1 + o], vData[t * 3 + 2 + o]);
			}
		}
		o += numPoints * 3;
		if (type & 2)
		{
			for (int t = 0; t < numPoints; t++)
			{
				geometry->m_vertices[t].texcoord = glm::vec2(vData[t * 2 + o], vData[t * 2 + 1 + o]);
			}
		}
		for (int t = 0; t < numIndices; t++)
		{
			geometry->m_indices.push_back(iData[t]);
		}

		for (int i = 0; i < numIndices - 2; i += 3)
		{
			KIRK::Mesh::face face;
			face.vertex_index_a = geometry->m_indices[i];
			face.vertex_index_b = geometry->m_indices[i + 1];
			face.vertex_index_c = geometry->m_indices[i + 2];
			face.material_index = mData[i / 3]; //old binary file => every material_index = 0
			geometry->m_faces.push_back(face);
		}
		delete[] vData;
		delete[] iData;
	}
	else //=> new Binary File
	{
		std::vector<float> vData(datNum);
		memcpy(vData.data(), (float*)(&uncompressed_data.data()[3 * sizeof(int)]), datNum * sizeof(float));

		std::vector<unsigned int> iData(numIndices);
		memcpy(iData.data(), (unsigned int*)(&uncompressed_data.data()[3 * sizeof(int) + datNum * sizeof(float)]), numIndices * sizeof(unsigned int));

		std::vector<unsigned int> materialIndex(numFaces);
		memcpy(materialIndex.data(), (unsigned int*)(&uncompressed_data.data()[3 * sizeof(int) + datNum * sizeof(float) + numIndices * sizeof(unsigned int)]), numFaces * sizeof(unsigned int));

		file.close();

		std::vector<unsigned int> matIndices;

		for (int t = 0; t < numPoints; t++)
		{
			vertex.position = glm::vec3(vData[t * 3], vData[t * 3 + 1], vData[t * 3 + 2]);
			geometry->m_vertices.push_back(vertex);
		}
		int o = numPoints * 3;
		if (type & 1)
		{
			for (int t = 0; t < numPoints; t++)
			{
				geometry->m_vertices[t].normal = glm::vec3(vData[t * 3 + o], vData[t * 3 + 1 + o], vData[t * 3 + 2 + o]);
			}
		}
		o += numPoints * 3;
		if (type & 2)
		{
			for (int t = 0; t < numPoints; t++)
			{
				geometry->m_vertices[t].texcoord = glm::vec2(vData[t * 2 + o], vData[t * 2 + 1 + o]);
			}
		}
		for (int t = 0; t < numIndices; t++)
		{
			geometry->m_indices.push_back(iData[t]);
		}
		for (int t = 0; t < numFaces; t++)
		{
			matIndices.push_back(materialIndex[t]);
		}

		for (int i = 0; i < numIndices - 2; i += 3)
		{
			KIRK::Mesh::face face;
			face.vertex_index_a = geometry->m_indices[i];
			face.vertex_index_b = geometry->m_indices[i + 1];
			face.vertex_index_c = geometry->m_indices[i + 2];
			face.material_index = matIndices[i / 3];
			geometry->m_faces.push_back(face);
		}
	}

	geometry->m_materials.emplace_back(new KIRK::Material("default"));

	LOG_INFO("Finished reading binary 3d file of '%' ", path.c_str());

	return geometry;
}

int KIRK::writeBinaryGeometry(std::string path, const KIRK::Mesh *geom, bool invertNormals, bool compress)
{
	std::size_t point = path.find_last_of(".");
	path = path.erase(point, path.size());
	path.append(".b3df");

	LOG_INFO("Starting to write binary 3d file...");

	int numPoints = geom->m_vertices.size();
	int numIndices = geom->m_indices.size();
	int numFaces = geom->m_faces.size();

	//fill vectors in order to access normals and texcoords for writing and to make iteration in for-loops possible
	std::vector<glm::vec3> positions, normals;
	std::vector<glm::vec2> texcoords;
	std::vector<unsigned int> matIndices;
	for (int t = 0; t < numPoints; t++)
	{
		positions.push_back(geom->m_vertices[t].position);
		normals.push_back(geom->m_vertices[t].normal);
		texcoords.push_back(geom->m_vertices[t].texcoord);
	}
	for (int t = 0; t < numFaces; t++)
	{
		matIndices.push_back(geom->m_faces[t].material_index);
	}

	//do write
	int type = 0;
	if (normals.size() == numPoints) //check if amount of normals == numPoints
		type |= 1;
	if (texcoords.size() == numPoints) //check if amount of UV-Coord == numPoints
		type |= 2;

	float normalInverter = 1.0f;
	if (invertNormals)
		normalInverter = -1.0f;

	std::vector<float> vData;
	int datNum = numPoints * 3;
	for (glm::vec3 v : positions)
	{
		vData.push_back(v.x);
		vData.push_back(v.y);
		vData.push_back(v.z);
	}
	if (type & 1)
	{
		datNum += numPoints * 3;
		for (glm::vec3 v : normals)
		{
			vData.push_back(v.x * normalInverter);
			vData.push_back(v.y * normalInverter);
			vData.push_back(v.z * normalInverter);
		}
	}
	if (type & 2)
	{
		datNum += numPoints * 2;
		for (glm::vec2 v : texcoords)
		{
			vData.push_back(v.x);
			vData.push_back(v.y);
		}
	}
	std::vector<unsigned int> iData;
	for (unsigned int i : geom->m_indices)
	{
		iData.push_back(i);
	}

	std::ofstream file(path, std::ios::out | std::ios::binary);
	if (compress)
	{
		std::vector<unsigned char> data;
		//data.reserve(3 * sizeof(int) + datNum * sizeof(float) + numIndices * sizeof(unsigned int));
		data.insert(data.end(), (unsigned char*)&type, (unsigned char*)&type + sizeof(int));
		data.insert(data.end(), (unsigned char*)&numPoints, (unsigned char*)&numPoints + sizeof(int));
		data.insert(data.end(), (unsigned char*)&numIndices, (unsigned char*)&numIndices + sizeof(int));
		data.insert(data.end(), (unsigned char*)vData.data(), (unsigned char*)vData.data() + datNum * sizeof(float));
		data.insert(data.end(), (unsigned char*)iData.data(), (unsigned char*)iData.data() + numIndices * sizeof(unsigned int));

		data.insert(data.end(), (unsigned char*)matIndices.data(), (unsigned char*)matIndices.data() + numFaces * sizeof(unsigned int));

		std::vector<unsigned char> compressed_data;
		lodepng::compress(compressed_data, data);

		file.write((char*)compressed_data.data(), compressed_data.size());

		LOG_INFO("Compressed data size: %", compressed_data.size());
	}
	else
	{
		file.write((char *)&type, sizeof(int));
		file.write((char *)&numPoints, sizeof(int));
		file.write((char *)&numIndices, sizeof(int));

		file.write((char *)vData.data(), datNum * sizeof(float));
		file.write((char *)iData.data(), numIndices * sizeof(unsigned int));

		file.write((char *)matIndices.data(), numFaces * sizeof(unsigned int));
	}

	file.close();

	LOG_INFO("Wrote % Points.", numPoints);
	LOG_INFO("Wrote % Indices.", numIndices);
	LOG_INFO("Wrote % VertexData objects.", vData.size());
	LOG_INFO("Wrote % IndexData objects.", iData.size());
	LOG_INFO("Wrote % MaterialIndices.", matIndices.size());
	LOG_INFO("Finished writing binary 3d file to '%' ", path.c_str());

	return 0;
}
