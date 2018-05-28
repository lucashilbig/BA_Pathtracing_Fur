/*
 * RT
 * ShaderFactory.h
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the ShaderFactory class
 *
 * Copyright (c) 2017 Hendrik Schwanekamp
 *
 */

#ifndef RT_SHADERFACTORY_H
#define RT_SHADERFACTORY_H

// includes
//--------------------
#include <map>
#include <memory>
#include <functional>

//--------------------

// namespace
//--------------------
namespace KIRK {
//--------------------

class Shader;

//-------------------------------------------------------------------
/**
 * class ShaderFactory
 * use this to create shaders from there string name
 *
 *
 */
class ShaderFactory
{
public:
    static ShaderFactory& getInstance();

    ShaderFactory(const ShaderFactory& other)=delete;
    ShaderFactory(ShaderFactory&& other)=delete;
    ShaderFactory operator=(ShaderFactory& other)= delete;
    ShaderFactory operator=(ShaderFactory&& other)= delete;

    void registerShader(const std::string& name, const std::function <std::shared_ptr<Shader>()>& creator);
    std::shared_ptr<Shader> getShader(const std::string& name);
    void doForAllShaders(std::function<void(const std::string &)> func);

private:
    ShaderFactory()= default;

    typedef std::pair<std::function<std::shared_ptr<Shader>()>, std::weak_ptr<Shader>> ShaderMapContent;
    std::map< std::string, ShaderMapContent> shaderLookupTable;
};

/**
 * Helper to register shaders at compile time
 */
template <class shader_type>
class ShaderRegistrator
{
public:
    ShaderRegistrator(const std::string& sName)
    {
        ShaderFactory::getInstance().registerShader(sName, [sName](){
            return std::dynamic_pointer_cast<Shader>(std::make_shared<shader_type>(sName));
        });
    }
};

}

#endif //RT_SHADERFACTORY_H
