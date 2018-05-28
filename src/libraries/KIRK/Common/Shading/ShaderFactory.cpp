/*
 * RT
 * ShaderFactory.cpp
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the ShaderFactory class
 *
 * Copyright (c) 2017 Hendrik Schwanekamp
 *
 */

// includes
//--------------------
#include <KIRK/Utils/Log.h>
#include "ShaderFactory.h"
#include "Shader.h"
//--------------------

// namespace
//--------------------
namespace KIRK {
//--------------------

// function definitions of the ShaderFactory class
//-------------------------------------------------------------------

ShaderFactory& ShaderFactory::getInstance()
{
    static ShaderFactory instance;
    return instance;
}

void ShaderFactory::registerShader(const std::string& name, const std::function<std::shared_ptr<Shader>()>& creator)
{
    if(shaderLookupTable.count(name) != 0)
        throw std::invalid_argument("You have multiple Shaders with the same name! Don't do that!");

    shaderLookupTable.emplace(name, ShaderMapContent(creator, std::shared_ptr<Shader>()));
}

std::shared_ptr<Shader> ShaderFactory::getShader(const std::string& name)
{
    if(shaderLookupTable.count(name) != 1)
    {
        LOG_ERRORs() << "The Shader " << name << " does not exist. Please fix your scene or implement the Shader (See Tutorialson the wiki page)";
        throw std::invalid_argument("Shader "+name+" existiert nicht");
    }
    ShaderMapContent& pair = shaderLookupTable.at(name);
    if ( pair.second.expired())
    {
        std::shared_ptr<Shader> ptr = pair.first();
        pair.second = ptr;
        return ptr;
    }
    else
        return pair.second.lock();
}

void ShaderFactory::doForAllShaders(std::function<void(const std::string &)> func)
{
    for(auto &&item :shaderLookupTable)
    {
        func(item.first);
    }
}

}