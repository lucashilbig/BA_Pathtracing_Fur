/*
 * RT
 * BsdfFactory.cpp
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the BsdfFactory class
 *
 * Copyright (c) 2017 Hendrik Schwanekamp
 *
 */

// includes
//--------------------
#include "BsdfFactory.h"
#include "KIRK/Utils/Log.h"
//--------------------

// namespace
//--------------------
namespace KIRK {
//--------------------

// function definitions of the BsdfFactory class
//-------------------------------------------------------------------

BsdfFactory &BsdfFactory::getInstance()
{
    static BsdfFactory instance;
    return instance;
}

void BsdfFactory::registerBSDF(const std::string &name, const std::shared_ptr<BSDF>& bsdf)
{
    bsdfLookup.emplace(name, bsdf);
}

std::shared_ptr<BSDF> BsdfFactory::getBsdf(const std::string &name)
{
    if(bsdfLookup.count(name) != 1)
    {
        LOG_ERRORs() << "The Bsdf " << name << " does not exist. Please fix your scene or implement the BSDF (See Tutorialson the wiki page)";
        throw std::invalid_argument("BSDF "+name+" existiert nicht");
    }
    return bsdfLookup.at(name);
}

void BsdfFactory::doForAllBsdfs(std::function<void(const std::string &)> func)
{
    for(auto &&bsdf : bsdfLookup)
    {
        func(bsdf.first);
    }
}

}