/*
 * RT
 * BsdfFactory.h
 *
 * @author: Hendrik Schwanekamp
 * @mail:   hendrik.schwanekamp@gmx.net
 *
 * Implements the BsdfFactory class
 *
 * Copyright (c) 2017 Hendrik Schwanekamp
 *
 */

#ifndef RT_BSDFFACTORY_H
#define RT_BSDFFACTORY_H

// includes
//--------------------
#include <string>
#include <functional>
#include <memory>
#include <map>
//--------------------

// namespace
//--------------------
namespace KIRK {
//--------------------

class BSDF;

//-------------------------------------------------------------------
/**
 * class BsdfFactory
 *
 * usage:
 *
 */
class BsdfFactory
{
public:
    static BsdfFactory &getInstance();

    BsdfFactory(const BsdfFactory &other) = delete;
    BsdfFactory(BsdfFactory &&other) = delete;
    BsdfFactory operator=(BsdfFactory &other)= delete;
    BsdfFactory operator=(BsdfFactory &&other)= delete;

    void registerBSDF(const std::string &name, const std::shared_ptr<BSDF>& bsdf);
    std::shared_ptr<BSDF> getBsdf(const std::string& name);
    void doForAllBsdfs(std::function<void(const std::string &)> func);

private:
    BsdfFactory()= default;
    std::map<std::string, std::shared_ptr<BSDF>> bsdfLookup;

};

template <class bsdf_type>
class BsdfRegistrator
{
public:
    BsdfRegistrator(const std::string& sName)
    {
        BsdfFactory::getInstance().registerBSDF(sName, std::make_shared<BSDF>(sName, bsdf_type::localSample, bsdf_type::evaluateLight));
    }
};

}
#endif //RT_BSDFFACTORY_H
