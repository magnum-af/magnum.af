#pragma once
#include "zee.hpp"
#include "arrayfire.h"

namespace magnumafcpp
{

class AtomisticExternalField : public ExternalField
{
public:
    using ExternalField::ExternalField;
    double E(const State &state) override;                     //Energy contribution
    double E(const State &state, const af::array &h) override; ///< Calculating the micromagnetic energy for a already calculated h field

private:
};

} // namespace magnumafcpp
