#pragma once

namespace downstreamConsumerProg {

enum class Precision
{
    Double,
    DoubleDouble
};

double solve_with_runtime_dispatch(Precision precision);

}  // namespace downstreamConsumerProg
