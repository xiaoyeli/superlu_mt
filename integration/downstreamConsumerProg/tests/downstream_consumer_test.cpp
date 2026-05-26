#include "downstream_consumer.h"

#include <cmath>
#include <stdexcept>

int main()
{
    const double x = downstreamConsumerProg::solve_with_runtime_dispatch(
        downstreamConsumerProg::Precision::Double);
    const double y = downstreamConsumerProg::solve_with_runtime_dispatch(
        downstreamConsumerProg::Precision::DoubleDouble);

    if (std::fabs(x - 8.0) > 1e-12) {
        throw std::runtime_error("unexpected double result");
    }
    if (std::fabs(y - 8.0) > 1e-12) {
        throw std::runtime_error("unexpected double-double result");
    }
    return 0;
}
