#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include "OMtools.hpp"
#include "researchlib.hpp"
#include "programs.hpp"

int main() 
{
    return programs::always_weakly_abstractly_reducible<3, 8>();
}