#include <iostream>

#include "particleutils.h"

int main()
{
    std::cout << terminalSettlingVelocity( 0.09e-3*meter, 2650 * kilogrammes_per_cubic_metre,
                                          1000 * kilogrammes_per_cubic_metre) << std::endl;
    std::cout << terminalSettlingVelocity( 0.5e-3*meter, 2650 * kilogrammes_per_cubic_metre,
                                          1000 * kilogrammes_per_cubic_metre) << std::endl;
    std::cout << terminalSettlingVelocity( 2.e-3*meter, 2650 * kilogrammes_per_cubic_metre,
                                          1000 * kilogrammes_per_cubic_metre) << std::endl;
    
    return 0;
}
