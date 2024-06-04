#include "params.hpp"
#include "hgsw.hpp"

using namespace hgsw;
using Params = hgsw_params::prm;

/* Driver program */
int main() { 
    
    HGSW<Params>(); // call Half-GSW encryption scheme
 
    return 0;
}
