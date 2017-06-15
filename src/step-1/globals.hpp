#ifndef H_GLOBALS__
#define H_GLOBALS__

//#include <grvy.h>

#define PRINT(x) {std::cout << #x <<": "<<x << std::endl;}

namespace heat{

typedef double real;  //Interested in changing from double to floats  

enum myBoundaryID {dir_Minus, dir_Plus, neu_Minus, neu_Plus};
}

#endif
