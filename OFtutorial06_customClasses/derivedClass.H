#include "fvCFD.H"

// use a standard C++ stringstream object
#include <sstream>

// This class is derived from the IOdictionary, declared in
// $FOAM_SRC/OpenFOAM/db/IOobjects/IOdictionary .
// We derive all the methods by setting maximum access level as public, thus allowing
// all the inherited methods to be accessed in the same way as they were for
// the base class.
class myDict : public IOdictionary
{
public:
    // Main constructor which accepts an IOobject (a file)
    myDict(const IOobject& ioObj);
    ~myDict();

    // custom method which prints the tokens (variables) defined in the file
    void printTokensInTheDict() const;
};
