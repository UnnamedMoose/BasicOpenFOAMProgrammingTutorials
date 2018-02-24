#include "customClass.H"

customClass::customClass()
{
    myInt_= 0;
}

customClass::~customClass()
{}

label customClass::basicFunction() const
{
    Info << "Calling customClass::basicFunction()" << endl;
    return myInt_*2;
}

void customClass::meshOpFunction(fvMesh& mesh)
{
    Info << "Custom class got a mesh with " << mesh.C().size() << " cells" << endl;
    myInt_ = mesh.C().size();
}
