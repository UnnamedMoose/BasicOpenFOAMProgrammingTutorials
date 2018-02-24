#include "fvCFD.H"

// This method simply implements the calculation of the distance of each
// cell centre from x0; it accepts the scalarField r as a reference to avoid
// passing large amounts of information by value, as this is expensive.
// It returns the maximum value found.
scalar computeR(const fvMesh& mesh, scalarField& r, dimensionedVector x0);

// This computes the velocity field. The reference to the pressure is obtained
// through the mesh object, using the name of the p field only. This assmes a
// default value which may be redefined, if necessary.
void computeU(const fvMesh& mesh, volVectorField& U, word pName = "p");
