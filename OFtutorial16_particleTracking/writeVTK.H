// writing vtk file for the particles.
// for VTK file format, plz refer:
// https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

Info << nl << "writing particle's path as VTK file...";

// making a separate directory for storing VTK file
fileName VTK_dir = mesh.time().path()/"VTK";
mkDir(VTK_dir);

// creating file pointer
autoPtr<OFstream> vtkFilePtr;
vtkFilePtr.reset(new OFstream(VTK_dir/"particle_path.vtk"));

// writing vtk header
vtkFilePtr() << "# vtk DataFile Version 2.0" << endl << "particle_path" << endl << "ASCII" << endl << "DATASET POLYDATA" << endl;
vtkFilePtr() << nl;

// writing points header
vtkFilePtr() << "POINTS " << particlePositions.size() << " DOUBLE" << endl;

// writing point coordinates
forAll(particlePositions,pt)
{
  vtkFilePtr() << particlePositions[pt].x() << " " << particlePositions[pt].y() << " " << particlePositions[pt].z() << endl;
}
vtkFilePtr() << nl;

// writing lines header and list
vtkFilePtr() << "LINES 1 " << particlePositions.size() + 1 << endl; // it is LINES <no of lines> <size of list>
vtkFilePtr() << particlePositions.size() << endl;           // <no of points> <point list .....>
forAll(particlePositions,pt)
{
  vtkFilePtr() << pt << endl;
}

Info << tab << " done." << endl;
