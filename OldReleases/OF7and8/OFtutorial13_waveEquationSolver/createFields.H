Info << nl << "Reading transportProperties\n" << endl;

// reading transportProperties dictionary file from constant/ directory.
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// reading the wave speed "C" value from transportProperties dictionary
dimensionedScalar C
(
    "C",
    dimensionSet(0,1,-1,0,0),
    transportProperties.lookup("C")
);

Info << nl << "Reading field h" << endl;

// reading the amplitude field "h" values.
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
