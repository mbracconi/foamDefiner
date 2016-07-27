Info<< "Reading foamDefinerDict options\n" << endl;
IOdictionary foamOptions
(
	IOobject
	(
		"foamDefinerDict",
		runTime.constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

// Surface patch
word  surfacePatch_(foamOptions.lookup("foamSurfacePatch"));

// Switch -> off solid fraction | on void fraction
Switch  fractionType_(foamOptions.lookup("fractionType"));

