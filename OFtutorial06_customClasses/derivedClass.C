#include "derivedClass.H"

myDict::myDict(const IOobject& ioObj)
:
    // Call the base class constructor to make sure all inherited members get
    // initialised as required
    IOdictionary(ioObj)
{
// do nothing since we do not have any bespoke fields
}

myDict::~myDict()
{}

void myDict::printTokensInTheDict() const
{
    // retrieve the list of non-space characters in the file using the
    // method defined in dictionary.H, from which the IOdictionary object itself
    // is derived.
    List<token> characters(this->tokens());

    // Create a stream which will hold the message to be printed out.
    // Important to remember about the namespace.
    std::stringstream ss;
    ss << "Tokens in the file:";

    // go over each token in the file
    forAll(characters,i)
        // if the entry is a word, add it to the message
        if (characters[i].isWord())
            ss << "\n" << tab << characters[i].wordToken();

    // print the message - convert to a C-style char array to make sure the
    // printout looks good
    Info << ss.str().c_str() << endl;
}
