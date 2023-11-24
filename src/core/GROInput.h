//
// Created by james on 25/08/15.
//

#ifndef CGTOOL_GROINPUT_H
#define CGTOOL_GROINPUT_H

#include "TrjInput.h"

#include <fstream>
#include <set>

class GROInput : public TrjInput
{
private:
    /** \brief Input file stream. */
    std::ifstream file_;

    /** \brief Set of amino acid residue names. */
    const std::set<std::string> aminoAcids_ = {
        "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN",
        "GLX", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
        "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

    /** \brief Open and prepare input file. */
    int openFile(const std::string &filename);
    /** \brief Close input file. */
    int closeFile();

public:
    /** \brief Constructor.  Calls openFile(). */
    GROInput(const std::string &filename);
    /** \brief Destructor.  Calls closeFile(). */
    ~GROInput();

    /** \brief Read a Frame from input file. */
    int readFrame(Frame &frame);

    void readResidues(std::vector<Residue> &residues);

    friend class Frame;
};

#endif // CGTOOL_GROINPUT_H
