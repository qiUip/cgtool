#ifndef FRAME_H_
#define FRAME_H_

#include "xdrfile.h"

#include <vector>
#include <string>
#include <map>

#include <string.h>

#include "bondset.h"

/**
* \brief Struct to hold atom data
*/
struct Atom{
    /** A serial number. no longer needed */
    int atom_num;
    /** Residue name in the itp file */
    std::string resname;
    /** Residue number */
    int resnum;
    /** Atomtype as a string.  I don't want to be dealing with *char */
    std::string atom_type;
    /** Atomic coordinates in x, y, z */
    double coords[3];
    /** Atomic charge from the force field */
    float charge = 0.f;
    /** Atomic mass */
    float mass = 0.f;
    /** \brief Create an Atom, set its number and zero its coordinates */
    Atom(int num){atom_num = num; coords[0] = 0.f; coords[1] = 0.f; coords[2] = 0.f;};
    /** Create a blank Atom instance */
    Atom(){coords[0] = 0.f; coords[1] = 0.f; coords[2] = 0.f;};
};


enum class BoxType{CUBIC, TRICLINIC};

/**
* \brief Class to hold a single frame of an XTC file
*
* Holds a std::vector<Atom> and contains member functions to operate on this
*/
class Frame{
protected:
    /** Has the XTC output been setup yet? */
    bool outputSetup_ = false;
    /** GROMACS xtc file to import frames */
    XDRFILE *xtcInput_ = nullptr;
    /** GROMACS xtc file to export frames */
    XDRFILE *xtcOutput_ = nullptr;
    /** XTC precision; not used internally, just for XTC input/output */
    float prec_ = 0.f;
    /** Holds atomic coordinates for GROMACS */
    rvec *x_ = NULL;
    /** Name of the Frame; taken from comment in the GRO file */
    std::string name_;
    /** Size of the simulation box */
    float box_[3][3];
    /** What box shape do we have?  Currently must be cubic */
    BoxType boxType_ = BoxType::CUBIC;
    /** What is the resname of the molecule we want to map - column 4 of the itp */
    std::string mapResname_;

public:
    /** Has the Frame been properly setup yet? */
    bool isSetup_ = false;
    /** The number of atoms stored in this frame that we find interesting */
    int numAtomsTrack_ = 0;
    /** Vector of Atoms; Each Atom contains position and type data */
    std::vector<Atom> atoms_;
    /** The number of atoms stored in this frame */
    int numAtoms_ = 0;
    /** Is frame invalid for some reason - molecule lies on pbc */
    bool invalid_ = false;
    /** The simulation time of this frame, in picoseconds */
    float time_ = 0.f;
    /** The number of this Frame, starting at 0 */
    int num_ = 0;
    /** The simulation step corresponding to this frame */
    int step_ = 0;
    /** Map mapping atom names to numbers for each residue */
    std::map<std::string, int> nameToNum_;
//    /** Vector of maps mapping atom names to numbers for each residue */
//    std::vector<std::map<std::string, int>> nameToNum_;
    /** How many atoms are in this residue? */
    int numAtomsPerResidue_;
    /** How many of this residue are there? */
    int numResidues_;


    /** \brief Create Frame passing frame number, number of atoms to store and the frame name
    * If we don't know the number of atoms at creation
    * this can be set later using Frame::allocateAtoms() */
    Frame(const int num, const int natoms, const std::string name);

    /** \brief Create Frame passing config files.
    * Replaces calls to the function Frame::setupFrame() */
    Frame(const std::string topname, const std::string xtcname);

    /** \brief Create Frame by copying data from another Frame
    * Intended for creating a CG Frame from an atomistic one.  Atoms are not copied. */
    Frame(const Frame &frame);

    /** \brief Destructor to free memory allocated by GROMACS functions */
    ~Frame();

    /**
    * \brief Create Frame, allocate atoms and read in data from start of XTC file
    * \throws logic_error if Frame has already been setup
    *
    * Uses libxdrfile to get number of atoms and allocate storage.
    * This function uses this data to create a Frame object to process this data.
    */
    bool setupFrame(const std::string &topname, const std::string &xtcname);

    /**
    * \brief Read next frame from the open XTC file
    */
    bool readNext();

    /**
    * \brief Prepare to write XTC output.
    * \throws std::runtime_error if memory cannot be allocated for atom array
    * \throws std::runtime_error if output TOP file cannot be opened
    * Allocate necessary atom array and create TOP file.
    */
    void setupOutput(const std::string &xtcnameout, const std::string &topnameout);

    /** \brief Write Frame to XTC output file */
    bool writeToXtc();

    /** \brief Recentre simulation box on an atom
    * Avoids problems where a residue is split by the periodic boundary,
    * causing bond lengths to be calculated incorrectly */
    void recentreBox(const int atom_num);

    /** Print info for all atoms up to n.  Default print all. */
    void printAtoms(int natoms=-1);

    /** Print all atoms up to n to GRO file.  Default print all. */
    void printGRO(const std::string &filename, int natoms=-1);

    /** \brief Calculate distance between two atoms */
    double bondLength(const int a, const int b);

    /**
    * \brief Calculate distance between two atoms in a BondStruct object
    * Wrapper around float bondLength(int, int)
    */
    double bondLength(BondStruct &bond);

    /**
    * \brief Calculate angle between vectors a->b and c->d
    * To be used for bond angles (b=c) and dihedrals (b=/=c)
    */
    double bondAngle(const int a, const int b, const int c, const int d);


    /**
    * \brief Calculate angle or dihedral between atoms in a BondStruct object
    * Wrapper around float bondAngle(int, int, int, int)
    */
    double bondAngle(BondStruct &bond);
};

#endif
