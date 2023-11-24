#ifndef CGMAP_H_
#define CGMAP_H_

#include <map>
#include <string>
#include <vector>

#include "frame.h"
#include "residue.h"

/**
 * \brief Struct to hold the mapping for a single CG bead
 */
struct BeadMap
{
    /** The name of this CG bead */
    std::string name;
    /** The number of the CG bead */
    int num;
    /** The particle type of the bead - e.g. MARTINI P4, C2, etc. */
    std::string type;
    /** The number of atoms this bead contains */
    int num_atoms;
    /** The atoms which should be mapped into this bead */
    std::vector<std::string> atoms;
    /** The atoms which should be mapped into this bead, by order in Frame */
    std::vector<int> atom_nums;
    /** Total mass of bead */
    double mass;
    /** Total charge on bead */
    double charge;
    /** Lennard-Jones C6 parameter */
    double c06 = 0.;
    /** Lennard-Jones C12 parameter */
    double c12 = 0.;
};

enum class MapType
{
    CM,
    GC,
    ATOM
};

/**
 * \brief Contains data and functions related to the CG mapping
 *
 * Has functions to read in a CG mapping from file and apply it to an atomistic
 * Frame Mostly just a wrapper around a BeadMap vector
 */
class CGMap
{
private:
    /** \brief What type of mapping are we going to apply?  CM, GC, or atom
     * centred Default is geometric centre of component atoms. */
    MapType mapType_ = MapType::GC;

    const std::vector<Residue> &aaRes_;
    std::vector<Residue> &cgRes_;

    /** \brief Correct LJ parameters for CG */
    double calcLJ(const std::vector<int> &ljs);

    /** Number of beads defined */
    int numBeads_;

    /** Vector of BeadMap; holds the mappings for every bead */
    std::vector<BeadMap> mapping_;

public:
    /**
     * \brief Constructor to create an instance from the mapping file provided
     */
    CGMap(const std::vector<Residue> &aa_res, std::vector<Residue> &cg_res,
          const std::string &filename = "")
        : aaRes_(aa_res), cgRes_(cg_res)
    {
        if (filename != "")
            fromFile(filename);
    };

    /**
     * \brief Read in CG mapping from file
     *
     * \throws std::runtime_error if file cannot be opened
     */
    void fromFile(const std::string &filename);

    /**
     * \brief Setup a CG Frame object that has already been declared
     *
     * Allocates space for each bead and copies over constant data from the
     * atomistic Frame
     */
    void initFrame(const Frame &aa_frame, Frame &cg_frame);

    /**
     * \brief Apply CG mapping to an atomistic Frame
     *
     * \throws std::runtime_error if Frame hasn't been setup.
     * Requires that initFrame has already been called to setup the CG Frame.
     */
    bool apply(const Frame &aa_frame, Frame &cg_frame);

    /** \brief Calculate dipoles from atomistic frame */
    void calcDipoles(const Frame &aa_frame, Frame &cg_frame);

    std::vector<BeadMap> getMapping() const
    {
        return mapping_;
    }
};

#endif
