#ifndef BONDSTRUCT_H_
#define BONDSTRUCT_H_

#include <string>
#include <vector>

#include "frame.h"

enum class BondType
{
    LENGTH   = 2,
    ANGLE    = 3,
    DIHEDRAL = 4
};
enum class FunctionalForm
{
    HARMONIC,
    COS,
    COSHARMONIC
};

/**
 * \brief Class to hold atoms in bonds, angles and dihedrals.
 */
class BondStruct
{
private:
    /** \brief What type of bond measure is it?  Length, angle or dihedral */
    const BondType type_;
    /** \brief The values of the bond parameter (length, angle, dih) for each
     * Frame */
    std::vector<double> values_;
    /** \brief Vector of atom numbers for this bond property; For a bond length
     * will contain two names; three for angle; four for dihedral */
    std::vector<int> atomNums_;
    /** \brief R^2 of fitting gaussian to bond distribution
     * A low value indicates that the bond is probably bimodal */
    double rsqr_ = 0.;
    /** \brief Force constant */
    double forceConstant_ = 0.;
    /** \brief Average of bond parameter.  Double avoids overflow on summing
     * over large numbers of Frames */
    double avg_ = 0.;

public:
    /** Constructor to set size (bond/angle/dihedral) */
    BondStruct(const BondType type);
    /** Copy constructor - required to push into vectors */
    BondStruct(const BondStruct &other)
        : type_(other.type_), atomNums_(other.atomNums_){};

    /**
     * \brief Calculate distance between two atoms in a BondStruct object
     * Wrapper around float bondLength(int, int)
     */
    double bondLength(const Frame &frame, const int offset) const;

    /**
     * \brief Calculate angle or dihedral between atoms in a BondStruct object
     * Wrapper around float bondAngle(int, int, int, int)
     */
    double bondAngle(const Frame &frame, const int offset) const;

    double bondDihedral(const Frame &frame, const int offset) const;

    BondType getType()
    {
        return type_;
    }

    std::vector<double> &getValues()
    {
        return values_;
    }

    double getValue(int val) const
    {
        return values_[val];
    }

    std::vector<int> &getAtomNums()
    {
        return atomNums_;
    }

    int getAtomNum(int val) const
    {
        return atomNums_[val];
    }

    void setAtomNum(int i, int val)
    {
        atomNums_[i] = val;
    }

    double getRsqr() const
    {
        return rsqr_;
    }

    void setRsqr(double rsqr)
    {
        rsqr_ = rsqr;
    }

    double getForceConstant() const
    {
        return forceConstant_;
    }

    void setForceConstant(double forceConstant)
    {
        forceConstant_ = forceConstant;
    }

    double getAvg() const
    {
        return avg_;
    }

    void setAvg(double avg)
    {
        avg_ = avg;
    }
};

#endif
