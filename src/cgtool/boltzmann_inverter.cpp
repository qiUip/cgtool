#include "boltzmann_inverter.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

BoltzmannInverter::BoltzmannInverter(const double temp, const int bins)
    : temp_(temp), bins_(bins)
{
    histogram_.init(bins_);
    gaussian_.alloc(bins_);
    harmonic_.alloc(bins_);
}

void BoltzmannInverter::calculate(BondStruct &bond)
{
    histogram_.zero();
    gaussian_.zero();
    harmonic_.zero();
    n_ = bond.values_.size();
    statisticalMoments(bond.values_);
    bond.avg_ = mean_;
    binHistogram(bond.values_);
    bond.rsqr_          = gaussianRSquared();
    type_               = bond.type_;
    bond.forceConstant_ = invertGaussianSimple();
}

double BoltzmannInverter::invertGaussianSimple()
{
    const double RT = 8.314 * temp_ / 1000.;

    switch (type_)
    {
        case BondType::LENGTH:
            return RT / (sdev_ * sdev_);
        case BondType::ANGLE:
        {
            const double sinmean = sin(mean_ * M_PI / 180.);
            const double sdevrad = sdev_ * M_PI / 180.;
            return RT / (sinmean * sinmean * sdevrad * sdevrad);
        }
        case BondType::DIHEDRAL:
        {
            // Assumes multiplicity 1 - CG tends to be
            // TODO try FFT to account for other multiplicities
            const double sdevrad = sdev_ * M_PI / 180.;
            return RT / (sdevrad * sdevrad);
        }
    }

    // This code can never be reached but GCC still complains without it
    return -1.;
}

void BoltzmannInverter::binHistogram(const vector<double> &vec)
{
    step_    = (max_ - min_) / (bins_ - 1);
    meanBin_ = static_cast<int>((mean_ - min_) / step_);

    int loc = 0;
    for (const double val : vec)
    {
        loc = static_cast<int>((val - min_) / step_);
        histogram_.increment(loc);
    }
}

double BoltzmannInverter::gaussianRSquared()
{
    const double prefactor  = 1. / (sdev_ * sqrt(2. * M_PI));
    const double postfactor = 1. / (2. * var_);

    double y_bar = 0.;
    // First pass to calculate mean and gaussian integral
    for (int i = 0; i < bins_; i++)
    {
        const double x = min_ + (i + 0.5) * step_;
        const double gau =
            prefactor * exp(-(x - mean_) * (x - mean_) * postfactor);
        gaussian_(i) = gau;
        y_bar += histogram_.at(i);
    }

    y_bar /= bins_;
    integral_ = n_ / gaussian_.sum();

    double ss_res = 0., ss_tot = 0.;
    // double ss_reg = 0.;
    // double sse = 0.;
    // Second pass to calculate R^2
    for (int i = 0; i < bins_; i++)
    {
        int actual = histogram_.at(i);
        gaussian_(i) *= integral_;
        const double gau = gaussian_.at(i);
        // ss_reg += (gau - y_bar) * (gau - y_bar);
        ss_res += (gau - actual) * (gau - actual);
        ss_tot += (actual - y_bar) * (actual - y_bar);
        // sse += (actual - gau) * (actual - gau);
    }
    const double r_sqr = 1 - ss_res / ss_tot;
    return r_sqr;
}

double BoltzmannInverter::statisticalMoments(const vector<double> &vec)
{
    double sum = 0.;
    if (n_ == 0)
        n_ = vec.size();
    // Calculate mean with first pass
    for (const double val : vec)
        sum += val;
    mean_ = sum / n_;
    max_  = mean_;
    min_  = mean_;

    // double ep = 0.;
    var_      = 0.;
    adev_     = 0.;
    // Calculate deviations with second pass
    for (const double val : vec)
    {
        const double dev = val - mean_;
        // ep += dev;
        adev_ += fabs(dev);
        var_ += dev * dev;
        min_ = std::min(min_, val);
        max_ = std::max(max_, val);
    }

    adev_ /= n_;
    var_  = var_ / (n_ - 1);
    sdev_ = sqrt(var_);
    return mean_;
}
