//
// Created by james on 28/04/15.
//

#include "membrane.h"

#include "small_functions.h"

using std::abs;
using std::array;
using std::map;
using std::set;
using std::string;
using std::vector;

Membrane::Membrane(const vector<Residue> &residues, const Frame &frame,
                   const int resolution, const int blocks, const bool header)
    : residues_(residues), header_(header)
{
    setResolution(resolution);
    sortBilayer(frame, blocks);
    prepCSVAreaPerLipid();
    prepCSVAvgThickness();
}

Membrane::~Membrane()
{
    fclose(aplFile_);
    fclose(avgFile_);
};

void Membrane::sortBilayer(const Frame &frame, const int blocks)
{
    // Reset running values
    numLipids_ = 0;
    upperHeads_.clear();
    lowerHeads_.clear();
    protAtoms_.clear();
    lowerNumRes_.clear();
    upperNumRes_.clear();

    // Copy box from Frame - assume orthorhombic
    box_[0] = frame.box_[0][0];
    box_[1] = frame.box_[1][1];
    box_[2] = frame.box_[2][2];

    // Find middle of membrane in z coord
    // Do this in blocks to account for curvature - more curvature needs more
    // blocks
    LightArray<double> block_avg_z(blocks, blocks);
    LightArray<double> block_tot_residues(blocks, blocks);

    for (const Residue &res : residues_)
    {
        if (res.ref_atom < 0)
            continue;
        numLipids_ += res.num_residues;
        for (int i = 0; i < res.num_residues; i++)
        {
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x   = wrap(static_cast<int>(frame.atoms_[num].coords[0] *
                                                blocks / box_[0]),
                                 0, blocks);
            const int y   = wrap(static_cast<int>(frame.atoms_[num].coords[1] *
                                                blocks / box_[1]),
                                 0, blocks);
            block_avg_z(x, y) += frame.atoms_[num].coords[2];
            block_tot_residues(x, y)++;
        }
    }

    block_avg_z /= block_tot_residues;

    double minz = block_avg_z.mean();
    double maxz = minz;

    // Separate membrane into upper and lower
    for (const Residue &res : residues_)
    {
        if (res.ref_atom < 0)
            continue;
        int num_in_leaflet[2] = {0, 0};
        for (int i = 0; i < res.num_residues; i++)
        {
            const int num = res.ref_atom + i * res.num_atoms + res.start;
            const int x =
                int(frame.atoms_[num].coords[0] * blocks / box_[0]) % blocks;
            const int y =
                int(frame.atoms_[num].coords[1] * blocks / box_[1]) % blocks;
            const double z = frame.atoms_[num].coords[2];

            minz = std::min(minz, z);
            maxz = std::max(maxz, z);

            if (z < block_avg_z(x, y))
            {
                lowerHeads_.push_back(num);
                num_in_leaflet[0]++;
                lowerNumRes_[res.resname]++;
            }
            else
            {
                upperHeads_.push_back(num);
                num_in_leaflet[1]++;
                upperNumRes_[res.resname]++;
            }
        }
        printf("%5s: %'4d lower, %'4d upper\n", res.resname.c_str(),
               num_in_leaflet[0], num_in_leaflet[1]);
    }

    for (const Residue &res : residues_)
    {
        if (res.resname != "PROT")
            continue;
        if (res.ref_atom_name == "ALL")
        {
            for (int i = res.start; i < res.end; i++)
            {
                const double z = frame.atoms_[i].coords[2];
                if (minz < z && z < maxz)
                    protAtoms_.push_back(i);
            }
            protein_ = true;
        }
    }
}

double Membrane::thickness(const Frame &frame, const bool with_reset)
{
    if (with_reset)
        reset();

    // Copy box from Frame - assume orthorhombic
    box_[0]  = frame.box_[0][0];
    box_[1]  = frame.box_[1][1];
    box_[2]  = frame.box_[2][2];
    step_[0] = box_[0] / grid_;
    step_[1] = box_[1] / grid_;

    double avg_thickness = 0;
    //    sortBilayer(frame, 4);

#pragma omp parallel sections reduction(+ : avg_thickness) default(shared)
    {
#pragma omp section
        {
            makePairs(frame, upperHeads_, lowerHeads_, upperPair_);
            avg_thickness += closestLipid(frame, upperHeads_, upperPair_,
                                          upperResPPL_, closestUpper_);
        }
#pragma omp section
        {
            makePairs(frame, lowerHeads_, upperHeads_, lowerPair_);
            avg_thickness += closestLipid(frame, lowerHeads_, lowerPair_,
                                          lowerResPPL_, closestLower_);
        }
    }

    avg_thickness /= 2;
    fprintf(avgFile_, "%8.3f%8.3f\n", frame.time_, avg_thickness);

    numFrames_++;
    return avg_thickness;
}

void Membrane::makePairs(const Frame &frame, const vector<int> &ref,
                         const vector<int> &other, map<int, double> &pairs)
{
    // For each reference particle in the ref leaflet
    for (const int i : ref)
    {
        double min_dist_2    = box_[0] * box_[1];
        array<double, 3> r_i = frame.atoms_[i].coords;

        // Find the closest reference particle in the other leaflet
        array<double, 3> r_j = {{0., 0., 0.}};
        for (const int j : other)
        {
            r_j[0] = frame.atoms_[j].coords[0];
            r_j[1] = frame.atoms_[j].coords[1];

            const double dist_2 = distSqrPlane(r_i, r_j, frame.boxDiag_);
            if (dist_2 < min_dist_2)
            {
                min_dist_2 = dist_2;
                r_j[2]     = frame.atoms_[j].coords[2];
            }
        }

        pairs[i] = abs(r_i[2] - r_j[2]);
    }
}

// TODO optimise this more - it takes ~80% of runtime
double Membrane::closestLipid(const Frame &frame, const vector<int> &ref,
                              const map<int, double> &pairs,
                              map<string, int> &resPPL,
                              LightArray<int> &closest)
{
    const double box_diag2 = box_[0] * box_[1];

    const size_t ref_len = ref.size();
    vector<array<double, 3>> ref_cache(ref.size());
    vector<int> ref_lookup(ref.size());
    {
        int it = 0;
        for (const int r : ref)
        {
            ref_cache[it]  = frame.atoms_[r].coords;
            ref_lookup[it] = r;
            it++;
        }
    }

    const size_t prot_len = protAtoms_.size();
    vector<array<double, 3>> prot_cache(prot_len);
    if (protein_)
    {
        int it = 0;
        for (const int r : protAtoms_)
        {
            prot_cache[it] = frame.atoms_[r].coords;
            it++;
        }
    }

    double sum = 0;
    int n_vals = 0;
    //make a copy of frame.boxDiag_ to avoid having to share the whole class object between threads.
    std::array<double, 3> boxDiag = frame.boxDiag_;

#pragma omp parallel for default(none)                                         \
    shared(boxDiag, ref, pairs, closest, ref_cache, ref_lookup, prot_cache,      \
               resPPL, box_diag2, ref_len, prot_len)                           \
    reduction(+ : sum, n_vals)
    for (int i = 0; i < grid_; i++)
    {
        array<double, 3> grid_coords;
        grid_coords[0] = (i + 0.5) * step_[0];

        for (int j = 0; j < grid_; j++)
        {
            grid_coords[1]   = (j + 0.5) * step_[1];
            double min_dist2 = box_diag2;

            // Find closest lipid in reference leaflet
            int closest_int = -1;
            for (int k = 0; k < ref_len; k++)
            {
                const double dist2 =
                    distSqrPlane(grid_coords, ref_cache[k], boxDiag);
                if (dist2 < min_dist2)
                {
                    closest_int = k;
                    min_dist2   = dist2;
                }
            }

            bool is_protein = false;
            if (protein_)
            {
                for (int k = 0; k < prot_len; k++)
                {
                    const double dist2 = distSqrPlane(
                        grid_coords, prot_cache[k], boxDiag);
                    if (dist2 < min_dist2)
                    {
                        is_protein = true;
                        break;
                    }
                }
            }

            // if (is_protein)
            // {
            //     // #pragma omp atomic update
            //     //                 resPPL["PROT"]++;
            // }
            // else
            if (!is_protein){
                const int close_ref = ref_lookup[closest_int];
                closest(i, j)       = close_ref;
                const double tmp    = pairs.at(close_ref);
                for (const Residue &res : residues_)
                {
                    if (close_ref >= res.start && close_ref < res.end)
                    {
#pragma omp atomic update
                        resPPL[res.resname]++;
                    }
                }

                n_vals++;
                sum += tmp;
#pragma omp atomic update
                thickness_(i, j) += tmp;
            }
        }
    }

    return sum / n_vals;
}

void Membrane::curvature(const Frame &frame)
{
    LightArray<double> avg_z(grid_, grid_);

    LightArray<double> respect_to_x(grid_, grid_);
    LightArray<double> respect_to_y(grid_, grid_);

    const double inv_h2_x = 1. / (step_[0] * step_[0]);
    const double inv_h2_y = 1. / (step_[1] * step_[1]);

    // double curv_x_avg = 0.;
    // double curv_y_avg = 0.;

    // Calculate average z coord on grid
    for (int i = 0; i < grid_; i++)
    {
        for (int j = 0; j < grid_; j++)
        {
            avg_z(i, j) = (frame.atoms_[closestUpper_.at(i, j)].coords[2] +
                           frame.atoms_[closestLower_.at(i, j)].coords[2]) /
                          2.;
        }
    }

    // Do finite differences wrt x and y
    for (int i = 1; i < grid_ - 1; i++)
    {
        for (int j = 1; j < grid_ - 1; j++)
        {
            respect_to_x(i, j) = inv_h2_x * (avg_z(i + 1, j) + avg_z(i - 1, j) -
                                             2 * avg_z(i, j));
            respect_to_y(i, j) = inv_h2_y * (avg_z(i, j + 1) + avg_z(i, j - 1) -
                                             2 * avg_z(i, j));

            // curv_x_avg += respect_to_x(i, j);
            // curv_y_avg += respect_to_y(i, j);
        }
    }

    // curv_x_avg /= grid_ * grid_;
    // curv_y_avg /= grid_ * grid_;

    for (int i = 0; i < grid_; i++)
    {
        for (int j = 0; j < grid_; j++)
        {
            curvMean_(i, j) =
                (respect_to_x.at(i, j) + respect_to_y.at(i, j)) / 2.;
            curvGaussian_(i, j) = respect_to_x.at(i, j) * respect_to_y.at(i, j);
        }
    }

    curvMean_.smooth(5);
    curvGaussian_.smooth(5);
}

void Membrane::printCSVCurvature(const std::string &filename) const
{
    const string file = filename + ".dat";
    // Backup using small_functions.h
    backup_old_file(file);
    FILE *f = fopen(file.c_str(), "w");
    if (f == nullptr)
        throw std::runtime_error("Could not open output file.");

    if (header_)
    {
        fprintf(f, "@legend Membrane mean curvature\n");
        fprintf(f, "@xlabel X (nm)\n");
        fprintf(f, "@ylabel Y (nm)\n");
        fprintf(f, "@xwidth %f\n", box_[0]);
        fprintf(f, "@ywidth %f\n", box_[1]);
    }

    //    printf("Printing\n");
    for (int i = 0; i < grid_; i++)
    {
        for (int j = 0; j < grid_; j++)
        {
            fprintf(f, "%8.3f", curvMean_.at(i, j));
            //            printf("%8.3f", curvMean_.at(i, j));
        }
        fprintf(f, "\n");
        //        printf("\n");
    }

    fclose(f);
}

void Membrane::prepCSVAreaPerLipid()
{
    const string file = "APL.dat";
    // Backup using small_functions.h
    backup_old_file(file);
    aplFile_ = fopen(file.c_str(), "w");
    if (!aplFile_)
        throw std::runtime_error("Could not open output file");

    if (header_)
    {
        fprintf(aplFile_, "@legend Area Per Lipid\n");
        fprintf(aplFile_, "@xlabel Time (ps)\n");
        fprintf(aplFile_, "@ylabel APL (nm^2)\n");
    }

    for (const auto &ppl : upperResPPL_)
        fprintf(aplFile_, "%5s", ppl.first.c_str());
    for (const auto &ppl : lowerResPPL_)
        fprintf(aplFile_, "%5s", ppl.first.c_str());
    fprintf(aplFile_, "\n");
}

void Membrane::printCSVAreaPerLipid(const float time) const
{
    fprintf(aplFile_, "%12.3f", time);
    for (const auto &ppl : upperResPPL_)
    {
        const int num     = ppl.second;
        const double area = num * step_[0] * step_[1];
        const double APL  = area / (numFrames_ * upperNumRes_.at(ppl.first));
        fprintf(aplFile_, "%12.3f", APL);
    }
    for (const auto &ppl : lowerResPPL_)
    {
        const int num     = ppl.second;
        const double area = num * step_[0] * step_[1];
        const double APL  = area / (numFrames_ * lowerNumRes_.at(ppl.first));
        fprintf(aplFile_, "%12.3f", APL);
    }
    fprintf(aplFile_, "\n");
}

void Membrane::prepCSVAvgThickness()
{
    const string file = "avg_thickness.dat";
    // Backup using small_functions.h
    backup_old_file(file);
    avgFile_ = fopen(file.c_str(), "w");
    if (!avgFile_)
        throw std::runtime_error("Could not open output file");

    if (header_)
    {
        fprintf(aplFile_, "@legend Average Thickness\n");
        fprintf(aplFile_, "@xlabel Time (ps)\n");
        fprintf(aplFile_, "@ylabel Thickness (nm)\n");
    }
}

double Membrane::mean() const
{
    return thickness_.mean();
}

void Membrane::normalize(const int smooth_iter)
{
    // Apply iterations of Gauss-Seidel smoother
    thickness_.smooth(smooth_iter);
    thickness_ /= 2 * numFrames_;
}

void Membrane::printCSV(const std::string &filename) const
{
    if (header_)
    {
        const string file = filename + ".dat";
        // Backup using small_functions.h
        backup_old_file(file);
        FILE *f = fopen(file.c_str(), "w");

        fprintf(f, "@legend Membrane thickness\n");
        fprintf(f, "@xlabel X (nm)\n");
        fprintf(f, "@ylabel Y (nm)\n");
        fprintf(f, "@xwidth %f\n", box_[0]);
        fprintf(f, "@ywidth %f\n", box_[1]);
        fclose(f);

        // Print CSV - true suppresses backup - file has been opened already
        thickness_.printCSV(filename, true);
    }
    else
    {
        thickness_.printCSV(filename);
    }
}

void Membrane::setResolution(const int n)
{
    grid_ = n;
    thickness_.alloc(grid_, grid_);

    closestUpper_.alloc(n, n);
    closestLower_.alloc(n, n);
    curvMean_.alloc(grid_, grid_);
    curvGaussian_.alloc(grid_, grid_);
}

void Membrane::reset()
{
    thickness_.zero();
    for (auto &ppl : upperResPPL_)
        ppl.second = 0;
    for (auto &ppl : lowerResPPL_)
        ppl.second = 0;
    numFrames_ = 0;
}
