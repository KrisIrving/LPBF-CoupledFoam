// SPDX-License-Identifier: GPL-3.0-or-later
// Demo contract: two persistent, consecutive IDs; global gathers are not
// a production distributed particle mapping implementation.
#ifndef LPBF_PARTICLE_BACKEND_H
#define LPBF_PARTICLE_BACKEND_H
#include "library.h"
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

namespace lpbf
{
using Vec = std::array<double, 3>;
struct Particle
{
    int id, owner;
    Vec x, v, omega;
    double mass, radius;
};
class ParticleBackend
{
public:
    virtual ~ParticleBackend() = default;
    virtual std::vector<Particle> state() = 0;
    virtual void forces(const std::vector<Vec>&) = 0;
    virtual void advance(int) = 0;
    virtual double time() = 0;
};

class LiggghtsBackend final : public ParticleBackend
{
    MPI_Comm comm_ = MPI_COMM_NULL;
    void* handle_ = nullptr;
    bool forceInstalled_ = false;
    void command(const std::string& s) { lammps_command(handle_, s.c_str()); }
public:
    explicit LiggghtsBackend(MPI_Comm comm, const std::string& input)
    {
        MPI_Comm_dup(comm, &comm_);
        char name[] = "embedded-liggghts";
        char logOpt[] = "-log", logFile[] = "log.dem";
        char screenOpt[] = "-screen", screenFile[] = "none";
        char* argv[] = {name, logOpt, logFile, screenOpt, screenFile};
        lammps_open(5, argv, comm_, &handle_);
        if (!handle_) throw std::runtime_error("LIGGGHTS handle is null");
        lammps_file(handle_, input.c_str());
        if (lammps_get_natoms(handle_) != 2)
            throw std::runtime_error("Demo requires exactly two particles");
    }
    ~LiggghtsBackend() override
    {
        if (handle_) lammps_close(handle_);
        if (comm_ != MPI_COMM_NULL) MPI_Comm_free(&comm_);
    }
    std::vector<Particle> state() override
    {
        if (lammps_get_natoms(handle_) != 2)
            throw std::runtime_error("Particle count changed; demo IDs must persist");
        int ids[2] = {0, 0};
        double x[6] = {}, v[6] = {}, omega[6] = {}, mass[2] = {}, radius[2] = {};
        lammps_gather_atoms(handle_, "id", 0, 1, ids);
        if (ids[0] != 1 || ids[1] != 2)
            throw std::runtime_error("Consecutive demo IDs unavailable");
        for (const char* field : {"x", "v", "omega", "rmass", "radius"})
            if (!lammps_extract_atom(handle_, field))
                throw std::runtime_error(std::string("Missing DEM field: ") + field);
        lammps_gather_atoms(handle_, "x", 1, 3, x);
        lammps_gather_atoms(handle_, "v", 1, 3, v);
        lammps_gather_atoms(handle_, "omega", 1, 3, omega);
        lammps_gather_atoms(handle_, "rmass", 1, 1, mass);
        lammps_gather_atoms(handle_, "radius", 1, 1, radius);
        auto* nlocal = static_cast<int*>(lammps_extract_global(handle_, "nlocal"));
        auto* tags = static_cast<int*>(lammps_extract_atom(handle_, "id"));
        if (!nlocal || (*nlocal && !tags))
            throw std::runtime_error("Cannot inspect local particle ownership");
        int rank; MPI_Comm_rank(comm_, &rank);
        int owners[2] = {}, counts[2] = {}, allOwners[2], allCounts[2];
        for (int j=0; j<*nlocal; ++j)
        {
            if (tags[j]<1 || tags[j]>2) throw std::runtime_error("Unexpected local ID");
            owners[tags[j]-1] = rank+1;
            counts[tags[j]-1] += 1;
        }
        MPI_Allreduce(owners, allOwners, 2, MPI_INT, MPI_SUM, comm_);
        MPI_Allreduce(counts, allCounts, 2, MPI_INT, MPI_SUM, comm_);
        std::vector<Particle> result;
        for (int i=0; i<2; ++i)
        {
            if (allCounts[i]!=1 || !std::isfinite(mass[i]) || mass[i]<=0)
                throw std::runtime_error("Invalid ownership or mass");
            result.push_back({ids[i], allOwners[i]-1,
                {x[3*i], x[3*i+1], x[3*i+2]},
                {v[3*i], v[3*i+1], v[3*i+2]},
                {omega[3*i], omega[3*i+1], omega[3*i+2]}, mass[i], radius[i]});
        }
        return result;
    }
    void forces(const std::vector<Vec>& values) override
    {
        if (values.size()!=2) throw std::runtime_error("Force count mismatch");
        for (int i=0; i<2; ++i)
        {
            const std::string fix = "hostForce" + std::to_string(i+1);
            if (forceInstalled_) command("unfix " + fix);
            std::ostringstream s; s << std::setprecision(17)
                << "fix " << fix << " particle" << i+1 << " addforce";
            for (double f : values[i])
            {
                if (!std::isfinite(f)) throw std::runtime_error("Nonfinite force");
                s << ' ' << f;
            }
            command(s.str());
        }
        forceInstalled_ = true;
    }
    void advance(int n) override
    {
        if (n<=0) throw std::runtime_error("Invalid DEM substeps");
        // Each run reinitialises force evaluation after changing host force.
        command("run " + std::to_string(n) + " pre yes post no");
    }
    double time() override
    {
        char name[] = "demoTime";
        auto* t = static_cast<double*>(lammps_extract_variable(handle_, name, nullptr));
        if (!t) throw std::runtime_error("DEM clock unavailable");
        const double result = *t; lammps_free(t); return result;
    }
};
}
#endif
