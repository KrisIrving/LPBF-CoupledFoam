// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LPBF_MECHANICAL_BACKEND_H
#define LPBF_MECHANICAL_BACKEND_H
#include <mpi.h>
#include "library.h"
#include <array>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <string>

namespace lpbfM2
{
struct Particle
{
    int id, owner;
    double x[3], v[3], omega[3], mass, radius;
};
class MechanicalBackend
{
    MPI_Comm comm_;
    void* handle_ = nullptr;
    // Apply hydrodynamic half kicks on owned particles only. Sphere inertia is
    // isotropic; DEM NVE/sphere advances the intervening contact/drift step.
    // This first backend is restricted to noncontact spherical verification.
    void kick(const std::vector<Particle>& particles,
              const std::vector<std::array<double,6>>& loads, double halfDt)
    {
        auto* n = static_cast<int*>(lammps_extract_global(handle_, "nlocal"));
        if (!n) throw std::runtime_error("Missing nlocal");
        if (!*n) return;
        auto* ids = static_cast<int*>(lammps_extract_atom(handle_, "id"));
        auto** v = static_cast<double**>(lammps_extract_atom(handle_, "v"));
        auto** w = static_cast<double**>(lammps_extract_atom(handle_, "omega"));
        if (!ids || !v || !w) throw std::runtime_error("Missing owned DEM state");
        for (int j=0; j<*n; ++j)
        {
            auto it=std::find_if(particles.begin(),particles.end(),
                [&](const Particle& p){return p.id==ids[j];});
            if (it==particles.end()) throw std::runtime_error("Unmapped owned ID");
            const auto i=static_cast<std::size_t>(it-particles.begin());
            const double inertia=0.4*it->mass*it->radius*it->radius;
            for (int k=0; k<3; ++k)
            {
                v[j][k] += halfDt*loads[i][k]/it->mass;
                w[j][k] += halfDt*loads[i][k+3]/inertia;
            }
        }
    }
public:
    explicit MechanicalBackend(const std::string& input)
    {
        MPI_Comm_dup(MPI_COMM_WORLD,&comm_);
        char a[]="embedded-dem", b[]="-log", c[]="log.dem", d[]="-screen", e[]="none";
        char* argv[]={a,b,c,d,e};
        lammps_open(5,argv,comm_,&handle_);
        if (!handle_) throw std::runtime_error("Null DEM handle");
        lammps_file(handle_,input.c_str());
    }
    ~MechanicalBackend()
    {
        if (handle_) lammps_close(handle_);
        MPI_Comm_free(&comm_);
    }
    MechanicalBackend(const MechanicalBackend&)=delete;
    MechanicalBackend& operator=(const MechanicalBackend&)=delete;
    std::vector<Particle> state()
    {
        int rank, size; MPI_Comm_rank(comm_,&rank); MPI_Comm_size(comm_,&size);
        auto* n=static_cast<int*>(lammps_extract_global(handle_,"nlocal"));
        if (!n) throw std::runtime_error("Missing DEM ownership count");
        std::vector<Particle> local;
        if (*n)
        {
            auto* id=static_cast<int*>(lammps_extract_atom(handle_,"id"));
            auto** x=static_cast<double**>(lammps_extract_atom(handle_,"x"));
            auto** v=static_cast<double**>(lammps_extract_atom(handle_,"v"));
            auto** w=static_cast<double**>(lammps_extract_atom(handle_,"omega"));
            auto* m=static_cast<double*>(lammps_extract_atom(handle_,"rmass"));
            auto* r=static_cast<double*>(lammps_extract_atom(handle_,"radius"));
            if (!id||!x||!v||!w||!m||!r) throw std::runtime_error("Missing DEM sphere fields");
            for (int j=0; j<*n; ++j)
            {
                Particle p{}; p.id=id[j]; p.owner=rank; p.mass=m[j]; p.radius=r[j];
                for (int k=0;k<3;++k) {p.x[k]=x[j][k];p.v[k]=v[j][k];p.omega[k]=w[j][k];}
                local.push_back(p);
            }
        }
        const int bytes=static_cast<int>(local.size()*sizeof(Particle));
        std::vector<int> counts(size), offsets(size);
        MPI_Allgather(&bytes,1,MPI_INT,counts.data(),1,MPI_INT,comm_);
        int total=0;
        for(int i=0;i<size;++i) {offsets[i]=total;total+=counts[i];}
        // Homogeneous ABI only, as in the current x86_64 fixture. Global
        // replication supports sparse IDs but is not a production scaling claim.
        std::vector<Particle> all(total/sizeof(Particle));
        MPI_Allgatherv(local.data(),bytes,MPI_BYTE,all.data(),counts.data(),offsets.data(),MPI_BYTE,comm_);
        std::sort(all.begin(),all.end(),[](const Particle&a,const Particle&b){return a.id<b.id;});
        if (all.empty()) throw std::runtime_error("No DEM particles");
        for(std::size_t i=0;i<all.size();++i)
        {
            if (all[i].id<=0 || (i && all[i].id==all[i-1].id)
                || !std::isfinite(all[i].mass) || all[i].mass<=0
                || !std::isfinite(all[i].radius) || all[i].radius<=0)
                throw std::runtime_error("Invalid ID, ownership or sphere data");
            for(int k=0;k<3;++k)
                if(!std::isfinite(all[i].x[k])||!std::isfinite(all[i].v[k])||!std::isfinite(all[i].omega[k]))
                    throw std::runtime_error("Nonfinite DEM motion");
        }
        return all;
    }
    double time()
    {
        char key[]="demoTime";
        auto* ptr=static_cast<double*>(lammps_extract_variable(handle_,key,nullptr));
        if(!ptr) throw std::runtime_error("Missing DEM clock");
        double t=*ptr; lammps_free(ptr); return t;
    }
    void advance(const std::vector<Particle>& p,
                 const std::vector<std::array<double,6>>& loads, int steps, double dt)
    {
        if(steps<1 || loads.size()!=p.size()) throw std::runtime_error("Invalid coupling window");
        for(const auto& f:loads) for(double v:f)
            if(!std::isfinite(v)) throw std::runtime_error("Invalid hydrodynamic load");
        for(int s=0;s<steps;++s)
        {
            kick(p,loads,0.5*dt);
            lammps_command(handle_,"run 1 pre yes post no");
            kick(p,loads,0.5*dt);
        }
    }
    void writeRestart(const std::string& path)
    {
        lammps_command(handle_,("write_restart "+path).c_str());
    }
};
}
#endif
