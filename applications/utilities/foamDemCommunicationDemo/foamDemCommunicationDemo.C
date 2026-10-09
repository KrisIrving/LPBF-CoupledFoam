// SPDX-License-Identifier: GPL-3.0-or-later
#include "argList.H"
#include "Time.H"
#include "Pstream.H"
#include "particleBackend.H"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace Foam;
int main(int argc, char* argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    int initialized=0; MPI_Initialized(&initialized);
    const bool ownMPI = !initialized;
    if (ownMPI) MPI_Init(&argc, &argv);
    int rank=0, ranks=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        if (ranks != Pstream::nProcs() || (ranks!=1 && ranks!=2))
            throw std::runtime_error("Use serial or mpirun -np 2 with -parallel");
        if (mag(runTime.deltaTValue()-0.01)>1e-14
            || mag(runTime.value())>1e-14 || mag(runTime.endTime().value()-0.2)>1e-14)
            throw std::runtime_error("Demo requires start=0, dt=0.01, end=0.2");
        {
            lpbf::LiggghtsBackend backend(MPI_COMM_WORLD, "input.dem");
            auto particles = backend.state();
            std::ofstream history;
            if (rank==0)
            {
                history.open("communication-history.csv");
                if (!history) throw std::runtime_error("Cannot write history");
                history << std::setprecision(17)
                    << "window,id,owner,time_of,time_dem,x,y,z,vx,vy,vz,fx,fy,fz,mass,radius,x_error,v_error,impulse_error,clock_error\n";
            }
            const double expectedMass = 2500.0 * 4.0/3.0 * std::acos(-1.0)*std::pow(0.005,3);
            for (int i=0; i<2; ++i)
                if (std::abs(particles[i].mass/expectedMass-1)>1e-12
                    || std::abs(particles[i].radius-0.005)>1e-14
                    || std::abs(particles[i].x[0]-(i==0 ? 0.45 : 0.55))>1e-12
                    || std::abs(particles[i].v[0]-(i==0 ? 0.4 : -0.4))>1e-12)
                    throw std::runtime_error("Initial particle fixture differs");
            double maxX=0, maxV=0, maxImpulse=0, maxClock=0;
            int migrations=0;
            for (int window=1; window<=20; ++window)
            {
                const auto before = particles;
                std::vector<lpbf::Vec> forces(2);
                for (int i=0; i<2; ++i)
                {
                    const double target = i==0 ? 0.6 : -0.6;
                    forces[i] = {before[i].mass*2.0*(target-before[i].v[0]), 0, 0};
                }
                backend.forces(forces);
                backend.advance(10); // DEM dt=0.001; window dt=0.01.
                ++runTime;
                particles = backend.state();
                const double demTime = backend.time();
                const double clock = std::abs(demTime-runTime.value());
                if (!std::isfinite(clock)) throw std::runtime_error("Invalid clock");
                maxClock = std::max(maxClock, clock);
                for (int i=0; i<2; ++i)
                {
                    double xe=0, ve=0, ie=0;
                    if (particles[i].owner!=before[i].owner) ++migrations;
                    for (int k=0; k<3; ++k)
                    {
                        const double a=forces[i][k]/before[i].mass, dt=0.01;
                        const double ex=before[i].x[k]+before[i].v[k]*dt+0.5*a*dt*dt;
                        const double ev=before[i].v[k]+a*dt;
                        const double reaction=-forces[i][k]*dt;
                        xe=std::max(xe,std::abs(particles[i].x[k]-ex));
                        ve=std::max(ve,std::abs(particles[i].v[k]-ev));
                        ie=std::max(ie,std::abs(before[i].mass*(particles[i].v[k]-before[i].v[k])+reaction));
                        if (!std::isfinite(particles[i].x[k]) || !std::isfinite(particles[i].v[k])
                            || std::abs(particles[i].omega[k])>1e-12)
                            throw std::runtime_error("Invalid particle state");
                    }
                    maxX=std::max(maxX,xe); maxV=std::max(maxV,ve); maxImpulse=std::max(maxImpulse,ie);
                    if (rank==0)
                    {
                        history << window << ',' << particles[i].id << ',' << particles[i].owner
                            << ',' << runTime.value() << ',' << demTime;
                        for (double q:particles[i].x) history << ',' << q;
                        for (double q:particles[i].v) history << ',' << q;
                        for (double q:forces[i]) history << ',' << q;
                        history << ',' << particles[i].mass << ',' << particles[i].radius
                            << ',' << xe << ',' << ve << ',' << ie << ',' << clock << '\n';
                    }
                }
            }
            // Reaction is an accounting ledger only; no fluid equation is solved.
            if (maxX>1e-10 || maxV>1e-10 || maxImpulse>1e-12 || maxClock>1e-12
                || (ranks==2 && migrations<2))
                throw std::runtime_error("Communication gates failed; inspect history");
            if (rank==0)
                std::cout << "DEMO_PASS ranks=" << ranks << " windows=20 migrations=" << migrations
                    << " max_x_error=" << maxX << " max_v_error=" << maxV
                    << " max_impulse_error=" << maxImpulse << " max_clock_error=" << maxClock << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "DEMO_FAIL rank=" << rank << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (ownMPI) MPI_Finalize();
    return 0;
}
