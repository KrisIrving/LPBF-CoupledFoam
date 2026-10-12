// SPDX-License-Identifier: GPL-3.0-or-later
// Opt-in diagnostic; no pressure/transfer algorithm or MPI transport changes.
#ifndef LPBF_JOINT_COMMUNICATION_AUDIT_H
#define LPBF_JOINT_COMMUNICATION_AUDIT_H
#include <mpi.h>
#include <fstream>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
namespace lpbfM2
{
class JointCommunicationAudit
{
    MPI_Comm comm_=MPI_COMM_NULL;
    int rank_=0;
    unsigned long long sequence_=0;
    [[noreturn]] static void fail(const char* message)
    {
        std::cerr<<message<<std::endl;
        // A local I/O or MPI error must not unwind into a collective destructor
        // while other ranks remain inside pressure communication.
        MPI_Abort(MPI_COMM_WORLD,1);
        throw std::runtime_error(message);
    }
public:
    explicit JointCommunicationAudit()
    {
        if(MPI_Comm_dup(MPI_COMM_WORLD,&comm_)!=MPI_SUCCESS)
            fail("Communication audit communicator duplication failed");
        MPI_Comm_rank(comm_,&rank_);
    }
    ~JointCommunicationAudit(){if(comm_!=MPI_COMM_NULL) MPI_Comm_free(&comm_);}
    JointCommunicationAudit(const JointCommunicationAudit&)=delete;
    JointCommunicationAudit& operator=(const JointCommunicationAudit&)=delete;
    void stage(unsigned long long code,double time,int correction,
               const std::vector<double>& input,bool cacheMiss)
    {
        static_assert(sizeof(double)==sizeof(unsigned long long),"Audit requires 64-bit double");
        unsigned long long hash=14695981039346656037ULL;
        for(double value:input)
        {
            unsigned long long bits=0;std::memcpy(&bits,&value,sizeof(bits));
            for(int i=0;i<8;++i){hash^=(bits>>(8*i))&255ULL;hash*=1099511628211ULL;}
        }
        const unsigned long long local[]={++sequence_,code,input.size(),hash,cacheMiss ? 1ULL:0ULL};
        // Rewrite one per-rank breadcrumb before the collective, so a later
        // crash retains the last entered stage without an unbounded log.
        std::ofstream log("joint-communication-rank"+std::to_string(rank_)+".txt",std::ios::trunc);
        log.precision(17);
        log<<"sequence="<<sequence_<<" stage="<<code<<" time="<<time
           <<" correction="<<correction<<" size="<<input.size()<<" hash="<<hash
           <<" cacheMiss="<<cacheMiss<<" status=entered\n";
        log.flush();
        if(!log) fail("Communication audit breadcrumb write failed");
        unsigned long long minimum[5],maximum[5];
        if(MPI_Allreduce(local,minimum,5,MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm_)!=MPI_SUCCESS
           ||MPI_Allreduce(local,maximum,5,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm_)!=MPI_SUCCESS)
            fail("Communication audit collective failed");
        for(int i=0;i<5;++i) if(minimum[i]!=maximum[i])
            fail("Joint communication stage/input/cache differs across ranks");
        log<<"status=agreed\n";log.flush();
        if(!log) fail("Communication audit breadcrumb flush failed");
    }
};
}
#endif
