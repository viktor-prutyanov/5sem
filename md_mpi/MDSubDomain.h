#include <list>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cassert>
#include <mpi.h>

#include "Particle.h"

#define CRYSTAL
#define NDEBUG
#define TAU 0.1
#define RND_PARAM 0.05

class MDSubDomain
{
    public:
        MDSubDomain(unsigned long int n, double l, int mSize, int mRank, double targetT);
        void MacroDump(FILE *file, double t);
        void MicroDump(FILE *file);
        void Step(double dt, bool isThOn);
        double K;
        double U;
        unsigned long int N_real;
        ~MDSubDomain();

    private:
        std::list<Particle> particles; 
        unsigned long int N;
        double L;
        double d;
        int mRank;
        int mSize;
        double lambda;
        double targetK;
        int leftRank;
        int rghtRank;
        std::vector<Vector3> sendBuf;
        std::vector<Vector3> recvLeftBuf;
        std::vector<Vector3> recvRghtBuf;
        std::vector<Vector3> sendLeftBuf;
        std::vector<Vector3> sendRghtBuf;
};

MDSubDomain::MDSubDomain(unsigned long int N, double L, int mSize, int mRank, double targetT)
    :K(0.), U(0.), N_real(0), 
    particles(std::list<Particle>()), N(N), L(L), d(L / mSize), 
    mRank(mRank), mSize(mSize), 
    lambda(1.), targetK(1.5 * targetT * N),
    sendBuf(std::vector<Vector3>()), 
    recvLeftBuf(std::vector<Vector3>()), 
    recvRghtBuf(std::vector<Vector3>()),
    sendLeftBuf(std::vector<Vector3>()), 
    sendRghtBuf(std::vector<Vector3>())
{
    unsigned int p = (unsigned int)ceil(pow(N / (mSize*mSize*mSize), 1./3.));
    unsigned int q = mSize * p;
    double step = L / q;
    unsigned long int rem = (q*q*q - N) / mSize + ((((q*q*q - N) % mSize) > mRank) ? 1 : 0);

    for (unsigned long int i = 0; i < q; ++i)
        for (unsigned long int j = 0; j < q; ++j)
            for (unsigned long int k = 0; k < p; ++k)
            {
                double x = step * (0.5 + i + RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX);
                double y = step * (0.5 + j + RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX);
                double z = step * (0.5 + k + RND_PARAM * (rand() * 1. - RAND_MAX/2) / RAND_MAX) + mRank * d;
                particles.emplace_back(x, y, z);  
            }
    
    std::vector<unsigned long int> v;
    for (unsigned long int i = 0; i < q*q*p; ++i)
        v.push_back(i);

    std::random_shuffle(v.begin(), v.end());

    unsigned long int victim, k;
    for (unsigned long int i = 0; i < rem; ++i)
    {
        victim = v[i];
        k = 0;
        for (auto j = particles.begin(); j != particles.end(); ++j)
        {
            if (k == victim)
                (*j).R = Vector3();
            ++k;
        }
    }
    
    particles.remove_if([](Particle p) -> bool { return ((p.R.X == 0.) && (p.R.Y == 0.) && (p.R.Z == 0.)); });

    if (mRank == 0)
    {
        rghtRank = mRank + 1;
        leftRank = mSize - 1;
    }
    else if (mRank == mSize - 1)
    {
        rghtRank = 0;
        leftRank = mRank - 1;
    }
    else
    {
        rghtRank = mRank + 1;
        leftRank = mRank - 1;
    }

    sendBuf.reserve(N / mSize);
}

MDSubDomain::~MDSubDomain()
{}

void MDSubDomain::Step(double dt, bool isThOn)
{
    /*
    double u = 0.;
    sendBuf.clear();
    
    for (auto i = particles.begin(); i != particles.end(); ++i)
    {
        (*i).F = Vector3();
        for (auto j = particles.begin(); j != particles.end(); ++j)
            if (i != j)
            {
                (*i).F -= (*i).LJForce((*j).R, L);
                u += (*i).LJPotential((*j).R, L) / 2;
            }
    }
    
    for (auto p : particles)
        sendBuf.push_back(p.R);
    
    particles.remove_if([d, mSize](Particle p) -> bool { return (p.GetRank(d, mSize) != mRank); });

    unsigned long int recvLeftNum, recvRghtNum, sendNum = sendBuf.size();
    
    MPI_Send(&sendNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD);
    MPI_Send(&sendNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD);
    MPI_Recv(&recvLeftNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&recvRghtNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD, &status);

    recvLeftBuf.resize(recvLeftNum);
    recvRghtBuf.resize(recvRghtNum);
    
    MPI_Sendrecv(sendBuf.data(), sendBuf.size() * 3, MPI_DOUBLE, leftRank, 0, recvRghtBuf.data(), recvRghtBuf.size() * 3, MPI_DOUBLE, rghtRank, 0, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(sendBuf.data(), sendBuf.size() * 3, MPI_DOUBLE, rghtRank, 0, recvLeftBuf.data(), recvLeftBuf.size() * 3, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD, &status);
    
    for (auto& p : particles)
    {
        for (auto r : recvLeftBuf)
        {
            p.F -= p.LJForce(r, L);
            u += p.LJPotential(r, L) / 2;
        }
        for (auto r : recvRghtBuf)
        {
            p.F -= p.LJForce(r, L);
            u += p.LJPotential(r, L) / 2;
        }
    }
*/
    double u = 0.;
    for (auto i = particles.begin(); i != particles.end(); ++i)
    {
        (*i).F = Vector3();
        for (auto j = particles.begin(); j != particles.end(); ++j)
            if (i != j)
            {
                (*i).F -= (*i).LJForce((*j).R, L);
                u += (*i).LJPotential((*j).R, L) / 2;
            }
    }

    sendBuf.clear();
    for (auto p : particles)
        sendBuf.push_back(p.R);

    unsigned long int recvLeftNum, recvRghtNum, sendNum = sendBuf.size();
    MPI_Status status;

    MPI_Send(&sendNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD);
    MPI_Send(&sendNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD);
    MPI_Recv(&recvLeftNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&recvRghtNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD, &status);
    
    recvLeftBuf.resize(recvLeftNum);
    recvRghtBuf.resize(recvRghtNum);

    //printf("sendNum = %lu, recvLeftNum = %lu, recvRghtNum = %lu\n", sendBuf.size(), recvLeftBuf.size(), recvRghtBuf.size());
    MPI_Sendrecv(sendBuf.data(), sendBuf.size() * 3, MPI_DOUBLE, leftRank, 0, recvRghtBuf.data(), recvRghtBuf.size() * 3, MPI_DOUBLE, rghtRank, 0, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(sendBuf.data(), sendBuf.size() * 3, MPI_DOUBLE, rghtRank, 0, recvLeftBuf.data(), recvLeftBuf.size() * 3, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD, &status);

    sendLeftBuf.clear();
    sendRghtBuf.clear();

    for (auto& p : particles)
    {
        for (auto r : recvLeftBuf)
        {
            p.F -= p.LJForce(r, L);
            u += p.LJPotential(r, L) / 2;
        }
        for (auto r : recvRghtBuf)
        {
            p.F -= p.LJForce(r, L);
            u += p.LJPotential(r, L) / 2;
        }
        p.V += p.F * dt;
        p.R += p.V * dt;
        p.ReduceR(L);
        
        if (p.GetRank(d, mSize) == leftRank)
        {
            sendLeftBuf.push_back(p.R);
            sendLeftBuf.push_back(p.V);
        }
        else if (p.GetRank(d, mSize) == rghtRank)
        {
            sendRghtBuf.push_back(p.R);
            sendRghtBuf.push_back(p.V);
        }
        else if (p.GetRank(d, mSize) != mRank)
        {
            printf("rank = %d, mRank = %d", p.GetRank(d, mSize), mRank);
            assert(("Particle goes further then neigbour rank.", false));
        }
    }
    
    particles.remove_if([this](Particle p) -> bool { return (p.GetRank(d, mSize) != mRank); });

    unsigned long int sendLeftNum = sendLeftBuf.size();
    unsigned long int sendRghtNum = sendRghtBuf.size();
   
    MPI_Send(&sendLeftNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD);
    MPI_Send(&sendRghtNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD);
    MPI_Recv(&recvLeftNum, 1, MPI_LONG_INT, leftRank, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&recvRghtNum, 1, MPI_LONG_INT, rghtRank, 0, MPI_COMM_WORLD, &status);

    recvLeftBuf.resize(recvLeftNum);
    recvRghtBuf.resize(recvRghtNum);
    
    MPI_Send(sendLeftBuf.data(), sendLeftNum * 3, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD);
    MPI_Send(sendRghtBuf.data(), sendRghtNum * 3, MPI_DOUBLE, rghtRank, 0, MPI_COMM_WORLD);
    MPI_Recv(recvLeftBuf.data(), recvLeftNum * 3, MPI_DOUBLE, leftRank, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(recvRghtBuf.data(), recvRghtNum * 3, MPI_DOUBLE, rghtRank, 0, MPI_COMM_WORLD, &status);

    for (unsigned long int i = 0; i < recvLeftNum; i += 2)
        particles.push_back(Particle(recvLeftBuf[i], recvLeftBuf[i + 1]));
    
    for (unsigned long int i = 0; i < recvRghtNum; i += 2)
        particles.push_back(Particle(recvRghtBuf[i], recvRghtBuf[i + 1]));

    double k = 0.;
    for (auto p : particles)
        k += p.GetK(); 
   
    unsigned long int n = particles.size();

    MPI_Allreduce(&n, &N_real, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&k, &K, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&u, &U, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    lambda = sqrt(1 + dt/TAU * (targetK/K - 1));
    if (isThOn)
        for (auto& p : particles)
            p.V *= lambda;
}

void MDSubDomain::MacroDump(FILE *file, double t)
{
    fprintf(file, "%lg,%lg,%lg,%lg,%lu\n", t, K, U, K + U, N_real);
}

void MDSubDomain::MicroDump(FILE *file)
{
    long unsigned int i = 0;
    for (auto p: particles)
    {
        printf("\t(%d/%d) [%lu/%lu] r={%lg %lg %lg} v={%lg %lg %lg} f={%lg %lg %lg}\n",
            mRank, mSize, i, particles.size(), 
            p.R.X, p.R.Y, p.R.Z, p.V.X, p.V.Y, p.V.Z, p.F.X, p.F.Y, p.F.Z);
        ++i;
    }
}
