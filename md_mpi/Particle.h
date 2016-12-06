#include "Vector3.h"

#define Rc 4
#define Ucs 0.00097632 //0.00547944174

class Particle
{
    public:
        Particle();
        Particle(double x, double y, double z);
        Particle(Vector3 r, Vector3 v);
        ~Particle();
        Vector3 EffDist(Vector3 v, double L);
        Vector3 LJForce(Vector3 v, double L);
        double LJPotential(Vector3 v, double L);
        void ReduceR(double L);
        int GetRank(double d, int size);
        double GetK();
        Vector3 R;
        Vector3 V;
        Vector3 F;
};

Particle::Particle(double x, double y, double z)
{
    R = Vector3(x, y, z);
    V = Vector3();
    F = Vector3();
}

Particle::Particle(Vector3 r, Vector3 v)
    :R(r), V(v)
{
    F = Vector3();
}

Particle::Particle()
{
    R = Vector3(); 
    V = Vector3();
    F = Vector3();
}

Particle::~Particle()
{}

double Particle::GetK()
{
    return V * V / 2;
}

Vector3 Particle::EffDist(Vector3 v, double L)
{
    Vector3 dist = v - R;
 
    if (dist.X > L/2)
        dist.X -= L;
    else if (dist.X < -L/2)
        dist.X += L;

    if (dist.Y > L/2)
        dist.Y -= L;
    else if (dist.Y < -L/2)
        dist.Y += L;
    
    if (dist.Z > L/2)
        dist.Z -= L;
    else if (dist.Z < -L/2)
        dist.Z += L;
    
    return dist;
}

int Particle::GetRank(double d, int size)
{
    int rank = ((int)floor(R.Z / d)) % size;
    return ((rank >= 0) ? rank : (rank + size));
}

void Particle::ReduceR(double L)
{
    R.X = (R.X > 0) ? (fmod(R.X - L/2, L) + L/2) : (fmod(R.X + L/2, L) - L/2);
    R.Y = (R.Y > 0) ? (fmod(R.Y - L/2, L) + L/2) : (fmod(R.Y + L/2, L) - L/2);
}

double Particle::LJPotential(Vector3 v, double L)
{
    Vector3 r = EffDist(v, L);
    double r2 = r * r;
    if (r2 > Rc * Rc)
    {
        return 0;
    }
    double r6 = r2 * r2 * r2;
    return 4 * (1 / r6 - 1) / r6 + Ucs;
}

Vector3 Particle::LJForce(Vector3 v, double L)
{
    Vector3 r = EffDist(v, L);
    double r2 = r * r;
    if (r2 > Rc * Rc)
        return Vector3();
    double r6 = r2 * r2 * r2;
    return r * (24 * (2 / r6 - 1) / (r6 * r2));
}
