class Vector3
{
    public:
        Vector3();
        Vector3(double x, double y, double z);
        ~Vector3();
        double X, Y, Z;
};

Vector3::Vector3()
    :X(0.), Y(0.), Z(0.)
{}

Vector3::Vector3(double x, double y, double z) 
    : X(x), Y(y), Z(z)
{}

Vector3::~Vector3()
{}

double operator*(const Vector3 &a, const Vector3 &b)
{
    return a.X * b.X + a.Y * b.Y + a.Z * b.Z;
}

Vector3 operator-(const Vector3 &a, const Vector3 &b)
{
    return Vector3(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
}

Vector3 operator+(const Vector3 &a, const Vector3 &b)
{
    return Vector3(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
}

Vector3& operator-=(Vector3 &a, const Vector3 &b)
{
    a.X -= b.X;
    a.Y -= b.Y;
    a.Z -= b.Z;
    return a;
}

Vector3& operator*=(Vector3 &a, double k)
{
    a.X *= k;
    a.Y *= k;
    a.Z *= k;
    return a;
}

Vector3& operator+=(Vector3 &a, const Vector3 &b)
{
    a.X += b.X;
    a.Y += b.Y;
    a.Z += b.Z;
    return a;
}

Vector3 operator*(const Vector3 &a, const double k)
{
    return Vector3(a.X * k, a.Y * k, a.Z * k);
}
