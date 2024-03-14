#include "triple.h"

#include <cmath> // fmin
#include <exception>
#include <iostream>

#include "json/json.h"

using namespace std;
using json = nlohmann::json;

// --- Constructors ------------------------------------------------------------

Triple::Triple(double X, double Y, double Z) : x(X), y(Y), z(Z) {}

Triple::Triple(json const &node)
{
  if (!node.is_array())
    throw runtime_error("Triple(): JSON node is not an array");

  if (node.size() != 3)
    throw runtime_error("Triple(): JSON node.size() != 3");

  if (!node[0].is_number())
    throw runtime_error("Triple(): JSON node is not a number");

  set(node[0], node[1], node[2]);
}

// --- Operators ---------------------------------------------------------------

Triple Triple::operator+(Triple const &t) const { return Triple(x + t.x, y + t.y, z + t.z); }

Triple Triple::operator+(double f) const { return Triple(x + f, y + f, z + f); }

Triple Triple::operator-() const { return Triple(-x, -y, -z); }

Triple Triple::operator-(Triple const &t) const { return Triple(x - t.x, y - t.y, z - t.z); }

Triple Triple::operator-(double f) const { return Triple(x - f, y - f, z - f); }

Triple Triple::operator*(Triple const &t) const { return Triple(x * t.x, y * t.y, z * t.z); }

Triple Triple::operator*(double f) const { return Triple(x * f, y * f, z * f); }

Triple Triple::operator/(double f) const
{
  double invf = 1.0 / f;
  return Triple(x * invf, y * invf, z * invf);
}

// --- Compound operators ------------------------------------------------------

Triple &Triple::operator+=(Triple const &t)
{
  x += t.x;
  y += t.y;
  z += t.z;
  return *this;
}

Triple &Triple::operator+=(double f)
{
  x += f;
  y += f;
  z += f;
  return *this;
}

Triple &Triple::operator-=(Triple const &t)
{
  x -= t.x;
  y -= t.y;
  z -= t.z;
  return *this;
}

Triple &Triple::operator-=(double f)
{
  x -= f;
  y -= f;
  z -= f;
  return *this;
}

Triple &Triple::operator*=(double f)
{
  x *= f;
  y *= f;
  z *= f;
  return *this;
}

Triple &Triple::operator/=(double f)
{
  double invf = 1.0 / f;
  x *= invf;
  y *= invf;
  z *= invf;
  return *this;
}

// --- Vector Operators --------------------------------------------------------

double Triple::dot(Triple const &t) const { return x * t.x + y * t.y + z * t.z; }

Triple Triple::cross(Triple const &t) const { return Triple(y * t.z - z * t.y, z * t.x - x * t.z, x * t.y - y * t.x); }

double Triple::length() const { return sqrt(length_2()); }

double Triple::length_2() const { return x * x + y * y + z * z; }

Triple Triple::normalized() const { return (*this) / length(); }

void Triple::normalize()
{
  double len = length();
  double invlen = 1.0 / len;
  x *= invlen;
  y *= invlen;
  z *= invlen;
}

// --- Color functions ---------------------------------------------------------

void Triple::set(double f)
{
  r = f;
  g = f;
  b = f;
}

void Triple::set(double f, double maxValue) { set(f / maxValue); }
void Triple::set(double red, double green, double blue)
{
  r = red;
  g = green;
  b = blue;
}

void Triple::set(double red, double green, double blue, double maxValue)
{
  set(red / maxValue, green / maxValue, blue / maxValue);
}

Triple &Triple::clamp(double maxValue)
{
  r = fmin(r, maxValue);
  g = fmin(g, maxValue);
  b = fmin(b, maxValue);
  return *this;
}

// --- Free Operators ----------------------------------------------------------

Triple operator+(double f, Triple const &t) { return Triple(f + t.x, f + t.y, f + t.z); }

Triple operator-(double f, Triple const &t) { return Triple(f - t.x, f - t.y, f - t.z); }

Triple operator*(double f, Triple const &t) { return Triple(f * t.x, f * t.y, f * t.z); }

Triple reflect(Triple const &incident, Triple const &normal) { return incident - 2.0 * normal.dot(incident) * normal; }

Triple refract(Triple const &incident, Triple const &normal, double const &nt)
{
  Vector N = normal;
  double cosi = -incident.dot(N);
  double n_i = 1.0;
  double n_t = nt;

  if (cosi > 0) {
    // ray travels from inside to outside
    N = -N;
    cosi = -cosi;
    double temp = n_t;
    n_t = n_i;
    n_i = temp;
  }
  double cosT2 = 1.0 - n_i * n_i * (1.0 - cosi * cosi) / (n_t * n_t);

  double eta = n_i / n_t;
  // If cosT2 is negative, it means total internal reflection occurs
  if (cosT2 < 0)
    return Vector(0, 0, 0); // Return zero vector
  Triple T = eta * (incident - cosi * N) - N * sqrt(cosT2);
  return T;
}

double schlick(double n1, double n2, const Triple &normal, const Triple &incident)
{
  double r0 = pow((n1 - n2) / (n1 + n2), 2);
  double cosX = -normal.x * incident.x - normal.y * incident.y - normal.z * incident.z;

  // If n1 > n2, apply Schlick's approximation
  if (n1 > n2)
  {
    double n = n1 / n2;
    double sinT2 = n * n * (1.0 - cosX * cosX);

    // Total internal reflection occurs
    if (sinT2 > 1.0)
      return 1.0;

    // Compute reflection coefficient using Schlick's approximation
    double cosT = sqrt(1.0 - sinT2);
    double x = (1.0 - cosT) / (1.0 + cosT);
    return r0 + (1.0 - r0) * pow(x, 5);
  }
  else
  {
    // Otherwise, use Fresnel equations directly
    double n = n2 / n1;
    double sinT2 = n * n * (1.0 - cosX * cosX);

    // Total internal reflection occurs
    if (sinT2 > 1.0)
      return 1.0;

    // Compute reflection coefficient using Fresnel equations
    double cosT = sqrt(1.0 - sinT2);
    double Rs = pow((n1 * cosX - n2 * cosT) / (n1 * cosX + n2 * cosT), 2);
    double Rp = pow((n1 * cosT - n2 * cosX) / (n1 * cosT + n2 * cosX), 2);
    return 0.5 * (Rs + Rp);
  }
}

// --- IO Operators ------------------------------------------------------------

istream &operator>>(istream &is, Triple &t)
{
  double x, y, z;
  //  is >> x >> y >> z;      // is not guaranteed to work pre C++17
  is >> x;
  is >> y;
  is >> z;
  t.set(x, y, z); // only assign if everything is extracted
  return is;
}

ostream &operator<<(ostream &os, Triple const &t)
{
  // format: [x, y, z] (no newline)
  os << '[' << t.x << ", " << t.y << ", " << t.z << ']';
  return os;
}
