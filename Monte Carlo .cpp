#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
// Modified version of the Numerical Recipes ran2 generator.
class Ran
{
private:
 uint64_t u, v, w;
 inline uint64_t int64()
 {
 u = u * 2862933555777941757UL + 704602954386353087UL;
 v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
 w = 4294957665U * (w & 0xffffffff) + (w >> 32);
 uint64_t x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
 return (x + v) ^ w;
 }
public:
 inline explicit Ran(uint64_t seed = 123456789UL) :
 u(0), v(4101842887655102017UL), w(1)
 {
 u = seed ^ v; int64();
 v = u; int64();
 w = v; int64();
 }
 inline operator double() { return 5.42101086242752217e-20 * int64(); }
};
Ran rng(4949939);
class Spins
{
private:
 std::vector<int> data_;
 int n_;
public:
 Spins(int n) : n_(n), data_((n+1) * (n+1)) {}
 int& operator() (int x, int y) {
 x = (x + n_) % n_; // Wrap periodic boundary conditions
 y = (y + n_) % n_;
 return data_[x + (n_ + 1) * y];
 }
 int operator() (int x, int y) const {
 x = (x + n_) % n_; // Wrap periodic boundary conditions
 y = (y + n_) % n_;
 return data_[x + (n_ + 1) * y];
 }
 int n() const { return n_; }
 };
int deltaKronecker(int i, int j) {
 if (i == j) {
 return 1;
 }
 return 0; // i and j different
}
void update(double beta, Spins& sigma, int q)
{
 int n = sigma.n();
 // Proposal : replace σ at a single site with a randomly chosen value
 for (int x = 0; x < n; x++)
 {
 for (int y = 0; y < n; y++)
 {
 int newSpin = floor(rng * q);
 int oldSpin = sigma(x, y);
 // Compute current local actionS
 int localActionS = (1 - deltaKronecker(sigma(x, y), sigma(x + 1, y))) + (1 - deltaKronecker(sigma(x, y), sigma(x - 1, y))) +
 (1 - deltaKronecker(sigma(x, y), sigma(x, y + 1))) + (1 - deltaKronecker(sigma(x, y), sigma(x, y - 1)));
 //cout << "Action S(sigma): " << actionS << "\n";
 sigma(x, y) = newSpin;
 // Compute new action
 int newLocalActionS = (1 - deltaKronecker(sigma(x, y), sigma(x + 1, y))) + (1 - deltaKronecker(sigma(x, y), sigma(x - 1, y))) +
 (1 - deltaKronecker(sigma(x, y), sigma(x, y + 1))) + (1 - deltaKronecker(sigma(x, y), sigma(x, y - 1)));
 // then compute min {1, exp{ −β(S(σ') − S(σ)) }
 double ds = beta * (newLocalActionS - localActionS);
 if (newLocalActionS < localActionS)
 {
 sigma(x, y) = newSpin;
 //cout << "Case 1 exp (-beta * (S(sigma') -S(sigma)) : " << exp(-ds) << " action/new Action " << localActionS << "/" <<
newLocalActionS << " [i,j] = " << x << "," << y << " New spin:" << newSpin << "\n";
 }
 else
 {
 if (exp(-ds) > rng)
 {
 sigma(x, y) = newSpin;
 //cout << "Case 2 exp (-beta * (S(sigma') -S(sigma)) : " << exp(-ds) << " action/new Action " << localActionS << "/" <<
newLocalActionS << " [i,j] = " << x << "," << y << " New spin:" << newSpin << "\n";
 }
 else
 {
 sigma(x, y) = oldSpin;
 //cout << "Case 3 exp (-beta * (S(sigma') -S(sigma)) : " << exp(-ds) << " action/new Action " << localActionS << "/" <<
newLocalActionS << " [i,j] = " << x << "," << y << " New spin:" << newSpin << "\n";
 }
 }
 }
 }
}
void initializeSpins(Spins& sigma, const int q)
{
 for (int x = 0; x < sigma.n(); x++)
 {
 for (int y = 0; y < sigma.n(); y++)
 {
 sigma(x, y) = floor(rng*q);
 }
 }
}
void printSpins(Spins& sigma)
{
 for (int x = 0; x < sigma.n(); x++)
 {
 for (int y = 0; y < sigma.n(); y++)
 {
 cout << sigma(x, y) << "\t";
 }
 cout << "\n";
 }
}
double fractionOfMostFrequentlyOccuringSpin(const Spins& sigma, const int q)
{
 std::vector<int> occurenceOfSpin(q,0);
 double max = 0.0;
 int L = sigma.n();
 for (int x = 0; x < sigma.n(); x++)
 {
 for (int y = 0; y < sigma.n(); y++)
 {
 occurenceOfSpin[sigma(x, y)]++;
 }
 }
 for (int k = 0; k < q; k++)
 {
 if (occurenceOfSpin[k] > max)
 {
 max = occurenceOfSpin[k];
 }
 }
/*
 for (int i = 0; i < q; i++)
 {
 cout << "Occurence Of Spin: [" << i << "] = " << occurenceOfSpin[i] << "\n";
 }
*/
 return (double) max/(L*L);
}
int main()
{
 // Regular LxL grid
 const int L = 20;
 // q-state Potts model
 int q = 3;
 double f;
 double mavg;
 Spins sigma(L);
 cout << "Grid size is " << sigma.n() << " L " << "\n";;
 cout << q << "-state Potts model" << "\n";
 initializeSpins(sigma, q);
 for (double beta = 0.5; beta <= 1.501; beta += 0.01)
 {
 mavg = 0;
 f = fractionOfMostFrequentlyOccuringSpin(sigma, q);
 //cout << "Initial configuration" << "\n";
 //printSpins(sigma);
 //cout << "Fraction of the Most Frequently Occuring Spin:" << fractionOfMostFrequentlyOccuringSpin(sigma, q) << "\n";
 //cout << "Fractional magnetisation of spin configuration M(sigma):" << (double)(q * f - 1) / (q - 1) << "\n";
 for (int t = 0; t < 10000; t++) {
 update(beta, sigma, q);
 f = fractionOfMostFrequentlyOccuringSpin(sigma, q);
 //cout << "Final configuration at " << beta << "\n";
 //printSpins(sigma);
 //cout << "Fraction of the Most Frequently Occuring Spin:" << f << " at " << beta << "\n";
 mavg += (double)(q * f - 1) / (q - 1);
 }
 cout << mavg/10000 << "\n";
 }
 }
