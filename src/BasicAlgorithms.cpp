//
//  BasicAlgorithms.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 16/01/2025.
//

#include "BigInt.hpp"
#include <sstream>
#include <thread>

std::pair<BigInt, BigInt> BruteForce(const BigInt &n)
{
    for (BigInt i = 2; i * i <= n; i++)
        if (n % i == 0ll)
            return std::pair(i, n / i);
    std::stringstream ErrorString;
    ErrorString << "Failed to factorize number " << n << " with brute force method";
    throw std::runtime_error(ErrorString.str());
}

std::pair<BigInt, BigInt> Fermat(const BigInt &n)
{
    for (BigInt a = n.SquareRoot(); a < n; a++)
    {
        BigInt b = a * a - n;
        if (b < 0ll) continue;
        BigInt c = b.SquareRoot();
        if (c * c == b) return std::pair(a - c, a + c);
    }
    std::stringstream ErrorString;
    ErrorString << "Failed to factorize number " << n << " with Fermat's method";
    throw std::runtime_error(ErrorString.str());
}

static std::pair<BigInt, BigInt> Result;

void Rho(const BigInt &n)
{
    while (true)
    {
        std::cout << "Starting Pollard's Rho\n";
        BigInt x = BigInt::Random(n - 2) + 2;
        BigInt y = x;
        BigInt c = BigInt::Random(n - 1) + 1;
        BigInt Divisor = 1;
        while (Divisor == 1)
        {
            if (Result.first != 0ll) return;
            x = ((x * x) % n + c) % n;
            y = ((y * y) % n + c) % n;
            y = ((y * y) % n + c) % n;
            Divisor = gcd(n, (x - y).abs());
        }
        if (Divisor != n && Result.first == 0ll) { Result = std::pair(Divisor, n / Divisor); return; }
    }
}

std::pair<BigInt, BigInt> PollardRho(const BigInt &n, uint32_t NumThreads = 1)
{
    NumThreads--;
    BigInt::Random(1);
    Result = { 0ll, 0ll };
    std::thread *Threads = static_cast<std::thread *>(::operator new(sizeof(std::thread) * NumThreads));
    for (int i = 0; i < NumThreads; i++) Threads[i] = std::thread(Rho, n);
    Rho(n);
    for (int i = 0; i < NumThreads; i++) Threads[i].join();
    for (int i = 0; i < NumThreads; i++) Threads[i].~thread();
    ::operator delete[](Threads);
    return Result;
}
