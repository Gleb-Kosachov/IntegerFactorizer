//
//  QuadraticSieve.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 17/01/2025.
//

#include "BigInt.hpp"
#include "blanczos.h"
#include <algorithm>
#include <fstream>
#include <thread>
#include <string>
#include <cmath>
#include <set>
#include <map>

int64_t BinPow(int64_t Base, uint64_t Power, int64_t Mod)
{
    if (Power == 0ll) return 1ll;
    if (Power == 1ll) return Base;
    int64_t Result = BinPow(Base, Power / 2, Mod);
    Result = (Result * Result) % Mod;
    if (Power % 2 == 1) Result = (Result * Base) % Mod;
    return Result;
}

BigInt SquareRootModP(const BigInt &a, int64_t p)
{
    if (a % p == 0ll) return 0ll;
    if (p == 2) return a % p;
    int64_t Mod8 = p % 8ll;
    if (Mod8 == 1)
    {
        int64_t t = p - 1, s = 0;
        while (t % 2 == 0)
        {
            s++;
            t /= 2;
        }
        int64_t d;
        for (d = 2;; d++)
            if (BinPow(d, (p - 1) / 2, p) == p - 1)
                break;
        int64_t A = BinPow((a % p).Int(), t, p);
        int64_t D = BinPow(d, t, p);
        int64_t m = 0;
        for (int64_t i = 0; i < s; i++)
            if (BinPow((A * BinPow(D, m, p)) % p, BinPow(2ll, s - 1 - i, p), p) == p - 1) m += (1 << i);
        return (BinPow((a % p).Int(), (t + 1) / 2, p) * BinPow(D, m / 2, p)) % p;
    }
    else if (Mod8 == 3 || Mod8 == 7)
        return BinPow((a % p).Int(), (p + 1) / 4, p);
    else if (Mod8 == 5)
    {
        BigInt Result = BinPow((a % p).Int(), (p + 3) / 8, p);
        if ((Result * Result) % p != a % p) Result = (Result * BinPow(2ll, (p - 1) / 4ll, p)) % p;
        return Result;
    }
    throw std::runtime_error(std::string("Modular square root of ") + a.String(10) + " does not exist");
}

BigInt ChooseMultiplier(const BigInt &n)
{
    static constexpr uint8_t MultiplierList[] = {
        1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19,
        21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38,
        39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59,
        61, 62, 65, 66, 67, 69, 70, 71, 73};
    double Scores[sizeof(MultiplierList)];
    for (int i = 0; i < sizeof(MultiplierList); i++)
    {
        Scores[i] = 0.5 * std::log(MultiplierList[i]);
        switch (((n * MultiplierList[i]) % 8).Int())
        {
            case 1:
                Scores[i] -= 2 * std::log(2);
                break;
            case 5:
                Scores[i] -= std::log(2);
                break;
            case 3:
            case 7:
                Scores[i] -= 0.5 * std::log(2);
                break;
        }
    }
    for (uint32_t Current = 3; Current < 2000; Current++)
    {
        bool prime = true;
        for (uint32_t Divisor = 2; Divisor * Divisor <= Current; Divisor++)
            if (Current % Divisor)
                prime = false;
        if (!prime) continue;
        double Contribution = std::log(static_cast<double>(Current)) / static_cast<double>(Current - 1);
        for (int i = 0; i < sizeof(MultiplierList); i++)
        {
            if (((n * MultiplierList[i]) % Current) == 0ll)
                Scores[i] -= Contribution;
            else if (BinPow(((n * MultiplierList[i]) % Current).Int(), (Current - 1) / 2, Current) == 1)
                Scores[i] -= 2.0 * Contribution;
        }
    }
    double BestScore = std::numeric_limits<double>::max();
    uint8_t BestMultiplier;
    for (int i = 0; i < sizeof(MultiplierList); i++)
    {
        if (Scores[i] < BestScore)
        {
            BestScore = Scores[i];
            BestMultiplier = MultiplierList[i];
        }
    }
    return n * BestMultiplier;
}

int64_t ChoosePrimeBound(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {20, 1200}, {30, 4000}, {40, 20000}, {60, 160000}, {80, 1300000}, {100, 2700000} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 6; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select B");
}

int64_t ChooseSieveSize(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {20, 65536}, {30, 65536}, {40, 65536}, {60, 65536}, {80, 3 * 65536}, {100, 9 * 65536} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 6; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select M");
}

uint32_t ChooseLargePrimeCutoff(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {20, 0}, {30, 10}, {40, 15}, {60, 30}, {80, 30}, {100, 50} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 6; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select large prime cutoff");
}

std::vector<int64_t> CalculateFactorBase(const BigInt &n, int64_t B)
{
    std::vector<bool> Sieve(B);
    std::vector<int64_t> Result = { 2 };
    for (uint64_t i = 1; i < B; i++) Sieve[i] = true;
    for (uint64_t i = 2; i <= B; i++)
    {
        if (!Sieve[i - 1]) continue;
        for (uint64_t j = i; i + j <= B; j += i)
            Sieve[i + j - 1] = false;
        if (i == 2) continue;
        if (BinPow((n % i).Int(), (i - 1) / 2, i) == 1) Result.push_back(i);
    }
    return Result;
}

std::vector<uint32_t> GetFactorOfASizes(const BigInt &OptimalValueOfA, const std::vector<int64_t> &FactorBase, uint32_t NumPrimesInFactorBase, int64_t M, std::vector<uint32_t> &FactorBounds)
{
    uint32_t OptimalSizeOfA = OptimalValueOfA.SizeInBase(2);
    FactorBounds.resize(std::log2(FactorBase.back()) + 3);
    uint32_t StartBits = 0;
    int i;

    for (i = 10; i < NumPrimesInFactorBase; i++)
    {
        uint32_t Length = std::log2(static_cast<double>(FactorBase[i])) + 0.5;
        if (Length > 15) break;
        if (Length > StartBits)
        {
            StartBits = Length;
            FactorBounds[StartBits] = i;
        }
    }
    FactorBounds[StartBits + 1] = i;

    if (OptimalSizeOfA > 210) StartBits = 15;
    else if (OptimalSizeOfA > 190) StartBits = 13;
    else if (OptimalSizeOfA > 180) StartBits = 12;
    else StartBits = 11;

    uint32_t NumFactorsOfA, Remainder;
    for (i = StartBits; i >= 7; i--)
    {
        NumFactorsOfA = OptimalSizeOfA / i;
        Remainder = OptimalSizeOfA % i;
        if (FactorBounds[i] == 0 || NumFactorsOfA == 1) continue;
        if (Remainder == 0 && NumFactorsOfA > 2 && FactorBounds[i + 1] > 0) break;
        if (Remainder <= NumFactorsOfA && NumFactorsOfA > 2 && FactorBounds[i + 1] > 0 && FactorBounds[i + 2] > 0) break;
        if ((i - Remainder) <= NumFactorsOfA && FactorBounds[i + 1] > 0 && FactorBounds[i - 1] > 0) break;
    }
    if (i < 7 || NumFactorsOfA < 2 || NumFactorsOfA > 20) throw std::runtime_error("Failed to select the coefficient a for polynomial");

    std::vector<uint32_t> FactorBits(NumFactorsOfA);
    for (int j = 0; j < NumFactorsOfA; j++) FactorBits[j] = i;
    if (Remainder <= NumFactorsOfA)
        for (int j = 0; j < Remainder; j++)
            FactorBits[j]++;
    else if ((i - Remainder) <= NumFactorsOfA)
    {
        FactorBits.push_back(i);
        NumFactorsOfA++;
        for (int j = 0; j < (i - Remainder); j++)
            FactorBits[j]--;
    }

    i = 0;
    for (int j = 0; j < NumFactorsOfA; j++)
        i += FactorBits[j];
    if (i != OptimalSizeOfA) throw std::runtime_error("Failed to select the coefficient a for polynomial");

    if (NumFactorsOfA >= 8 && NumFactorsOfA < 15)
    {
        if (FactorBits.front() > FactorBits.back())
        {
            switch (NumFactorsOfA) {
                default:
                    FactorBits[3]--;
                    FactorBits[2]--;
                case 9:
                    FactorBits[1]--;
                case 8:
                    FactorBits[0]--;
            }
        }
        else
        {
            switch (NumFactorsOfA) {
                default:
                    FactorBits[NumFactorsOfA - 4]--;
                    FactorBits[NumFactorsOfA - 3]--;
                case 9:
                    FactorBits[NumFactorsOfA - 2]--;
                case 8:
                    FactorBits[NumFactorsOfA - 1]--;
            }
        }
    }
    return FactorBits;
}

BigInt GetA(const BigInt &OptimalValueOfA, const std::vector<int64_t> &FactorBase, const std::vector<uint32_t> &FactorBounds, const std::vector<uint32_t> &FactorBits, std::vector<uint32_t> &FactorsOfA)
{
    BigInt a = 1;
    for (int i = 0; i < FactorBits.size() - 1;)
    {
        uint32_t Bits = FactorBits[i];
        uint32_t Range = FactorBounds[Bits + 1] - FactorBounds[Bits];
        FactorsOfA[i] = FactorBounds[Bits] + std::rand() % Range;
        int j;
        for (j = 0; j < i; j++)
            if (FactorsOfA[i] == FactorsOfA[j])
                break;
        if (j == i) a *= FactorBase[FactorsOfA[i++]];
    }
    int64_t TargetLastDivisor = (OptimalValueOfA / a).Int();
    uint32_t Left = std::lower_bound(FactorBase.begin(), FactorBase.end(), TargetLastDivisor) - FactorBase.begin();
    while (FactorBase[Left] > TargetLastDivisor) Left--;
    uint32_t Right = Left + 1;
    while (true)
    {
        if (TargetLastDivisor - FactorBase[Left] <= FactorBase[Right] - TargetLastDivisor)
        {
            if (std::find(FactorsOfA.begin(), FactorsOfA.end() - 1, Left) == FactorsOfA.end() - 1)
            {
                FactorsOfA.back() = Left;
                a *= FactorBase[Left];
                break;
            }
            else Left--;
        }
        else
        {
            if (std::find(FactorsOfA.begin(), FactorsOfA.end() - 1, Right) == FactorsOfA.end() - 1)
            {
                FactorsOfA.back() = Right;
                a *= FactorBase[Right];
                break;
            }
            else Right++;
        }
    }
    std::sort(FactorsOfA.begin(), FactorsOfA.end());
    return a;
}

void TrialDivision(const BigInt &n, const BigInt &a, const std::vector<uint32_t> &CurrentFactorIndices, const BigInt &b, const int64_t *RootsModP, int64_t M, uint32_t LargePrimeCutoff, const std::vector<int64_t> &FactorBase, uint32_t NumPrimesInFactorBase, uint32_t *Sieve, std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> *Matrix, std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> *P, std::mutex *Mutex, int64_t i, uint64_t &Counter)
{
    BigInt CongruentNum = a * i + b;
    CongruentNum *= CongruentNum;
    BigInt CurrentNum = CongruentNum - n;
    std::vector<uint32_t> Divisors;
    if (CurrentNum < 0ll) { Divisors.push_back(NumPrimesInFactorBase); CurrentNum = -CurrentNum; }
    std::vector<uint32_t>::const_iterator NextFactorIndex = CurrentFactorIndices.begin();
    for (int j = 0; j < NumPrimesInFactorBase; j++)
    {
        if (NextFactorIndex != CurrentFactorIndices.end() && *NextFactorIndex == j)
            NextFactorIndex++;
        else
        {
            int64_t jmodp = i % FactorBase[j];
            if (jmodp < 0) jmodp += FactorBase[j];
            if (jmodp != RootsModP[j * 2] && jmodp != RootsModP[j * 2 + 1]) continue;
        }
        while (CurrentNum % FactorBase[j] == 0ll)
        {
            CurrentNum /= FactorBase[j];
            if (!Divisors.empty() && Divisors.back() == j) Divisors.pop_back();
            else Divisors.push_back(j);
        }
    }
    if (CurrentNum == 1)
    {
        Mutex->lock();
        Matrix->emplace_back(std::make_pair(CongruentNum - n, std::move(CongruentNum)), std::move(Divisors));
        Mutex->unlock();
        Counter++;
    }
    else
    {
        Mutex->lock();
        bool FoundPair = P->contains(CurrentNum);
        if (FoundPair)
        {
            std::pair<BigInt, std::vector<uint32_t>> Pair = P->at(CurrentNum);
            Mutex->unlock();
            std::pair<BigInt, BigInt> Relation;
            Relation.first = (CongruentNum - n) * (Pair.first - n);
            Relation.second = CongruentNum * Pair.first;
            Mutex->lock();
            Matrix->emplace_back(std::move(Relation), std::vector<uint32_t>());
            int i = 0, j = 0;
            if (Divisors[i] == NumPrimesInFactorBase)
            {
                if (Pair.second[j] == NumPrimesInFactorBase) { i++; j++; }
                else Matrix->back().second.push_back(Divisors[i++]);
            }
            else if (Pair.second[j] == NumPrimesInFactorBase) Matrix->back().second.push_back(Pair.second[j++]);
            while (i < Divisors.size() || j < Pair.second.size())
            {
                if (i < Divisors.size() && j < Pair.second.size())
                {
                    if (Divisors[i] < Pair.second[j])
                        Matrix->back().second.push_back(Divisors[i++]);
                    else if (Divisors[i] > Pair.second[j])
                        Matrix->back().second.push_back(Pair.second[j++]);
                    else { i++; j++; }
                }
                else if (i < Divisors.size())
                    Matrix->back().second.push_back(Divisors[i++]);
                else Matrix->back().second.push_back(Pair.second[j++]);
            }
            if (Matrix->back().second.empty())
            {
                Matrix->pop_back();
                Counter--;
            }
            Mutex->unlock();
            Counter++;
        }
        else
        {
            P->emplace(CurrentNum, std::make_pair(std::move(CongruentNum), std::move(Divisors)));
            Mutex->unlock();
        }
    }
}

void AccumulateLogs(const BigInt &n, const BigInt &a, const std::vector<uint32_t> &CurrentFactorIndices, const BigInt &b, const int64_t *RootsModP, int64_t M, uint32_t LargePrimeCutoff, const std::vector<int64_t> &FactorBase, uint32_t NumPrimesInFactorBase, const std::vector<uint32_t> &FactorBaseLogs, uint32_t *Sieve, std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> *Matrix, std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> *P, std::mutex *Mutex)
{
    int64_t TwoM = 2 * M;
    Sieve[0] = a.SizeInBase(2);
    for (int64_t i = 1; i < 2 * M; i++) Sieve[i] = Sieve[i - 1];
    std::vector<uint32_t>::const_iterator NextFactorOfA = CurrentFactorIndices.begin();
    uint32_t const *SieveEnd = Sieve + TwoM;
    for (int i = 0; i < NumPrimesInFactorBase; i++)
    {
        if (NextFactorOfA != CurrentFactorIndices.end() && *NextFactorOfA == i) { NextFactorOfA++; continue; }
        int64_t CurrentPrime = FactorBase[i];
        int64_t j = -((M / CurrentPrime + 1) * CurrentPrime); j -= CurrentPrime;
        int64_t k = j + RootsModP[i * 2] + M;
        j += RootsModP[i * 2 + 1] + M;
        while (j < 0) j += CurrentPrime;
        while (k < 0) k += CurrentPrime;
        uint32_t log = FactorBaseLogs[i];
        if (k != j)
        {
            uint32_t *kptr = Sieve + k;
            while (kptr < SieveEnd)
            {
                *kptr += log;
                kptr += CurrentPrime;
            }
        }
        uint32_t *jptr = Sieve + j;
        while (jptr < SieveEnd)
        {
            *jptr += log;
            jptr += CurrentPrime;
        }
    }
}

void Sieve(const BigInt &n, const BigInt &a, const std::vector<uint32_t> &CurrentFactorIndices, const BigInt &b, const int64_t *RootsModP, int64_t M, uint32_t LargePrimeCutoff, const std::vector<int64_t> &FactorBase, uint32_t NumPrimesInFactorBase, const std::vector<uint32_t> &FactorBaseLogs, uint32_t *Sieve, std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> *Matrix, std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> *P, std::mutex *Mutex)
{
    AccumulateLogs(n, a, CurrentFactorIndices, b, RootsModP, M, LargePrimeCutoff, FactorBase, NumPrimesInFactorBase, FactorBaseLogs, Sieve, Matrix, P, Mutex);
    uint64_t Counter = 0;
    constexpr uint32_t Interval = 4;
    BigInt Num = a * (-M - Interval) + b;
    Num *= Num; Num -= n;
    BigInt SummandIncrease = 2 * a * a * Interval * Interval;
    BigInt Summand = 2 * a * a * (-M - 2 * Interval) * Interval + a * a * Interval * Interval + 2 * a * b * Interval;
    uint32_t CurrentSize = Num.SizeInBase(2);
    uint32_t *CurrentSievePtr = Sieve;
    for (int64_t i = -M, j = 0; i < M; i++, j++)
    {
        if (j % Interval == 0)
        {
            Summand += SummandIncrease;
            Num += Summand;
            CurrentSize = Num.SizeInBase(2);
        }
        if (CurrentSize > *(CurrentSievePtr++) + LargePrimeCutoff) continue;
        TrialDivision(n, a, CurrentFactorIndices, b, RootsModP, M, LargePrimeCutoff, FactorBase, NumPrimesInFactorBase, Sieve, Matrix, P, Mutex, i, Counter);
    }
}

void TraverseB(const BigInt &n, const BigInt &a, const std::vector<uint32_t> &CurrentFactorIndices, BigInt *B, int64_t M, uint32_t LargePrimeCutoff, const std::vector<int64_t> &FactorBase, uint32_t NumPrimesInFactorBase, const std::vector<uint32_t> &FactorBaseLogs, const std::vector<int64_t> &SqrtNModP, int64_t *InverseOfAModP, uint32_t NumFactorsOfA, int64_t *BTimesInverseA, int64_t *RootsModP, uint32_t *SieveArray, std::mutex *Mutex, std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> *Matrix, std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> *P)
{
    for (int i = 0; i < NumPrimesInFactorBase; i++)
        InverseOfAModP[i] = BinPow((a % FactorBase[i]).Int(), FactorBase[i] - 2, FactorBase[i]);
    for (int i = 0; i < NumFactorsOfA; i++)
        B[i] = (SqrtNModP[CurrentFactorIndices[i]] * (a / FactorBase[CurrentFactorIndices[i]]) * BinPow(((a / FactorBase[CurrentFactorIndices[i]]) % FactorBase[CurrentFactorIndices[i]]).Int(), FactorBase[CurrentFactorIndices[i]] - 2, FactorBase[CurrentFactorIndices[i]])) % a;

    uint32_t Mask = 0;
    BigInt b = 0ll;
    for (int i = 0; i < NumFactorsOfA; i++)
        b += B[i];
    for (int i = 0; i < NumPrimesInFactorBase; i++)
        for (int j = 0; j < NumFactorsOfA; j++)
            BTimesInverseA[i * NumFactorsOfA + j] = ((B[j] * InverseOfAModP[i]) % FactorBase[i]).Int();
    for (int i = 0; i < NumFactorsOfA; i++)
        B[i] *= 2ll;
    for (int i = 0; i < NumPrimesInFactorBase; i++)
    {
        RootsModP[i * 2] = (SqrtNModP[i] * InverseOfAModP[i]) % FactorBase[i];
        for (int j = 0; j < NumFactorsOfA; j++) RootsModP[i * 2] = (RootsModP[i * 2] + BTimesInverseA[i * NumFactorsOfA + j]) % FactorBase[i];
        RootsModP[i * 2] %= FactorBase[i];
        RootsModP[i * 2] = (FactorBase[i] - RootsModP[i * 2]) % FactorBase[i];
        RootsModP[i * 2 + 1] = (RootsModP[i * 2] + (2 * SqrtNModP[i] * InverseOfAModP[i]) % FactorBase[i]) % FactorBase[i];
    }
    Sieve(n, a, CurrentFactorIndices, b, RootsModP, M, LargePrimeCutoff, FactorBase, NumPrimesInFactorBase, FactorBaseLogs, SieveArray, Matrix, P, Mutex);
    for (int i = 1; i < (1 << (NumFactorsOfA - 1)); i++)
    {
        uint32_t Index = 0;
        while (((1 << Index) & i) == 0) Index++;
        bool Add = Mask & (1 << Index);
        Mask ^= (1 << Index);
        if (!Add) b -= B[Index];
        else b += B[Index];
        for (int j = 0; j < NumPrimesInFactorBase; j++)
        {
            RootsModP[j * 2] += FactorBase[j] * 2ll;
            RootsModP[j * 2 + 1] += FactorBase[j] * 2ll;
            if (!Add)
            {
                RootsModP[j * 2] += BTimesInverseA[j * NumFactorsOfA + Index] * 2ll;
                RootsModP[j * 2] %= FactorBase[j];
                RootsModP[j * 2 + 1] += BTimesInverseA[j * NumFactorsOfA + Index] * 2ll;
                RootsModP[j * 2 + 1] %= FactorBase[j];
            }
            else
            {
                RootsModP[j * 2] -= BTimesInverseA[j * NumFactorsOfA + Index] * 2ll;
                RootsModP[j * 2] %= FactorBase[j];
                RootsModP[j * 2 + 1] -= BTimesInverseA[j * NumFactorsOfA + Index] * 2ll;
                RootsModP[j * 2 + 1] %= FactorBase[j];
            }
        }
        Sieve(n, a, CurrentFactorIndices, b, RootsModP, M, LargePrimeCutoff, FactorBase, NumPrimesInFactorBase, FactorBaseLogs, SieveArray, Matrix, P, Mutex);
    }
}

std::pair<BigInt, BigInt> QuadraticSieve(const BigInt &n, uint32_t NumThreads = 1)
{
    std::chrono::high_resolution_clock::time_point beg = std::chrono::high_resolution_clock::now();
    std::cout << "Size of n is " << n.SizeInBase(2) << " bits\n";
    BigInt kn = ChooseMultiplier(n);
    std::vector<std::thread> Threads(NumThreads);
    std::mutex Mutex;
    int64_t PrimeBound = ChoosePrimeBound(kn);
    std::cout << "Factor base upper bound: " << PrimeBound << "\n";
    std::vector<int64_t> FactorBase = CalculateFactorBase(kn, PrimeBound);
    uint32_t NumPrimesInFactorBase = FactorBase.size();
    std::cout << "Factor base size: " << NumPrimesInFactorBase << "\n";
    std::vector<uint32_t> FactorBaseLogs(NumPrimesInFactorBase);
    for (int i = 0; i < NumPrimesInFactorBase; i++)
        FactorBaseLogs[i] = std::log2(static_cast<double>(FactorBase[i])) + 0.5;
    int64_t M = ChooseSieveSize(n);
    std::cout << "M: " << M << "\n";
    uint32_t LargePrimeCutoff = ChooseLargePrimeCutoff(n);
    std::cout << "Large prime cutoff: " << LargePrimeCutoff << "\n";
    std::vector<int64_t> SqrtNModP(NumPrimesInFactorBase);
    for (int i = 0; i < NumPrimesInFactorBase; i++)
        SqrtNModP[i] = SquareRootModP(kn, FactorBase[i]).Int();
    std::vector<uint32_t> FactorBaseBitBounds;
    BigInt OptimalValueOfA = ((kn * 2ll) / (M * M)).SquareRoot();
    std::vector<uint32_t> FactorOfABits = GetFactorOfASizes(OptimalValueOfA, FactorBase, NumPrimesInFactorBase, M, FactorBaseBitBounds);
    std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> Matrix;

    std::vector<std::vector<uint32_t>> CurrentFactorIndices(NumThreads);
    std::vector<std::vector<int64_t>> InverseOfAModP(NumThreads);
    std::vector<std::vector<BigInt>> B(NumThreads);
    std::vector<std::vector<int64_t>> BTimesInverseA(NumThreads);
    std::vector<std::vector<int64_t>> RootsModP(NumThreads);
    std::vector<std::vector<uint32_t>> Sieves(NumThreads);
    for (int i = 0; i < NumThreads; i++)
    {
        CurrentFactorIndices[i] = std::vector<uint32_t>(FactorOfABits.size());
        InverseOfAModP[i] = std::vector<int64_t>(NumPrimesInFactorBase);
        B[i] = std::vector<BigInt>(FactorOfABits.size());
        BTimesInverseA[i] = std::vector<int64_t>(FactorOfABits.size() * NumPrimesInFactorBase);
        RootsModP[i] = std::vector<int64_t>(NumPrimesInFactorBase * 2);
        Sieves[i] = std::vector<uint32_t>(M * 2);
    }
    std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> P;
    std::set<BigInt> UsedA;
    std::vector<BigInt> a(NumThreads);
    while (Matrix.size() < NumPrimesInFactorBase + 100)
    {
        for (int i = 0; i < NumThreads; i++)
        {
            a[i] = GetA(OptimalValueOfA, FactorBase, FactorBaseBitBounds, FactorOfABits, CurrentFactorIndices[i]);
            while (UsedA.contains(a[i]))
                a[i] = GetA(OptimalValueOfA, FactorBase, FactorBaseBitBounds, FactorOfABits, CurrentFactorIndices[i]);
            UsedA.insert(a[i]);
        }

        std::chrono::high_resolution_clock::time_point beg = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < NumThreads; i++)
            Threads[i] = std::thread(TraverseB, kn, a[i], CurrentFactorIndices[i], B[i].data(), M, LargePrimeCutoff, FactorBase, NumPrimesInFactorBase, FactorBaseLogs, SqrtNModP, InverseOfAModP[i].data(), static_cast<uint32_t>(FactorOfABits.size()), BTimesInverseA[i].data(), RootsModP[i].data(), Sieves[i].data(), &Mutex, &Matrix, &P);
        for (int i = 0; i < NumThreads; i++)
            Threads[i].join();
        std::cout << "Sieving completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - beg).count() << " milliseconds\nMatrix size: " << Matrix.size() << std::endl;

        if (Matrix.size() >= NumPrimesInFactorBase + 100)
        {
            std::sort(Matrix.begin(), Matrix.end());
            Matrix.resize(std::unique(Matrix.begin(), Matrix.end()) - Matrix.begin());
        }
    }

    FactorBase.push_back(-1);
    NumPrimesInFactorBase++;
    std::vector<uint32_t> NonZeroEntries;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].second.size(); j++)
        {
            NonZeroEntries.push_back(Matrix[i].second[j]);
            NonZeroEntries.push_back(i);
        }
    }
    std::vector<uint64_t> NullSpaces(Matrix.size());
    uint32_t NumNullSpaces = blanczos(NonZeroEntries.data(), NonZeroEntries.size() / 2, NumPrimesInFactorBase, Matrix.size(), NullSpaces.data());
    for (int i = 0; i < NumNullSpaces; i++)
    {
        BigInt A = 1, B = 1;
        for (int j = 0; j < Matrix.size(); j++)
        {
            if ((1ull << i) & NullSpaces[j])
            {
                A *= Matrix[j].first.first;
                B *= Matrix[j].first.second;
            }
        }
        A = A.SquareRoot(); B = B.SquareRoot();
        BigInt gcd1 = gcd(n, (A - B).abs()), gcd2 = gcd(n, (A + B).abs());
        if (gcd1 != 1 && gcd1 != n) return std::make_pair(gcd1, n / gcd1);
        if (gcd2 != 1 && gcd2 != n) return std::make_pair(gcd2, n / gcd2);
    }
    return std::pair(0, 0);
}
