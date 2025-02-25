//
//  QuadraticSieve.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 17/01/2025.
//

#include "BigInt.hpp"
#include "blanczos.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <memory>
#include <thread>
#include <chrono>
#include <mutex>
#include <string>
#include <vector>
#include <barrier>
#include <cmath>
#include <set>
#include <map>


struct QuadraticSieveData
{
    BigInt n;
    BigInt OptimalValueOfA;
    
    int64_t PrimeBound; // The upper bound of primes in the factor base (normally denoted as B in literature)
    int64_t M; // The sieve size (actual size is 2 * M)
    int64_t Delta; // The value used to calculate the trial division threshold
    uint32_t SmallPrimeBound; // The number of small primes to be skipped when sieving
    
    std::vector<int64_t> FactorBase;
    std::vector<int64_t> SqrtNModP; // Modular square roots of n for each prime in the factor base
    std::vector<uint8_t> FactorBaseLogs; // The logarithm (base 2) of each prime in the factor base
    
    std::vector<uint32_t> FactorOfABits;
    std::vector<uint32_t> FactorBaseBitBounds;
    
    std::vector<std::vector<uint32_t>> CurrentFactorIndices; // Indices of factors of the coefficient a in the factor base
    std::vector<std::vector<int64_t>> InverseOfAModP; // Modular inverses
    std::vector<std::vector<BigInt>> B; // The values used to compute the modular square root of n via the chinese remainder theorem
    std::vector<std::vector<int64_t>> BTimesInverseA;
    std::vector<std::vector<int64_t>> RootsModP; // The roots of the polynomial for each prime in the factor base
    std::vector<std::vector<uint8_t>> Sieves;
    std::vector<std::vector<uint8_t *>> SievePointers;
    uint8_t SmallThreshold, LargeThreshold;
    
    std::vector<BigInt> a; // A coefficient of the polynomial
    std::set<BigInt> UsedA;
    std::mutex AMutex;
    
    std::vector<std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>>> Relations;
    std::mutex RelationMutex;
    
    std::map<BigInt, std::pair<BigInt, std::vector<uint32_t>>> P; // A set of partial relations used for the large prime variation
    std::mutex PMutex;
    
    bool Running;
};


// Binary exponentiation algorithm
int64_t BinPow(int64_t Base, uint64_t Power, int64_t Mod)
{
    if (Power == 0ll) return 1ll;
    if (Power == 1ll) return Base;
    int64_t Result = BinPow(Base, Power / 2, Mod);
    Result = (Result * Result) % Mod;
    if (Power % 2 == 1) Result = (Result * Base) % Mod;
    return Result;
}

// Tonelli-Shanks algorithm
BigInt SquareRootModP(const BigInt &a, int64_t p)
{
    if (a % p == 0ll) return static_cast<int64_t>(0);
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

int64_t ChoosePrimeBound(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {30, 4000}, {40, 20000}, {60, 160000}, {80, 1300000}, {100, 2700000}, {110, 6300000} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 7; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select prime bound");
}

int64_t ChooseSieveSize(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {30, 65536}, {40, 65536}, {60, 65536}, {80, 3 * 65536}, {100, 9 * 65536}, {110, 13 * 65536} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 7; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select M");
}

uint32_t ChooseDelta(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {30, 10}, {40, 15}, {60, 30}, {80, 30}, {100, 50}, {110, 50} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 7; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select delta");
}

uint32_t ChooseSmallPrimeBound(const BigInt &n)
{
    std::pair<uint32_t, int64_t> Values[] = { {30, 0}, {40, 10}, {60, 100}, {80, 200}, {100, 200}, {110, 200} };
    uint32_t Size = n.SizeInBase(10);
    for (int i = 0; i < 7; i++)
        if (std::max(Size, Values[i].first) - std::min(Size, Values[i].first) <= 3)
            return Values[i].second;
    throw std::runtime_error("Failed to select small prime bound");
}

void CalculateFactorBase(QuadraticSieveData *Data)
{
    std::vector<bool> Sieve(Data->PrimeBound);
    Data->FactorBase.push_back(2);
    for (uint64_t i = 1; i < Data->PrimeBound; i++) Sieve[i] = true;
    for (uint64_t i = 2; i <= Data->PrimeBound; i++)
    {
        if (!Sieve[i - 1]) continue;
        for (uint64_t j = i; i + j <= Data->PrimeBound; j += i)
            Sieve[i + j - 1] = false;
        if (i == 2) continue;
        if (BinPow((Data->n % i).Int(), (i - 1) / 2, i) == 1) Data->FactorBase.push_back(i);
    }
    
    Data->FactorBaseLogs.resize(Data->FactorBase.size());
    for (int i = 0; i < Data->FactorBase.size(); i++)
        Data->FactorBaseLogs[i] = std::log2(static_cast<double>(Data->FactorBase[i])) + 0.5;
}

// The following three functions (selection of the multiplier k for the number n and the coefficient a for the polynomial f(x)) are taken from msieve
BigInt ChooseMultiplier(const BigInt &n)
{
    static constexpr uint8_t MultiplierList[] = {
        1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19,
        21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38,
        39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59,
        61, 62, 65, 66, 67, 69, 70, 71, 73};
    std::vector<BigInt> kn(sizeof(MultiplierList));
    for (int i = 0; i < sizeof(MultiplierList); i++)
        kn[i] = n * MultiplierList[i];
    double Scores[sizeof(MultiplierList)];
    for (int i = 0; i < sizeof(MultiplierList); i++)
    {
        Scores[i] = 0.5 * std::log(MultiplierList[i]);
        switch ((kn[i] % 8).Int())
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
        
        // In msieve, the contribution of a prime p is log(p) / (p - 1). However, unlinke msieve,
        // my implementation does not deal with squares of prime numbers in the factor base. For
        // for this reason, the contribution is log(p) / p
        double Contribution = std::log2(static_cast<double>(Current)) / static_cast<double>(Current);
        for (int i = 0; i < sizeof(MultiplierList); i++)
        {
            if ((kn[i] % Current) == 0ll)
                Scores[i] -= Contribution;
            else if (BinPow((kn[i] % Current).Int(), (Current - 1) / 2, Current) == 1)
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
    std::cout << "Using multiplier " << static_cast<int>(BestMultiplier) << "\n";
    return n * BestMultiplier;
}

void GetFactorOfASizes(QuadraticSieveData *Data)
{
    uint32_t OptimalSizeOfA = Data->OptimalValueOfA.SizeInBase(2);
    Data->FactorBaseBitBounds.resize(std::log2(Data->FactorBase.back()) + 3);
    uint32_t StartBits = 0;
    int i;

    for (i = 10; i < Data->FactorBase.size(); i++)
    {
        uint32_t Length = std::log2(static_cast<double>(Data->FactorBase[i])) + 0.5;
        if (Length > 15) break;
        if (Length > StartBits)
        {
            StartBits = Length;
            Data->FactorBaseBitBounds[StartBits] = i;
        }
    }
    Data->FactorBaseBitBounds[StartBits + 1] = i;

    if (OptimalSizeOfA > 210) StartBits = 15;
    else if (OptimalSizeOfA > 190) StartBits = 13;
    else if (OptimalSizeOfA > 180) StartBits = 12;
    else StartBits = 11;

    uint32_t NumFactorsOfA, Remainder;
    for (i = StartBits; i >= 7; i--)
    {
        NumFactorsOfA = OptimalSizeOfA / i;
        Remainder = OptimalSizeOfA % i;
        if (Data->FactorBaseBitBounds[i] == 0 || NumFactorsOfA == 1) continue;
        if (Remainder == 0 && NumFactorsOfA > 2 && Data->FactorBaseBitBounds[i + 1] > 0) break;
        if (Remainder <= NumFactorsOfA && NumFactorsOfA > 2 && Data->FactorBaseBitBounds[i + 1] > 0 && Data->FactorBaseBitBounds[i + 2] > 0) break;
        if ((i - Remainder) <= NumFactorsOfA && Data->FactorBaseBitBounds[i + 1] > 0 && Data->FactorBaseBitBounds[i - 1] > 0) break;
    }
    if (i < 7 || NumFactorsOfA < 2 || NumFactorsOfA > 20) throw std::runtime_error("Failed to select the coefficient a for polynomial");

    Data->FactorOfABits.resize(NumFactorsOfA);
    for (int j = 0; j < NumFactorsOfA; j++) Data->FactorOfABits[j] = i;
    if (Remainder <= NumFactorsOfA)
        for (int j = 0; j < Remainder; j++)
            Data->FactorOfABits[j]++;
    else if ((i - Remainder) <= NumFactorsOfA)
    {
        Data->FactorOfABits.push_back(i);
        NumFactorsOfA++;
        for (int j = 0; j < (i - Remainder); j++)
            Data->FactorOfABits[j]--;
    }

    i = 0;
    for (int j = 0; j < NumFactorsOfA; j++)
        i += Data->FactorOfABits[j];
    if (i != OptimalSizeOfA) throw std::runtime_error("Failed to select the coefficient a for polynomial");

    if (NumFactorsOfA >= 8 && NumFactorsOfA < 15)
    {
        if (Data->FactorOfABits.front() > Data->FactorOfABits.back())
        {
            switch (NumFactorsOfA) {
                default:
                    Data->FactorOfABits[3]--;
                    Data->FactorOfABits[2]--;
                case 9:
                    Data->FactorOfABits[1]--;
                case 8:
                    Data->FactorOfABits[0]--;
            }
        }
        else
        {
            switch (NumFactorsOfA) {
                default:
                    Data->FactorOfABits[NumFactorsOfA - 4]--;
                    Data->FactorOfABits[NumFactorsOfA - 3]--;
                case 9:
                    Data->FactorOfABits[NumFactorsOfA - 2]--;
                case 8:
                    Data->FactorOfABits[NumFactorsOfA - 1]--;
            }
        }
    }
}

BigInt GetA(QuadraticSieveData *Data, std::vector<uint32_t> &FactorsOfA)
{
    BigInt a = 1;
    for (int i = 0; i < Data->FactorOfABits.size() - 1;)
    {
        uint32_t Bits = Data->FactorOfABits[i];
        uint32_t Range = Data->FactorBaseBitBounds[Bits + 1] - Data->FactorBaseBitBounds[Bits];
        FactorsOfA[i] = Data->FactorBaseBitBounds[Bits] + std::rand() % Range;
        int j;
        for (j = 0; j < i; j++)
            if (FactorsOfA[i] == FactorsOfA[j])
                break;
        if (j == i) a *= Data->FactorBase[FactorsOfA[i++]];
    }
    int64_t TargetLastDivisor = (Data->OptimalValueOfA / a).Int();
    uint32_t Left = std::lower_bound(Data->FactorBase.begin(), Data->FactorBase.end(), TargetLastDivisor) - Data->FactorBase.begin();
    while (Data->FactorBase[Left] > TargetLastDivisor) Left--;
    uint32_t Right = Left + 1;
    while (true)
    {
        if (TargetLastDivisor - Data->FactorBase[Left] <= Data->FactorBase[Right] - TargetLastDivisor)
        {
            if (std::find(FactorsOfA.begin(), FactorsOfA.end() - 1, Left) == FactorsOfA.end() - 1)
            {
                FactorsOfA.back() = Left;
                a *= Data->FactorBase[Left];
                break;
            }
            else Left--;
        }
        else
        {
            if (std::find(FactorsOfA.begin(), FactorsOfA.end() - 1, Right) == FactorsOfA.end() - 1)
            {
                FactorsOfA.back() = Right;
                a *= Data->FactorBase[Right];
                break;
            }
            else Right++;
        }
    }
    std::sort(FactorsOfA.begin(), FactorsOfA.end());
    return a;
}

void TrialDivision(QuadraticSieveData *Data, const BigInt &a, const BigInt &b, const int64_t *RootsModP, int64_t i)
{
    uint32_t NumPrimesInFactorBase = Data->FactorBase.size();
    BigInt CongruentNum = a * i;
    CongruentNum += b;
    CongruentNum *= CongruentNum;
    BigInt CurrentNum = CongruentNum - Data->n;
    
    std::vector<uint32_t> Divisors;
    if (CurrentNum < 0ll) { Divisors.push_back(NumPrimesInFactorBase); CurrentNum = -CurrentNum; }
    
    for (int j = 0; j < NumPrimesInFactorBase; j++, RootsModP += 2)
    {
        if (RootsModP[0] != INT64_MAX)
        {
            int64_t imodp = i % Data->FactorBase[j];
            if (imodp < 0) imodp += Data->FactorBase[j];
            if (imodp != RootsModP[0] && imodp != RootsModP[1]) continue;
        }
        while (CurrentNum % Data->FactorBase[j] == 0ll)
        {
            CurrentNum /= Data->FactorBase[j];
            if (!Divisors.empty() && Divisors.back() == j) Divisors.pop_back();
            else Divisors.push_back(j);
        }
    }
    
    if (CurrentNum == 1)
    {
        Data->RelationMutex.lock();
        Data->Relations.emplace_back(std::make_pair(CongruentNum - Data->n, std::move(CongruentNum)), std::move(Divisors));
        Data->RelationMutex.unlock();
    }
    else
    {
        Data->PMutex.lock();
        bool FoundPair = Data->P.find(CurrentNum) != Data->P.end();
        if (FoundPair)
        {
            std::pair<BigInt, std::vector<uint32_t>> Pair = Data->P.at(CurrentNum);
            Data->PMutex.unlock();
            
            std::pair<std::pair<BigInt, BigInt>, std::vector<uint32_t>> MatrixElement;
            MatrixElement.first.first = (CongruentNum - Data->n) * (Pair.first - Data->n);
            MatrixElement.first.second = CongruentNum * Pair.first;
            
            int i = 0, j = 0;
            if (Divisors[i] == NumPrimesInFactorBase)
            {
                if (Pair.second[j] == NumPrimesInFactorBase) { i++; j++; }
                else MatrixElement.second.push_back(Divisors[i++]);
            }
            else if (Pair.second[j] == NumPrimesInFactorBase) MatrixElement.second.push_back(Pair.second[j++]);
            while (i < Divisors.size() || j < Pair.second.size())
            {
                if (i < Divisors.size() && j < Pair.second.size())
                {
                    if (Divisors[i] < Pair.second[j])
                        MatrixElement.second.push_back(Divisors[i++]);
                    else if (Divisors[i] > Pair.second[j])
                        MatrixElement.second.push_back(Pair.second[j++]);
                    else { i++; j++; }
                }
                else if (i < Divisors.size())
                    MatrixElement.second.push_back(Divisors[i++]);
                else MatrixElement.second.push_back(Pair.second[j++]);
            }
            
            if (!MatrixElement.second.empty())
            {
                Data->RelationMutex.lock();
                Data->Relations.emplace_back(std::move(MatrixElement));
                Data->RelationMutex.unlock();
            }
        }
        else
        {
            Data->P.emplace(CurrentNum, std::make_pair(std::move(CongruentNum), std::move(Divisors)));
            Data->PMutex.unlock();
        }
    }
}

void GetSievePointers(QuadraticSieveData *Data, uint8_t **SievePointers, const int64_t *RootsModP, uint8_t *Sieve)
{
    uint32_t NumPrimesInFactorBase = Data->FactorBase.size();
    int64_t *FactorBase = Data->FactorBase.data() + Data->SmallPrimeBound;
    RootsModP += Data->SmallPrimeBound * 2;
    SievePointers += Data->SmallPrimeBound * 2;
    for (int i = Data->SmallPrimeBound; i < NumPrimesInFactorBase; i++, SievePointers += 2, RootsModP += 2, FactorBase++)
    {
        if (RootsModP[0] == INT64_MAX) continue;
        int64_t CurrentPrime = *FactorBase;
        int64_t j = (Data->M + RootsModP[0]) % CurrentPrime;
        int64_t k = (Data->M + RootsModP[1]) % CurrentPrime;
        if (k == j) SievePointers[1] = nullptr;
        else SievePointers[1] = Sieve + k;
        *SievePointers = Sieve + j;
    }
}

void AddLogs(QuadraticSieveData *Data, uint8_t *Sieve, uint8_t **SievePointers)
{
    int64_t *FactorBase = Data->FactorBase.data();
    uint32_t NumPrimesInFactorBase = Data->FactorBase.size();
    uint8_t *FactorBaseLogs = Data->FactorBaseLogs.data();
    
    uint8_t const *SieveEnd = Sieve + 2 * Data->M;
    SievePointers += Data->SmallPrimeBound * 2;
    for (int i = Data->SmallPrimeBound; i < NumPrimesInFactorBase; i++, SievePointers += 2)
    {
        if (!SievePointers[0]) continue;
        int64_t CurrentPrime = FactorBase[i];
        uint8_t log = FactorBaseLogs[i];
        uint8_t *ptr = SievePointers[0];
        while (ptr < SieveEnd)
        {
            *ptr += log;
            ptr += CurrentPrime;
        }
        ptr = SievePointers[1];
        if (ptr)
        {
            while (ptr < SieveEnd)
            {
                *ptr += log;
                ptr += CurrentPrime;
            }
        }
    }
}

void Sieve(QuadraticSieveData *Data, const BigInt &a, const BigInt &b, const int64_t *RootsModP, uint8_t *Sieve, uint8_t **SievePointers)
{
    GetSievePointers(Data, SievePointers, RootsModP, Sieve);
    AddLogs(Data, Sieve, SievePointers);
    
    uint8_t *CurrentSievePtr = Sieve;
    for (int64_t i = -Data->M; i < Data->M; i++, CurrentSievePtr++)
    {
        if (*CurrentSievePtr < Data->SmallThreshold) continue;
        const int64_t *CurrentRoots = RootsModP;
        for (int j = 0; j < Data->SmallPrimeBound; j++, CurrentRoots += 2)
        {
            if (CurrentRoots[0] == INT64_MAX) continue;
            int64_t imodp = i % Data->FactorBase[j];
            if (imodp < 0) imodp += Data->FactorBase[j];
            if (imodp == CurrentRoots[0] || imodp == CurrentRoots[1]) *CurrentSievePtr += Data->FactorBaseLogs[j];
        }
        if (*CurrentSievePtr < Data->LargeThreshold) continue;
        TrialDivision(Data, a, b, RootsModP, i);
    }
}

void FindRelations(QuadraticSieveData *Data, uint32_t ThreadIndex)
{
    while (Data->Running)
    {
        Data->AMutex.lock();
        std::cout << "\rFound " << Data->Relations.size() << " relations" << std::flush;
        do Data->a[ThreadIndex] = GetA(Data, Data->CurrentFactorIndices[ThreadIndex]);
        while (Data->UsedA.find(Data->a[ThreadIndex]) != Data->UsedA.end());
        Data->AMutex.unlock();
        
        uint32_t NumPrimesInFactorBase = Data->FactorBase.size();
        uint32_t NumFactorsOfA = Data->FactorOfABits.size();
        std::vector<uint32_t> &CurrentFactorIndices = Data->CurrentFactorIndices[ThreadIndex];
        std::vector<int64_t> &InverseOfAModP = Data->InverseOfAModP[ThreadIndex];
        std::vector<BigInt> &B = Data->B[ThreadIndex];
        std::vector<int64_t> &BTimesInverseA = Data->BTimesInverseA[ThreadIndex];
        std::vector<int64_t> &RootsModP = Data->RootsModP[ThreadIndex];
        BigInt &a = Data->a[ThreadIndex];
        
        for (int i = 0; i < NumPrimesInFactorBase; i++)
            InverseOfAModP[i] = BinPow((a % Data->FactorBase[i]).Int(), Data->FactorBase[i] - 2, Data->FactorBase[i]);
        
        for (int i = 0; i < NumFactorsOfA; i++)
        {
            int64_t Factor = Data->FactorBase[CurrentFactorIndices[i]];
            BigInt AOverFactor = a / Factor;
            B[i] = Data->SqrtNModP[CurrentFactorIndices[i]] * AOverFactor;
            AOverFactor %= Factor;
            B[i] *= BinPow(AOverFactor.Int(), Factor - 2, Factor);
            B[i] %= a;
        }
        
        BigInt b = static_cast<int64_t>(0);
        for (int i = 0; i < NumFactorsOfA; i++)
            b += B[i];
        
        int64_t *ptr = BTimesInverseA.data();
        for (int i = 0; i < NumFactorsOfA; i++)
        {
            for (int j = 0; j < NumPrimesInFactorBase; j++, ptr++)
            {
                BigInt Current = B[i] * InverseOfAModP[j];
                Current %= Data->FactorBase[j];
                *ptr = Current.Int();
            }
        }
        
        for (int i = 0; i < NumFactorsOfA; i++)
            B[i] *= 2ll;
        
        std::memset(Data->Sieves[ThreadIndex].data(), 0, Data->M * 2);
        uint32_t Mask = 0;
        int64_t *CurrentRoots = RootsModP.data();
        uint32_t *NextFactorOfA = CurrentFactorIndices.data();
        uint32_t *FactorsOfAEnd = NextFactorOfA + NumFactorsOfA;
        for (uint32_t i = 0; i < NumPrimesInFactorBase; i++, CurrentRoots += 2)
        {
            if (NextFactorOfA != FactorsOfAEnd && *NextFactorOfA == i)
            {
                CurrentRoots[0] = INT64_MAX;
                NextFactorOfA++;
                continue;
            }
            
            CurrentRoots[0] = Data->SqrtNModP[i] * InverseOfAModP[i];
            CurrentRoots[0] %= Data->FactorBase[i];
            for (int j = 0; j < NumFactorsOfA; j++)
            {
                CurrentRoots[0] += BTimesInverseA[j * NumPrimesInFactorBase + i];
                CurrentRoots[0] %= Data->FactorBase[i];
            }
            CurrentRoots[0] %= Data->FactorBase[i];
            CurrentRoots[0] = Data->FactorBase[i] - CurrentRoots[0];
            CurrentRoots[0] %= Data->FactorBase[i];
            
            CurrentRoots[1] = Data->SqrtNModP[i] * InverseOfAModP[i];
            CurrentRoots[1] *= 2;
            CurrentRoots[1] += CurrentRoots[0];
            CurrentRoots[1] %= Data->FactorBase[i];
        }
        Sieve(Data, a, b, RootsModP.data(), Data->Sieves[ThreadIndex].data(), Data->SievePointers[ThreadIndex].data());
        
        for (int i = 1; i < (1 << (NumFactorsOfA - 1)); i++)
        {
            uint32_t Index = 0;
            while (((1 << Index) & i) == 0) Index++;
            bool Add = Mask & (1 << Index);
            Mask ^= (1 << Index);
            
            if (!Add) b -= B[Index];
            else b += B[Index];
            
            std::memset(Data->Sieves[ThreadIndex].data(), 0, Data->M * 2);
            CurrentRoots = RootsModP.data();
            int64_t *CurrentBTimesInverseOfA = BTimesInverseA.data() + Index * NumPrimesInFactorBase;
            for (uint32_t j = 0; j < NumPrimesInFactorBase; j++, CurrentRoots += 2, CurrentBTimesInverseOfA++)
            {
                if (CurrentRoots[0] == INT64_MAX) continue;
                
                CurrentRoots[0] += Data->FactorBase[j] * 2ll;
                CurrentRoots[1] += Data->FactorBase[j] * 2ll;
                if (!Add)
                {
                    CurrentRoots[0] += *CurrentBTimesInverseOfA * 2ll;
                    CurrentRoots[0] %= Data->FactorBase[j];
                    CurrentRoots[1] += *CurrentBTimesInverseOfA * 2ll;
                    CurrentRoots[1] %= Data->FactorBase[j];
                }
                else
                {
                    CurrentRoots[0] -= *CurrentBTimesInverseOfA * 2ll;
                    CurrentRoots[0] %= Data->FactorBase[j];
                    CurrentRoots[1] -= *CurrentBTimesInverseOfA * 2ll;
                    CurrentRoots[1] %= Data->FactorBase[j];
                }
            }
            Sieve(Data, a, b, RootsModP.data(), Data->Sieves[ThreadIndex].data(), Data->SievePointers[ThreadIndex].data());
        }
        
        if (Data->Running && Data->Relations.size() >= Data->FactorBase.size())
        {
            Data->RelationMutex.lock();
            std::sort(Data->Relations.begin(), Data->Relations.end());
            Data->Relations.resize(std::unique(Data->Relations.begin(), Data->Relations.end()) - Data->Relations.begin());
            if (Data->Relations.size() >= Data->FactorBase.size() + 100)
                Data->Running = false;
            Data->RelationMutex.unlock();
        }
    }
}

std::pair<BigInt, BigInt> LinearAlgebra(const BigInt &n, QuadraticSieveData *Data, uint32_t NumThreads, std::vector<std::thread> &Threads)
{
    Data->FactorBase.push_back(-1);
    std::vector<uint32_t> NonZeroEntries;
    for (int i = 0; i < Data->Relations.size(); i++)
    {
        for (int j = 0; j < Data->Relations[i].second.size(); j++)
        {
            NonZeroEntries.push_back(Data->Relations[i].second[j]);
            NonZeroEntries.push_back(i);
        }
    }
    
    std::vector<uint64_t> NullSpaces(Data->Relations.size());
    uint32_t NumNullSpaces = blanczos(NonZeroEntries.data(), NonZeroEntries.size() / 2, Data->FactorBase.size(), Data->Relations.size(), NullSpaces.data());
    bool Finished = false;
    std::pair<BigInt, BigInt> Result = std::make_pair(0, 0);
    std::mutex Mutex;
    auto ThreadLambda = [&](uint32_t Index) {
        while (Index < NumNullSpaces && !Finished)
        {
            BigInt A = 1, B = 1;
            for (int j = 0; j < Data->Relations.size(); j++)
            {
                if (Finished) return;
                if ((1ull << Index) & NullSpaces[j])
                {
                    A *= Data->Relations[j].first.first;
                    B *= Data->Relations[j].first.second;
                }
            }
            if (Finished) return;
            A = A.SquareRoot(); B = B.SquareRoot();
            BigInt gcd1 = gcd(n, (A - B).abs()), gcd2 = gcd(n, (A + B).abs());
            if (gcd1 != 1 && gcd1 != n)
            {
                Mutex.lock();
                if (Finished)
                {
                    Mutex.unlock();
                    return;
                }
                Result = std::make_pair(gcd1, n / gcd1);
                Finished = true;
                Mutex.unlock();
            }
            if (gcd2 != 1 && gcd2 != n)
            {
                Mutex.lock();
                if (Finished)
                {
                    Mutex.unlock();
                    return;
                }
                Result = std::make_pair(gcd2, n / gcd2);
                Finished = true;
                Mutex.unlock();
            }
            Index += NumThreads;
        }
    };
    
    for (uint32_t i = 0; i < NumThreads - 1; i++)
        Threads[i] = std::thread(ThreadLambda, i);
    ThreadLambda(NumThreads - 1);
    
    for (int i = 0; i < NumThreads - 1; i++)
        Threads[i].join();
    
    return Result;
}

void Init(QuadraticSieveData *Data, const BigInt &n, uint32_t NumThreads)
{
    Data->n = ChooseMultiplier(n);
    Data->PrimeBound = ChoosePrimeBound(Data->n);
    std::cout << "Factor base upper bound: " << Data->PrimeBound << "\n";
    CalculateFactorBase(Data);
    std::cout << "Factor base size: " << Data->FactorBase.size() << "\n";
    Data->M = ChooseSieveSize(n);
    std::cout << "M: " << Data->M << "\n";
    Data->Delta = ChooseDelta(n);
    std::cout << "Delta: " << Data->Delta << "\n";
    Data->SmallPrimeBound = ChooseSmallPrimeBound(n);
    std::cout << "Small prime bound: " << Data->SmallPrimeBound << std::endl;
    
    Data->SqrtNModP.resize(Data->FactorBase.size());
    for (int i = 0; i < Data->FactorBase.size(); i++)
        Data->SqrtNModP[i] = SquareRootModP(Data->n, Data->FactorBase[i]).Int();
    Data->OptimalValueOfA = ((Data->n * 2ll) / (Data->M * Data->M)).SquareRoot();
    GetFactorOfASizes(Data);
    Data->LargeThreshold = ((Data->n * Data->M * Data->M) / 2ll).SquareRoot().SizeInBase(2) - Data->Delta;
    Data->SmallThreshold = Data->LargeThreshold - Data->Delta;

    Data->CurrentFactorIndices.resize(NumThreads);
    Data->InverseOfAModP.resize(NumThreads);
    Data->B.resize(NumThreads);
    Data->BTimesInverseA.resize(NumThreads);
    Data->RootsModP.resize(NumThreads);
    Data->Sieves.resize(NumThreads);
    Data->SievePointers.resize(NumThreads);
    for (int i = 0; i < NumThreads; i++)
    {
        Data->CurrentFactorIndices[i] = std::vector<uint32_t>(Data->FactorOfABits.size());
        Data->InverseOfAModP[i] = std::vector<int64_t>(Data->FactorBase.size());
        Data->B[i] = std::vector<BigInt>(Data->FactorOfABits.size());
        Data->BTimesInverseA[i] = std::vector<int64_t>(Data->FactorOfABits.size() * Data->FactorBase.size());
        Data->RootsModP[i] = std::vector<int64_t>(Data->FactorBase.size() * 2);
        Data->Sieves[i] = std::vector<uint8_t>(Data->M * 2);
        Data->SievePointers[i] = std::vector<uint8_t *>(Data->FactorBase.size() * 2);
    }
    Data->a.resize(NumThreads);
}

std::pair<BigInt, BigInt> QuadraticSieve(const BigInt &n, uint32_t NumThreads = 1)
{
    std::cout << "Size of n is " << n.SizeInBase(2) << " bits\n";
    
    std::unique_ptr<QuadraticSieveData> Data = std::make_unique<QuadraticSieveData>();
    Init(Data.get(), n, NumThreads);
    
    Data->Running = true;
    std::vector<std::thread> Threads;
    Threads.reserve(NumThreads - 1);
    for (uint32_t i = 0; i < NumThreads - 1; i++)
        Threads.emplace_back(FindRelations, Data.get(), i);
    FindRelations(Data.get(), NumThreads - 1);
    for (int i = 0; i < Threads.size(); i++)
        Threads[i].join();
    
    return LinearAlgebra(n, Data.get(), NumThreads, Threads);
}
