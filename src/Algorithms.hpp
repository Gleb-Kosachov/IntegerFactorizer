//
//  BasicAlgorithms.hpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 16/01/2025.
//

#ifndef BasicAlgorithms_hpp
#define BasicAlgorithms_hpp

#include "BigInt.hpp"

std::pair<BigInt, BigInt> BruteForce(const BigInt &);
std::pair<BigInt, BigInt> Fermat(const BigInt &);
std::pair<BigInt, BigInt> PollardRho(const BigInt &, uint32_t = 1);
std::pair<BigInt, BigInt> QuadraticSieve(const BigInt &, uint32_t = 1);

#endif /* BasicAlgorithms_hpp */
