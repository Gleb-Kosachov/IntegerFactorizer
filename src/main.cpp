//
//  main.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 26/12/2024.
//

#include <iostream>
#include <algorithm>
#include <set>
#include <ctime>
#include <fstream>
#include "Algorithms.hpp"
#include "BigInt.hpp"

int main()
{
//    std::srand(std::time(NULL));
    try
    {
//        std::string a, b; std::cin >> a >> b;
//        BigInt n = BigInt(a.data()) * BigInt(b.data());
//        std::string str = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
        std::string str = "24305374772286341905584488566453089543418508609421058188889925765889633292380493";
//        std::string str = "497806691290494324233539625829262221580368990850811330950699";
//        std::string str; std::cin >> str;
        BigInt n(str.data());
        std::pair<BigInt, BigInt> Factorization;
        std::chrono::high_resolution_clock::time_point Begin = std::chrono::high_resolution_clock::now();
//        for (int i = 0; i < 10; i++)
            Factorization = QuadraticSieve(n, 4);
        std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - Begin).count() << "\n";
        std::cout << Factorization.first << "\n" << Factorization.second << "\n";
        if (Factorization.first * Factorization.second != n)
            std::cout << "Something went wrong\n";
    }
    catch (std::exception *a)
    {
        std::cout << a->what() << "\n";
        std::abort();
    }
}
