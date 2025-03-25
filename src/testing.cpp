//
//  main.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 26/12/2024.
//

#include <iostream>
#include <algorithm>
#include <chrono>
#include <set>
#include <ctime>
#include <fstream>
#include "Algorithms.hpp"
#include "BigInt.hpp"

int main()
{
    std::srand(std::time(NULL));
    try
    {
        const std::string sizes[7] = { "3", "5", "10", "15", "20", "30", "40" };
        std::ofstream out;
        std::ifstream in;
        
        out.open("src/Output/TrialDivision.txt");
        for (int i = 0; i < 7; i++)
        {
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            uint64_t Time = 0;
            bool LessThanTenMinutes = true;
            for (int i = 0; i < 500; i++)
            {
                std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
                BruteForce(Numbers[std::rand() % 500].data());
                uint64_t TestTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin).count();
                if (TestTime > 600000)
                {
                    LessThanTenMinutes = false;
                    break;
                }
                Time += TestTime;
            }
            if (!LessThanTenMinutes) break;
            Time /= 500;
            out << sizes[i] << ": " << Time << "\n";
        }
        out.close();
        
        out.open("src/Output/Fermat.txt");
        for (int i = 0; i < 7; i++)
        {
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            uint64_t Time = 0;
            bool LessThanTenMinutes = true;
            for (int i = 0; i < 500; i++)
            {
                std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
                Fermat(Numbers[std::rand() % 500].data());
                uint64_t TestTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin).count();
                if (TestTime > 600000)
                {
                    LessThanTenMinutes = false;
                    break;
                }
                Time += TestTime;
            }
            if (!LessThanTenMinutes) break;
            Time /= 500;
            out << sizes[i] << ": " << Time << "\n";
        }
        out.close();
        
        out.open("src/Output/PollardRho.txt");
        for (int i = 0; i < 7; i++)
        {
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            uint64_t Time = 0;
            bool LessThanTenMinutes = true;
            for (int i = 0; i < 500; i++)
            {
                std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
                PollardRho(Numbers[std::rand() % 500].data());
                uint64_t TestTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin).count();
                if (TestTime > 600000)
                {
                    LessThanTenMinutes = false;
                    break;
                }
                Time += TestTime;
            }
            if (!LessThanTenMinutes) break;
            Time /= 500;
            out << sizes[i] << ": " << Time << "\n";
        }
        out.close();
        
        out.open("src/Output/QuadraticSieve.txt");
        for (int i = 3; i < 7; i++)
        {
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            uint64_t Time = 0;
            bool LessThanTenMinutes = true;
            for (int i = 0; i < 500; i++)
            {
                std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
                QuadraticSieve(Numbers[std::rand() % 500].data());
                uint64_t TestTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin).count();
                if (TestTime > 600000)
                {
                    LessThanTenMinutes = false;
                    break;
                }
                Time += TestTime;
            }
            if (!LessThanTenMinutes) break;
            Time /= 500;
            out << sizes[i] << ": " << Time << "\n";
        }
        out.close();
    }
    catch (std::exception *a)
    {
        std::cout << a->what() << "\n";
        std::abort();
    }
}
