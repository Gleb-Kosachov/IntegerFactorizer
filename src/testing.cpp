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
        
        for (int i = 0; i < 7; i++)
        {
            out.open("src/Output/TrialDivision.txt", std::ios::app);
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 500; i++)
                BruteForce(Numbers[std::rand() % 500].data());
            uint64_t Time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 500000000; 
            out << sizes[i] << ": " << Time << "\n";
            out.close();
	    if (Time >= 60000) break;
        }
        
        for (int i = 0; i < 7; i++)
        {
            out.open("src/Output/Fermat.txt", std::ios::app);
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 500; i++)
                Fermat(Numbers[std::rand() % 500].data());
            uint64_t Time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 500000000; 
            out << sizes[i] << ": " << Time << "\n";
            out.close();
	    if (Time >= 60000) break;
        }

        for (int i = 0; i < 7; i++)
        {
            out.open("src/Output/PollardRho.txt", std::ios::app);
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 500; i++)
                PollardRho(Numbers[std::rand() % 500].data());
            uint64_t Time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 500000000; 
            out << sizes[i] << ": " << Time << "\n";
            out.close();
	    if (Time >= 60000) break;
        }

        for (int i = 3; i < 7; i++)
        {
            out.open("src/Output/QuadraticSieve.txt", std::ios::app);
            in.open("src/Semiprimes/" + sizes[i] + ".txt");
            if (!in.is_open()) throw std::runtime_error("File not found!");
            std::string Numbers[500];
            for (int i = 0; i < 500; i++) in >> Numbers[i];
            in.close();
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            for (int i = 0; i < 500; i++)
                QuadraticSieve(Numbers[std::rand() % 500].data());
            uint64_t Time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 500000000; 
            out << sizes[i] << ": " << Time << "\n";
            out.close();
	    if (Time >= 60000) break;
        }
    }
    catch (std::exception *a)
    {
        std::cout << a->what() << "\n";
        std::abort();
    }
}
