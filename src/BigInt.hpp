//
//  BigInt.hpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 26/12/2024.
//

#ifndef BigInt_hpp
#define BigInt_hpp

#include <iostream>
#include <string>
#include <gmp.h>

class BigInt
{
public:
    BigInt();
    ~BigInt();
    BigInt(int64_t);
    BigInt(const char *);
    BigInt(const BigInt &) noexcept;
    BigInt(BigInt &&) noexcept;
    
    friend BigInt operator+(const BigInt &, const BigInt &);
    friend BigInt operator-(const BigInt &, const BigInt &);
    friend BigInt operator*(const BigInt &, const BigInt &);
    friend BigInt operator/(const BigInt &, const BigInt &);
    friend BigInt operator%(const BigInt &, const BigInt &);
    
    friend BigInt operator+(int64_t, const BigInt &);
    friend BigInt operator-(int64_t, const BigInt &);
    friend BigInt operator*(int64_t, const BigInt &);
    
    friend BigInt operator+(const BigInt &, int64_t);
    friend BigInt operator-(const BigInt &, int64_t);
    friend BigInt operator*(const BigInt &, int64_t);
    friend BigInt operator/(const BigInt &, int64_t);
    friend BigInt operator%(const BigInt &, int64_t);
    
    void operator+=(const BigInt &);
    void operator-=(const BigInt &);
    void operator*=(const BigInt &);
    void operator/=(const BigInt &);
    void operator%=(const BigInt &);
    void operator+=(int64_t);
    void operator-=(int64_t);
    void operator*=(int64_t);
    void operator/=(int64_t);
    void operator%=(int64_t);
    BigInt operator++(int);
    
    bool operator>(const BigInt &) const;
    bool operator<(const BigInt &) const;
    bool operator==(const BigInt &) const;
    bool operator!=(const BigInt &) const;
    bool operator>=(const BigInt &) const;
    bool operator<=(const BigInt &) const;
    
    bool operator>(int64_t) const;
    bool operator<(int64_t) const;
    bool operator==(int64_t) const;
    bool operator!=(int64_t) const;
    bool operator>=(int64_t) const;
    bool operator<=(int64_t) const;
    
    BigInt operator=(BigInt &&) noexcept;
    BigInt operator=(const BigInt &) noexcept;
    BigInt operator-() const;
    
    friend std::ostream &operator<<(std::ostream &, const BigInt &);
    friend BigInt gcd(const BigInt &a, const BigInt &b);
    
    BigInt SquareRoot() const;
    BigInt Root(uint32_t) const;
    BigInt abs() const;

    uint32_t SizeInBase(uint32_t) const;
    std::string String(uint32_t) const;
    int64_t Int() const;
    
    static BigInt Random(const BigInt &);
    
private:
    mpz_t n;
};

BigInt operator+(const BigInt &, const BigInt &);
BigInt operator-(const BigInt &, const BigInt &);
BigInt operator*(const BigInt &, const BigInt &);
BigInt operator/(const BigInt &, const BigInt &);
BigInt operator%(const BigInt &, const BigInt &);

BigInt operator+(int64_t, const BigInt &);
BigInt operator-(int64_t, const BigInt &);
BigInt operator*(int64_t, const BigInt &);

BigInt operator+(const BigInt &, int64_t);
BigInt operator-(const BigInt &, int64_t);
BigInt operator*(const BigInt &, int64_t);
BigInt operator/(const BigInt &, int64_t);
BigInt operator%(const BigInt &, int64_t);

std::ostream &operator<<(std::ostream &, const BigInt &);
BigInt gcd(const BigInt &a, const BigInt &b);

#endif /* BigInt_hpp */
