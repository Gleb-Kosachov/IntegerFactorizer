//
//  BigInt.cpp
//  FactorizationOfLargeSemiprimes
//
//  Created by Глеб Косачёв on 26/12/2024.
//

#include "BigInt.hpp"
#include <cmath>
#include <cstring>
#include <iomanip>
#include <complex>
#include <exception>
#include <mutex>

BigInt::BigInt()
{
    mpz_init(n);
    mpz_set_ui(n, 0);
}

BigInt::~BigInt()
{
    mpz_clear(n);
}

BigInt::BigInt(int64_t Num)
{
    mpz_init(n);
    mpz_set_ui(n, std::abs(Num));
    if (Num < 0) mpz_neg(n, n);
}

BigInt::BigInt(const char *Num)
{
    mpz_init(n);
    mpz_set_str(n, Num, 10);
}

BigInt::BigInt(const BigInt &other) noexcept
{
    mpz_init(n);
    mpz_set(n, other.n);
}

BigInt::BigInt(BigInt &&other) noexcept
{
    mpz_init(n);
    mpz_set_ui(n, 0);
    std::swap(n, other.n);
}

bool BigInt::operator>(const BigInt &other) const
{
    return mpz_cmp(n, other.n) > 0;
}

bool BigInt::operator<(const BigInt &other) const
{
    return mpz_cmp(n, other.n) < 0;
}

bool BigInt::operator==(const BigInt &other) const
{
    return mpz_cmp(n, other.n) == 0;
}

bool BigInt::operator!=(const BigInt &other) const
{
    return !this->operator==(other);
}

bool BigInt::operator>=(const BigInt &other) const
{
    return this->operator>(other) || this->operator==(other);
}

bool BigInt::operator<=(const BigInt &other) const
{
    return this->operator<(other) || this->operator==(other);
}

bool BigInt::operator>(int64_t other) const
{
    return mpz_cmp_si(n, other) > 0;
}

bool BigInt::operator<(int64_t other) const
{
    return mpz_cmp_si(n, other) < 0;
}

bool BigInt::operator==(int64_t other) const
{
    return mpz_cmp_si(n, other) == 0;
}

bool BigInt::operator!=(int64_t other) const
{
    return !this->operator==(other);
}

bool BigInt::operator>=(int64_t other) const
{
    return this->operator>(other) || this->operator==(other);
}

bool BigInt::operator<=(int64_t other) const
{
    return this->operator<(other) || this->operator==(other);
}

BigInt operator+(const BigInt &a, const BigInt &b)
{
    BigInt Result(a);
    Result += b;
    return Result;
}

BigInt operator-(const BigInt &a, const BigInt &b)
{
    BigInt Result(a);
    Result -= b;
    return Result;
}

BigInt operator*(const BigInt &a, const BigInt &b)
{
    BigInt Result(a);
    Result *= b;
    return Result;
}

BigInt operator/(const BigInt &a, const BigInt &b)
{
    BigInt Result(a);
    Result /= b;
    return Result;
}

BigInt operator%(const BigInt &a, const BigInt &b)
{
    BigInt Result(a);
    Result %= b;
    return Result;
}

BigInt operator+(int64_t a, const BigInt &b)
{
    if (a < 0) return b - (-a);
    BigInt Result;
    mpz_add_ui(Result.n, b.n, a);
    return Result;
}

BigInt operator-(int64_t a, const BigInt &b)
{
    BigInt Result = b - a;
    mpz_neg(Result.n, Result.n);
    return Result;
}

BigInt operator*(int64_t a, const BigInt &b)
{
    BigInt Result;
    mpz_mul_ui(Result.n, b.n, std::abs(a));
    if (a < 0) mpz_neg(Result.n, Result.n);
    return Result;
}

BigInt operator+(const BigInt &a, int64_t b)
{
    if (b < 0) return a - (-b);
    BigInt Result;
    mpz_add_ui(Result.n, a.n, b);
    return Result;
}

BigInt operator-(const BigInt &a, int64_t b)
{
    if (b < 0) return a + (-b);
    BigInt Result(a);
    mpz_sub_ui(Result.n, a.n, b);
    return Result;
}

bool Debug = false;

BigInt operator*(const BigInt &a, int64_t b)
{
    BigInt Result;
    mpz_mul_ui(Result.n, a.n, std::abs(b));
    if (b < 0) mpz_neg(Result.n, Result.n);
    return Result;
}

BigInt operator/(const BigInt &a, int64_t b)
{
    BigInt Result;
    mpz_div_ui(Result.n, a.n, std::abs(b));
    if (b < 0) mpz_neg(Result.n, Result.n);
    return Result;
}

BigInt operator%(const BigInt &a, int64_t b)
{
    BigInt Result;
    mpz_mod_ui(Result.n, a.n, b);
    return Result;
}

void BigInt::operator+=(const BigInt &other)
{
    mpz_add(n, n, other.n);
}

void BigInt::operator-=(const BigInt &other)
{
    mpz_sub(n, n, other.n);
}

void BigInt::operator*=(const BigInt &other)
{
    mpz_mul(n, n, other.n);
}

void BigInt::operator/=(const BigInt &other)
{
    mpz_div(n, n, other.n);
}

void BigInt::operator%=(const BigInt &other)
{
    mpz_mod(n, n, other.n);
}

void BigInt::operator+=(int64_t other)
{
    if (other < 0) this->operator-=(-other);
    else mpz_add_ui(n, n, other);
}

void BigInt::operator-=(int64_t other)
{
    if (other < 0) this->operator+=(-other);
    else mpz_sub_ui(n, n, other);
}

void BigInt::operator*=(int64_t other)
{
    mpz_mul_ui(n, n, std::abs(other));
    if (other < 0) mpz_neg(n, n);
}

void BigInt::operator/=(int64_t other)
{
    mpz_div_ui(n, n, std::abs(other));
    if (other < 0) mpz_neg(n, n);
}

void BigInt::operator%=(int64_t other)
{
    mpz_mod_ui(n, n, other);
}

BigInt BigInt::operator++(int)
{
    BigInt ReturnValue = *this;
    this->operator+=(1);
    return ReturnValue;
}

BigInt BigInt::operator=(BigInt &&other) noexcept
{
    std::swap(n, other.n);
    return *this;
}

BigInt BigInt::operator=(const BigInt &other) noexcept
{
    mpz_set(n, other.n);
    return *this;
}

BigInt BigInt::operator-() const
{
    BigInt Result;
    mpz_neg(Result.n, n);
    return Result;
}

std::ostream &operator<<(std::ostream &out, const BigInt &Num)
{
    out << Num.n;
    return out;
}

BigInt gcd(const BigInt &a, const BigInt &b)
{
    BigInt Result;
    mpz_gcd(Result.n, a.n, b.n);
    return Result;
}

BigInt BigInt::SquareRoot() const
{
    BigInt Result;
    mpz_sqrt(Result.n, n);
    return Result;
}

BigInt BigInt::Root(uint32_t n) const
{
    BigInt Result;
    mpz_root(Result.n, this->n, n);
    return Result;
}

BigInt BigInt::abs() const
{
    BigInt Result;
    mpz_abs(Result.n, n);
    return Result;
}

uint32_t BigInt::SizeInBase(uint32_t Base) const
{
    return mpz_sizeinbase(n, Base);
}

std::string BigInt::String(uint32_t Base) const
{
    char *str = mpz_get_str(nullptr, Base, n);
    std::string Result = str;
    void (*freefunc)(void *, size_t);
    mp_get_memory_functions (NULL, NULL, &freefunc);
    freefunc(str, std::strlen(str) + 1);
    return Result;
}

int64_t BigInt::Int() const
{
    return mpz_get_ui(n);
}

BigInt BigInt::Random(const BigInt &Max)
{
    static gmp_randstate_t state;
    static bool Initialized = false;
    if (!Initialized)
    {
        Initialized = true;
        gmp_randinit_default(state);
    }
    BigInt Result;
    mpz_urandomm(Result.n, state, Max.n);
    return Result;
}
