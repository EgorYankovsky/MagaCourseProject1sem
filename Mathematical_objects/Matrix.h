#pragma once

class Matrix {
protected:
    virtual double operator() (size_t i, size_t j) const = 0;
    virtual double& operator() (size_t i, size_t j) = 0;
    virtual size_t getSize() const = 0;
};