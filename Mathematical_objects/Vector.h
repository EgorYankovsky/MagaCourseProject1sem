#pragma once
class Vector {
    virtual double operator() (size_t i) const = 0;
    virtual double& operator() (size_t i) = 0;
    virtual size_t getSize() const = 0;
};

