#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include <random>
#include <vector>

template <typename T>
class RandomizingConstVector {
        const boost::random::uniform_int_distribution<size_t> uniform;
        const std::vector<T> myVector;

        public:
        RandomizingConstVector(const std::vector<T>& aVector):
                uniform(0, aVector.size() - 1), myVector(aVector) {}

        const T& getRandomElementByReference(std::mt19937_64& rng) const {
                return myVector[uniform(rng)];
        }
};