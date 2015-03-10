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

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const RandomizingConstVector&);
};

template <typename T>
class RandomizingVector {
        boost::random::uniform_int_distribution<size_t> uniform;
        std::vector<T> myVector;

        public:
        RandomizingVector(const std::vector<T>& aVector):
                uniform(0, aVector.size() - 1), myVector(aVector) {}

        const T& getRandomElementByReference(std::mt19937_64& rng) {
                return myVector[uniform(rng)];
        }

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const RandomizingVector&);
};



template <typename T>
std::ostream& operator<< (std::ostream& out, const RandomizingConstVector<T>& aCV) {
	out << aCV.myVector;
	return out;
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const RandomizingVector<T>& aV) {
	out << aV.myVector;
	return out;
}