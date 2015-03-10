#pragma once

#include <boost/random/uniform_int_distribution.hpp>
#include <iterator> // std::back_inserter
#include <random>
#include <vector>

template <typename T>
class RandomizingConstVector {
	const boost::random::uniform_int_distribution<size_t> uniform;
	const std::vector<T> myVector;

	public:
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef typename std::vector<T>::size_type size_type;
	typedef typename std::vector<T>::const_iterator const_iterator;


	RandomizingConstVector(const std::vector<T>& aVector):
		uniform(0, aVector.size() - 1), myVector(aVector) {}

	const T& getRandomElementByReference(std::mt19937_64& rng) const {
		return myVector[uniform(rng)];
	}

	size_type size() const {
		return myVector.size();
	}
	const_iterator begin() const {
		return myVector.begin();
	}
	const_iterator end() const {
		return myVector.end();
	}
	const_reference operator[](size_type pos) const {
		return myVector[pos];
	};

	/** Output */
	friend std::ostream& operator<< (std::ostream&, const RandomizingConstVector&);
};

template <typename T>
class RandomizingVector {
	boost::random::uniform_int_distribution<size_t> uniform;
	std::vector<T> myVector;

	public:
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef typename std::vector<T>::size_type size_type;
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;
	RandomizingVector(const std::vector<T>& aVector):
		uniform(0, aVector.size() - 1), myVector(aVector) {}

	const T& getRandomElementByReference(std::mt19937_64& rng) {
		return myVector[uniform(rng)];
	}

	size_type size() const {
		return myVector.size();
	}
	iterator begin() {
		return myVector.begin();
	}
	iterator end() {
		return myVector.end();
	}
	const_iterator begin() const {
		return myVector.begin();
	}
	const_iterator end() const {
		return myVector.end();
	}
	reference operator[](size_type pos) {
		return myVector[pos];
	}
	const_reference operator[](size_type pos) const {
		return myVector[pos];
	};

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