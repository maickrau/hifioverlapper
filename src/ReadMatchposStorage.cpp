#include <limits>
#include <cassert>
#include <cstring>
#include "ReadMatchposStorage.h"

ReadMatchposStorage::ReadMatchposStorageIterator::ReadMatchposStorageIterator(const ReadMatchposStorage& storage, size_t index) :
	storage(storage),
	index(index)
{
}

ReadMatchposStorage::ReadMatchposStorageIterator& ReadMatchposStorage::ReadMatchposStorageIterator::operator++()
{
	index += 1;
	return *this;
}

bool ReadMatchposStorage::ReadMatchposStorageIterator::operator==(const ReadMatchposStorageIterator& other) const
{
	return (&storage == &other.storage && index == other.index);
}

bool ReadMatchposStorage::ReadMatchposStorageIterator::operator!=(const ReadMatchposStorageIterator& other) const
{
	return !(*this == other);
}

std::tuple<uint32_t, uint32_t, uint32_t> ReadMatchposStorage::ReadMatchposStorageIterator::operator*() const
{
	return storage.getValue(index);
}

ReadMatchposStorage::ReadMatchposStorage() :
	values(nullptr),
	realsize(0),
	capacity(0)
{
}

ReadMatchposStorage::~ReadMatchposStorage()
{
	if (values != nullptr) delete [] values;
}

ReadMatchposStorage::ReadMatchposStorage(ReadMatchposStorage&& other) :
	values(nullptr),
	realsize(0),
	capacity(0)
{
	*this = std::move(other);
}

ReadMatchposStorage& ReadMatchposStorage::operator=(ReadMatchposStorage&& other)
{
	if (values != nullptr) delete [] values;
	values = other.values;
	realsize = other.realsize;
	capacity = other.capacity;
	other.values = nullptr;
	other.realsize = 0;
	other.capacity = 0;
	return *this;
}

std::tuple<uint32_t, uint32_t, uint32_t> ReadMatchposStorage::getValue(size_t index) const
{
	std::tuple<uint32_t, uint32_t, uint32_t> result;
	size_t startIndex = index*10;
	std::get<0>(result) = ((uint32_t)values[startIndex] << (uint32_t)24) + ((uint32_t)values[startIndex+1] << (uint32_t)16) + ((uint32_t)values[startIndex+2] << (uint32_t)8) + (uint32_t)values[startIndex+3];
	std::get<1>(result) = ((uint32_t)(values[startIndex+4] & 0x7F) << (uint32_t)16) + ((uint32_t)values[startIndex+5] << (uint32_t)8) + (uint32_t)values[startIndex+6];
	std::get<2>(result) = ((uint32_t)values[startIndex+7] << (uint32_t)16) + ((uint32_t)values[startIndex+8] << (uint32_t)8) + (uint32_t)values[startIndex+9];
	assert(std::get<2>(result) > std::get<1>(result));
	std::get<1>(result) += ((uint32_t)values[startIndex+4] & 0x80) << 24; // reverse orientation bit
	std::get<2>(result) += ((uint32_t)values[startIndex+4] & 0x80) << 24; // reverse orientation bit
	assert(std::get<2>(result) > std::get<1>(result));
	return result;
}

void ReadMatchposStorage::emplace_back(uint32_t read, uint32_t startpos, uint32_t endpos)
{
	assert(endpos > startpos);
	assert((startpos & 0x7FFFFFFF) < (1ull<<24ull));
	assert((endpos & 0x7FFFFFFF) < (1ull<<24ull));
	if (values == nullptr) resize(10);
	if (realsize == capacity) resize(2*capacity);
	size_t pos = realsize*10;
	values[pos] = read >> 24;
	values[pos+1] = read >> 16;
	values[pos+2] = read >> 8;
	values[pos+3] = read;
	values[pos+4] = startpos >> 16;
	values[pos+5] = startpos >> 8;
	values[pos+6] = startpos;
	values[pos+7] = endpos >> 16;
	values[pos+8] = endpos >> 8;
	values[pos+9] = endpos;
	values[pos+4] += (startpos & 0x80000000) >> 24; // reverse orientation bit
	realsize += 1;
}

void ReadMatchposStorage::emplace_back(__uint128_t value)
{
	uint32_t read = value >> 64;
	uint32_t startpos = value >> 32;
	uint32_t endpos = value;
	emplace_back(read, startpos, endpos);
}

size_t ReadMatchposStorage::size() const
{
	return realsize;
}

ReadMatchposStorage::ReadMatchposStorageIterator ReadMatchposStorage::begin() const
{
	return ReadMatchposStorage::ReadMatchposStorageIterator { *this, 0 };
}

ReadMatchposStorage::ReadMatchposStorageIterator ReadMatchposStorage::end() const
{
	return ReadMatchposStorage::ReadMatchposStorageIterator { *this, size() };
}

void ReadMatchposStorage::compact()
{
	resize(realsize);
}

void ReadMatchposStorage::resize(size_t newCapacity)
{
	uint8_t* newValues = new uint8_t[newCapacity*10];
	if (values != nullptr) memcpy(newValues, values, 10*realsize);
	std::memset(newValues+10*realsize, 0, (newCapacity-realsize)*10);
	capacity = newCapacity;
	std::swap(values, newValues);
	if (newValues != nullptr) delete [] newValues;
}

void ReadMatchposStorage::reserve(size_t newCapacity)
{
	resize(newCapacity);
}

void swap(ReadMatchposStorage& left, ReadMatchposStorage& right)
{
	std::swap(left.values, right.values);
	std::swap(left.capacity, right.capacity);
	std::swap(left.realsize, right.realsize);
}
