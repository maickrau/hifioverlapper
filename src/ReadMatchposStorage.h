#ifndef ReadMatchposStorage_h
#define ReadMatchposStorage_h

#include <vector>
#include <tuple>
#include <cstdint>
#include <cstdio>

class ReadMatchposStorage
{
public:
	class ReadMatchposStorageIterator
	{
	public:
		ReadMatchposStorageIterator(const ReadMatchposStorage& storage, size_t index);
		ReadMatchposStorageIterator& operator++();
		bool operator==(const ReadMatchposStorageIterator& other) const;
		bool operator!=(const ReadMatchposStorageIterator& other) const;
		std::tuple<uint32_t, uint32_t, uint32_t> operator*() const;
	private:
		const ReadMatchposStorage& storage;
		size_t index;
	};
	ReadMatchposStorage();
	~ReadMatchposStorage();
	ReadMatchposStorage(const ReadMatchposStorage& other) = delete; // not needed for now
	ReadMatchposStorage(ReadMatchposStorage&& other);
	ReadMatchposStorage& operator=(const ReadMatchposStorage& other) = delete; // not needed for now
	ReadMatchposStorage& operator=(ReadMatchposStorage&& other);
	std::tuple<uint32_t, uint32_t, uint32_t> getValue(size_t index) const;
	void emplace_back(uint32_t read, uint32_t startpos, uint32_t endpos);
	void emplace_back(__uint128_t value);
	size_t size() const;
	void compact();
	ReadMatchposStorageIterator begin() const;
	ReadMatchposStorageIterator end() const;
	void reserve(size_t newCapacity);
private:
	void resize(size_t newCapacity);
	uint16_t* values;
	uint32_t realsize; // uint32_t because even this saves a few Gb in large datasets
	uint32_t capacity; // uint32_t because even this saves a few Gb in large datasets
	friend void swap(ReadMatchposStorage& left, ReadMatchposStorage& right);
};

void swap(ReadMatchposStorage& left, ReadMatchposStorage& right);

#endif
