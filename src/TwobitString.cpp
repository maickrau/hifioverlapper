#include <cassert>
#include "TwobitString.h"

TwobitString::TwobitString(const std::string& str)
{
	resize(str.size());
	for (size_t i = 0; i < str.size(); i++)
	{
		switch(str[i])
		{
			case 'a':
			case 'A':
				set(i, 0);
				break;
			case 'c':
			case 'C':
				set(i, 1);
				break;
			case 'g':
			case 'G':
				set(i, 2);
				break;
			case 't':
			case 'T':
				set(i, 3);
				break;
			default:
				assert(false);
				break;
		}
	}
}

std::string TwobitString::toString() const
{
	std::string result;
	result.resize(size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = "ACGT"[get(i)];
	}
	return result;
}

void TwobitString::resize(size_t size)
{
	realSize = size;
	bits.resize((realSize+3)/4);
}

uint8_t TwobitString::get(size_t i) const
{
	size_t index = i / 4;
	size_t offset = (i % 4) * 2;
	return (bits[index] >> offset) & 3;
}

void TwobitString::set(size_t i, uint8_t v)
{
	assert(v < 4);
	size_t index = i / 4;
	size_t offset = (i % 4) * 2;
	uint8_t removeMask = ~(3 << offset);
	bits[index] &= removeMask;
	uint8_t addMask = (uint8_t)v << offset;
	bits[index] |= addMask;
}

void TwobitString::emplace_back(uint8_t v)
{
	realSize += 1;
	if (realSize % 4 == 0) bits.emplace_back(0);
	set(realSize-1, v);
}

size_t TwobitString::size() const
{
	return realSize;
}
