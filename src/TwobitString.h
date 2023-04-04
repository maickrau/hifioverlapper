#ifndef TwobitString_h
#define TwobitString_h

#include <cstdint>
#include <string>
#include <vector>

class TwobitString
{
public:
	TwobitString() = default;
	TwobitString(const std::string&);
	std::string toString() const;
	std::string substr(size_t start, size_t size) const;
	void resize(size_t size);
	uint8_t get(size_t i) const;
	void set(size_t i, uint8_t v);
	void emplace_back(uint8_t v);
	size_t size() const;
private:
	size_t realSize;
	std::vector<uint8_t> bits;
};

#endif
