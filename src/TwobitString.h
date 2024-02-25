#ifndef TwobitString_h
#define TwobitString_h

#include <cstdint>
#include <string>
#include <vector>

class TwobitString
{
public:
	TwobitString();
	TwobitString(const std::string&);
	std::string toString() const; // output chars are ATCG
	std::string substr(size_t start, size_t size) const; // output chars are ATCG
	void resize(size_t size);
	uint8_t get(size_t i) const; // output char is 0123
	void set(size_t i, uint8_t v); // input char is 0123
	void emplace_back(uint8_t v); // input char is 0123
	size_t size() const;
	bool operator==(const TwobitString& other) const;
	bool operator<(const TwobitString& other) const;
private:
	size_t realSize;
	std::vector<uint8_t> bits;
};

#endif
