#pragma once
#include "FRAGHIST.h"
#include <unordered_map>
#include <vector>

class FRAGHIST_EXTERNAL {
public:
	std::unordered_map<std::string, std::vector<std::vector<int>>> data{};

	// ----- functions ------
	FRAGHIST_EXTERNAL() {};
	int getHistFromFile(std::string fileName, int startIdx);
	FRAGHIST convertToFRAGHIST();
};