#include "FRAGHIST_EXTERNAL.h"
#include "MOB.h"
#include <fstream>
#include <sstream>
#include <algorithm>


/* =========================================== */
// merge histogram in the file into the one in the running program,
// and return how many histogram items existed in the file
int FRAGHIST_EXTERNAL::getHistFromFile(std::string fileName, int startIdx) {
	std::string line;
	std::ifstream histFile(fileName);
	while (std::getline(histFile, line)) {
		if (line == "$") {
			break;
		}
	}
	int count = 0;
	while (std::getline(histFile, line)) {
		std::string inchi = line;
		std::vector<std::vector<int>> dataItem{};

		std::getline(histFile, line, ':');
		std::vector<int> y{ std::stoi(line) };
		dataItem.push_back(y);
		std::getline(histFile, line);

		std::getline(histFile, line);
		std::istringstream lineStream(line);
		std::string temp;
		std::vector<int> ID;
		while (std::getline(lineStream, temp, ')')) {
			size_t i = temp.find("(");
			if (i != std::string::npos) {
				ID = retrieveBondsID(temp.substr(i + 1));
				for (int& j : ID) {
					j += startIdx;
				}
				dataItem.push_back(ID);
			}
		}
		std::getline(histFile, line); // read the empty line

		// if this hist item is already exist, then append
		auto thisOne = data.insert({inchi, dataItem});
		if (!thisOne.second) {
			std::vector<std::vector<int>>& it = thisOne.first->second;
			it[0][0] += dataItem[0][0];
			for (size_t j = 1; j < dataItem.size(); ++j) {
				it.push_back(dataItem[j]);
			}
		}

		count++;
	}
	histFile.close();
	return count;
}


/* =========================================== */
FRAGHIST FRAGHIST_EXTERNAL::convertToFRAGHIST() {
	FRAGHIST hist0{};
	std::vector<std::string> x_orig{}; // vector of identifier
	std::vector<float> y_mean_orig{};

	for (auto& dataItem : data) {
		x_orig.push_back(dataItem.first);
		y_mean_orig.push_back(static_cast<float>(dataItem.second[0][0]));
		std::set<std::string> IDs{};
		for (size_t i = 1; i < dataItem.second.size(); ++i) {
			IDs.insert(getIdentifierStr(dataItem.second[i]));
		}
		hist0.xID_orig.insert({ dataItem.first, IDs});
	}


	// sort them by meanVal (descending)
	std::vector<std::pair<size_t, float> > order(y_mean_orig.size());
	size_t n = 0;
	for (auto it = y_mean_orig.begin(); it != y_mean_orig.end(); ++it, ++n)
		order[n] = std::make_pair(n, *it);
	std::sort(order.begin(), order.end(),
		[](auto& a, auto& b) { return a.second > b.second; });

	// filter items where mean-value != xID.size()
	for (auto& i : order) {
		hist0.x_orig.push_back(x_orig[i.first]);
		hist0.y_mean_orig.push_back(y_mean_orig[i.first]);
	}


	// filer non-duplications
	for (size_t i = 0; i < hist0.x_orig.size(); ++i) {
		if (int(std::round(hist0.y_mean_orig[i])) == 1) {
			continue;
		}
		hist0.x.push_back(hist0.x_orig[i]);
		hist0.y.push_back(int(std::round(hist0.y_mean_orig[i])));
		hist0.xID.insert({ hist0.x_orig[i], hist0.xID_orig[hist0.x_orig[i]] });
	}


	return hist0;
}