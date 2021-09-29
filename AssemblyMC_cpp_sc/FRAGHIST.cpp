#include "FRAGHIST.h"
#include <cmath>
#include <algorithm>


/* =========================================== */
// add one cutting scheme to histogram
void FRAGHIST::addOnce(std::vector<std::tuple<
	std::string, std::string, MYINT>>& histItem) {

	// new item in the most specific sense {count, M0}
	// ID: identifierStr2: then we have newItemSpec
	std::vector<MYINT> newItemSpec{1, 0};

	for (auto& thisItem : histItem) {
		std::string& ID = std::get<0>(thisItem);
		std::string& identifierStr2 = std::get<1>(thisItem);
		MYINT& M0 = std::get<2>(thisItem);

		auto it = histOrig.find(ID);
		if (it == histOrig.end()) { // didn't find
			newItemSpec[1] = M0;
			histOrig[ID] = std::unordered_map<std::string, std::vector<MYINT>>
				{ {identifierStr2, newItemSpec} };
		}
		else { // found
			auto ItemSpec = it->second.find(identifierStr2);
			if (ItemSpec == it->second.end()) {
				newItemSpec[1] = M0;
				it->second[identifierStr2] = newItemSpec;
			}
			else {
				ItemSpec->second[0] += 1;
			}
		}
	}
	totalMexp += 1;
}


/* =========================================== */
// calculate the final histogram that will be used in MC step 2.
// zValue is the quantile of a standard normal distribution: 95% -> 1.96
void FRAGHIST::getHistOverall(MYINT totalM, long double z) {
	std::vector<std::string> x_orig0{};
	std::vector<float> y_mean_orig0{}, y_err_orig0{};
	xID_orig.clear();

	long double z2 = z * z, z2d2 = z2 / 2.0, z2d4 = z2 / 4.0;
	long double PdP = (totalMexp + z2) / totalM;
	long double count, meanVal, errVal;

	for (auto& it : histOrig) {// go through every fragment
		meanVal = 0;
		errVal = 0;
		std::set<std::string> idTemp{};
		for (auto& identifierItem : it.second) { // identifierItem.second = (count, M0)
			count = static_cast<long double>(identifierItem.second[0]);
			meanVal += (count + z2d2) / identifierItem.second[1] / PdP;
			errVal += pow( z / PdP * (
				sqrt(count*(totalMexp-count)/totalMexp + z2d4)
				/ identifierItem.second[1]
				), 2);
			idTemp.insert(identifierItem.first);
		}
		errVal = sqrt(errVal) / meanVal;
		
		xID_orig[it.first] = idTemp;
		x_orig0.push_back(it.first);
		y_mean_orig0.push_back(float(meanVal));
		y_err_orig0.push_back(float(errVal));
	}
	
	// sort them by meanVal (descending)
	std::vector<std::pair<size_t, float> > order(y_mean_orig0.size());
	size_t n = 0;
	for (auto it = y_mean_orig0.begin(); it != y_mean_orig0.end(); ++it, ++n)
		order[n] = std::make_pair(n, *it);
	std::sort(order.begin(), order.end(),
		[](auto& a, auto& b) { return a.second > b.second; });

	x_orig.clear();
	y_mean_orig.clear();
	y_err_orig.clear();
	for (auto& i : order) {
		if (int(xID_orig[x_orig0[i.first]].size())
			== int(std::round(y_mean_orig0[i.first]))) {
			x_orig.push_back(x_orig0[i.first]);
			y_mean_orig.push_back(y_mean_orig0[i.first]);
			y_err_orig.push_back(y_err_orig0[i.first]);
		}
	}
}


/* =========================================== */
// get the final hist we will use for Monte Carlo step 2.
void FRAGHIST::getFinalHist(bool filterLargeErr, double ErrThreshold) {
	x.clear();
	y.clear();
	xID.clear();
	for (size_t i = 0; i < x_orig.size(); ++i) {
		if (filterLargeErr && y_err_orig[i] > ErrThreshold) {
			continue;
		}
		if (int(std::round(y_mean_orig[i])) == 1) {
			continue;
		}
		x.push_back(x_orig[i]);
		y.push_back(int( std::round(y_mean_orig[i]) ));
		auto it = xID.find(x_orig[i]);
		if (it == xID.end()) {
			xID[x_orig[i]] = xID_orig[x_orig[i]];
		}
	}
}
