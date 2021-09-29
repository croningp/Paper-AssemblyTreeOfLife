#pragma once
#include <unordered_map>
#include <set>
#include <vector>
typedef unsigned long long int MYINT;
const long double z_quantile_Wilson_score = 1.96;
const double Hist_Err_Threshold = 0.05;


class FRAGHIST {
public:
	std::vector<std::string> x_orig{}; // vector of identifier
	std::vector<float> y_mean_orig{}, y_err_orig{};
	std::vector<std::string> x{};// fragment inchi string
	std::vector<size_t> y{};
	std::unordered_map<std::string, std::set<std::string>> 
		xID_orig{}, xID{}; // the inchi corresponding bonds list

	MYINT totalMexp = 0;// how many cuts went through
	std::unordered_map<std::string, 
		std::unordered_map<std::string, std::vector<MYINT>> > histOrig{};

	// --------- functions -----------
	FRAGHIST() {};
	void addOnce(std::vector<std::tuple<
		std::string, std::string, MYINT>>& histItem);
	void getHistOverall(MYINT totalM, long double z = z_quantile_Wilson_score);

	void getFinalHist(bool filterLargeErr = true, double ErrThreshold = Hist_Err_Threshold);
};