#include "common.h"

std::ostream& operator<<(std::ostream& os, const Path& path)
{
	for (const auto& state : path)
	{
		os << state.location << ", "; // << "(" << state.is_single() << "),";
	}
	return os;
}


bool isSamePath(const Path& p1, const Path& p2)
{
	if (p1.size() != p2.size())
		return false;
	for (unsigned i = 0; i < p1.size(); i++)
	{
		if (p1[i].location != p2[i].location)
			return false;
	}
	return true;
}

pair<vector<int>, int> myfindMaxValuePair(
    unordered_map<vector<int>, int, container_hash<vector<int>>> const &x)
{
    return *std::max_element(x.begin(), x.end(),
                             [](const pair<vector<int>, int> &p1,
                                const pair<vector<int>, int> &p2)
                             {
                                return p1.second < p2.second;
                             });
}