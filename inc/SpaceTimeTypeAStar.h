#pragma once
#include "SingleAgentSolver.h"
#include "SpaceTimeAStar.h"


class SpaceTimeTypeAStar: public SpaceTimeAStar
{
public:
	// find path by time-space A* search
	// Returns a shortest path that satisfies the constraints of the give node  while
	// minimizing the number of internal conflicts (that is conflicts with known_paths for other agents found so far).
	// lowerbound is an underestimation of the length of the path in order to speed up the search.
	Path findOptimalPath(const HLNode& node, const ConstraintTable& initial_constraints,
						const vector<Path*>& paths, int agent, int lower_bound);
	
	// return the path and the lowerbound
	pair<Path, int> findSuboptimalPath(const HLNode& node, const ConstraintTable& initial_constraints,
		const vector<Path*>& paths, int agent, int lowerbound, double w, int other_sum_lb=0, int other_sum_cost=0,
		int outer_sum_lb=0, double single_flex=0.0, int hl_h_val=0);

	string getName() const { return "TypeAStar"; }

	SpaceTimeTypeAStar(const Instance& instance, int agent):
		SpaceTimeAStar(instance, agent) {}

private:
	// define typedefs and handles for heap
	unordered_map<int, heap_focal_t> type_lists;

	// Updates the path datamember
	void updateFocalList(int lowerbound, int other_sum_lb, int other_sum_cost, int outer_sum_lb, double sinfle_flex, int hl_h_val);
	inline AStarNode* popNode();
	inline void pushNode(AStarNode* node);
	inline void pushNodeToTypes(AStarNode* node);
	void releaseNodes();
	int getNumInTypes(void);
	bool isTypeValid(void);

	double upperbound;  // For FEECBS
	uint64_t iter_counter = 0;
};
