#pragma once
#include "CBS.h"
#include "ECBSNode.h"


class ECBS : public CBS
{
public:
	ECBS(const Instance& instance, bool sipp, int screen) : CBS(instance, sipp, screen) {}
	ECBS(vector<SingleAgentSolver*>& search_engines, 
		const vector<ConstraintTable>& init_constraints,
		vector<Path>& init_paths, 
		vector<int>& init_min_f_vals, int screen) : CBS(search_engines, init_constraints, init_paths, screen){
		paths_found_initially.resize(num_of_agents);
		for (int _ag_ = 0; _ag_ < num_of_agents; _ag_ ++)
		{
			setInitialPath(_ag_, init_paths[_ag_], init_min_f_vals[_ag_]);
			setMinFVal(_ag_, init_min_f_vals[_ag_]);
		}
	}
	void setInitialPath(int agent, Path _path, int _min_f_val=0)
	{ 
		if (paths_found_initially.empty())
			paths_found_initially.resize(num_of_agents);
		paths_found_initially[agent].first = _path;
		paths_found_initially[agent].second = _min_f_val;
	}
	int getInitialPathLength(int agent) const override {return (int) paths_found_initially[agent].first.size() - 1; }

	////////////////////////////////////////////////////////////////////////////////////////////
	// Runs the algorithm until the problem is solved or time is exhausted 
	bool solve(double time_limit, int _cost_lowerbound = 0, int _cost_upperbound = MAX_COST) override;
    void clear() override; // used for rapid random  restart

private:
	vector< pair<Path, int> > paths_found_initially;  // contain initial paths found
	pairing_heap< ECBSNode*, compare<ECBSNode::compare_node_by_f> > cleanup_list; // it is called open list in ECBS
	pairing_heap< ECBSNode*, compare<ECBSNode::compare_node_by_inadmissible_f> > open_list; // this is used for EES
	pairing_heap< ECBSNode*, compare<ECBSNode::compare_node_by_d> > focal_list; // this is ued for both ECBS and EES

	void adoptBypass(ECBSNode* curr, ECBSNode* child, const vector<int>& fmin_copy, const vector<Path*>& path_copy);

	// node operators
	void pushNode(ECBSNode* node);
	ECBSNode* selectNode();
	bool reinsertNode(ECBSNode* node);

	// high level search
	bool generateChild(ECBSNode* child, ECBSNode* curr, int child_idx=0);
	bool generateRoot();
	bool findPathForSingleAgent(ECBSNode*  node, int ag);
	bool findPathForMetaAgent(ECBSNode* node, const vector<int>& ma1, const vector<int>& ma2=vector<int>());
	void classifyConflicts(ECBSNode &node);
	void computeConflictPriority(shared_ptr<Conflict>& con, ECBSNode& node);

	// For NFECBS and NFEECBS
	void getFlex(const vector<int>& agent);

	//update information
	void updatePaths(ECBSNode* curr);
	void printPaths() const;
	void printAgentPath(int ag, Path* path_ptr=nullptr) const;

	// Visualization results
	bool validatePaths(ECBSNode* node);
	bool validatePathswithCurrConflicts(ECBSNode* node);
	void initializeIterAnalysis(void);
	const vector<pair<Path*, int>> collectPaths(ECBSNode* node);
};