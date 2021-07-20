#pragma once
#include "CBSHeuristic.h"
#include "RectangleReasoning.h"
#include "CorridorReasoning.h"
#include "MutexReasoning.h"

enum high_level_solver_type { ASTAR, ASTAREPS, NEW, EES };

class CBS
{
public:
	bool randomRoot = false; // randomize the order of the agents in the root CT node

	/////////////////////////////////////////////////////////////////////////////////////
	// stats
	double runtime = 0;
	double runtime_generate_child = 0; // runtime of generating child nodes
	double runtime_build_CT = 0; // runtime of building constraint table
	double runtime_build_CAT = 0; // runtime of building conflict avoidance table
	double runtime_path_finding = 0; // runtime of finding paths for single agents
	double runtime_detect_conflicts = 0;
	double runtime_preprocessing = 0; // runtime of building heuristic table for the low level
	double runtime_solver = 0;  // runtime for using inner solver for meta-agents
	double runtime_sort_ma = 0;  // runtime for computing priority heuristic

	uint64_t num_cardinal_conflicts = 0;
	uint64_t num_corridor_conflicts = 0;
	uint64_t num_rectangle_conflicts = 0;
	uint64_t num_target_conflicts = 0;
	uint64_t num_mutex_conflicts = 0;
	uint64_t num_standard_conflicts = 0;

	uint64_t num_adopt_bypass = 0; // number of times when adopting bypasses
	uint64_t num_push_focal = 0;  // number of child CT node being pushed to FOCAL

	uint64_t num_HL_expanded = 0;
	uint64_t num_HL_generated = 0;
	uint64_t num_LL_expanded = 0;
	uint64_t num_LL_generated = 0;

	uint64_t solver_counter;
	uint64_t solver_num_HL_expanded = 0;
	uint64_t solver_num_HL_generated = 0;
	uint64_t solver_num_LL_expanded = 0;
	uint64_t solver_num_LL_generated = 0;

	uint64_t num_cleanup = 0; // number of expanded nodes chosen from cleanup list
	uint64_t num_open = 0; // number of expanded nodes chosen from open list
	uint64_t num_focal = 0; // number of expanded nodes chsoen from focal list

	// statistics for branch and every iteration
	int cleanup_head_lb;
    std::shared_ptr<vector<int>> iter_sum_lb;
    std::shared_ptr<vector<int>> br_sum_lb;
    std::shared_ptr<vector<int>> all_sum_lb;
    std::shared_ptr<vector<int>> open_sum_lb;

    std::shared_ptr<vector<int>> iter_sum_cost;
    std::shared_ptr<vector<int>> br_sum_cost;
    std::shared_ptr<vector<int>> all_sum_cost;
    std::shared_ptr<vector<int>> open_sum_cost;
    
    std::shared_ptr<vector<int>> iter_num_conflicts;
    std::shared_ptr<vector<int>> br_num_conflicts;
    std::shared_ptr<vector<int>> all_num_conflicts;
    std::shared_ptr<vector<int>> open_num_conflicts;

    std::shared_ptr<vector<double>> iter_remained_flex;
    std::shared_ptr<vector<double>> br_remained_flex;
    std::shared_ptr<vector<double>> all_remained_flex;
    std::shared_ptr<vector<double>> open_remained_flex;

	std::shared_ptr<vector<double>> iter_subopt;
    std::shared_ptr<vector<double>> br_subopt;
    std::shared_ptr<vector<double>> all_subopt;

    std::shared_ptr<vector<uint64_t>> iter_sum_ll_generate;
    std::shared_ptr<vector<uint64_t>> br_sum_ll_generate;
    std::shared_ptr<vector<uint64_t>> all_sum_ll_generate;

    std::shared_ptr<vector<int>> iter_node_idx;
    std::shared_ptr<vector<int>> br_node_idx;
    std::shared_ptr<vector<int>> open_node_idx;
    std::shared_ptr<vector<int>> all_node_idx;

	std::shared_ptr<vector<uint64_t>> iter_num_focal;
	std::shared_ptr<vector<uint64_t>> iter_num_open;
	std::shared_ptr<vector<uint64_t>> iter_num_cleanup;
	std::shared_ptr<vector<int>> iter_node_type;
	// end of statistics for branch and every iteration

	// CBSNode* dummy_start = nullptr;
	// CBSNode* goal_node = nullptr;
	HLNode* dummy_start = nullptr;
	HLNode* goal_node = nullptr;



	bool solution_found = false;
	int solution_cost = -2;

	/////////////////////////////////////////////////////////////////////////////////////////
	// set params
	void setHeuristicType(heuristics_type h, heuristics_type h_hat)
	{
	    heuristic_helper.type = h;
	    heuristic_helper.setInadmissibleHeuristics(h_hat);
	}
	void setPrioritizeConflicts(bool p) {PC = p;	heuristic_helper.PC = p; }
	void setRectangleReasoning(bool r) {rectangle_reasoning = r; heuristic_helper.rectangle_reasoning = r; }
	void setCorridorReasoning(bool c) {corridor_reasoning = c; heuristic_helper.corridor_reasoning = c; }
	void setTargetReasoning(bool t) {target_reasoning = t; heuristic_helper.target_reasoning = t; }
	void setMutexReasoning(bool m) {mutex_reasoning = m; heuristic_helper.mutex_reasoning = m; }
	void setDisjointSplitting(bool d) {disjoint_splitting = d; heuristic_helper.disjoint_splitting = d; }
	void setBypass(bool b) { bypass = b; } // 2-agent solver for heuristic calculation does not need bypass strategy.
	void setConflictSelectionRule(conflict_selection c) { conflict_selection_rule = c; heuristic_helper.conflict_selection_rule = c; }
	void setNodeSelectionRule(node_selection n) { node_selection_rule = n; heuristic_helper.node_selection_rule = n; }
	void setSavingStats(bool s) { save_stats = s; heuristic_helper.save_stats = s; }
	void setHighLevelSolver(high_level_solver_type s, double w)
	{
		solver_type = s;
		suboptimality = w;
	}
	void setNodeLimit(int n) { node_limit = n; }
	void setMergeThreshold(int b) { merge_th = b; }
	void setMergeRestart(bool mr) { mr_active = mr; }
	void setMASolver(shared_ptr<CBS> in_solver) { inner_solver = in_solver; }
	void setMetaAgents(vector<vector<int>> in_meta_agents) { meta_agents = in_meta_agents; }
	void setMAVector(vector<bool> in_ma_vec) { ma_vec = in_ma_vec; }
	void setInitialPath(int agent, Path _path) { paths_found_initially[agent] = _path; }
	void setInitConstraints(int agent, ConstraintTable _table) {initial_constraints[agent].init(_table);}
	void setInitSumLB (int _sum_lb) { init_sum_lb = _sum_lb; }
	void setTimeLimit(double tl) { time_limit = tl; }
	void setIsStart(bool _st) { is_start = _st; }

	////////////////////////////////////////////////////////////////////////////////////////////
	// Runs the algorithm until the problem is solved or time is exhausted 
	bool solve(double time_limit, int cost_lowerbound = 0, int cost_upperbound = MAX_COST);

	int getLowerBound() const { return cost_lowerbound; }
	Path getPath(int agent) const { return *paths[agent]; }
	vector<Path> getPath(void) const;

	CBS(const Instance& instance, bool sipp, int screen);
	CBS(vector<SingleAgentSolver*>& search_engines,
		const vector<ConstraintTable>& constraints,
		vector<Path>& paths_found_initially, int screen);
	void clearSearchEngines();
	~CBS();

	// Save results
	void saveResults(const string &fileName, const string &instanceName) const;
	void saveStats(const string &fileName, const string &instanceName);
	void saveCT(const string &fileName) const; // write the CT to a file
    void savePaths(const string &fileName) const; // write the paths to a file
	void clear(); // used for rapid random  restart

	int getInitialPathLength(int agent) const {return (int) paths_found_initially[agent].size() - 1; }
	int getminFVal(int agent) const { return min_f_vals[agent]; }

protected:
    bool rectangle_reasoning;  // using rectangle reasoning
	bool corridor_reasoning;  // using corridor reasoning
	bool target_reasoning;  // using target reasoning
	bool disjoint_splitting;  // disjoint splitting
	bool mutex_reasoning;  // using mutex reasoning
	bool bypass; // using Bypass1
	bool PC; // prioritize conflicts
	bool save_stats;
	high_level_solver_type solver_type; // the solver for the high-level search
	conflict_selection conflict_selection_rule;
	node_selection node_selection_rule;

	MDDTable mdd_helper;	
	RectangleReasoning rectangle_helper;
	CorridorReasoning corridor_helper;
	MutexReasoning mutex_helper;
	CBSHeuristic heuristic_helper;

	list<HLNode*> allNodes_table; // this is ued for both ECBS and EES


	string getSolverName() const;

	int screen;
	
	double time_limit;
	double suboptimality = 1.0;
	int cost_lowerbound = 0;
	int inadmissible_cost_lowerbound;
	int node_limit = MAX_NODES;
	int cost_upperbound = MAX_COST;

	vector<ConstraintTable> initial_constraints;
	clock_t start;

	int num_of_agents;

	vector<Path*> paths;
	vector<Path> paths_found_initially;  // contain initial paths found
	// vector<MDD*> mdds_initially;  // contain initial paths found
	vector < SingleAgentSolver* > search_engines;  // used to find (single) agents' paths and mdd

	// For nested framework
	shared_ptr<CBS> inner_solver;  // inner (E)CBS for solving meta-agents
	vector<vector<int>> meta_agents;
	vector<bool> ma_vec;
	vector<vector<int>> conflict_matrix;
	bool is_solver = false;
	bool mr_active = false;
	bool is_start = false;  // This is for Merge and Restart
	int merge_th = INT_MAX;
	int init_sum_lb = 0;  // Obtain from outer (E)CBS, may be useless
	double flex = 0.0;  // flex for the meta-agent
	vector<int> min_f_vals; // lower bounds of the cost of the shortest path

	vector<int> findMetaAgent(int __ag__) const;
	bool shouldMerge(const vector<int>& __ma1__, const vector<int>& __ma2__, int mode=0) const;

	template <typename T>
    void sortMetaAgents(const vector<T>& sort_based, bool is_ascending);
	// End nested framework

	void addConstraints(const HLNode* curr, HLNode* child1, HLNode* child2) const;
	set<int> getInvalidAgents(const list<Constraint>& constraints); // return agents that violates the constraints
	//conflicts
	void findConflicts(HLNode& curr);
	void findConflicts(HLNode& curr, int a1, int a2);
	shared_ptr<Conflict> chooseConflict(const HLNode &node) const;
	static void copyConflicts(const list<shared_ptr<Conflict>>& conflicts,
		list<shared_ptr<Conflict>>& copy, const list<int>& excluded_agent) ;
	void removeLowPriorityConflicts(list<shared_ptr<Conflict>>& conflicts) const;
	void computeSecondPriorityForConflict(Conflict& conflict, const HLNode& node);

	inline void releaseNodes();

	// print and save
	void printResults() const;
	static void printConflicts(const HLNode &curr) ;

	bool validateSolution() const;
	inline int getAgentLocation(int agent_id, size_t timestep) const;

	vector<int> shuffleAgents() const;  //generate random permuattion of agent indices
	bool terminate(HLNode* curr); // check the stop condition and return true if it meets
	void computeConflictPriority(shared_ptr<Conflict>& con, CBSNode& node); // check the conflict is cardinal, semi-cardinal or non-cardinal

	void getBranchEval(HLNode* __node__, int open_head_lb);
	void saveEval(void);
	void saveNumNodesInLists(void);

private: // CBS only, cannot be used by ECBS
	pairing_heap< CBSNode*, compare<CBSNode::compare_node_by_f> > cleanup_list; // it is called open list in ECBS
	pairing_heap< CBSNode*, compare<CBSNode::compare_node_by_inadmissible_f> > open_list; // this is used for EES
	pairing_heap< CBSNode*, compare<CBSNode::compare_node_by_d> > focal_list; // this is ued for both ECBS and EES

	// node operators
	inline void pushNode(CBSNode* node);
	CBSNode* selectNode();
	inline bool reinsertNode(CBSNode* node);

		 // high level search
	bool generateChild(CBSNode* child, CBSNode* curr);
	bool generateRoot();
	bool findPathForSingleAgent(CBSNode*  node, int ag, int lower_bound = 0);
	void classifyConflicts(CBSNode &parent);
		 //update information
	inline void updatePaths(CBSNode* curr);
	void printPaths() const;
};
