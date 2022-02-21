#pragma once
#include "MDD.h"
#include "RectangleReasoning.h"
#include "CorridorReasoning.h"


enum heuristics_type { ZERO, CG, DG, WDG, GLOBAL, PATH, LOCAL, CONFLICT, STRATEGY_COUNT }; //  GREEDY,

struct DoubleConstraintsHasher // Hash a CT node by constraints on two agents (look-up table entry)
{
	int a1{};
	int a2{};
	HLNode* n{};

	DoubleConstraintsHasher() = default;
	DoubleConstraintsHasher(int a1, int a2, HLNode* n) : a1(a1), a2(a2), n(n) {};

	struct EqNode
	{
		bool operator() (const DoubleConstraintsHasher& h1, const DoubleConstraintsHasher& h2) const
		{
			if (h1.a1 != h2.a1 || h1.a2 != h2.a2)
				return false;
			std::set<Constraint> cons1[2], cons2[2];
			HLNode* curr = h1.n;
			while (true)
			{
				for (const auto& constraint : curr->constraints)
				{
					if (get<4>(constraint) == constraint_type::LEQLENGTH ||
						get<4>(constraint) == constraint_type::POSITIVE_VERTEX ||
						get<4>(constraint) == constraint_type::POSITIVE_EDGE)
					{
						cons1[0].insert(constraint);
						cons2[0].insert(constraint);
					}
					else if (get<0>(constraint) == h1.a1)
						cons1[0].insert(constraint);
					else if (get<0>(constraint) == h1.a2)
						cons2[0].insert(constraint);
				}
				if (curr->parent != nullptr)
					curr = curr->parent;
				else  // The root CT node may containt constraints (for inner solver)
					break;
			}
			curr = h2.n;
			while (true)
			{
				for (const auto& constraint : curr->constraints)
				{
					if (get<4>(constraint) == constraint_type::LEQLENGTH ||
						get<4>(constraint) == constraint_type::POSITIVE_VERTEX ||
						get<4>(constraint) == constraint_type::POSITIVE_EDGE)
					{
						cons1[1].insert(constraint);
						cons2[1].insert(constraint);
					}
					else if (get<0>(constraint) == h2.a1)
						cons1[1].insert(constraint);
					else if (get<0>(constraint) == h2.a2)
						cons2[1].insert(constraint);
				}
				if (curr->parent != nullptr)
					curr = curr->parent;
				else  // The root CT node may containt constraints (for inner solver)
					break;
			}
			if (cons1[0].size() != cons1[1].size() || cons2[0].size() != cons2[1].size())
				return false;

			if (!equal(cons1[0].begin(), cons1[0].end(), cons1[1].begin()))
				return false;
			return equal(cons2[0].begin(), cons2[0].end(), cons2[1].begin());
		}
	};


	struct Hasher
	{
		size_t operator()(const DoubleConstraintsHasher& entry) const
		{
			HLNode* curr = entry.n;
			size_t cons1_hash = 0, cons2_hash = 0;
			while (true)
			{
				for (const auto& constraint : curr->constraints)
				{
					if (get<0>(constraint) == entry.a1)
					{
						cons1_hash += 3 * std::hash<int>()(std::get<0>(constraint)) +
							5 * std::hash<int>()(std::get<1>(constraint)) +
							7 * std::hash<int>()(std::get<2>(constraint)) +
							11 * std::hash<int>()(std::get<3>(constraint));
					}
					else if (get<0>(constraint) == entry.a2)
					{
						cons2_hash += 3 * std::hash<int>()(std::get<0>(constraint)) +
							5 * std::hash<int>()(std::get<1>(constraint)) +
							7 * std::hash<int>()(std::get<2>(constraint)) +
							11 * std::hash<int>()(std::get<3>(constraint));
					}
					else if (get<4>(constraint) == constraint_type::LEQLENGTH ||
						get<4>(constraint) == constraint_type::POSITIVE_VERTEX ||
						get<4>(constraint) == constraint_type::POSITIVE_EDGE)
					{
						cons1_hash += 3 * std::hash<int>()(std::get<0>(constraint)) +
							5 * std::hash<int>()(std::get<1>(constraint)) +
							7 * std::hash<int>()(std::get<2>(constraint)) +
							11 * std::hash<int>()(std::get<3>(constraint));
						cons2_hash += 3 * std::hash<int>()(std::get<0>(constraint)) +
							5 * std::hash<int>()(std::get<1>(constraint)) +
							7 * std::hash<int>()(std::get<2>(constraint)) +
							11 * std::hash<int>()(std::get<3>(constraint));
					}
				}
				if (curr->parent != nullptr)
					curr = curr->parent;
				else  // The root CT node may containt constraints (for inner solver)
					break;
			}
			return cons1_hash ^ (cons2_hash << 1);
		}
	};
};

// <h value, num of CT nodes, 0> for CBS
// <h value, a1 f at root, a2 f at root> for ECBS
typedef unordered_map<DoubleConstraintsHasher, tuple<int, int, int>, DoubleConstraintsHasher::Hasher, DoubleConstraintsHasher::EqNode> HTable;


class CBSHeuristic
{
public:
	heuristics_type type;
	bool rectangle_reasoning; // using rectangle reasoning
	bool corridor_reasoning; // using corridor reasoning
	bool target_reasoning; // using target reasoning
	bool mutex_reasoning; // using mutex reasoning
	bool disjoint_splitting; // disjoint splitting
	bool PC; // prioritize conflicts
	bool is_solver;  // check if we are calculating the heuristic of outer ECBS or Inner ECBS

	bool save_stats;
	conflict_selection conflict_selection_rule;
	node_selection node_selection_rule;

	double runtime_build_dependency_graph = 0;
	double runtime_solve_MVC = 0;
	uint64_t num_solve_MVC = 0;
	uint64_t num_merge_MDDs = 0;
	uint64_t num_solve_2agent_problems = 0;
	uint64_t num_memoization = 0; // number of times when memorization helps

	// stats
	list<tuple<int, int, const HLNode*, uint64_t, int> > sub_instances; 	// <agent 1, agent 2, node, number of expanded CT nodes, h value> 

	CBSHeuristic(int num_of_agents,
				vector<Path*>& paths,
				vector<SingleAgentSolver*>& search_engines,
				vector<ConstraintTable>& initial_constraints,
				MDDTable& mdd_helper) : num_of_agents(num_of_agents),
		paths(paths), search_engines(search_engines), initial_constraints(initial_constraints), mdd_helper(mdd_helper) {}
	
	void init()
	{
		if (type == heuristics_type::DG || type == heuristics_type::WDG)
		{
			lookupTable.resize(num_of_agents);
			for (int i = 0; i < num_of_agents; i++)
			{
				lookupTable[i].resize(num_of_agents);
			}
		}
	}

	inline void setNumOfAgents(int num_of_ags)
	{
		num_of_agents = num_of_ags;
	}

	inline void setMetaAgents(vector<int> in_ags)
	{
		meta_agents.clear();
		ma_vec = vector<bool>(num_of_agents, false);
		for (const int& _ag_ : in_ags)
		{
			meta_agents.push_back(vector<int>({_ag_}));
			ma_vec[_ag_] = true;
		}
	}

	inline void setMetaAgents(vector<vector<int>> in_ma)
	{
		meta_agents = in_ma;
		for (const vector<int>& _ma_ : meta_agents)
			for (const int& _ag_ : _ma_)
				ma_vec[_ag_] = true;
	}

	inline void setInitConstraints(vector<ConstraintTable> in_cons_table)
	{
		initial_constraints = in_cons_table;
	}

	inline void setInitConstraints(ConstraintTable in_cons, int agent)
	{
		initial_constraints[agent] = in_cons;
	}

	void setInadmissibleHeuristics(heuristics_type h)
    {
        inadmissible_heuristic = h;
        if (h == heuristics_type::CONFLICT)
        {
            sum_distance_errors.assign(int(conflict_type::TYPE_COUNT), 0);  //(int(conflict_priority::PRIORITY_COUNT) * int(conflict_type::TYPE_COUNT), 0);
            sum_cost_errors.assign(int(conflict_type::TYPE_COUNT), 0);
            num_of_errors.assign(int(conflict_type::TYPE_COUNT), 1);
        }
        else
        {
            sum_distance_errors.assign(1, 0);
            sum_cost_errors.assign(1, 0);
            num_of_errors.assign(1, 0);
        }
    }
	bool computeInformedHeuristics(CBSNode& curr, double time_limit);
	bool computeInformedHeuristics(ECBSNode& curr, const vector<int>& min_f_vals, double time_limit);
	void computeQuickHeuristics(HLNode& curr);
	void updateOnlineHeuristicErrors(CBSNode& curr);
	void updateOnlineHeuristicErrors(ECBSNode& curr);
    void updateInadmissibleHeuristics(HLNode& curr);

	// EES heuristics
	// int MVConAllConflicts(HLNode& curr);
	// int greedyWDG(CBSNode& curr, double time_limit);
	double getCostError(int i = 0) const { return (num_of_errors[i] == 0)? 0 : sum_cost_errors[i] / num_of_errors[i]; }
	double getDistanceError(int i = 0) const { return (num_of_errors[i] == 0)? 0 : sum_distance_errors[i]  / num_of_errors[i]; }
	heuristics_type getInadmissibleHeuristics(void) {return inadmissible_heuristic;};
	// void copyConflictGraph(HLNode& child, const HLNode& parent);
	void clear()
	{
		lookupTable.clear(); 
		meta_agents.clear();
		ma_vec.clear();
	}

private:
    heuristics_type inadmissible_heuristic;

	int screen = 0;
	int num_of_agents;
	vector<vector<HTable> > lookupTable;

	vector<vector<int>> meta_agents;
	vector<bool> ma_vec;

	// double sum_distance_error = 0;
	// double sum_cost_error = 0;
	// int num_of_errors = 0;
	// int initialize_online_learning = 0;
    vector<double> sum_distance_errors;
    vector<double> sum_cost_errors;
    vector<int> num_of_errors;

	double time_limit;
	int node_limit = 4;  // terminate the sub CBS solver if the number of its expanded nodes exceeds the node limit.
	double start_time;
	int ILP_node_threshold = 5; // when #nodes >= ILP_node_threshold, use ILP solver; otherwise, use DP solver
	int ILP_edge_threshold = 10; // when #edges >= ILP_edge_threshold, use ILP solver; otherwise, use DP solver
	int ILP_value_threshold = 32; // when value >= ILP_value_threshold, use ILP solver; otherwise, use DP solver
	// TODO: run some experiments to pick a good ILP_node_threshold
	vector<Path*>& paths;
	const vector<SingleAgentSolver*>& search_engines;
	vector<ConstraintTable>& initial_constraints;  // This is for solve2agent function
	MDDTable& mdd_helper;

	void buildConflictGraph(vector<bool>& HG, const HLNode& curr);
	void buildCardinalConflictGraph(CBSNode& curr, vector<int>& CG, int& num_of_CGedges);
	bool buildDependenceGraph(CBSNode& node, vector<int>& CG, int& num_of_CGedges);
	bool buildWeightedDependencyGraph(CBSNode& curr, vector<int>& CG);
	bool buildWeightedDependencyGraph(ECBSNode& node, AdjEdges& WDG, vector<int>& ag_fmin, vector<int>& ag_fmax);
	bool dependent(int a1, int a2, HLNode& node); // return true if the two agents are dependent
	pair<int, int> solve2Agents(int a1, int a2, const CBSNode& node, bool cardinal); // return h value and num of CT nodes
    tuple<int, int, int> solve2Agents(int a1, int a2, const ECBSNode& node); // return h value and num of CT nodes
	static bool SyncMDDs(const MDD &mdd1, const MDD& mdd2); 	// Match and prune MDD according to another MDD.
	// void setUpSubSolver(CBS& cbs) const;
	int minimumVertexCover(const vector<int>& CG); // mvc on disjoint components
	int minimumVertexCover(const vector<int>& CG, int old_mvc, int cols, int num_of_edges); // incremental mvc
	bool KVertexCover(const vector<int>& CG, int num_of_CGnodes, int num_of_CGedges, int k, int cols);
	int greedyMatching(const vector<bool>& CG, int cols);
    static int greedyMatching(const std::vector<int>& CG,  int cols);
    int greedyWeightedMatching(const vector<int>& CG, int cols);
	int minimumWeightedVertexCover(const vector<int>& HG);
	// int minimumConstrainedWeightedVertexCover(const vector<int>& CG);
	int weightedVertexCover(const vector<int>& CG);
	int DPForWMVC(vector<int>& x, int i, int sum, const vector<int>& CG, const vector<int>& range, int& best_so_far); // dynamic programming
	// int ILPForWMVC(const vector<int>& CG, const vector<int>& range) const; // Integer linear programming
	int ILPForEWMVC(const AdjEdges& WDG, const vector<int>& ag_fmin, const vector<int>& ag_fmax);  // Integer linear programming for hyper-edge WMVC
	int ILPForConstrainedWMVC(const std::vector<int>& CG, const std::vector<int>& range);
	int DPForConstrainedWMVC(vector<bool>& x, int i, int sum, const vector<int>& CG, const vector<int>& range, int& best_so_far);
};




