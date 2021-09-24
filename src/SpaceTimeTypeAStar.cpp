#include "SpaceTimeTypeAStar.h"


Path SpaceTimeTypeAStar::findOptimalPath(const HLNode& node, const ConstraintTable& initial_constraints,
	const vector<Path*>& paths, int agent, int lowerbound)
{
	return findSuboptimalPath(node, initial_constraints, paths, agent, lowerbound, 1).first;
}


// find path by time-space A* search
// Returns a bounded-suboptimal path that satisfies the constraints of the give node  while
// minimizing the number of internal conflicts (that is conflicts with known_paths for other agents found so far).
// lowerbound is an underestimation of the length of the path in order to speed up the search.
pair<Path, int> SpaceTimeTypeAStar::findSuboptimalPath(const HLNode& node, const ConstraintTable& initial_constraints,
	const vector<Path*>& paths, int agent, int lowerbound, double w, int other_sum_lb, int other_sum_cost,
	int outer_sum_lb, double single_flex, int hl_h_val)
{
	this->w = w;
	this->use_focal = true;
	Path path;
	num_expanded = 0;
	num_generated = 0;

	// build constraint table
	auto t = clock();
	constraint_table.init(initial_constraints);
	constraint_table.build(node, agent);
	runtime_build_CT = (double)(clock() - t) / CLOCKS_PER_SEC;
	if (constraint_table.constrained(start_location, 0))
	{
		return {path, 0};
	}

	int holding_time = constraint_table.getHoldingTime(); // the earliest timestep that the agent can hold its goal location. The length_min is considered here.
	t = clock();
	constraint_table.buildCAT(agent, paths, node.makespan + 1);
	runtime_build_CAT = (double)(clock() - t) / CLOCKS_PER_SEC;

    lowerbound =  max(holding_time, lowerbound);

	// generate start and add it to the OPEN & FOCAL list
	auto start = new AStarNode(start_location, 0, max(lowerbound, my_heuristic[start_location]), nullptr, 0, 0, false);

	num_generated++;
	start->open_handle = open_list.push(start);
	start->focal_handle = focal_list.push(start);
	start->in_openlist = true;
	allNodes_table.insert(start);
	min_f_val = (int) start->getFVal();

	pushNodeToTypes(start);

	// upperbound = w * (max(min_f_val, lowerbound + hl_h_val) + other_sum_lb) - other_sum_cost;
	upperbound = w * max(outer_sum_lb, min_f_val + other_sum_lb) - other_sum_cost + single_flex;
	assert(min_f_val <= upperbound);
	assert(my_heuristic[start_location] <= upperbound);

	while (!open_list.empty())
	{
		// update FOCAL if min f-val increases
		updateFocalList(lowerbound, other_sum_lb, other_sum_cost, outer_sum_lb, single_flex, hl_h_val);
		auto* curr = popNode();
        assert(curr->location >= 0);
		// check if the popped node is a goal
		if (curr->location == goal_location && // arrive at the goal location
			!curr->wait_at_goal && // not wait at the goal location
			curr->timestep >= holding_time) // the agent can hold the goal location afterward
		{
			updatePath(curr, path);
			break;
		}

		// this is for range constraint: agent should arrive at or before timestep length_max
		if (curr->timestep >= constraint_table.length_max)
			continue;

		auto next_locations = instance.getNeighbors(curr->location);
		next_locations.emplace_back(curr->location);
		for (int next_location : next_locations)
		{
			int next_timestep = curr->timestep + 1;
			if (max((int)node.makespan, constraint_table.latest_timestep) + 1 < curr->timestep)
			{ // now everything is static, so switch to space A* where we always use the same timestep
				if (next_location == curr->location)
				{
					continue;
				}
				next_timestep--;
			}

			if (constraint_table.constrained(next_location, next_timestep) ||
				constraint_table.constrained(curr->location, next_location, next_timestep))
				continue;

			// compute cost to next_id via curr node
			int next_g_val = curr->g_val + 1;
			int next_h_val = max(lowerbound - next_g_val, my_heuristic[next_location]);
			if (next_g_val + next_h_val > constraint_table.length_max)
				continue;
			int next_internal_conflicts = curr->num_of_conflicts +
				constraint_table.getNumOfConflictsForStep(curr->location, next_location, next_timestep);

			// generate (maybe temporary) node
			auto next = new AStarNode(next_location, next_g_val, next_h_val,
				curr, next_timestep, next_internal_conflicts, false);
			if (next_location == goal_location && curr->location == goal_location)
				next->wait_at_goal = true;

			// try to retrieve it from the hash table
			auto it = allNodes_table.find(next);
			if (it == allNodes_table.end())
			{
				pushNode(next);
				allNodes_table.insert(next);
				continue;
			}
			// update existing node's if needed (only in the open_list)

			auto existing_next = *it;
			if (existing_next->getFVal() > next->getFVal() || // if f-val decreased through this new path
				(existing_next->getFVal() == next->getFVal() &&
					existing_next->num_of_conflicts > next->num_of_conflicts)) // or it remains the same but there's fewer conflicts
			{
				if (!existing_next->in_openlist) // if its in the closed list (reopen)
				{
					existing_next->copy(*next);
					pushNode(existing_next);
				}
				else
				{
					bool add_to_focal = false;  // check if it was above the focal bound before and now below (thus need to be inserted)
					bool update_in_focal = false;  // check if it was inside the focal and needs to be updated (because f-val changed)
					bool update_open = false;
					if ((next_g_val + next_h_val) <= upperbound && use_focal)
					{  // if the new f-val qualify to be in FOCAL
						if (existing_next->getFVal() > upperbound)
							add_to_focal = true;  // and the previous f-val did not qualify to be in FOCAL then add
						else
							update_in_focal = true;  // and the previous f-val did qualify to be in FOCAL then update
					}
					if (existing_next->getFVal() > next_g_val + next_h_val)
						update_open = true;

					existing_next->copy(*next);	// update existing node

					if (update_open)
						open_list.increase(existing_next->open_handle);  // increase because f-val improved
					if (add_to_focal)
					{
						existing_next->focal_handle = focal_list.push(existing_next);
						pushNodeToTypes(existing_next);
						assert(getNumInTypes() == focal_list.size());
						assert(isTypeValid());
					}
					if (update_in_focal)
					{
						focal_list.update(existing_next->focal_handle);  // should we do update? yes, because number of conflicts may go up or down
						type_lists[(int) existing_next->getFVal()].update(existing_next->type_handle);
						assert(getNumInTypes() == focal_list.size());
						assert(isTypeValid());
					}
				}
			}
			delete(next);  // not needed anymore -- we already generated it before
		}  // end for loop that generates successors
		iter_counter ++;
	}  // end while loop

	releaseNodes();
	constraint_table.clear();
	return {path, min_f_val};
}


inline AStarNode* SpaceTimeTypeAStar::popNode()
{
	AStarNode* node;
	if (use_focal)
	{
		if (iter_counter % 2 == 1)  // Focal search steop
		{
			cout << "====================================================" << endl;
			cout << "Before: " << endl;
			cout << "type_lists.size(): " << type_lists.size() << endl;
			for (auto itr = type_lists.begin(); itr != type_lists.end(); itr++)
			{
				cout << "\ntype_lists[" << itr->first << "]: " << endl;
				for (auto tmp_n : itr->second)
				{
					cout << "\t" << tmp_n << tmp_n->getFVal() << ", " << tmp_n->location << ", " << tmp_n->num_of_conflicts << endl;
				}
				cout << "---------------------------------" << endl;
			}

			AStarNode *node = focal_list.top(); focal_list.pop();
			type_lists[(int) node->getFVal()].erase(node->type_handle);

			cout << "After: " << endl;
			cout << "type_lists.size(): " << type_lists.size() << endl;
			for (auto itr = type_lists.begin(); itr != type_lists.end(); itr++)
			{
				cout << "\ntype_lists[" << itr->first << "]: " << endl;
				for (auto tmp_n : itr->second)
				{
					cout << "\t" << tmp_n->getFVal() << ", " << tmp_n->location << ", " << tmp_n->num_of_conflicts << endl;
				}
				cout << "---------------------------------" << endl;
			}

			assert(getNumInTypes() == focal_list.size());
			assert(isTypeValid());
		}
		else  // Type-based focal search
		{
			cout << "type_lists.size(): " << type_lists.size() << endl;
			for (auto itr = type_lists.begin(); itr != type_lists.end(); itr++)
			{
				cout << "\ntype_lists[" << itr->first << "]: " << endl;
				for (auto tmp_n : itr->second)
				{
					cout << "\t" << tmp_n->getFVal() << ", " << tmp_n->location << ", " << tmp_n->num_of_conflicts << endl;
				}
				cout << "---------------------------------" << endl;
			}

			// Select the node with minimum number of conflicts in the given random type
			cout << "Before focal_list.size(): " << focal_list.size() << endl;
			assert(getNumInTypes() == focal_list.size());
			assert(isTypeValid());

			auto rand_it = type_lists.begin();
			advance(rand_it, rand() % type_lists.size());

			// Debug
			cout << "rand_it: (" << rand_it->first << ", " << rand_it->second.top() << ")" << endl;
			cout << "rand_it->second.size(): " << rand_it->second.size() << endl;

			// Select node
			node = rand_it->second.top();
			rand_it->second.pop();

			cout << "node->getFVal()   " << node->getFVal() << endl;
			cout << "After rand_it->second.size(): " << rand_it->second.size() << endl;

			// Erase node pointer from type_lists and focal_list
			if (rand_it->second.empty())
				type_lists.erase(rand_it->first);
			focal_list.erase(node->focal_handle);

			cout << "After focal_list.size(): " << focal_list.size() << endl;

			// Debug
			assert(getNumInTypes() == focal_list.size());
			assert(isTypeValid());
		}

		open_list.erase(node->open_handle);
	}
	else
	{
		node = open_list.top(); open_list.pop();
	}
	node->in_openlist = false;
	num_expanded++;
	return node;
}


inline void SpaceTimeTypeAStar::pushNode(AStarNode* node)
{
	node->open_handle = open_list.push(node);
	node->in_openlist = true;
	num_generated++;

	if (num_generated > nl_ratio*node_limit && nl_ratio > 0 && w > 1 && upperbound > w * min_f_val)
		use_focal = false;
	else if (node->getFVal() <= upperbound)
	{
		node->focal_handle = focal_list.push(node);		
		pushNodeToTypes(node);
	}

	// Debug
	cout << "++++++++++ push node ++++++++++" << endl;
	cout << "focal_list.size(): " << focal_list.size() << endl;
	for (auto t_n : focal_list)
	{
		cout << t_n << endl;
	}
	for (auto tmp_type : type_lists)
	{
		cout << "tmp_type.second.size(): " << tmp_type.second.size();
		for (auto t_n : tmp_type.second)
			cout << ", " << tmp_type.second.top();
		cout << endl;
	}
	assert(getNumInTypes() == focal_list.size());
	assert(isTypeValid());
}

inline void SpaceTimeTypeAStar::pushNodeToTypes(AStarNode* node)
{
	int fval = (int) node->getFVal();
	if (type_lists.find(fval) == type_lists.end())  // Type of the node is not found
	{
		heap_focal_t tmp_focal;
		node->type_handle = tmp_focal.push(node);
		type_lists.insert(make_pair(fval, tmp_focal));
	}
	else  // The key already exists
	{
		node->type_handle = type_lists[fval].push(node);
	}
}


int SpaceTimeTypeAStar::getNumInTypes(void)
{
	int number_nodes_in_type_lists = 0;
	for (auto tmp_type : type_lists)
		number_nodes_in_type_lists += tmp_type.second.size();
	return number_nodes_in_type_lists;
}


bool SpaceTimeTypeAStar::isTypeValid(void)
{
	for (auto tmp_type : type_lists)
	{
		for (AStarNode *tmp_node : tmp_type.second)
		{
			if (find(focal_list.begin(), focal_list.end(), tmp_node) == focal_list.end())
			{
				cout << tmp_node << " not found in FOCAL" << endl;
				return false;
			}
		}
	}
	return true;
}


void SpaceTimeTypeAStar::updateFocalList(int lowerbound, int other_sum_lb, int other_sum_cost, int outer_sum_lb, double single_flex, int hl_h_val)
{
	auto open_head = open_list.top();
	if (open_head->getFVal() > min_f_val)
	{
		int new_min_f_val = (int)open_head->getFVal();

		double new_upper_bound;  // Get new_upper_bound
		new_upper_bound = w * max(outer_sum_lb, new_min_f_val + other_sum_lb) - other_sum_cost + single_flex;
		assert(new_min_f_val <= new_upper_bound);
		if (use_focal)  // Update focal list
		{
			for (auto n : open_list)
			{
				if (n->getFVal() >  upperbound && n->getFVal() <= new_upper_bound)
				{
					n->focal_handle = focal_list.push(n);
					pushNodeToTypes(n);
				}
			}
		}
		min_f_val = new_min_f_val;
		upperbound = new_upper_bound;
		assert((double) min_f_val <= upperbound);
	}
}

void SpaceTimeTypeAStar::releaseNodes()
{
	open_list.clear();
	focal_list.clear();
	type_lists.clear();
	for (auto node: allNodes_table)
		delete node;
	allNodes_table.clear();
}