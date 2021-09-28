#include "ECBS.h"


bool ECBS::solve(double time_limit, int _cost_lowerbound, int _cost_upperbound)
{
	this->cost_lowerbound = _cost_lowerbound;
	this->inadmissible_cost_lowerbound = 0;
	this->time_limit = time_limit;

	if (screen > 0) // 1 or 2
	{
		string name = getSolverName();
		name.resize(35, ' ');
		cout << name << ": ";
	}

	if (!is_start)  // set timer
	{
		is_start = true;
		start = clock();
	}

	generateRoot();

	while (!cleanup_list.empty() && !solution_found)
	{
		// Debug
		auto cleanup_head = cleanup_list.top();
		cleanup_head_lb = cleanup_list.top()->g_val;
		if (screen > 3)
		{
			open_node_idx->push_back(cleanup_head->time_generated);
			open_sum_lb->push_back(cleanup_head_lb);
			open_sum_fval->push_back(cleanup_head->getFVal());
			open_sum_cost->push_back(cleanup_head->sum_of_costs);
			open_num_conflicts->push_back(cleanup_head->conflicts.size() + cleanup_head->unknownConf.size());
			open_remained_flex->push_back(suboptimality * cleanup_head_lb - cleanup_head->sum_of_costs);
			
			if (screen == 5)  // Check the number of CT nodes in FOCAL, OPEN, and CLEANUP
			{
				iter_num_focal->push_back(focal_list.size());
				iter_num_open->push_back(open_list.size());
				iter_num_cleanup->push_back(cleanup_list.size());
			}
		}
		// End debug

		auto curr = selectNode();  // update focal list and select the CT node

		// Debug
		assert((double) curr->sum_of_costs <= suboptimality * curr->getFVal());
		assert((double) curr->sum_of_costs <= suboptimality * cleanup_head->getFVal());

		if (screen > 3)
		{
			iter_node_idx->push_back(curr->time_generated);
			iter_sum_lb->push_back(curr->g_val);
			iter_sum_fval->push_back(curr->getFVal());
			iter_sum_cost->push_back(curr->sum_of_costs);
			iter_num_conflicts->push_back(curr->conflicts.size() + curr->unknownConf.size());
			iter_remained_flex->push_back(suboptimality * curr->g_val - curr->sum_of_costs);
			iter_subopt->push_back((double) curr->sum_of_costs / (double) cleanup_head_lb);
			iter_sum_ll_generate->push_back(curr->ll_generated);

			for (const vector<int>& ma : meta_agents)
			{
				for (const int& ag : ma)
				{
					iter_ag_lb->at(ag).push_back(min_f_vals[ag]);
					iter_ag_cost->at(ag).push_back(paths[ag]->size()-1);
				}
			}
			
			if (screen == 5)  // Debug
			{
				iter_use_flex->push_back(curr->use_flex);
				iter_no_more_flex->push_back(curr->no_more_flex);
				iter_cannot_use_flex->push_back(curr->cannot_use_flex);
				if (curr->chosen_from == "cleanup")
					iter_node_type->push_back(0);
				else if (curr->chosen_from == "open")
					iter_node_type->push_back(1);
				else if (curr->chosen_from == "focal")
					iter_node_type->push_back(2);
			}
		}
		// End debug

		if (terminate(curr))
			return solution_found;

		if (use_flex && curr->chosen_from == "cleanup")  // Early replanning for FEECBS
		{
			restart_cnt ++;
			if (restart_cnt > restart_th && restart_th > -1)
			{
				restart = true;
				break;
			}
		}

		if ((curr == dummy_start || curr->chosen_from == "cleanup") &&
		     !curr->h_computed) // heuristics has not been computed yet
		{
            runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
            bool succ = heuristic_helper.computeInformedHeuristics(*curr, min_f_vals, time_limit - runtime);
            runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
            if (!succ) // no solution, so prune this node
            {
                if (screen > 1)
                    cout << "	Prune " << *curr << endl;
                curr->clear();
                continue;
            }

            if (reinsertNode(curr))
                continue;
		}

        classifyConflicts(*curr);

		//Expand the node
		num_HL_expanded++;
		curr->time_expanded = num_HL_expanded;
		if (bypass && curr->chosen_from != "cleanup")
		{
			bool foundBypass = true;
			while (foundBypass)
			{
				if (terminate(curr))
					return solution_found;
				foundBypass = false;
				ECBSNode* child[2] = { new ECBSNode() , new ECBSNode() };

				// if (screen == 3)
				// 	debugChooseConflict(*curr);
				curr->conflict = chooseConflict(*curr);  // TODO: choose conflict that is between two meta_agents

				// update conflict_matrix and joint meta_agent
				vector<int> ma1 = findMetaAgent(curr->conflict->a1);
				vector<int> ma2 = findMetaAgent(curr->conflict->a2);
				assert(ma1 != ma2);

				for (const int& a1 : ma1)
				{
					for (const int& a2 : ma2)
					{
						conflict_matrix[a1][a2] += 1;
						conflict_matrix[a2][a1] += 1;
					}
				}

				addConstraints(curr, child[0], child[1]);
				if (screen > 1)
					cout << "	Expand " << *curr << endl << 	"	on " << *(curr->conflict) << endl;

				bool solved[2] = { false, false };
				vector<Path*> path_copy(paths);
				vector<int> fmin_copy(min_f_vals);
				for (int i = 0; i < 2; i++)
				{
					if (i > 0)
					{
						paths = path_copy;
						min_f_vals = fmin_copy;
					}
					solved[i] = generateChild(child[i], curr, i);
					if (!solved[i])
					{
						delete (child[i]);
						continue;
					}
					else if (i == 1 && !solved[0])
						continue;
					else if (bypass &&
						child[i]->sum_of_costs <= suboptimality * cost_lowerbound &&
						child[i]->distance_to_go < curr->distance_to_go) // Bypass1
					{
						foundBypass = true;

						assert(curr->g_val <= child[i]->g_val);

						if (!use_flex)  // bypass condition checking for ECBS and EECBS
						{
							for (const auto& path : child[i]->paths)
							{
								//// path.first: agent id
								//// path.second.first: path
								//// path.second.second: lower bound of the path
								if ((double)path.second.first.size()-1 > suboptimality * fmin_copy[path.first])
								{
									// bypassing for EECBS
									foundBypass = false;
									break;
								}
							}
						}
						else  // bypass condition checking for FECBS and FEECBS
						{
							// if (curr->g_val < child[i]->g_val)  // do not bypass if the sum of lower bound increases
							// 	foundBypass = false;
							if (child[i]->sum_of_costs > suboptimality * curr->g_val)
								foundBypass = false;
						}

						if (foundBypass)
						{
							adoptBypass(curr, child[i], fmin_copy, path_copy);
							if (screen > 1)
								cout << "	Update " << *curr << endl;
							break;  // do not generate another child
						}
					}
				}
				if (foundBypass)
				{
					for (auto & i : child)
					{
						delete i;
					}
                    classifyConflicts(*curr); // classify the new-detected conflicts
				}

				else if (shouldMerge(ma1, ma2)) // Should merge
				{
					clock_t t1 = clock();
					// Update new meta-agent to global variable and curr CT node
					meta_agents.erase(std::remove(meta_agents.begin(), meta_agents.end(), ma1), meta_agents.end());
					meta_agents.erase(std::remove(meta_agents.begin(), meta_agents.end(), ma2), meta_agents.end());
					vector<int> joint_ma;
					std::merge(ma1.begin(), ma1.end(), ma2.begin(), ma2.end(), joint_ma.begin());
					meta_agents.push_back(joint_ma);

					if (screen > 1)
					{
						cout << "Merging agents: ";
						for (int _tmp_ag_ : joint_ma)
						{
							cout << _tmp_ag_ << ", ";
						}
						cout << endl;
					}

					bool foundPaths = findPathForMetaAgent(curr, joint_ma);
					joint_ma.clear();

					if (foundPaths)
					{
						// Copy new meta_agents to the current node
						curr->meta_agents = meta_agents;
						curr->ma_vec =ma_vec;
						curr->is_merged = true;

						findConflicts(*curr);
						heuristic_helper.computeQuickHeuristics(*curr);
						pushNode(curr);

						// // Debug
						// cout << "allNnodes_table: ";
						// for (auto tmp_node : allNodes_table)
						// {
						// 	cout << *tmp_node << ", ";
						// }
						// cout << endl;

						assert(curr->g_val <= curr->sum_of_costs);
						assert(curr->sum_of_costs <= suboptimality * curr->g_val);

						runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
					}

					else if (inner_solver->solution_cost == -1)  // Timeout
					{
						// solution_cost = -1;
						runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
						// runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
						continue;
					}

					else
					{
						runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
						if (screen == 3)
							cout << "Failed: No path for Meta-Agent!!!" << endl;
						continue;
					}
				}  // Should merge end

				else  // expansion
				{
					for (int i = 0; i < 2; i++)
					{
						if (solved[i])
						{
							pushNode(child[i]);
							curr->children.push_back(child[i]);
							if (screen > 1)
							{
								cout << "		Generate " << *child[i] << endl;
							}
							if (screen > 3)
							{
								all_node_idx->push_back(child[i]->time_generated);
								all_sum_lb->push_back(child[i]->g_val);
								all_sum_fval->push_back(child[i]->getFVal());
								all_sum_cost->push_back(child[i]->sum_of_costs);
								all_num_conflicts->push_back(child[i]->conflicts.size() + child[i]->unknownConf.size());
								all_remained_flex->push_back(suboptimality * child[i]->getFVal() - child[i]->sum_of_costs);
								all_subopt->push_back((double) child[i]->sum_of_costs / (double) child[i]->getFVal());
								all_sum_ll_generate->push_back(child[i]->ll_generated);
							}
						}
					}
				}
			}
		}
		else // no bypass
		{
			ECBSNode* child[2] = { new ECBSNode() , new ECBSNode() };
			curr->conflict = chooseConflict(*curr);

			// update conflict_matrix and joint meta_agent
			vector<int> ma1 = findMetaAgent(curr->conflict->a1);
			vector<int> ma2 = findMetaAgent(curr->conflict->a2);
			// if (ma1 == ma2)
			// {
			// 	for (int tmp_ag : ma1)
			// 		printAgentPath(tmp_ag);
			// }
			assert(ma1 != ma2);

			for (const int& a1 : ma1)
			{
				for (const int& a2 : ma2)
				{
					conflict_matrix[a1][a2] += 1;
					conflict_matrix[a2][a1] += 1;
				}
			}

			if (shouldMerge(ma1, ma2)) // Should merge
			{
				clock_t t1 = clock();
				// Update new meta-agent to global variable and curr CT node
				meta_agents.erase(std::remove(meta_agents.begin(), meta_agents.end(), ma1), meta_agents.end());
				meta_agents.erase(std::remove(meta_agents.begin(), meta_agents.end(), ma2), meta_agents.end());
				vector<int> joint_ma;
				std::merge(ma1.begin(), ma1.end(), ma2.begin(), ma2.end(), std::back_inserter(joint_ma));
				meta_agents.push_back(joint_ma);

				// if (screen > 1)
				// {
				// 	cout << "Merging agents: ";
				// 	for (int _tmp_ag_ : joint_ma)
				// 		cout << _tmp_ag_ << ", ";
				// 	cout << endl;
				// 	cout << "collected constraints:" << endl;
				// 	for (Constraint tmp_constraint : curr->constraints)
				// 	{
				// 		cout << "< a:" << get<0>(tmp_constraint) << ", v1:" << get<1>(tmp_constraint);
				// 		cout << ", v2:" << get<2>(tmp_constraint) << ", t:" << get<3>(tmp_constraint);
				// 		cout << ", " << get<4>(tmp_constraint) << ">" << endl;
				// 	}
				// }

				bool foundPaths = findPathForMetaAgent(curr, joint_ma);
				// cout << "Paths" << endl;
				// for (int tmp_ag : joint_ma)
				// {
				// 	printAgentPath(tmp_ag);
				// }
				// cout << "----------------------------------" << endl;
				joint_ma.clear();

				if (foundPaths)
				{
					// Copy new meta_agents to the current node
					curr->meta_agents = meta_agents;
					curr->ma_vec = ma_vec;
					curr->is_merged = true;

					findConflicts(*curr);
					// if (screen > 1)
					// {
					// 	cout << "After replanning the meta-agent" << endl;
					// 	printConflicts(*curr);
					// 	cout << endl;
					// }
					heuristic_helper.computeQuickHeuristics(*curr);
					pushNode(curr);

					// // Debug
					// cout << "push node done!" << endl;
					// cout << "allNnodes_table:" << endl;
					// for (auto tmp_node : allNodes_table)
					// {
					// 	cout << "Address: " << tmp_node << ", " << *tmp_node << endl;
					// }
					// cout << endl;

					assert(curr->g_val <= curr->sum_of_costs);
					assert(curr->sum_of_costs <= suboptimality * curr->g_val);

					runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
				}

				else if (inner_solver->solution_cost == -1)  // Timeout
				{
					runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
					// runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
					continue;
				}

				else
				{
					runtime_merge += (double)(clock() - t1) / CLOCKS_PER_SEC;
					if (screen == 3)
						cout << "Failed: No path for Meta-Agent!!!" << endl;
					continue;
				}
			}  // Should merge end

			else
			{
				addConstraints(curr, child[0], child[1]);

				if (screen > 1)
					cout << "	Expand " << *curr << endl << "	on " << *(curr->conflict) << endl;

				bool solved[2] = { false, false };
				vector<vector<PathEntry>*> path_copy(paths);
				vector<int> fmin_copy(min_f_vals);
				for (int i = 0; i < 2; i++)
				{
					if (i > 0)
					{
						paths = path_copy;
						min_f_vals = fmin_copy;
					}
					solved[i] = generateChild(child[i], curr, i);
					if (!solved[i])
					{
						delete (child[i]);
						continue;
					}
					pushNode(child[i]);
					curr->children.push_back(child[i]);
					if (screen > 1)
						cout << "		Generate " << *child[i] << endl;
					if (screen > 3)
					{
						all_node_idx->push_back(child[i]->time_generated);
						all_sum_lb->push_back(child[i]->g_val);
						all_sum_fval->push_back(child[i]->getFVal());
						all_sum_cost->push_back(child[i]->sum_of_costs);
						all_num_conflicts->push_back(child[i]->conflicts.size() + child[i]->unknownConf.size());
						all_remained_flex->push_back(suboptimality * child[i]->getFVal() - child[i]->sum_of_costs);
						all_subopt->push_back((double) child[i]->sum_of_costs / (double) child[i]->getFVal());
						all_sum_ll_generate->push_back(child[i]->ll_generated);
					}
				}
			}
		}
		switch (curr->conflict->type)
		{
		case conflict_type::RECTANGLE:
			num_rectangle_conflicts++;
			break;
		case conflict_type::CORRIDOR:
			num_corridor_conflicts++;
			break;
		case  conflict_type::TARGET:
			num_target_conflicts++;
			break;
		case conflict_type::STANDARD:
			num_standard_conflicts++;
			break;
		case conflict_type::MUTEX:
			num_mutex_conflicts++;
			break;
		default:
			break;
		}
		if (curr->chosen_from == "cleanup")
			num_cleanup++;
		else if (curr->chosen_from == "open")
			num_open++;
		else if (curr->chosen_from == "focal")
			num_focal++;

		if (curr->conflict->priority == conflict_priority::CARDINAL)
			num_cardinal_conflicts++;
        if (!curr->children.empty())
            heuristic_helper.updateOnlineHeuristicErrors(*curr); // update online heuristic errors

		if (!curr->is_merged)  // Only clear conflicts if curr is not reinserted to lists (due to merging)
			curr->clear();
	}  // end of while loop

	// Restart from FEECBS to EECBS
	clock_t curr_t = clock();
	if (restart && runtime < time_limit)
	{
		vector<pair<Path,int>> backup_initial_paths = paths_found_initially;  // We do not want to clear initial paths
		clear();

		assert(use_flex);
		use_flex = false;
		srand(0);
		bool debug_foundSol = solve(time_limit, 0);
		assert(debug_foundSol == solution_found);
	}

	return solution_found;
}

void ECBS::adoptBypass(ECBSNode* curr, ECBSNode* child, const vector<int>& fmin_copy, const vector<Path*>& path_copy)
{
	num_adopt_bypass++;
	curr->sum_of_costs = child->sum_of_costs;
	curr->cost_to_go = child->cost_to_go;
	curr->distance_to_go = child->distance_to_go;
	curr->conflicts = child->conflicts;
	curr->unknownConf = child->unknownConf;
	curr->conflict = nullptr;
	curr->makespan = child->makespan;
	for (const auto& path : child->paths) // update paths
	{
		auto p = curr->paths.begin();
		while (p != curr->paths.end())
		{
			if (path.first == p->first)  // if path is already stored in curr node
			{
				assert(fmin_copy[path.first] == p->second.second);
				assert(path_copy[path.first]->size() - 1 == p->second.first.size() - 1);
				assert(path_copy[path.first] == &p->second.first);

				p->second.first = path.second.first;  // replace the path in curr node with path in child CT node
				paths[p->first] = &p->second.first;  // update the new path to global paths
                min_f_vals[p->first] = p->second.second;  // update back the old fmin to global
				break;
			}
			++p;
		}
		if (p == curr->paths.end())  // if path is not stored in curr node (e.g. root)
		{
			curr->paths.emplace_back(path);
			curr->paths.back().second.second = fmin_copy[path.first];
			paths[path.first] = &curr->paths.back().second.first;
			min_f_vals[path.first] = fmin_copy[path.first];  // update back to curr fmin
		}
	}
}

// takes the paths_found_initially and UPDATE all (constrained) paths found for agents from curr to start
// also, do the same for ll_min_f_vals and paths_costs (since its already "on the way").
void ECBS::updatePaths(ECBSNode* curr)
{
	for (int i = 0; i < num_of_agents; i++)
	{
		paths[i] = &paths_found_initially[i].first;
		min_f_vals[i] = paths_found_initially[i].second;
	}
	vector<bool> updated(num_of_agents, false);  // initialized for false

	while (curr != nullptr)
	{
		for (auto & path : curr->paths)
		{
			int agent = path.first;
			if (!updated[agent])
			{
				paths[agent] = &path.second.first;
				min_f_vals[agent] = path.second.second;
				updated[agent] = true;
			}
		}
		curr = curr->parent;
	}
}


bool ECBS::generateRoot()
{
	if (meta_agents.empty())  // initialize for outer ECBS
	{
		for (int ag = 0; ag < num_of_agents; ag ++)
		{
			meta_agents.push_back(vector<int>({ag}));
			ma_vec[ag] = true;
		}
	}

	paths.resize(num_of_agents, nullptr);
	init_min_f_vals.resize(num_of_agents);
	mdd_helper.init(num_of_agents);
	heuristic_helper.init();
	if (!path_initialize)
	{
		// initialize paths_found_initially
		paths_found_initially.clear();
		paths_found_initially.resize(num_of_agents);
		min_f_vals.clear();
		min_f_vals.resize(num_of_agents);
	}

	int initial_g_val = 0;
	int initial_soc = 0;
	vector<vector<int>> init_meta_agent = meta_agents;
	for (const vector<int>& ma : meta_agents)
	{
		for (const int& ag : ma)
		{
			if (path_initialize)
			{
				// already initialize min_f_val at findPathForMetaAgent
				assert(search_engines[ag]->my_heuristic[search_engines[ag]->start_location] <= min_f_vals[ag]);
				initial_soc += paths_found_initially[ag].first.size() - 1;
			}
			else
			{ 
				min_f_vals[ag] = search_engines[ag]->my_heuristic[search_engines[ag]->start_location];
			}

			initial_g_val += min_f_vals[ag];
		}
	}

	vector<vector<int>> sort_based;
	vector<bool> sort_ascend;
	vector<int> fmin_based = vector<int>(meta_agents.size(), 0);
	vector<int> conf_based = vector<int>(meta_agents.size(), 0);
	for (size_t ma_id = 0; ma_id < meta_agents.size(); ma_id++)
		for (const int& tmp_ag : meta_agents[ma_id])
			fmin_based[ma_id] += min_f_vals[tmp_ag];

	if (random_init)
		std::random_shuffle(meta_agents.begin(), meta_agents.end());  // generate random permutation of agent indices
	else
		sortMetaAgents(fmin_based, fmin_ascend, conf_based, conf_ascend);  // sort according to the low-level heuristics

	// // debug
	// cout << endl;
	// for (const vector<int>& ma : meta_agents)
	// {
	// 	for (const int& ag : ma)
	// 	{
	// 		cout << "agent: " << ag << "-> " << fmin_based[ag] << endl;
	// 	}
	// }

	auto root = new ECBSNode();
	root->g_val = initial_g_val;
	root->sum_of_costs = (path_initialize)? initial_soc : initial_g_val;
	root->meta_agents = meta_agents;
	root->ma_vec = ma_vec;

	bool first_time = true;
	while (true)
	{
		if (path_initialize)  // Already initialize paths at outer solver
			break;

		auto candidate = new ECBSNode();
		candidate->g_val = initial_g_val;
		candidate->sum_of_costs = initial_g_val;
		candidate->meta_agents = meta_agents;
		candidate->ma_vec = ma_vec;

		vector<pair<Path, int>> candidate_paths;
		candidate_paths.resize(num_of_agents);

		for (const vector<int>& ma : meta_agents)
		{
			if (ma.size() > 1)  // this is a meta agent
			{
				// TODO: Add mr_active
				// shuffle agents in the meta-agent first
				assert(mr_active);
				continue;
			}
			else  // this is a single agent
			{
				int ag = ma.front();
				if (candidate_paths[ag].first.empty())
				{
					candidate_paths[ag] = search_engines[ag]->findSuboptimalPath(*candidate, initial_constraints[ag], paths, ag, 0, suboptimality);
					// if (use_flex)
					// {
					// 	candidate_paths[ag] = search_engines[ag]->findSuboptimalPath(*candidate, initial_constraints[ag], paths, ag, 0, suboptimality, 
					// 		candidate->g_val - min_f_vals[ag], candidate->sum_of_costs - min_f_vals[ag], init_sum_lb, flex);
					// }
					// else
					// {
					// 	candidate_paths[ag] = search_engines[ag]->findSuboptimalPath(*candidate, initial_constraints[ag], paths, ag, 0, suboptimality);
					// }

					num_LL_expanded += search_engines[ag]->num_expanded;
					num_LL_generated += search_engines[ag]->num_generated;
				}
								
				if (candidate_paths[ag].first.empty())
				{
					cout << "No path exists for agent " << ag << endl;
					return false;
				}

				paths[ag] = &candidate_paths[ag].first;
				candidate->makespan = max(candidate->makespan, candidate_paths[ag].first.size() - 1);
				candidate->g_val = candidate->g_val + max(candidate_paths[ag].second - min_f_vals[ag], 0);
				candidate->sum_of_costs = candidate->sum_of_costs + ((int)paths[ag]->size() - 1) - min_f_vals[ag];
				min_f_vals[ag] = max(candidate_paths[ag].second, min_f_vals[ag]);
				candidate_paths[ag].second = min_f_vals[ag];
			}
		}

		findConflicts(*candidate);

		if (first_time || candidate->unknownConf.size() < root->unknownConf.size())
		{
			first_time = (first_time)? false: first_time;
			root = candidate;
			paths_found_initially = candidate_paths;

			if (root_replan && root->unknownConf.size() > 0)
			{
				vector<int> conf_num = vector<int>(num_of_agents, 0);
				for (const auto& conf : candidate->unknownConf)
				{
					conf_num[conf->a1] ++;
					conf_num[conf->a2] ++;
				}

				fill(paths.begin(), paths.end(), nullptr);
				meta_agents = init_meta_agent;
				conf_based = vector<int>(meta_agents.size(), 0);
				for (size_t ma_id = 0; ma_id < meta_agents.size(); ma_id++)
					for (const int& tmp_ag : meta_agents[ma_id])
						conf_based[ma_id] += conf_num[tmp_ag];
				sortMetaAgents(fmin_based, fmin_ascend, conf_based, conf_ascend);

				// // debug
				// cout << endl;
				// for (const vector<int>& ma : meta_agents)
				// {
				// 	for (const int& ag : ma)
				// 	{
				// 		cout << "agent: " << ag << "-> " << conf_based[ag] << ", " << fmin_based[ag] << endl;
				// 	}
				// }
				// // end debug

				// Initialize min_f_vals
				for (const vector<int>& ma : meta_agents)
					for (const int& ag : ma)
						min_f_vals[ag] = search_engines[ag]->my_heuristic[search_engines[ag]->start_location];
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}

	if (!path_initialize)  // Outer
	{
		for (const vector<int>& ma : meta_agents)
		{
			for (const int& ag : ma)
			{
				paths[ag] = &paths_found_initially[ag].first;
				min_f_vals[ag] = paths_found_initially[ag].second;
				init_min_f_vals[ag] = paths_found_initially[ag].second;
			}
		}
	}
	else  // Inner or restart from root
	{
		for (int ag = 0; ag < num_of_agents; ag++)
		{
			paths[ag] = &paths_found_initially[ag].first;
			init_min_f_vals[ag] = min_f_vals[ag];
			paths_found_initially[ag].second = min_f_vals[ag];
			root->makespan = max(root->makespan, paths[ag]->size() - 1);
		}
		findConflicts(*root);
	}
	
	path_initialize = true;

	if (screen > 3)
		root->ll_generated = num_LL_generated;

	root->h_val = 0;
	root->depth = 0;
    heuristic_helper.computeQuickHeuristics(*root);
	pushNode(root);
	dummy_start = root;

	if (screen >= 2) // print start and goals
	{
		printPaths();
		// printConflicts(*dummy_start);
	}
	return true;
}


bool ECBS::generateChild(ECBSNode*  node, ECBSNode* parent, int child_idx)
{
	clock_t t1 = clock();
	node->parent = parent;
	node->HLNode::parent = parent;
	node->g_val = parent->g_val;
	node->h_val = parent->h_val;
	node->sum_of_costs = parent->sum_of_costs;
	node->makespan = parent->makespan;
	node->meta_agents = parent->meta_agents;
	node->ma_vec = parent->ma_vec;
	node->depth = parent->depth + 1;
	auto agents = getInvalidAgents(node->constraints);
	assert(!agents.empty());
	for (auto agent : agents)
	{
		vector<int> invalid_ma = findMetaAgent(agent);
		if (invalid_ma.size() == 1)
		{
			if (!findPathForSingleAgent(node, agent))
			{
				if (screen > 1)
					cout << "	No paths for agent " << agent << ". Node pruned." << endl;
				runtime_generate_child += (double)(clock() - t1) / CLOCKS_PER_SEC;
				return false;
			}
		}
		else  // Replan the meta-agent
		{
			if (!findPathForMetaAgent(node, invalid_ma))
			{
				if (screen > 1)
					cout << "	No paths for agent " << agent << ". Node pruned." << endl;
				runtime_generate_child += (double)(clock() - t1) / CLOCKS_PER_SEC;
				return false;
			}
		}
	}
	assert(node->sum_of_costs <= suboptimality * node->getFVal());
	assert(parent->getFVal() <= node->getFVal());
	assert(parent->g_val <= node->g_val);

	findConflicts(*node);
	heuristic_helper.computeQuickHeuristics(*node);

	updateConflictImpacts(*node, *parent, child_idx);

	runtime_generate_child += (double)(clock() - t1) / CLOCKS_PER_SEC;
	return true;
}


bool ECBS::findPathForSingleAgent(ECBSNode*  node, int ag)
{
	clock_t t = clock();
	pair<Path, int> new_path;
	if (use_flex)
	{
		if (screen > 1)
		{
			assert((double) node->sum_of_costs <= suboptimality * node->getFVal());
			int tmp_cost = 0;
			int tmp_lb = 0;
			for (const vector<int>& ma : meta_agents)
			{
				for (const int& ag : ma)
				{
					tmp_cost += (int) paths[ag]->size() - 1;
					tmp_lb += min_f_vals[ag];
				}
			}
			assert(tmp_cost == node->sum_of_costs);
			assert(tmp_lb == node->g_val);
		}

		int other_sum_lb = node->g_val - min_f_vals[ag];
		int other_sum_cost = node->sum_of_costs - (int) paths[ag]->size() + 1;

		new_path = search_engines[ag]->findSuboptimalPath(*node, initial_constraints[ag], paths, ag, min_f_vals[ag],
			suboptimality, other_sum_lb, other_sum_cost, init_sum_lb, flex, node->h_val);

		// // Flex restriction
		// bool not_use_flex;
		// if (node == dummy_start || node->parent == dummy_start)
		// {
		// 	not_use_flex = true;
		// }
		// else
		// {
		// 	not_use_flex = node->parent->chosen_from == "cleanup" || 
		// 		node->parent->conflict->priority == conflict_priority::CARDINAL || 
		// 		node->parent->g_val < node->g_val;
		// }

		// if (suboptimality* (double) other_sum_lb - (double) other_sum_cost >= 0 && not_use_flex)
		// {
		// 	// Not use flex if the CT node is from cleanup or the conflict is cardinal
		// 	new_path = search_engines[ag]->findSuboptimalPath(*node, initial_constraints[ag], paths, ag, min_f_vals[ag], suboptimality);

		// 	node->use_flex = false;
		// 	if (not_use_flex)
		// 		node->cannot_use_flex = true;
		// }

		// else
		// {
		// 	new_path = search_engines[ag]->findSuboptimalPath(*node, initial_constraints[ag], paths, ag, min_f_vals[ag],
		// 		suboptimality, other_sum_lb, other_sum_cost, init_sum_lb, flex, node->h_val);

		// 	node->use_flex = true;
		// 	if (suboptimality* (double) other_sum_lb - (double) other_sum_cost < 0)
		// 		node->no_more_flex = true;
		// }
		// // End Flex restrictions
	}
	else
	{
		new_path = search_engines[ag]->findSuboptimalPath(*node, initial_constraints[ag], paths, ag, min_f_vals[ag], suboptimality);
	}

	num_LL_expanded += search_engines[ag]->num_expanded;
	num_LL_generated += search_engines[ag]->num_generated;
	runtime_build_CT += search_engines[ag]->runtime_build_CT;
	runtime_build_CAT += search_engines[ag]->runtime_build_CAT;
	runtime_path_finding += (double)(clock() - t) / CLOCKS_PER_SEC;
	if (new_path.first.empty())
		return false;
	assert(!isSamePath(*paths[ag], new_path.first));

	node->g_val = node->g_val + max(new_path.second - min_f_vals[ag], 0);
	if (node->parent != nullptr)
	    node->h_val = max(0, node->parent->getFVal() - node->g_val); // pathmax

	node->sum_of_costs = node->sum_of_costs - (int) paths[ag]->size() + (int) new_path.first.size();
	min_f_vals[ag] = max(new_path.second, min_f_vals[ag]);  // make sure the recorded lower bound is always the maximum
	new_path = make_pair(new_path.first, min_f_vals[ag]);
	node->paths.emplace_back(ag, new_path);
	paths[ag] = &node->paths.back().second.first;
	node->makespan = max(node->makespan, new_path.first.size() - 1);

	assert(node->sum_of_costs <= suboptimality * node->getFVal());
	assert(node->getFVal() >= node->parent->getFVal());

	return true;
}

// Plan a path for the meta agent in a child node
// Collect constraints outside the function. 
bool ECBS::findPathForMetaAgent(ECBSNode*  node, const vector<int>& meta_ag)
{
	if (inner_solver == nullptr)
	{
		cerr << "Failed: No solver for MetaAgent!" << endl;
		exit(-1);
	}
	else
	{
		inner_solver->clear();
	}

	solver_counter ++;

	// Set constraint to inner solver
	vector<bool> _ma_vec_ = vector<bool>(num_of_agents, false);
	int outer_lb = 0;
	for (const int& ag : meta_ag)
	{
		_ma_vec_[ag] = true;
		outer_lb += min_f_vals[ag];  // Determine sum of fmin of the meta-agent
		ConstraintTable _constraint_table;
		_constraint_table.init(initial_constraints[ag]);
		_constraint_table.build(*node, ag);
		inner_solver->setInitConstraints(ag, _constraint_table);
	}

	inner_solver->setMetaAgents(meta_ag);  // Set the meta-agent for inner solver
	inner_solver->setMAVector(_ma_vec_);
	inner_solver->setInitSumLB(outer_lb);

	double outer_flex = 0.0;
	for (const vector<int>& _ma_ : meta_agents)  // Set paths and outer flex for inner solver
	{
		for (const int& _ag_ : _ma_)
		{
			if (paths[_ag_] != nullptr)
			{
				if (!_ma_vec_[_ag_] && use_flex)  // Add flex from other agents outside the meta-agent
					outer_flex += suboptimality * min_f_vals[_ag_] - (paths[_ag_]->size() - 1);
				inner_solver->setInitialPath(_ag_, *paths[_ag_]);  // Initialize paths for inner solver
				inner_solver->setMinFVal(_ag_, min_f_vals[_ag_]);  // Initialize min_f_val for inner solver
			}
			else
			{
				if (!_ma_vec_[_ag_] && use_flex)  // Add flex from other agents outside the meta-agent
					outer_flex += suboptimality * min_f_vals[_ag_] - min_f_vals[_ag_];
			}
			
		}
	}
	// if (screen == 2)
	// {
	// 	cout << "--------------------------------------" << endl;
	// 	cout << "conflict of original CT ndoe" << endl;
	// 	printConflicts(*node);
	// 	cout << "--------------------------------------" << endl;
	// 	cout << "chosen onflict: " << *node->conflict << endl;
	// }
	// if (screen > 3)
	// {
	// 	cout << inner_solver->getInitialPathLength(0) << endl;
	// 	cout << endl;
	// }
	inner_solver->setFlex(outer_flex);
	inner_solver->setIsStart(false);
	inner_solver->path_initialize = true;
	clock_t t = clock();
	double inner_time_limit = time_limit - ((t-start) / CLOCKS_PER_SEC);
	bool foundSol = inner_solver->solve(inner_time_limit, outer_lb);  // Run solver for meta-agent

	// Update Outer ECBS parameters from Inner ECBS
	// cout << "inner solver done!" << endl;
	solver_num_HL_expanded += inner_solver->num_HL_expanded;
	solver_num_HL_generated += inner_solver->num_HL_generated;
	solver_num_LL_expanded += inner_solver->num_LL_expanded;
	solver_num_LL_generated += inner_solver->num_LL_generated; 

	runtime_solver += inner_solver->runtime;

	if (foundSol)
	{
		int old_sum_of_costs = 0;
		int new_sum_lb = 0;
		for (const int& ag : meta_ag)
		{
			old_sum_of_costs += (int) paths[ag]->size() - 1;
			new_sum_lb += inner_solver->getminFVal(ag);

			for (list<pair<int, pair<Path, int>>>::iterator p_it = node->paths.begin(); p_it != node->paths.end();)
			{
				if (ag == p_it->first)  // Delete repeated paths in node
				{
					paths[ag] = nullptr;
					p_it = node->paths.erase(p_it ++);
				}
				else
				{
					++ p_it;
				}
			}
		}

		node->g_val = node->g_val + max(new_sum_lb - outer_lb, 0);
		node->sum_of_costs = node->sum_of_costs - old_sum_of_costs + inner_solver->solution_cost;

		for (const int& ag : meta_ag)  // Update the paths and the lower bounds of the meta-agent
		{
			Path inner_path = inner_solver->getPath(ag);
			int inner_lb = inner_solver->getminFVal(ag);
			int ag_lb = (new_sum_lb < outer_lb)? min_f_vals[ag] : inner_lb;
			node->paths.emplace_back(ag, make_pair(inner_path, ag_lb));
			paths[ag] = &node->paths.back().second.first;
			node->makespan = max(node->makespan, inner_path.size() - 1);
		}
		
		inner_solver->setInitSumLB(0);  // Reset sum of fmin
		inner_solver->setFlex(0);  // Reset flex

		runtime_path_finding += (double)(clock() - t) / CLOCKS_PER_SEC;
		return true;
	}

	else if (inner_solver->runtime > inner_time_limit)  // Timeout in solver!!!
	{
		assert(inner_solver->solution_cost == -1);
		runtime_path_finding += (double)(clock() - t) / CLOCKS_PER_SEC;
		return false;
	}

	else  // No solution for the meta-agent
	{
		delete node;
		runtime_path_finding += (double)(clock() - t) / CLOCKS_PER_SEC;
		return false;
	}
}


inline void ECBS::pushNode(ECBSNode* node)
{
	if (node->time_generated == 0)
	{
		num_HL_generated++;
		node->time_generated = num_HL_generated;
	}

	// update handles
    node->cleanup_handle = cleanup_list.push(node);
	switch (solver_type)
	{
	case high_level_solver_type::ASTAREPS:  // cleanup_list is called open_list in ECBS
        if (node->sum_of_costs <= suboptimality * cost_lowerbound)
		{
			num_push_focal++;  // debug
            node->focal_handle = focal_list.push(node);
		}
		break;
	case high_level_solver_type::NEW:
		if (node->getFHatVal() <= suboptimality * cost_lowerbound)
		{
			node->focal_handle = focal_list.push(node);
			num_push_focal++;  // debug
		}
		break;
	case high_level_solver_type::EES:
		node->open_handle = open_list.push(node);
		if (node->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
		{
			node->focal_handle = focal_list.push(node);
			num_push_focal++;  // debug
		}
		break;
	case high_level_solver_type::CLEANUP:
		node->open_handle = open_list.push(node);
		if (node->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
		{
			node->focal_handle = focal_list.push(node);
			num_push_focal++;  // debug
		}
		break;
	case high_level_solver_type::CLEANUP_ASTAREPS:
		if (node->sum_of_costs <= suboptimality * cost_lowerbound)
		{
			num_push_focal++;  // debug
            node->focal_handle = focal_list.push(node);
		}
		break;
	default:
		break;
	}
	if (!node->isin_allNodesTable)
	{
		node->isin_allNodesTable = true;
		allNodes_table.push_back(node);
	}
	if (node->is_merged)
	{
		mergedNodes_table.try_emplace(node->time_generated, node);
	}
}


inline bool ECBS::reinsertNode(ECBSNode* node)
{

	switch (solver_type)
	{
	case high_level_solver_type::ASTAREPS:
        if (node->sum_of_costs <= suboptimality * cost_lowerbound)
            return false;
        node->cleanup_handle = cleanup_list.push(node);
        break;
	case high_level_solver_type::NEW:
        if (node->getFHatVal() <= suboptimality * cost_lowerbound)
            return false;
        node->cleanup_handle = cleanup_list.push(node);
		break;
	case high_level_solver_type::EES:
        node->cleanup_handle = cleanup_list.push(node);
		node->open_handle = open_list.push(node);
		if (node->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
			node->focal_handle = focal_list.push(node);
		break;
	case high_level_solver_type::CLEANUP:
        node->cleanup_handle = cleanup_list.push(node);
		node->open_handle = open_list.push(node);
		if (node->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
			node->focal_handle = focal_list.push(node);
		break;
	case high_level_solver_type::CLEANUP_ASTAREPS:
		if (node->sum_of_costs <= suboptimality * cost_lowerbound)
            return false;
        node->cleanup_handle = cleanup_list.push(node);
        break;
	default:
		break;
	}
	if (screen == 2)
	{
		cout << "	Reinsert " << *node << endl;
	}
	return true;
}


ECBSNode* ECBS::selectNode()
{
	ECBSNode* curr = nullptr;
	assert(solver_type != high_level_solver_type::ASTAR);
	switch (solver_type)
	{
	case high_level_solver_type::EES:
		// update the focal list if necessary
		if (open_list.top()->getFHatVal() != inadmissible_cost_lowerbound)
		{
			inadmissible_cost_lowerbound = open_list.top()->getFHatVal();
			double focal_list_threshold = suboptimality * inadmissible_cost_lowerbound;
			focal_list.clear();
			for (auto n : open_list)
			{
				if (n->getFHatVal() <= focal_list_threshold)
					n->focal_handle = focal_list.push(n);
			}
		}

		// choose the best node
		if (screen > 1 && cleanup_list.top()->getFVal() > cost_lowerbound)
			cout << "Lowerbound increases from " << cost_lowerbound << " to " << cleanup_list.top()->getFVal() << endl;
		cost_lowerbound = max(cleanup_list.top()->getFVal(), cost_lowerbound);
		if (focal_list.top()->sum_of_costs <= suboptimality * cost_lowerbound)
		{ // return best d
			curr = focal_list.top();
			curr->chosen_from = "focal";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go; */
			focal_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
			open_list.erase(curr->open_handle);
		}
		else if (open_list.top()->sum_of_costs <= suboptimality * cost_lowerbound)
		{ // return best f_hat
			curr = open_list.top();
			curr->chosen_from = "open";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			open_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
			focal_list.erase(curr->focal_handle);
		}
		else
		{ // return best f
			curr = cleanup_list.top();
			curr->chosen_from = "cleanup";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			cleanup_list.pop();
			open_list.erase(curr->open_handle);
			if (curr->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
				focal_list.erase(curr->focal_handle);
		}
		break;
	case high_level_solver_type::ASTAREPS:
		// update the focal list if necessary
		if (cleanup_list.top()->getFVal() > cost_lowerbound)
		{
			if (screen == 3)
			{
				cout << "  Note -- FOCAL UPDATE!! from |FOCAL|=" << focal_list.size() << " with |OPEN|=" << cleanup_list.size() << " to |FOCAL|=";
			}
			double old_focal_list_threshold = suboptimality * cost_lowerbound;
			cost_lowerbound = max(cost_lowerbound, cleanup_list.top()->getFVal());
			double new_focal_list_threshold = suboptimality * cost_lowerbound;
			for (auto n : cleanup_list)
			{
				if (n->sum_of_costs > old_focal_list_threshold && n->sum_of_costs <= new_focal_list_threshold)
					n->focal_handle = focal_list.push(n);
			}
			if (screen == 3)
			{
				cout << focal_list.size() << endl;
			}
		}

		// choose best d in the focal list
		curr = focal_list.top();
		curr->chosen_from = "focal";
		/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
		curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
		curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
		curr->f_of_best_in_focal = focal_list.top()->getFVal();
		curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
		curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
		focal_list.pop();
		cleanup_list.erase(curr->cleanup_handle);
		break;
	case high_level_solver_type::NEW:
		// update the focal list if necessary
		if (cleanup_list.top()->getFVal() > cost_lowerbound)
		{
			if (screen == 3)
			{
				cout << "  Note -- FOCAL UPDATE!! from |FOCAL|=" << focal_list.size() << " with |OPEN|=" << cleanup_list.size() << " to |FOCAL|=";
			}
			double old_focal_list_threshold = suboptimality * cost_lowerbound;
			cost_lowerbound = max(cost_lowerbound, cleanup_list.top()->getFVal());
			double new_focal_list_threshold = suboptimality * cost_lowerbound;
			focal_list.clear();
			for (auto n : cleanup_list)
			{
				heuristic_helper.updateInadmissibleHeuristics(*n);
				if (n->getFHatVal() <= new_focal_list_threshold)
					n->focal_handle = focal_list.push(n);
			}
			if (screen == 3)
			{
				cout << focal_list.size() << endl;
			}
		}

		if (focal_list.empty()) // choose best f in the cleanup list (to improve the lower bound)
		{
			curr = cleanup_list.top();
			curr->chosen_from = "cleanup";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;*/
			cleanup_list.pop();
		}
		else // choose best d in the focal list
		{
			curr = focal_list.top();
			curr->chosen_from = "focal";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			focal_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
		}
		break;
	case high_level_solver_type::CLEANUP:
		// update the focal list if necessary
		if (open_list.top()->getFHatVal() != inadmissible_cost_lowerbound)
		{
			inadmissible_cost_lowerbound = open_list.top()->getFHatVal();
			double focal_list_threshold = suboptimality * inadmissible_cost_lowerbound;
			focal_list.clear();
			for (auto n : open_list)
			{
				if (n->getFHatVal() <= focal_list_threshold)
					n->focal_handle = focal_list.push(n);
			}
		}

		// choose the best node
		if (screen > 1 && cleanup_list.top()->getFVal() > cost_lowerbound)
			cout << "Lowerbound increases from " << cost_lowerbound << " to " << cleanup_list.top()->getFVal() << endl;
		
		// Update node counter node_cnt
		if (cleanup_list.top()->getFVal() > cost_lowerbound)
			node_cnt = 0;
		else
			node_cnt ++;
	
		cost_lowerbound = max(cleanup_list.top()->getFVal(), cost_lowerbound);

		if (node_cnt < cleanup_th)  // Keep expanding from CLEANUP to increase the min f value
		{
			curr = cleanup_list.top();
			curr->chosen_from = "cleanup";
			cleanup_list.pop();
			open_list.erase(curr->open_handle);
			if (curr->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
				focal_list.erase(curr->focal_handle);
		}
		else if (focal_list.top()->sum_of_costs <= suboptimality * cost_lowerbound)
		{ // return best d
			if (screen > 1)
			{
				cout << "node_cnt: " << node_cnt << endl;
				cout << "cleanup_list.top()->getFVal(): " << cleanup_list.top()->getFVal() << endl;
				cout << "cost_lowerbound: " << cost_lowerbound << endl;
			}
			curr = focal_list.top();
			curr->chosen_from = "focal";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go; */
			focal_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
			open_list.erase(curr->open_handle);
		}
		else if (open_list.top()->sum_of_costs <= suboptimality * cost_lowerbound)
		{ // return best f_hat
			curr = open_list.top();
			curr->chosen_from = "open";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			open_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
			focal_list.erase(curr->focal_handle);
		}
		else
		{ // return best f
			curr = cleanup_list.top();
			curr->chosen_from = "cleanup";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_open = open_list.top()->getFVal();
			curr->f_hat_of_best_in_open = open_list.top()->getFHatVal();
			curr->d_of_best_in_open = open_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			cleanup_list.pop();
			open_list.erase(curr->open_handle);
			if (curr->getFHatVal() <= suboptimality * inadmissible_cost_lowerbound)
				focal_list.erase(curr->focal_handle);
		}
		break;
	case high_level_solver_type::CLEANUP_ASTAREPS:
		// update the focal list if necessary
		if (cleanup_list.top()->getFVal() > cost_lowerbound)
		{
			if (screen == 3)
			{
				cout << "  Note -- FOCAL UPDATE!! from |FOCAL|=" << focal_list.size() << " with |OPEN|=" << cleanup_list.size() << " to |FOCAL|=";
			}
			double old_focal_list_threshold = suboptimality * cost_lowerbound;

			node_cnt = 0;  // Update node counter node_cnt
			cost_lowerbound = max(cost_lowerbound, cleanup_list.top()->getFVal());
			double new_focal_list_threshold = suboptimality * cost_lowerbound;
			for (auto n : cleanup_list)
			{
				if (n->sum_of_costs > old_focal_list_threshold && n->sum_of_costs <= new_focal_list_threshold)
					n->focal_handle = focal_list.push(n);
			}
			if (screen == 3)
			{
				cout << focal_list.size() << endl;
			}
		}
		else
		{
			node_cnt ++;
		}
		
	
		if (node_cnt < cleanup_th)  // Keep expanding from CLEANUP to increase the min f value
		{
			curr = cleanup_list.top();
			curr->chosen_from = "cleanup";
			cleanup_list.pop();
			if (curr->getFVal() <= suboptimality * cost_lowerbound)
				focal_list.erase(curr->focal_handle);
		}
		else
		{
			// choose best d in the focal list
			curr = focal_list.top();
			curr->chosen_from = "focal";
			/*curr->f_of_best_in_cleanup = cleanup_list.top()->getFVal();
			curr->f_hat_of_best_in_cleanup = cleanup_list.top()->getFHatVal();
			curr->d_of_best_in_cleanup = cleanup_list.top()->distance_to_go;
			curr->f_of_best_in_focal = focal_list.top()->getFVal();
			curr->f_hat_of_best_in_focal = focal_list.top()->getFHatVal();
			curr->d_of_best_in_focal = focal_list.top()->distance_to_go;*/
			focal_list.pop();
			cleanup_list.erase(curr->cleanup_handle);
		}
		
		break;
	default:
		break;
	}

	// takes the paths_found_initially and UPDATE all constrained paths found for agents from curr to dummy_start (and lower-bounds)
	updatePaths(curr);
	meta_agents = curr->meta_agents;

	if (screen > 1)
		cout << endl << "Pop " << *curr << endl;
	return curr;
}

void ECBS::printPaths() const
{
	for (int i = 0; i < num_of_agents; i++)
	{
		cout << "Agent " << i << " (" << paths_found_initially[i].first.size() - 1 << " -->" <<
			paths[i]->size() - 1 << "): ";
		for (const auto & t : *paths[i])
			cout << t.location << "->";
		cout << endl;
	}
}

void ECBS::printAgentPath(int ag) const
{
	cout << "Agent " << ag << " (" << paths_found_initially[ag].first.size() - 1 << " -->" <<
		paths[ag]->size() - 1 << "): ";
	for (const auto & t : *paths[ag])
		cout << t.location << "->";
	cout << endl;
}

void ECBS::classifyConflicts(ECBSNode &node)
{
    if (node.unknownConf.empty())
        return;
	// Classify all conflicts in unknownConf
	while (!node.unknownConf.empty())
	{
		shared_ptr<Conflict> con = node.unknownConf.front();
		int a1 = con->a1, a2 = con->a2;
		int timestep = get<3>(con->constraint1.back());
		constraint_type type = get<4>(con->constraint1.back());
		node.unknownConf.pop_front();

		if (PC)
		    if (node.chosen_from == "cleanup" ||
               // (min_f_vals[a1] * suboptimality >= min_f_vals[a1] + 1 &&
               //min_f_vals[a2] * suboptimality >= min_f_vals[a2] + 1))
               (int)paths[a1]->size() - 1 == min_f_vals[a1] ||
               (int)paths[a2]->size() - 1 == min_f_vals[a2]) // the path of at least one agent is its shortest path
			    computeConflictPriority(con, node);

		/*if (con->priority == conflict_priority::CARDINAL && heuristic_helper.type == heuristics_type::ZERO)
		{
			computeSecondPriorityForConflict(*con, node);
			node.conflicts.push_back(con);
			return;
		}*/

		// TODO: Mutex reasoning for ECBS
		/*if (mutex_reasoning)
		{
			// TODO mutex reasoning is per agent pair, don't do duplicated work...
			auto mdd1 = mdd_helper.getMDD(node, a1, paths[a1]->size());
			auto mdd2 = mdd_helper.getMDD(node, a2, paths[a2]->size());

			auto mutex_conflict = mutex_helper.run(a1, a2, node, mdd1, mdd2);

			if (mutex_conflict != nullptr)
			{
				computeSecondPriorityForConflict(*mutex_conflict, node);
				node.conflicts.push_back(mutex_conflict);
				continue;
			}
		}*/

		// Target Reasoning: highest secondary priority
		if (con->type == conflict_type::TARGET)
		{
			computeSecondPriorityForConflict(*con, node);
			node.conflicts.push_back(con);
			continue;
		}

		// Corridor reasoning: second secondary priority
		if (corridor_reasoning)
		{
			auto corridor = corridor_helper.run(con, paths, node);
			if (corridor != nullptr)
			{
				corridor->loc1 = con->loc1;
				corridor->loc2 = con->loc2;
				corridor->priority = con->priority;
				computeSecondPriorityForConflict(*corridor, node);
				node.conflicts.push_back(corridor);
				continue;
			}
		}


		// Rectangle reasoning: less secondary priority
		if (rectangle_reasoning &&
		    (int)paths[a1]->size() - 1 == min_f_vals[a1] && // the paths for both agents are their shortest paths
		    (int)paths[a2]->size() - 1 == min_f_vals[a2] &&
            min_f_vals[a1] > timestep &&  //conflict happens before both agents reach their goal locations
			min_f_vals[a2] > timestep &&
			// findMetaAgent(a1).size() == 1 &&  // maybe not needed if we change the constraint after finding out the agents are meta-agent
			// findMetaAgent(a2).size() == 1 &&
			type == constraint_type::VERTEX) // vertex conflict
		{
			// TODO: memorize the conflicts location for meta-agent to view it as vertex conflict
			auto mdd1 = mdd_helper.getMDD(node, a1, paths[a1]->size());
			auto mdd2 = mdd_helper.getMDD(node, a2, paths[a2]->size());
			auto rectangle = rectangle_helper.run(paths, timestep, a1, a2, mdd1, mdd2);
			if (rectangle != nullptr)
			{
                if (!PC)
                    rectangle->priority = conflict_priority::UNKNOWN;
				computeSecondPriorityForConflict(*rectangle, node);
				rectangle->loc1 = con->loc1;
				rectangle->loc2 = con->loc2;
				node.conflicts.push_back(rectangle);
				continue;
			}
		}

		computeSecondPriorityForConflict(*con, node);
		node.conflicts.push_back(con);
	}

	// remove conflicts that cannot be chosen, to save some memory
	removeLowPriorityConflicts(node.conflicts);
}


void ECBS::computeConflictPriority(shared_ptr<Conflict>& con, ECBSNode& node)
{
    int a1 = con->a1, a2 = con->a2;
	int timestep = get<3>(con->constraint1.back());
	constraint_type type = get<4>(con->constraint1.back());
	bool cardinal1 = false, cardinal2 = false;
	MDD *mdd1 = nullptr, *mdd2 = nullptr;
	if (timestep >= (int)paths[a1]->size())
		cardinal1 = true;
	else //if (!paths[a1]->at(0).is_single())
	{
		mdd1 = mdd_helper.getMDD(node, a1, paths[a1]->size());
	}
	if (timestep >= (int)paths[a2]->size())
		cardinal2 = true;
	else //if (!paths[a2]->at(0).is_single())
	{
		mdd2 = mdd_helper.getMDD(node, a2, paths[a2]->size());
	}

	if (type == constraint_type::EDGE) // Edge conflict
	{
		if (timestep < (int)mdd1->levels.size())
		{
			cardinal1 = mdd1->levels[timestep].size() == 1 &&
				mdd1->levels[timestep].front()->location == paths[a1]->at(timestep).location &&
				mdd1->levels[timestep - 1].size() == 1 &&
				mdd1->levels[timestep - 1].front()->location == paths[a1]->at(timestep - 1).location;
		}
		if (timestep < (int)mdd2->levels.size())
		{
			cardinal2 = mdd2->levels[timestep].size() == 1 &&
				mdd2->levels[timestep].front()->location == paths[a2]->at(timestep).location &&
				mdd2->levels[timestep - 1].size() == 1 &&
				mdd2->levels[timestep - 1].front()->location == paths[a2]->at(timestep - 1).location;
		}
	}
	else // vertex conflict or target conflict
	{
		if (!cardinal1 && timestep < (int)mdd1->levels.size())
		{
			cardinal1 = mdd1->levels[timestep].size() == 1 &&
				mdd1->levels[timestep].front()->location == paths[a1]->at(timestep).location;
		}
		if (!cardinal2 && timestep < (int)mdd2->levels.size())
		{
			cardinal2 = mdd2->levels[timestep].size() == 1 &&
				mdd2->levels[timestep].front()->location == paths[a2]->at(timestep).location;
		}
	}

	if (cardinal1 && cardinal2)
	{
		con->priority = conflict_priority::CARDINAL;
	}
	else if (cardinal1 || cardinal2)
	{
		con->priority = conflict_priority::SEMI;
	}
	else
	{
		con->priority = conflict_priority::NON;
	}
}


// used for rapid random  restart
void ECBS::clear()
{
    mdd_helper.clear();
    heuristic_helper.clear();
    paths.clear();
    paths_found_initially.clear();
    min_f_vals.clear();
	init_min_f_vals.clear();

    open_list.clear();
    cleanup_list.clear();
    focal_list.clear();
    for (auto& node : allNodes_table)
        delete node;
    allNodes_table.clear();

	path_initialize = false;
    dummy_start = nullptr;
    goal_node = nullptr;
    solution_found = false;
    solution_cost = -2;

	// for nested framework
	meta_agents.clear();
	ma_vec.clear();
	ma_vec.resize(num_of_agents, false);  // checking if need to solve agent

	if (is_solver)
	{
		num_HL_generated = 0;
		num_HL_expanded = 0;
		num_LL_generated = 0;
		num_LL_expanded = 0;
	}

	conflict_matrix.clear();
	conflict_matrix.resize(num_of_agents, vector<int>(num_of_agents, 0));
}