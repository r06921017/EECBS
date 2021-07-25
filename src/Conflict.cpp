#include "Conflict.h"
#include "RectangleReasoning.h"
#include "MDD.h"


std::ostream& operator<<(std::ostream& os, const Constraint& constraint)
{
	os << "<" << std::get<0>(constraint) << "," << std::get<1>(constraint) << "," <<
		std::get<2>(constraint) << "," << std::get<3>(constraint) << ",";
	switch (get<4>(constraint))
	{
		case constraint_type::VERTEX:
			os << "V";
			break;
		case constraint_type::POSITIVE_VERTEX:
			os << "V+";
			break;
		case constraint_type::EDGE:
			os << "E";
			break;
		case constraint_type::POSITIVE_EDGE:
			os <<"E+";
			break;
		case constraint_type::BARRIER:
			os << "B";
			break;
		case constraint_type::RANGE:
			os << "R";
			break;
		case constraint_type::GLENGTH:
			os << "G";
			break;
		case constraint_type::LEQLENGTH:
			os << "L";
			break;
	}
	os << ">";
	return os;
}


std::ostream& operator<<(std::ostream& os, const Conflict& conflict)
{
	switch (conflict.priority)
	{
		case conflict_priority::CARDINAL:
			os << "cardinal ";
			break;
		case conflict_priority::PSEUDO_CARDINAL:
			os << "pseudo-cardinal ";
			break;
		case conflict_priority::SEMI:
			os << "semi-cardinal ";
			break;
		case conflict_priority::NON:
			os << "non-cardinal ";
			break;
        case conflict_priority::PRIORITY_COUNT:
            break;
    }
	switch (conflict.type)
	{
		case conflict_type::STANDARD:
			os << "standard";
			break;
		case conflict_type::RECTANGLE:
			os << "rectangle";
			break;
		case conflict_type::CORRIDOR:
			os << "corrdior";
			break;
		case conflict_type::TARGET:
			os << "target";
			break;
		case conflict_type::MUTEX:
			os << "mutex";
			break;
	    case conflict_type::TYPE_COUNT:
            break;
    }
	os << " conflict:  " << conflict.a1 << " with ";
	for (auto con : conflict.constraint1)
		os << con << ",";		
	os << " and " << conflict.a2 << " with ";
	for (auto con : conflict.constraint2)
		os << con << ",";		
	return os;
}

bool operator < (const Conflict& conflict1, const Conflict& conflict2) // return true if conflict2 has higher priority
{
	if (conflict1.priority == conflict2.priority)
	{
		if (conflict1.type == conflict2.type)
		{
			if (conflict1.secondary_priority == conflict2.secondary_priority)
			{
				// Get minimum increased_flex
				double impact1 = max(conflict1.getImpactVal(0, impact_type::FLEX), conflict1.getImpactVal(1, impact_type::FLEX));
				double impact2 = max(conflict2.getImpactVal(0, impact_type::FLEX), conflict2.getImpactVal(1, impact_type::FLEX));

				if (impact1 == impact2)
				{
					impact1 = max(conflict1.getImpactVal(0, impact_type::LB), conflict1.getImpactVal(1, impact_type::LB));
					impact2 = max(conflict2.getImpactVal(0, impact_type::LB), conflict2.getImpactVal(1, impact_type::LB));
					if (impact1 == impact2)
					{
						impact1 = max(conflict1.getImpactVal(0, impact_type::REDUCED_CONFLICTS), conflict1.getImpactVal(1, impact_type::REDUCED_CONFLICTS));
						impact2 = max(conflict2.getImpactVal(0, impact_type::REDUCED_CONFLICTS), conflict2.getImpactVal(1, impact_type::REDUCED_CONFLICTS));
						if (impact1 == impact2)
						{
							if (max(conflict1.impacts[0].count, conflict1.impacts[1].count) == max(conflict2.impacts[0].count, conflict2.impacts[1].count))
							{
								return rand() % 2;
							}
							return max(conflict1.impacts[0].count, conflict1.impacts[1].count) > max(conflict2.impacts[0].count, conflict2.impacts[1].count);
						}
						return impact1 < impact2;
					}
					return impact1 < impact2;
				}
				return impact1 < impact2;
			}
			return conflict1.secondary_priority > conflict2.secondary_priority;
		}
		return conflict1.type > conflict2.type;
	}
	return conflict1.priority > conflict2.priority;
}

double Conflict::getImpactVal(int child_idx, int _mode_) const
{
	if (impacts[child_idx].count > 0)
	{
		switch (_mode_)
		{
		case impact_type::FLEX:
			return impacts[child_idx].increased_flex / (double) impacts[child_idx].count;
			break;
		case impact_type::LB:
			return (double) impacts[child_idx].increased_lb / (double) impacts[child_idx].count;
			break;
		case impact_type::REDUCED_CONFLICTS:
			return (double) impacts[child_idx].reduced_num_conflicts / (double) impacts[child_idx].count;
			break;
		default:
			break;
		}
	}
	else
	{
		assert(impacts[child_idx].count == 0);
		return 0;
	}
}