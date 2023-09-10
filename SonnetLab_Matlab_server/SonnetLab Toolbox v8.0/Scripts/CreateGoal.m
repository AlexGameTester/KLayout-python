function goal = CreateGoal(enabled, response_goal, relation, target_goal, weight)
%createGoal(enabled, response_goal, relation, target_goal, weight) 
%Creates a goal to be used in a goal set via GoalSetBlock.addGoal
%   response_goal and target_goal should be made with 
%   CreateResponse(...) or the SonnetGoal*(...)
%   constructors (which CreateResponse uses).
%		(SonnetGoalValue, SonnetGoalNetwork, SonnetGoalFile)
%   weight defaults to 1.0
%   Example usage:
%       response = CreateResponse('Network', 'DB', 'S', [1 1])
%       target = CreateResponse('Value', 20)
%       goal = CreateGoal('Y', response, '<', target, 1.0)
if nargin == 4
    weight = 1;
end
goal = SonnetGoal(enabled, response_goal, relation, target_goal, weight);
end

