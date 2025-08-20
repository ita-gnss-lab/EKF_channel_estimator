function [stateAPriori] = propagate_state(stateVector, modelTransitionMatrix)
%   PROPAGATE_STATE 
%   calculates the state a priori given the current state and its model
%   transition matrix

stateAPriori = modelTransitionMatrix * stateVector;

end

