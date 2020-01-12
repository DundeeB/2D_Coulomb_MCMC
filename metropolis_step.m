function [state, accept_flag] = metropolis_step(state, i, dr)
orig = state.spheres;
E0 = Energy(state, i);
new_state = state;
new_state.spheres(i,:) = cyclic(state.spheres(i,:) + dr, state.cyclic_boundary);
E1 = Energy(new_state, i);
accept_flag = 0;
p = min(1, exp(-state.beta*(E1-E0)));
if rand < p  % accept step if E1-E0<0, or randomally accept it if E1-E0>0
    state = new_state;
    accept_flag = 1;
end    
end

