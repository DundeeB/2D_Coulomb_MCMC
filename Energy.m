function [E] = Energy(state,i)
[N,~] = size(state.spheres);
E = 0;
p1 = state.spheres(i,1:2);
q = @(z) 2*(z>state.H/2)-1;
q1 = q(state.spheres(i,3));
for j=[1:i-1 i+1:N]
    p2 = state.spheres(j,1:2);
    q2 = q(state.spheres(j,3));
    for i_p=1:length(p1)
        dr = cyclic_dist(p1,p2,state.cyclic_boundary);
        if dr<2*state.rad
            E = inf;
            return
        else
            E = E + q1*q2/dr;
        end
    end
end
end