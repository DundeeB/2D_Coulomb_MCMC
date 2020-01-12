clear all; units;
cd('C:\Users\Daniel Abutbul\OneDrive - Technion\2D Coulomb Naive');
addpath ../'3D Metropolis Monte Carlo'/;
T = [1e4*Kelvin];  % beta = beta*K*q^2, where exp(-beta*H)=exp(-beta*Energy) and energy is 1/r
I = ones(1, length(T));
eta_arr = 0.9*I;  % eta = N*sig^2/A
n_row_arr = 10*I;
n_col_arr = 10*I;

N_real  = 5e4;  % 1e4*N/9*5;  % TBD!
N_save = 1e2;  % 1e4;  % TBD!
N_start = 2e2;  % 1e4;  % TBD!
f = 1;  % factor step size

code_dir = pwd;
simulations_dir = 'C:\Users\Daniel Abutbul\OneDrive - Technion\simulation-results-Coulomb\';

for j = 1:length(n_col_arr)
    tic;
    n_row = n_row_arr(j);
    n_col = n_col_arr(j);
    
    N = n_row * n_col;
    beta = e^2/(k_B*T(j));  % exp(-H/kT)=exp(-e^2/r^2*1/kT) = exp(-beta/r^2) for beta=e^2/kT
    state.rad = r_ion;  % hard sphere rejection still exist
    A = N*(2*state.rad)^2/eta_arr(j);
    
    state.H = 4*state.rad;  % z>H/2 -> charge = +1, z< H/2-> chage = -1
    state.beta = beta;
    state.cyclic_boundary = sqrt(A)*[1 1];
    state.spheres = antiferro_rect_starting_cond3D([n_col n_row 1],...
        [state.cyclic_boundary state.H],state.rad);
    sim_name = ['N=' num2str(N) '_eta=' num2str(eta_arr(j)) '_T=' num2str(T(j))]
    title_name = ['N=' num2str(N) ' \eta=' num2str(eta_arr(j)) ' T=' num2str(T(j)) 'K'];
    %%
    addpath('.');
    cd(simulations_dir); mkdir(sim_name); cd(sim_name);
    %%
    a = sqrt(A/N)-2*state.rad;
    assert(a>0, "Too many spheres!");
    step_size = f*a;
    %%
    save('Input_parameters');
    %%
    q = 0;
    for i=1:N_real
        t = rand*2*pi;
        i_p = randi(N);
        [state, q_] = metropolis_step(state, i_p, step_size*[cos(t) sin(t) 0]);
        q = q + q_;
        if mod(i,N_real/10)==0
            disp([num2str(i/N_real*100) '%']);
        end

        if (mod(i,N_save) == 0 && i > N_start) || i==N_real
            dlmwrite(num2str(i),state.spheres,'\t');
        end
    end
    %% write first and last sphere
    subplot(1,2,2); 
    plot_circles(state); 
    title([num2str(i) ' step. Acceptance rate: ' ...
        num2str(100*q/i) '%, steps per sphere: ' num2str(q/N)]);
    load('Input_parameters');
    subplot(1,2,1); 
    plot_circles(state); 
    title(title_name);
    savefig('Starting_and_Final_configuration.fig')
%     close all;
    %%
    cd(code_dir);
    toc;
    
    tic;
    post_process([simulations_dir sim_name],true, 'output_psi14');
    toc;
end