function [] = post_process(lib, plot_flag, output_file_name, N_realizations)
display(['Post processing for library: ' lib]);
switch nargin
    case 4
        psi14 = psi_post_process_for_lib(lib, plot_flag, N_realizations);
%         [b, M, N_sp] = M_frustration_post_proccess_for_lib(lib, plot_flag, ...
%             N_realizations);
    case 3
        psi14= psi_post_process_for_lib(lib, plot_flag);
%         [b, M, N_sp] = M_frustration_post_proccess_for_lib(lib, plot_flag);
end
save([lib '\' output_file_name],'psi14');
end