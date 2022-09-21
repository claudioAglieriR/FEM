clear all

%% Third homework - Claudio Aglieri Rinella - 2058086

%% --------------------------------------------------------------

format long

maxit= 5500;
n_mesh= 0:4;

% Creating the vectors associated with errors of ichol and jacobi pcg
epsilon_ichol= zeros(length(n_mesh),1);
eps_ichol_ratio= ones(length(n_mesh),1);

epsilon_jacob= zeros(length(n_mesh),1);
eps_jacob_ratio= ones(length(n_mesh),1);

% Computing all the necessary values for each mesh thanks to a for cycle

for i=1:1:length(n_mesh)

    n_file= n_mesh(i);

    % Loading the input data
    topol= load(sprintf('mesh%i/mesh%i.topol', n_file, n_file));
    coord= load(sprintf('mesh%i/mesh%i.coord', n_file, n_file));
    bound= load(sprintf('mesh%i/mesh%i.bound', n_file, n_file));

    % Storing the coordinates given by the matrix "coord" in two vectors
    x= coord(:,1);
    y= coord(:,2);

    % Computing the surface measure relative to each node
    [~, s_node]= surfaces(topol, coord);

    % Computing the values that the exact solution takes on the nodes
    u_exact= zeros(length(coord(:,1)), 1);

    for j=1:1:length(coord(:,1))

        u_exact(j)= x(j)^2 + y(j)^2 - (x(j)^2*y(j)^2) - 1;
    end

    % Computing the stiffness matrix
    H= stiffness_mat(topol, coord);
    n= size(H, 2);

    % Computing the right-hand side vector
    f= right_hand_side(coord, s_node);

    % Boundary conditions enforcement
    H= bound_condition(H, bound);

    % Computing the tolerance and assigning the initial solution
    tol= (1e-8)*(norm(f, 2));
    x0= zeros(n,1);

    % Computing the ichol preconditioner and solving the linear system with pcg
    L= ichol(H);
    [uh_ichol, flag_ichol(i), relres_ichol(i), iter_ichol(i), resvec_ichol{i}]= pcg(H, f, tol, maxit, L, L.', x0);

    % Computing the jacobi preconditioner and solving the linear system with pcg
    M= sparse(diag(diag(H)));
    [uh_jacob, flag_jacob(i), relres_jacob(i), iter_jacob(i), resvec_jacob{i}]= pcg(H, f, tol, maxit, M);

    % Computing the error for both Jacobi and Ichol preconditioner(CG)
    for j=1:1:length(u_exact)

        epsilon_ichol(i)= epsilon_ichol(i) + (((uh_ichol(j)-u_exact(j))^2)*(s_node(j)));
        epsilon_jacob(i)= epsilon_jacob(i) + (((uh_jacob(j)-u_exact(j))^2)*(s_node(j)));
    end

    epsilon_ichol(i)= sqrt(epsilon_ichol(i));
    epsilon_jacob(i)= sqrt(epsilon_jacob(i));

    % Computing the error ratio for both Jacobi and Ichol preconditioner(CG)
    if i>1

        eps_ichol_ratio(i)= epsilon_ichol(i)/epsilon_ichol(i-1);
        eps_jacob_ratio(i)= epsilon_jacob(i)/epsilon_jacob(i-1);
    end


    % Semi-logarithmic convergence plots
    figure('Name',sprintf('mesh%i', n_file))
    semilogy(0:iter_ichol(i), resvec_ichol{i},'--o', 0:iter_jacob(i), resvec_jacob{i},'--*');
    legend(sprintf('Ichol (mesh%i)', n_file), sprintf('Jacobi (mesh%i)', n_file));
    xlabel('Iteration');
    ylabel('Residual norm');
end



%% ____________________________________________________________________________________

%% Plot of the convergence history for all the refinement levels
%% using Incomplete Cholesky factorization

figure('Name', 'Ichol')

axes('YScale', 'log')
xlabel('Iteration');
ylabel('Residual norm');

hold on

for i=1:1:length(n_mesh)

    % Plotting the convergence profile for each run on the same figure
    plot(0:iter_ichol(i), resvec_ichol{i},'--*');
end

legend('mesh0 Ichol','mesh1 Ichol', 'mesh2 Ichol', 'mesh3 Ichol', 'mesh4 Ichol');

hold off


%% ____________________________________________________________________________________

%% Plot of the convergence history for all the refinement levels
%% using Jacobi preconditioner

figure('Name', 'Jacobi')

axes('YScale', 'log')
xlabel('Iteration');
ylabel('Residual norm');

hold on

for i=1:1:length(n_mesh)

    % Plotting the convergence profile for each run on the same figure
    plot(0:iter_jacob(i), resvec_jacob{i}, '--.');
end

legend('mesh0 Jacobi','mesh1 Jacobi', 'mesh2 Jacobi', 'mesh3 Jacobi', 'mesh4 Jacobi');

hold off

%% ____________________________________________________________________________________


%% Table showing the values of epsilon, epsilon's ratio and iter
 
figure('Name', 'Table values')

hold on

mesh= {'mesh0'; 'mesh1'; 'mesh2'; 'mesh3';'mesh4'};
epsilon_ichol= [epsilon_ichol(1); epsilon_ichol(2); epsilon_ichol(3); epsilon_ichol(4); epsilon_ichol(5)];
eps_ichol_ratio= [ eps_ichol_ratio(1); eps_ichol_ratio(2); eps_ichol_ratio(3); eps_ichol_ratio(4); eps_ichol_ratio(5) ];
iter_ichol= [ iter_ichol(1); iter_ichol(2); iter_ichol(3); iter_ichol(4); iter_ichol(5)];
epsilon_jacob= [epsilon_jacob(1); epsilon_jacob(2); epsilon_jacob(3); epsilon_jacob(4); epsilon_jacob(5)];
eps_jacob_ratio= [eps_jacob_ratio(1); eps_jacob_ratio(2); eps_jacob_ratio(3); eps_jacob_ratio(4); eps_jacob_ratio(5) ];
iter_jacob= [iter_jacob(1); iter_jacob(2); iter_jacob(3); iter_jacob(4); iter_jacob(5) ];

T= table(epsilon_ichol, eps_ichol_ratio, iter_ichol,  epsilon_jacob, eps_jacob_ratio, iter_jacob, 'RowNames', mesh);

uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

hold off



%% ____________________________________________________________________________________

% %% EXTRA 
% 
% last_resvec_ichol= zeros(length(n_mesh),1);
% last_resvec_jacob= zeros(length(n_mesh),1);
% 
% for i=1:1:length(n_mesh)
% 
%     temporary= resvec_ichol{i};
%     last_resvec_ichol(i)= temporary(end);
%     clear temporary;
%     temporary= resvec_jacob{i};
%     last_resvec_jacob(i)= temporary(end);
%     clear temporary;
% end
% 
% disp('Ichol')
% disp('flag_ichol')
% disp(flag_ichol.')
% disp('relres_ichol')
% disp(relres_ichol.')
% disp('last_resvec_ichol')
% disp(last_resvec_ichol)
% 
% disp('Jacobi')
% disp('flag_jacob')
% disp(flag_jacob.')
% disp('relres_jacob')
% disp(relres_jacob.')
% disp('last_resvec_jacob')
% disp(last_resvec_jacob)
