function H_structure= pattern_creation(topol)

%% 2 - Pattern creation for the stiffness matrix

    [m,n]= size(topol);
    row= zeros(m*n, 1);
    col= zeros(m*n, 1);

    j= 0;
    k= 1;

    while j<m

        % Assigning values to "row" vector. "row" has this form:
        % [1,1,1,2,2,2,3,3,3, ...]

        row(k)= k-(2*j);
        row(k+1)= k-(2*j);
        row(k+2)= k-(2*j);

        j= j+1;

        % Assigning values to "col" vector.
        % "col" contains all the values of topol matrix

        col(k)= topol(j,1);
        col(k+1)= topol(j,2);
        col(k+2)= topol(j,3);

        k= k+3;
    end

    row= row.';
    col= col.';

    % Creation of the data structure necessary to store H
    A= sparse(row,col,1);
    H_structure= A.'*(A);
    H_structure= H_structure*0;
end
