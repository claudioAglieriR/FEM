function H= stiffness_mat(topol, coord)

%% 3 - Computation of local stiffness matrices and their assembly

    % Creation of the matrix with the determined pattern, thanks to the
    % function "pattern_creation"
    H= pattern_creation(topol);

    % Computing the surface measure of each element of the mesh
    [s_triangles, ~]= surfaces(topol, coord);

    % x and y are the coordinates of each node
    x= coord(:,1);
    y= coord(:,2);

    for k=1:1:length(topol(:,1))

        ind= [1:3, 1:3];

        % Computation of coefficients necessary to compute the the local stiffness matrices
        b= zeros(3,1);
        c= zeros(3,1);

        % b(1,1) is "b_i", c(1,1) is "c_i"
        % b(2,1) is "b_j", c(2,1) is "c_j"
        % b(3,1) is "b_m", c(3,1) is "c_m"

        for i=1:1:3
            
            b(i,1)= ((y(topol(k,ind(1+i))))- y(topol(k,ind(2+i))));
            c(i,1)= ((x(topol(k,ind(2+i)))- x(topol(k,ind(1+i)))));
        end

        % ____________________________________________________________________________
        % a(1,1) would be "a_i", a(2,1) would be "a_j", a(3,1) would be "a_m"
        %
        %   a= zeros(3,1);
        %   a(1,1)= ((x(topol(k,2))*y(topol(k,3)))- (x(topol(k,3))*y(topol(k,2))));
        %   a(2,1)= ((x(topol(k,3))*y(topol(k,1))- (x(topol(k,1))*y(topol(k,3)))));
        %   a(3,1)= (((x(topol(k,1))*y(topol(k,2)))- (x(topol(k,2))*y(topol(k,1)))));
        % ____________________________________________________________________________


        % Computation of Hloc for each element, i.e. the local stiffness matrix
        bmat= zeros(3,3);
        cmat= zeros(3,3);

        for i=1:1:3
            for j=1:1:3

                bmat(i,j)= b(i)*b(j);
                cmat(i,j)= c(i)*c(j);
            end
        end

        Hloc= (1/(4*s_triangles(k)))*(bmat+cmat);

        % Assembly of the global stiffness matrix H

        for i=1:1:3
            
            row= topol(k,i);

            for j=1:1:3

                col= topol(k,j);
                H(row,col)= H(row,col)+ Hloc(i,j);
            end
        end
    end
end