function [H]= bound_condition(H, bound)

%% 5 - Boundary conditions enforcement

    m= size(H, 1);
    Rmax= 1e15;

    for i=1:1:m

        if (ismember(i,bound))

            H(i,i)= Rmax;
        end
    end
end