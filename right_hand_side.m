function f= right_hand_side(coord, s_node)

%% 4 - Computation of the right-hand side

    % x and y are the coordinates of each node
    x= coord(:,1);
    y= coord(:,2);

    for i=1:1:(length(coord(:,1)))
        
        f(i)= (-4 + 2*x(i)^2 + 2*y(i)^2)*s_node(i);
    end

    f= f.';
end