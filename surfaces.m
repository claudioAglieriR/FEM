function [s_triangles, s_node]= surfaces(topol, coord)

%% Computation of the surface measures associated with topol and coord files

    % s_triangles is the vector containing the elements (triangles) surface measure
    s_triangles= zeros(length(topol(:,1)), 1);

    % s_node is the vector containing the surface measure relative to each node
    s_node= zeros((length(coord(:,1))), 1);


    % Computation of each element surface measure
    for i=1:1:length(s_triangles)

        s_triangles(i)= 0.5* det( [1, coord(topol(i,1),:); 1, coord(topol(i,2),:); 1, coord(topol(i,3),:)] );
    end

    % Computation of the surface measure relative to each node
    for i=1:1:length(s_node)

        for j=1:length(s_triangles)

            if (ismember(i,topol(j,:)))

                s_node(i)= s_node(i)+s_triangles(j);
            end
        end
        
        % One third of the area of the triangle patch surrounding node i
        s_node(i)= s_node(i)/3;
    end
end
