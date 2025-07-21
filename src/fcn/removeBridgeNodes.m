function G_minimal = removeBridgeNodes(G_minimal)

allDeg = degree(G_minimal);

for i=1:numel(allDeg)
    % Check if node is a bridge, i.e., degree = 2.
    try % TODO replace try-catch statement
        if degree(G_minimal, string(i)) == 2
            % Identify the node's two neighbors
            nghbrs = neighbors(G_minimal, string(i));
            % Extract the weights of the bridge node's two edges.
            dists = distances(G_minimal, string(i), nghbrs);
            % Add a new edge between the two neighbors, the weight is the sum
            % of the two previous.
            G_minimal = addedge(G_minimal, nghbrs(1), nghbrs(2), (dists(1)+dists(2)));
            % Remove the bridge node from G_new
            G_minimal = rmnode(G_minimal, string(i));
        end
    catch
    end

end



end