%% Function 2: Helper function to plot a 3D grid
% This function visualizes a 3D mesh grid by drawing its lines.
function plotGrid(X, Y, Z)
    % Plot the mesh grid lines for each slice along the third dimension.
    for k = 1:size(Z, 3)
        mesh(X(:,:,k), Y(:,:,k), Z(:,:,k), 'EdgeColor', [0.7 0.7 0.7], 'FaceAlpha', 0);
    end

    % Plot the grid lines running parallel to the third dimension.
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            plot3(squeeze(X(i,j,:)), squeeze(Y(i,j,:)), squeeze(Z(i,j,:)), 'Color', [0.7 0.7 0.7]);
        end
    end

    % Plot the grid lines running parallel to the first dimension.
    for i = 1:size(X, 3)
        plot3(squeeze(X(:,:,i)), squeeze(Y(:,:,i)), squeeze(Z(:,:,i)), 'Color', [0.7 0.7 0.7]);
    end
end