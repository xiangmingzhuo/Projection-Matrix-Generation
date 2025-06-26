%% 辅助函数：绘制三维网格
function plotGrid(X, Y, Z)
    % 绘制网格线
    for k = 1:size(Z,3)
        mesh(X(:,:,k), Y(:,:,k), Z(:,:,k), 'EdgeColor', [0.7 0.7 0.7], 'FaceAlpha', 0);
    end
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            plot3(squeeze(X(i,j,:)), squeeze(Y(i,j,:)), squeeze(Z(i,j,:)), 'Color', [0.7 0.7 0.7]);
        end
    end
    for i = 1:size(X,3)
        plot3(squeeze(X(:,:,i)), squeeze(Y(:,:,i)), squeeze(Z(:,:,i)), 'Color', [0.7 0.7 0.7]);
    end
end