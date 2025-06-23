%单条射线穿过网格的可视化
function visualize_results(P1, P2, vertices, P_enter, P_exit)
    figure;
    hold on;
    grid on;
    axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('射线穿过网格可视化');
    
    % 绘制网格
    faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 
             2 4 8 6; 4 3 7 8; 3 1 5 7];
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', 'cyan', 'FaceAlpha', 0.2, 'EdgeColor', 'b');
    
    % 绘制射线
    plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'r-', 'LineWidth', 1.5);
    
    % 标记关键点
    plot3(P_enter(1), P_enter(2), P_enter(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot3(P_exit(1), P_exit(2), P_exit(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % 添加标签
    text(P_enter(1), P_enter(2), P_enter(3), ...
        sprintf(' 穿入点\n (%.1f,%.1f,%.1f)', P_enter), 'FontSize', 8);
    text(P_exit(1), P_exit(2), P_exit(3), ...
        sprintf(' 穿出点\n (%.1f,%.1f,%.1f)', P_exit), 'FontSize', 8);
    
    legend('网格', '射线', '穿入点', '穿出点');
    view(3);
    rotate3d on;
end