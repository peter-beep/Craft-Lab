% Pete Rigas, Craft Lab Summer 2019

function Voronoi_(vector)
    % Obtain the x and y coordinates for the center of each dot that we are
    % counting
    
    vector_1 = vector(1,:);
    vector_2 = vector(2,:);
    
    [vx , vy] = voronoi(vector_2 , vector_1);
    [vx_1 , vy_1] = voronoi(vector_1 , vector_2);
    
    % First Figure (THIS FIGURE HAS ALL COMBINATIONS OF VECTORS WITH VX AND
    % VY).  
    
    figure(1)
    subplot(2,1,1)
    plot(vector_2,vector_1, 'r+' , vx , vy , 'b-')
    xlim([0 200])
    ylim([0 200])
%     subplot(5,1,2)
%     plot(vector_1,vector_2, 'r+' , vx , vy , 'b-')
%     xlim([0 200])
%     ylim([0 200])
%     subplot(3,1,3)
%     plot(vector_2,vector_1, 'r+' , vx_1 , vy_1 , 'b-')
%     xlim([0 200])
%     ylim([0 200])
%     subplot(3,1,2)
%     plot(vector_1,vector_2, 'r+' , vx_1 , vy_1 , 'b-')
%     xlim([0 200])
%     ylim([0 200])
    subplot(2,1,2)
    imshow('voro-1.png')
    
    % Final First Figure Plot
    
%     figure(1)
%     subplot(2,1,1)
%     plot(vector_2,vector_1, 'r+' , vx , vy , 'b-')
%     xlim([0 200])
%     ylim([0 200])
%     subplot(2,1,1)
%     imshow('Image (99).png')
%     xlim([0 200])
%     ylim([0 200])
    
    % Second Figure
    
%     figure(2)
%     
%     % First Plot
%     subplot(5,1,1)
%     plot(vector_2,vector_1, 'r+' , vx , vy , 'b-')
%     xlim([min(vector_1)+1 max(vector_1)+1])
%     ylim([min(vector_2)+1 max(vector_2)+1])
%     title('First Voronoi Plot')
%     
%     % Second Plot
%     subplot(5,1,2)
%     plot(vector_2,vector_1, 'r+' , vx_1 , vy_1 , 'b-')
%     xlim([min(vector_2) max(vector_2)])
%     ylim([min(vector_1) max(vector_1)])
%     title('Second Voronoi Plot')
%     
%     % Original Image that we were trying to compute the dots of
%     subplot(5,1,3)
%     imshow('Image (93).png')
%     title('Original Image')
    
    % Second Figure
    
%     for j = 10:-1:2
%     figure(j)
%     
%     % Zoom in Plot of Voronoi Cells from the First Plot
%     subplot(5,1,1)
%     plot(vector_2,vector_1, 'r+' , vx , vy , 'b-')
%     xlim([min(vector_1)-75*j min(vector_1)+75 * j])
%     ylim([min(vector_2)-75*j min(vector_2)+75 * j]) 
%     
%     % Zoom in Plot of Voronoi Cells from the Second Plot
%     subplot(5,1,2)
%     plot(vector_2,vector_1, 'r+' , vx_1 , vy_1 , 'b-')
%     xlim([min(vector_2)-75*j min(vector_2)+75 * j])
%     ylim([min(vector_1)-75*j min(vector_1)+75 * j])
%     
%     % Zoom in of Black Dots from Original Image
%     subplot(5,1,3)
%     imshow('Image (93).png')
%     title('Zoom In of Original Image')
%     xlim([min(vector_1)-75*j max(vector_1)+75 * j])
%     ylim([min(vector_2)-75*j max(vector_2)+75 * j])
%     end 
%     
%     disp('Voronoi Area I')
%     XX_1 = polyarea(vx,vy)
%     disp(XX_1)
%     disp('Voronoi Area II')
%     XX_2 = polyarea(vx_1 , vy_1)
%     disp(XX_2)
    
end 