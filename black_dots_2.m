% Pete Rigas, Craft Lab, Summer 2019

function [Blob_Lengths, Total_Blobs, zero_positions, New_zero_positions, New_zero_positions_half, Centers_of_Dots_vector] = black_dots_2(Image, Important_value)

% Change output parameters from the function if needed
% [Image_short, C, Blob_Lengths, Blob_Counts, Total_Blobs, New_zero_positions, New_zero_positions_half, Centers_of_Dots_vector] = black_dots_2(Image)

% black_dots = [];
count_vec = zeros(size(Image,1),size(Image,2));
count_vec_2 = zeros(size(Image,1),size(Image,2));

Image = double(Image);
Image_short = Image./(max(max(Image)));

% Creat count_vec for the for loop below

                % COMMENTED THIS PART OUT, NO LONGER NECESSARY, butit does
                % give a way to look at pixels which could be useful for
                % random colorings of Voronoi domains

% for i = 1 : size(Image,2)
%     for j = 1: size(Image,1)
%         if Image_short(j,i) > 0 && Image_short(j,i) < 1
%             if (Image_short(j,i) >= 0.5 && Image_short(j,i) < 1)
%                 % Transitioning from black pixels
%                 count_vec(j,i) = 0;
%             elseif (Image_short(j,i) > 0 && Image_short(i,j) < 0.5 || isequal(Image_short(j,i),0))
%                 % Transitioning from white pixels
%                 count_vec(j,i) = 1;
%             
%             end
%         end
%     end
% end 

                % COMMENTED THIS PART OUT, NO LONGER NECESSARY, butit does
                % give a way to look at pixels which could be useful for
                % random colorings of Voronoi domains


% [row, col] = find(count_vec == 0);
% quantity = [row col];
% [row col] = find(Image_short==0);
% A = [row col];
% maximum_number_of_black_dots = length(A(:,1));

                %  THIS PART IS NO LONGER NECESSARY AS WELL. This part was
                % used to compute the total number of black and white
                % pixels in the image, which it can still accomplish.
            

% total_pixels = size(Image_short,1) * size(Image_short,2);
% sum_vec_1 = zeros(1,size(Image,2));
% sum_vec_2 = zeros(1,size(Image,2));
% 
% for i = 1 : size(Image_short,2)
%     x_white = sum(Image_short(:,i) ==0);
%     x_black = sum(Image_short(:,i) == 1);
%     sum_vec_1(i) = x_white;
%     sum_vec_2(i) = x_black;
% end 
% 
% total_white = sum(sum_vec_1);
% total_black = sum(sum_vec_2);
% 
% sum_vec_3 = size(Image_short,1)-sum_vec_2;
% sum_vec_transition = sum_vec_3 - sum_vec_1;
% 
% if ~(isequal(sum(sum_vec_1)+sum(sum_vec_2)+sum(sum_vec_transition), total_pixels))
%     disp('Error: total number of pixels in image is not equal to the total number of black, white and grey pixels')
% end 


C = { };

for z = 1: size(Image_short,2)
    [row col] = find(Image_short(:,z)==Important_value);
    vector_zeros = [row col];
    C{z} = vector_zeros;
end 

% blob_size = zeros(1, length(vector_zeros));

C_1=cellfun('length',C);

% C_1(C_1 == 0) = [ ];

Blob_Lengths = { };
Blob_Counts = zeros(1,size(Image,2));
% Initially had circle_count and circle_length variables here

for z = 2 : size(Image_short,2)
    % Declare a vector for the number of 0's in a fixed column of Image_short

        % Instantiate the vector blob_sizes, from the number of 0's that we count in a
        % fixed column of Image_short
        
        if C_1(z) == 0
            
            disp('No zeros found') 
            
        else
        
        X = C{z}(1: C_1(z)-1 , 1);
            
            
            for J = 2 : length(X)
                
                circle_count = 1;
                circle_length = 0;
                
                if isequal(abs(C{z}(J-1,1)- C{z}(J,1)),1)
                     Blob_Counts(z) = circle_count;
                     % epsilon = abs(C{z}(J-1,1) - C{z}(J,1));
                     circle_length = circle_length + 1;
                     Blob_Lengths{z}(J,1) = circle_length;
                     
                elseif ~(isequal(C{z}(J-1,1)-C{z}(J,1),1))
                    % circle_center = abs(C{z}(J-1,1)- C{z}(J,1));
                    % Blob_Lengths{z}(J,1) = circle_length;
                    circle_count = circle_count+1;
                    Blob_Counts(z) = circle_count;
                end

            end 
            
        end 
        
end
        

Total_Blobs = sum(Blob_Counts)/length(Blob_Counts);
Total_Blobs = Total_Blobs/2;

% Round Total_Blobs to an integer number
Total_Blobs = round(Total_Blobs);

% Calculate the length of each blob, as well as the total number of blobs

% Blob_Final_Lengths = { };
Array_2 = { };

zero_positions = { };

for P = 1: length(Blob_Lengths)
        x = find(Blob_Lengths{P}==Important_value);
        zero_positions{P} = x;
end

for AP = 1 : length(zero_positions)
    for q = 1 : length(zero_positions{AP})
        if zero_positions{AP}(end) < length(Blob_Lengths{AP})
            zero_positions{P}(end+1) = length(Blob_Lengths{AP});
        end
    end 
end

New_zero_positions = { };
New_zero_positions_half = { };

% Instantiate the New_zero_positions and New_zero_positions_half array.
% With the half array, we will be able to roughly determine where the 

for P = 1: length(Blob_Lengths)
    W = zero_positions{P};
    if (isequal(length(W),1))
            
            disp('Column of matrix only has ONE black pixel')
            New_zero_positions{P}(1) = zero_positions{P}(1)-1;
            
            New_zero_positions_half{P}(1) = New_zero_positions{P}(1)/2;
            
            New_zero_positions{P}(2) = length(Blob_Lengths{P}) - zero_positions{P}(1)-1;
            
            New_zero_positions_half{P}(2) = New_zero_positions{P}(2)/2;
            
    % Check if there are 2 0's next to each other in the Blob_Lengths Cell  
            
%     elseif (length(W)> 2)
%         for U = 2 : length(W)
%             if isequal(abs(W(U-1)-W(U)),1)
%                 New_zero_positions{P}(1) = 2;
%             end 
%         end 
            
    elseif (isequal(length(W),2))
        
       for Y = 2 : length(W)    
    
        previous_zero = zero_positions{P}(Y);

        New_zero_positions{P}(Y-1) = previous_zero - zero_positions{P}(Y-1);
        
        New_zero_positions_half{P}(Y-1) = New_zero_positions{P}(Y-1)/2;
        
        % New_zero_positions{P}(Y) = zero_positions{P}(end) - zero_positions{P}(Y-1);
    
       % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;
       
       end
       
       New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - New_zero_positions{P}(Y-1)-1;
       
       New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;

    elseif (isequal(length(W),3))
         for Y = 2 : length(W)
             
           previous_zero = zero_positions{P}(Y);
             
           New_zero_positions{P}(Y-1) = previous_zero - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y-1) = New_zero_positions{P}(Y-1)/2;
        
          %  New_zero_positions{P}(Y) = zero_positions{P}(end) - zero_positions{P}(Y-1);
    
          %  New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;  
           
         end
         
         New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
         
         New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),4))
         for Y = 2 : length(W)
             
           previous_zero = zero_positions{P}(Y);
           
           New_zero_positions{P}(Y-1) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y-1) = New_zero_positions{P}(Y-1)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end
         
            New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
            
            New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
          elseif (isequal(length(W),5))
         for Y = 2 : length(W)
           
           previous_zero = zero_positions{P}(Y);
             
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = length(Blob_Lengths{P}) - New_zero_positions{P}(Y+1);
           
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),6))
         for Y = 2 : length(W)
             
             previous_zero = zero_positions{P}(Y);
             
             
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
          New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),7))
         for Y = 2 : length(W)
             
           previous_zero = zero_positions{P}(Y);
             
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
          New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),8))
         for Y = 2 : length(W)
             
             previous_zero = zero_positions{P}(Y);
             
             
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
          % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
          New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;

         elseif (isequal(length(W),9))
         for Y = 2 : length(W)
             
             previous_zero = zero_positions{P}(Y);
 
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
          New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),10))
         for Y = 2 : length(W)
             
             previous_zero = zero_positions{P}(Y);
 
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
         
         elseif (isequal(length(W),11))
         for Y = 2 : length(W)
             
             previous_zero = zero_positions{P}(Y);
             
           New_zero_positions{P}(Y) = zero_positions{P}(Y) - zero_positions{P}(Y-1);
           
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;
        
           % New_zero_positions{P}(Y+1) = zero_positions{P}(end) - zero_positions{P}(Y-1) - New_zero_positions{P}(Y) + 1;
    
           % New_zero_positions{P}(Y+2) = abs(length(Blob_Lengths{P}(end)- (New_zero_positions{P}(Y) + New_zero_positions{P}(Y+1))))+3;     
         end 
         
          New_zero_positions{P}(Y) = length(Blob_Lengths{P}) - sum(New_zero_positions{P})-1;
          
           New_zero_positions_half{P}(Y) = New_zero_positions{P}(Y)/2;

    end 
  
end

% Separate last chunk of code to troubleshoot the iteration ending with an
% off by one error.

Centers_of_Dots_vector = [ ];
 
% Not completely sure why we have to subtract 3. This is the only logical
% possibility that prevents the code from crashing, from the dimensions of
% the other objects that I have created.

% From cases that did not work, subtracting by the difference in parenthese
% below prevents the error from occuring. Admittedly, I do not completely
% understand why it works.

for P = 1: size(Image,2) - (size(Image,2) - size(New_zero_positions,2))
    
    % x = find(Blob_Lengths{P}==0);
    
    % 
    
    x = zero_positions{P};
    
    % y = New_zero_positions_half{P};
    
    % It is important to keep in these alternative branches of the if
    % statement.
    
    % if isempty(New_zero_positions_half{P-3})
       % disp('There are no lowest pixel values in this column')
  
    % elseif ~(isequal(length(New_zero_positions_half{P}),0))
    if ~(isequal(length(New_zero_positions_half{P}),0))
    for A = 1 : length(Blob_Lengths{P})
        for AB = 2 : length(New_zero_positions_half{P})
        for B = 1 : length(x)
                if isequal(A, x(B))
                Y_1 = x(B);
                Y_2 = New_zero_positions_half{P}(AB-1);
                Y_3 = round(Y_1+Y_2);
               %  Centers_of_Dots_Cell{P}(A) = round(x(B) + New_zero_positions_half{P}(B));
                disp('_____________NEW DOT_____________')
                disp('Position of zero:')
                disp(Y_1)
                disp('Position of center of black/white dot, from zero position above:')
                disp(Y_2)
                disp('x coordinate center of dot:')
                disp(Y_3)
                disp('y coordinate center of dot:')
                disp(P)
                % Store variables for x coordinates and y coordinates in
                % Centers_of_Dots_vector
                Centers_of_Dots_vector(1,P) = Y_3;
                Centers_of_Dots_vector(2,P) = P;
                
                end 
        end 
        end
    end
    
    elseif isequal(length(New_zero_positions_half{P}),1)
        disp('Error')
    end  
    
    
end  

end 
                % if isequal(mod(Blob_Lengths{P}(1:length(Blob_Lengths{P}),1),2),1)
                
                    

%                         vector_A1 = Blob_Lengths{P}(x(A1-1)+1:x(A1)-1,1);
%                         Array_2{P}(:,J-1) = vector_A1;
                        
                        
                       %  Array_2{P}(:,J+1) = zeros(abs(x(A1) - x(A1+1))-1 ,1);
                    
                    
                 
                    
%                         vector_A2 = Blob_Lengths{P}(x(A2-1) + 1 : x(A2) -1 , 1);
%                         Array_2{P}(:,A2) = vector_A2;
                        
                        % Array_2{P}(:,A2) = zeros(abs(x(A2)-x(A2+1))-1,1);
            
            % elseif

            % Blob_Lengths_Numbers{P}(x(F),1) = abs(Blob_Lengths{P}(x(F-1),1)-Blob_Lengths{P}(X(F),1)); 
 

    
% Q = find(Blob_counts == 0);
% final_count = zeros(1, length(Blob_Counts)-length(Q));
% for I = 2 : length(Blob_counts)
%         for O = 1 : length(Q)
%            if isequal(I , Q(O))
%               disp('Entry is equal to zero')
%            elseif Blob_Counts(I) > 0
%                if isequal(abs(Blob_Counts(I-1), Blob_Counts(I)),1)
%                    % final_count
%                end 
%            end 
%         end 
% end 

   
                 
%                   ~(isequal(mod(length(vector_zeros),2),0))
%                      % Make vector_zeros have even length
%                        vector_zeros(end) = [];
%                        blob_size(y) = abs(Image_short(C{z}(1: length(C{z}),1),z)-Image_short(C{z}(1:length(C{z}),1),z));
%                  end 
%                 
%                   
% end 
         
 
 
 






%  for u = 2: 399
%      for t = 2: 393
%         %for b = 1 : size(quantity,2)
%             %for a = 1: size(quantity,1)
%                     if ~(isequal(Image_short(t,u),1))
%                         if Image_short(t,u) < 1
%                         if isequal(Image_short(t,u),0)
%                             if isequal(Image_short(t,u-1), Image_short(t,u))
%                                 % Check whether entries surrounding entry u,t of
%                                 % count_vec are also 1
%                     
%                             % update maximum number of black dots, in fact, decrease the
%                             % number of dots to reflect the fact that a
%                             % 'cluster' is surrounding the 1 value that we
%                             % observed
%                                 Image_short(t,u-1)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                         
%                             elseif isequal(Image_short(t,u+1), Image_short(t,u))
%                                 
%                                 Image_short(t,u+1)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                                 
%                             elseif isequal(Image_short(t-1,u), Image_short(t,u)) 
%                                 
%                                 Image_short(t-1,u)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                                 
%                             elseif isequal(Image_short(t+1,u) , Image_short(t,u))
%                                 
%                                 Image_short(t+1,u)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%    
%                             end
% 
%                             elseif (Image_short(t,u) > 0 && Image_short(t,u) < 1)
%  
%                                if isequal(Image_short(t,u-1), Image_short(t,u))
%                                 % Check whether entries surrounding entry u,t of
%                                 % count_vec are also 1
%                     
%                             % update maximum number of black dots, in fact, decrease the
%                             % number of dots to reflect the fact that a
%                             % 'cluster' is surrounding the 1 value that we
%                             % observed
%                             
%                                 Image_short(t,u-1) = 7;
%                         
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                         
%                             elseif isequal(Image_short(t,u+1), Image_short(t,u))
%                                 
%                                 Image_short(t,u+1) = 7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                                 
%                             elseif isequal(Image_short(t-1,u), Image_short(t,u))
%                                 
%                                 Image_short(t-1,u)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                                 
%                             elseif isequal(Image_short(t+1,u) , Image_short(t,u))
%                                 
%                                 Image_short(t+1,u)=7;
%                                 maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                                 
%                             elseif Image_short(t,u-1) , Image_short(t,u) > 0
%                                 
%                                 Image_short(t,u-1)=7;
%                                 
%                                 
%                             elseif Image_short(t,u-1) , Image_short(t,u) < 0
%                                 
%                                 Image_short(t,u-1)=7;
%                                 
                                
%                               end 
%                        
%                         end 
%                         
%                         end 
%                     end 
%                         
          % Check whether entries surrounding entry u,t of
                    % count_vec are also 1
                    
                        % update maximum number of black dots, in fact, decrease the
                        % number of dots to reflect the fact that a
                        % 'cluster' is surrounding the 1 value that we
                        % observed
                        
%      end
%  end 
                        
%                         if (isequal(count_vec(t-1,u) , count_vec(t,u)) && isequal(count_vec(t+1,u) , count_vec(t,u)))
%                        
%                             maximum_number_of_black_dots = maximum_number_of_black_dots-1;
%                         
%                         
% %                 elseif (isequal(count_vec(t+1,u+1) , count_vec(t,u)) || isequal(count_vec(t-1,u-1), count_vec(t,u)))
% %                         
% %                         maximum_number_of_black_dots = maximum_number_of_black_dots- 0.5;
%                         end 
%                     end
%      end
%  end 


%     [row col] = find(count_vec == 1)
%     x = [row col]
%     
%     max_number_dots = [ ];
%     
%     for a = 1: size(x,2)
%         for b = 1:size(x,1)
%             % Calculate maximum number of black dots, per column
%             
%         end 
%     end
            
            
            
%             if abs(Image_short(j+1,i) - Image_short(j,i)) > 0 && abs(Image_short(j+1,i) - Image_short(j,i)) <= 0.5
%                 
%                 
%                 count_vec(j,i) = 0;
%                     for j_0 = j : 393
%                              if abs(Image_short(j_0,i) - Image_short(j_0+1,i)) < 0.5
%                                  % Add NO more counts to the count_vector,
%                                  % because the change in pixel value is not
%                                  % big enough, and we are still in the same
%                                  % black dot!
%                                  
%                                  count_vec(j_0,i) = 0;
%                                  
%                              elseif abs(Image_short(j_0,i) - Image_short(j_0+1,i)) > 0.5
%                                  for a = j_0 : size(Image_short,1)-1
%                                      count_vec(a,i) = count_vec(a,i) + 1;
%                                  end 
%                              end 
% 
%                     end 
%             
%             end
%         end 
%                      
%             if abs(Image_short(j+1,i) - Image_short(j,i)) == 0
%                 black_dots(j,i) = 0;
%                 
%             elseif abs(Image_short(j+1,i) - Image_short(j,i)) > 0.5 && abs(Image_short(j+1,i) - Image_short(j,i)) < 1
%                 black_dots(j,i) = 0;
%             end 
%     end  
% end  
 

