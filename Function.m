% Separate last chunk of code to troubleshoot the iteration ending with an
% off by one error.

% Centers_of_Dots_vector = [ ];
 
% Not completely sure why we have to subtract 3. This is the only logical
% possibility that prevents the code from crashing, from the dimensions of
% the other objects that I have created.

% 

% for P = 1: size(Image,2) - 40
%     
%     % x = find(Blob_Lengths{P}==0);
%     
%     % 
%     
%     x = zero_positions{P};
%     
%     % y = New_zero_positions_half{P};
%     
%     % It is important to keep in these alternative branches of the if
%     % statement.
%     
%     % if isempty(New_zero_positions_half{P-3})
%        % disp('There are no lowest pixel values in this column')
%   
%     % elseif ~(isequal(length(New_zero_positions_half{P}),0))
%     if ~(isequal(length(New_zero_positions_half{P}),0))
%     for A = 1 : length(Blob_Lengths{P})
%         for AB = 2 : length(New_zero_positions_half{P})
%         for B = 1 : length(x)
%                 if isequal(A, x(B))
%                 Y_1 = x(B);
%                 Y_2 = New_zero_positions_half{P}(AB-1);
%                 Y_3 = round(Y_1+Y_2);
%                %  Centers_of_Dots_Cell{P}(A) = round(x(B) + New_zero_positions_half{P}(B));
%                 disp('_____________NEW DOT_____________')
%                 disp('Position of zero:')
%                 disp(Y_1)
%                 disp('Position of center of black/white dot, from zero position above:')
%                 disp(Y_2)
%                 disp('x coordinate center of dot:')
%                 disp(Y_3)
%                 disp('y coordinate center of dot:')
%                 disp(P)
%                 % Store variables for x coordinates and y coordinates in
%                 % Centers_of_Dots_vector
%                 Centers_of_Dots_vector(1,P) = Y_3;
%                 Centers_of_Dots_vector(2,P) = P;
%                 
%                 end 
%         end 
%         end
%     end
%     
%     elseif isequal(length(New_zero_positions_half{P}),1)
%         disp('Error')
%     end  
%     
%     
% end  