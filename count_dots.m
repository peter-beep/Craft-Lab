function count_dots(Image)
    
    % Instantiate the threshold vector, with all possible threshold values,
    % to determine which one is the most appropriate.
    threshold_vector = [];
    x = unique(Image);    
    
    for j = 1 : length(x)
        threshold_vector(j) = x(j);
    end

    new_vector = [];

    for k = 1 : length(threshold_vector)
        for l = 1: length(Image)
            [labeledImage, Black_Dots] = bwlabel(Image(:,l) > threshold_vector(k));
            % Make sure that labeled_image_i is a summation of the number of
            % black dots, taken over all entries of the labeled_image_i vector
            labeledImage = sum(labeledImage);
            new_vector(l) = labeledImage;
        end
    end
    
    disp(new_vector)
    disp(max(new_vector))
    disp(min(new_vector))

end 