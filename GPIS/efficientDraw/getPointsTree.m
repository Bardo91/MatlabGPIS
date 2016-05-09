function [ points, col ] = getPointsTree( points, col, cell )
    if(length(cell{2}) == 0)
          points = [points; cell{5}];
          if(cell{6}(1) > 0)
             col = [col; +1]; 
          else
             col = [col;-1];
          end
    else
        for i = 1:4
           [points, col] = getPointsTree(points, col, cell{2}{i});
        end

    end

end

