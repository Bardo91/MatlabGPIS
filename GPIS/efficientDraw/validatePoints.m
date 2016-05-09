function [ cell ] = validatePoints( cell, root )
    if(length(cell{2}) == 0) %If has not childrens validate
        incX = cell{3}(2) - cell{3}(1);
        incY = cell{4}(2) - cell{4}(1);
        xup = cell{5} - [0,incY];
        xbot = cell{5} + [0,incY];
        xleft = cell{5} - [incX,0];
        xright = cell{5} + [incX,0];
        
        if(isInBounds(xup, root{3}, root{4}))
            valup = sign(checkVal(root, xup));
        else
            valup = sign(cell{6}(1));
        end
        if(isInBounds(xbot, root{3}, root{4}))
            valbot = sign(checkVal(root, xbot));
        else
            valbot = sign(cell{6}(1));
        end
        if(isInBounds(xleft, root{3}, root{4}))
            valleft = sign(checkVal(root, xleft));
        else
            valleft = sign(cell{6}(1));
        end
        if(isInBounds(xright, root{3}, root{4}))
            valright = sign(checkVal(root, xright));
        else
            valright = sign(cell{6}(1));
        end
        
        if(abs(valup + valbot + valleft + valright + sign(cell{6}(1))) == 5)
           cell{1} = 0;
        else
           cell{1} = 2;
        end
        
    else %If has childrens go deeper
        for i = 1:4
           if(cell{2}{i}{1} ~= 0)
               cell{2}{i} = validatePoints( cell{2}{i}, root);
           end
        end
    end
end

