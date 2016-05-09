function [cell, root] = regenerateCell( cell, X , evalFun, root)
    if(cell{1} ~= 1)
        if(length(cell{2}) == 0) %If has not childrens check value
            if(cell{1} == 0)
               cell{1} = 2;
               cell = expandCell(cell, evalFun, 1);
            end
        else
            for i=1:4
               if(isInBounds(X, cell{2}{i}{3}, cell{2}{i}{4}))
                    [cellRef{2}{i}, root] = regenerateCell(cell{2}{i}, X, evalFun, root);
               end
            end
        end
    end
end

