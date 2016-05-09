function [ cellRef ] = expandCell( cellRef, evalFun)
    if(length(cellRef{2}) == 0)
        xLim = cellRef{3};
        yLim = cellRef{4};
        centroid = cellRef{5};
        cellRef{2} = {  {1,{},[xLim(1), centroid(1)],[yLim(1), centroid(2)],[sum([xLim(1), centroid(1)])/2, sum([yLim(1), centroid(2)])/2], evalFun([sum([xLim(1), centroid(1)])/2, sum([yLim(1), centroid(2)])/2]')}, 
                        {1,{},[xLim(1), centroid(1)],[centroid(2), yLim(1)],[sum([xLim(1), centroid(1)])/2, sum([centroid(2), yLim(1)])/2], evalFun([sum([xLim(1), centroid(1)])/2, sum([centroid(2), yLim(1)])/2]')}, 
                        {1,{},[centroid(1), xLim(2)],[yLim(1), centroid(2)],[sum([centroid(1), xLim(2)])/2, sum([yLim(1), centroid(2)])/2], evalFun([sum([centroid(1), xLim(2)])/2, sum([yLim(1), centroid(2)])/2]')}, 
                        {1,{},[centroid(1), xLim(2)],[centroid(2), yLim(1)],[sum([centroid(1), xLim(2)])/2, sum([centroid(2), yLim(1)])/2], evalFun([sum([centroid(1), xLim(2)])/2, sum([centroid(2), yLim(1)])/2]')}};
    else
        for i = 1:4
           if(cellRef{2}{i}{1} == 2)
               cellRef{2}{i} = expandCell( cellRef{2}{i}, evalFun);
           end
        end
    end
end

