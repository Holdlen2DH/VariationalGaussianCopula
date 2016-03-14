% Calculate the modified Z-score of a query point from a window set 
% Shaobo Han
% 08/13/2015
function score = mzscore_test(x, windowset)
    tmpX = [x, windowset]; 
    tmpm = median(tmpX); 
    mad = median(abs(tmpX-tmpm)); 
    score = 0.6745*(x-tmpm)/mad; 
end