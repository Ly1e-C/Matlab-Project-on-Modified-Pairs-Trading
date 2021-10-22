function sh = pairsFun2(x,data,selectedPairs,spreadH,spreadL,scaling,cost)
% define pairs to accept vectorized inputs and return only sharpe ratio
%%
% Copyright 2010, The MathWorks, Inc.
% All rights reserved.

[row,col] = size(x);
sh  = zeros(row,1);
x = round(x);

if ~exist('scaling','var')
    scaling = sqrt(252);
end

if ~exist('cost','var')
    cost = 0.01;
end

if ~exist('spreadH', 'var')
    spreadH = 1.3;
end

if ~exist('spreadL', 'var')
    spreadL = 0.8;
end

% run simulation
parfor i = 1:row
    if col == 2
        [~,~,sh(i)] = pairs5(data, selectedPairs, x(i,1), x(i,2), spreadH,spreadL,scaling,cost);
    else
        [~,~,sh(i)] = pairs5(data, selectedPairs, x(i,1), x(i,2), x(i,3), scaling, cost);
    end
end