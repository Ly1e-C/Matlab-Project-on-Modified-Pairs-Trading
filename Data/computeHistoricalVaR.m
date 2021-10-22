function VaR = computeHistoricalVaR(returns,confidence_level,plot_flag)
% Inputs:
% returns               Vector of returns
% confidence_level      Confidence level (default 0.95)
% plot_flag             if true, visualize result (default is true)
%
% Outputs:
% VaR                   Value at Risk

% handle inputs
if nargin < 3
    plot_flag = true;
end
if nargin < 2
    confidence_level = 0.95;
end

% Sort returns from smallest to largest
sorted_returns = sort(returns);

% Store the number of returns
num_returns = numel(returns);

% Calculate the index of the sorted return that will be VaR
VaR_index = ceil((1-confidence_level)*num_returns);

% Use the index to extract VaR from sorted returns
VaR = sorted_returns(VaR_index);

% Plot results if requested
if plot_flag
    % Histogram data
    [count,bins] = hist(returns,30);
    % Create 2nd data set that is zero above Var point
    count_cutoff = count.*(bins < VaR);
    % Scale bins
    scale = (bins(2)-bins(1))*num_returns;
    % Plot full data set
    bar(bins,count/scale,'b');
    hold on;
    % Plot cutoff data set
    bar(bins,count_cutoff/scale,'r');
    grid on;
    hold off;
    title(['Histogram of Returns. Red Indicates Returns Below VaR: ',num2str(VaR)],'FontWeight','bold');
end