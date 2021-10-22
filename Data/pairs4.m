function varargout = pairs4(series2, M, N, spreadH, spreadL, scaling, cost, capital)
% PAIRS returns a trading signal for a simple pairs trading strategy

%%
% Copyright 2010, The MathWorks, Inc.
% All rights reserved.

%% Process input args
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

if ~exist('capital', 'var')
    capital = 1000;
end

if nargin == 1
    % default values
    M = 420; % lookback period
    N = 60; % default rebalance period
elseif nargin == 2
    error('PAIRS:NoRebalancePeriodDefined',...
        'When defining a lookback window, the rebalancing period must also be defined')
end

% Very often, the pairs will be convincingly cointegrated, or convincingly
% NOT cointegrated.  In these cases, a warning is issued not to read too
% much into the test statistic.  Since we don't use the test statistic, we
% can suppress these warnings.
warning('off', 'econ:egcitest:LeftTailStatTooSmall')
warning('off', 'econ:egcitest:LeftTailStatTooBig')

%% Sweep across the entire time series
% Every N periods, we use the previous M periods' worth of information to
% estimate the cointegrating relationship (if it exists).
%
% We then use this estimated relationship to identify trading opportunities
% until the next rebalancing date.

s = zeros(size(series2));
indicate = zeros(length(series2),1);
res = zeros(length(series2),1);
reg0 = zeros(length(series2),3);
rebalanceDayCounter = zeros(length(series2),1);
portValue = capital*ones(length(series2),1);
tradeState = 0;
lastDay = 0;
pxChange = zeros(length(series2),2);
sChange = zeros(length(series2),2);
r0 = zeros(length(series2),2);
r  = zeros(length(series2),1);
percentReturn  = zeros(length(series2),1);

for i = M+1:length(s)
    if lastDay == 1
        rebalanceDayCounter(i) = 1;
        lastDay = 0;
    else
        rebalanceDayCounter(i) = rebalanceDayCounter(i-1) + 1;
    end
    if rebalanceDayCounter(i) == N
        lastDay = 1;
    end
    if rebalanceDayCounter(i) == 1
        [h,~,~,~,reg1] = egcitest(series2(i-M:i-1, :));
        s(i, 2) = 0;
        s(i, 1) = 0;
        tradeState = 0;
        if i > length(s) - N + 1
            h = 0;
        end
    end
    if h ~= 0
        res(i) = series2(i-1,1) - (reg1.coeff(1) + reg1.coeff(2).*series2(i-1,2));
        reg0(i, 1) = reg1.coeff(1);
        reg0(i, 2) = reg1.coeff(2);
        reg0(i, 3) = reg1.RMSE;
        indicate(i) = res(i)./reg1.RMSE;
        if tradeState == 0
            if indicate(i) > spreadH
                s(i, 2) = (portValue(i-1)/2)/series2(i-1,2);
                s(i, 1) = -1/reg1.coeff(2) .* s(i, 2);
                tradeState = 1;
                rebalanceDayCounter(i) = 1;
            elseif indicate(i) < -spreadH
                s(i, 2) = (-portValue(i-1)/2)/series2(i-1,2);
                s(i, 1) = -1/reg1.coeff(2) .* s(i, 2);
                tradeState = 1;
                rebalanceDayCounter(i) = 1;
            else
                s(i, 2) = 0;
                s(i, 1) = 0;
                tradeState = 0;
            end
        else
            if abs(indicate(i)) < spreadL
                s(i, 2) = 0;
                s(i, 1) = 0;
                tradeState = 0;
                lastDay = 1;
            else
                s(i, 2) = s(i-1, 2);
                s(i, 1) = s(i-1, 1);
                tradeState = 1;
            end
        end
    end
    pxChange(i,:) = series2(i,:) - series2(i-1,:);
    sChange(i,:) = s(i,:) - s(i-1,:);
    r0(i,:) = [s(i,:).* pxChange(i,:) - abs(sChange(i,:))*cost/2];
    r(i) = sum(r0(i,:), 2);
    percentReturn(i) = r(i)/portValue(i-1);
    portValue(i) = portValue(i-1) + r(i);
end

%% Calculate performance statistics

sumr = cumsum(r);
sh = scaling*sharpe(percentReturn,0);
var95 = computeHistoricalVaR(percentReturn,0.95,false);
var99 = computeHistoricalVaR(percentReturn,0.99,false);
cVar95 = mean(percentReturn(percentReturn < var95));
cVar99 = mean(percentReturn(percentReturn < var99));
meanRet = mean(percentReturn);
medianRet = median(percentReturn);
volRet = std(percentReturn);
skewnessRet = skewness(percentReturn);
kurtosisRet = kurtosis(percentReturn);
%maxDD = maxdrawdown(portValue);
MAR = 0;
sortinoRatio = (mean(percentReturn) - MAR) / sqrt(lpm(percentReturn, MAR, 2));
load SPXret19902000.mat SPXret19902000
beta = regress(percentReturn, SPXret19902000);
warning('off');
adjustedAlpha = portalpha(percentReturn, SPXret19902000);
result0 = [series2 res indicate s r0 r sumr percentReturn portValue];

if nargout == 0
    %% Plot results
    ax(1) = subplot(3,1,1);
    plot(series2), grid on
    legend('Stock1','Stock2')
    title(['Pairs trading results, Sharpe Ratio = ',num2str(sh,3)])
    ylabel('Price (USD)')
    
    ax(2) = subplot(3,1,2);
    plot([indicate,spreadH*ones(size(indicate)),-spreadH*ones(size(indicate)),spreadL*ones(size(indicate)),-spreadL*ones(size(indicate))])
    grid on
    legend(['Indicator'],'Stock1: Over bought','Stock1: Over sold','Go Flat','Go Flat',...
        'Location','NorthWest')
    title(['Pairs indicator: rebalance every ' num2str(N)...
        ' minutes with previous ' num2str(M) ' minutes'' prices.'])
    ylabel('Indicator')
    
    ax(3) = subplot(3,1,3);
    plot([s,cumsum(r)]), grid on
    legend('Position for Stock1','Position for Stock2','Cumulative Return',...
        'Location', 'NorthWest')
    title(['Final Return = ',num2str(sum(r),3),' (',num2str(sum(r)/capital*100,3),'%)'])
    ylabel('Return (USD)')
    xlabel('Serial time number')
    linkaxes(ax,'x')
else
    %% Return values
    for i = 1:nargout
        switch i
            case 1
                varargout{1} = s; % signal
            case 2
                varargout{2} = r0;
            case 3
                varargout{3} = adjustedAlpha; % return (pnl)
            case 4
                varargout{4} = cVar95; % sharpe ratio
            case 5
                varargout{5} = cVar99; % indicator
            case 6
                varargout{6} = result0;
            otherwise
                warning('PAIRS:OutputArg',...
                    'Too many output arguments requested, ignoring last ones');
        end 
    end
end