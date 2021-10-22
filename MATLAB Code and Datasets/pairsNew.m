function varargout = pairsNew(SPX218, selectedPairs, M, N, spreadH, spreadL, scaling, cost, trainSet, capital)
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

if ~exist('trainSet', 'var')
    trainSet = true;
end

if nargin == 2
    % default values
    M = 378; % lookback period
    N = 63; % default rebalance period
elseif nargin == 3
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

series2 = zeros(height(SPX218),2);
s = zeros(height(SPX218),2);
indicate = zeros(height(SPX218),1);
res = zeros(height(SPX218),1);
reg0 = zeros(height(SPX218),3);
rebalanceDayCounter = zeros(height(SPX218),1);
portValue = capital*ones(height(SPX218),1);
tradeState = false;
lastDay = false;
pxChange = zeros(height(SPX218),2);
sChange = zeros(height(SPX218),2);
r0 = zeros(height(SPX218),2);
r  = zeros(height(SPX218),1);
percentReturn  = zeros(height(SPX218),1);
stat_pre = zeros(size(selectedPairs,1),1);
coolDown = 0;
dailyLossThreshold = -0.05;
periodDDThreshold = 0.15;
coolDownPeriod = N;
goodTrade = 0;
badTrade = 0;

if trainSet == true
    p = size(SPX218(SPX218.Date <= datetime('12/31/1992','InputFormat','MM/dd/uuuu'), 1),1);
else
    p = size(SPX218(SPX218.Date <= datetime('12/31/1996','InputFormat','MM/dd/uuuu'), 1),1);
end
startPeriod = p+1;

for i = startPeriod:length(s)
    if coolDown <= 0
        if lastDay == true
            rebalanceDayCounter(i) = 1;
            lastDay = false;
        else
            rebalanceDayCounter(i) = rebalanceDayCounter(i-1) + 1;
        end
        if rebalanceDayCounter(i) == N
            lastDay = true;
        end
        if rebalanceDayCounter(i) == 1
            for c1 = 1:size(selectedPairs,1)
                stock1 = selectedPairs{c1,1};
                stock2 = selectedPairs{c1,2};
                series2_temp(:, 1) = SPX218{i-M:i-1,stock1+1};
                series2_temp(:, 2) = SPX218{i-M:i-1,stock2+1};
                [~,~,stat_pre(c1),~,~] = egcitest(series2_temp);
            end
            [~, bestPairIdx] = min(stat_pre);
            stock1 = selectedPairs{bestPairIdx,1};
            stock2 = selectedPairs{bestPairIdx,2};
            series2(:, 1) = SPX218{:,stock1+1};
            series2(:, 2) = SPX218{:,stock2+1};
            [h,~,~,~,reg1] = egcitest(series2(i-M:i-1, :));
            if reg1.coeff(2)*series2(i-1,2) < 0.5*series2(i-1,1) || reg1.coeff(2)*series2(i-1,2) > 2*series2(i-1,1)
                h = 0;
            end
            s(i, 2) = 0;
            s(i, 1) = 0;
            tradeState = false;
            if i > length(s) - N + 1
                h = 0;
            end
            if s(i-1, 2) ~= 0
                badTrade = badTrade + 1;
            end
        end
        if h ~= 0
            res(i) = series2(i-1,1) - (reg1.coeff(1) + reg1.coeff(2).*series2(i-1,2));
            reg0(i, 1) = reg1.coeff(1);
            reg0(i, 2) = reg1.coeff(2);
            reg0(i, 3) = reg1.RMSE;
            indicate(i) = res(i)./reg1.RMSE;
            if tradeState == false
                if indicate(i) > spreadH
                    s(i, 2) = portValue(i-1)/(series2(i-1,2) + 1/reg1.coeff(2)*series2(i-1,1));
                    s(i, 1) = -1/reg1.coeff(2) * s(i, 2);
                    tradeState = true;
                    rebalanceDayCounter(i) = 1;
                elseif indicate(i) < -spreadH
                    s(i, 2) = -portValue(i-1)/(series2(i-1,2) + 1/reg1.coeff(2)*series2(i-1,1));
                    s(i, 1) = -1/reg1.coeff(2) * s(i, 2);
                    tradeState = true;
                    rebalanceDayCounter(i) = 1;
                else
                    s(i, 2) = 0;
                    s(i, 1) = 0;
                    tradeState = false;
                end
            else
                if abs(indicate(i)) < spreadL
                    s(i, 2) = 0;
                    s(i, 1) = 0;
                    tradeState = false;
                    goodTrade = goodTrade + 1;
                    lastDay = true;
                else
                    s(i, 2) = s(i-1, 2);
                    s(i, 1) = s(i-1, 1);
                    tradeState = true;
                end
            end
        end
        pxChange(i,:) = series2(i,:) - series2(i-1,:);
        sChange(i,:) = s(i,:) - s(i-1,:);
        r0(i,:) = s(i,:).* pxChange(i,:) - abs(sChange(i,:))*cost/2;
        r(i) = sum(r0(i,:), 2);
        percentReturn(i) = r(i)/portValue(i-1);
        if percentReturn(i) < dailyLossThreshold
            coolDown = coolDownPeriod;
        end
        portValue(i) = portValue(i-1) + r(i);
        maxDD_temp = maxdrawdown(portValue(i-N:i-1));
        if maxDD_temp > periodDDThreshold
            coolDown = coolDownPeriod;
        end
    else
        coolDown = coolDown - 1;
    end
end

%% Calculate performance statistics

sumr = cumsum(r);
percentReturnL = percentReturn;
percentReturn = percentReturn(startPeriod:end);
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
maxDD = maxdrawdown(portValue);
MAR = 0;
sortinoRatio = (mean(percentReturn) - MAR) / sqrt(lpm(percentReturn, MAR, 2));

if trainSet == true
    load 'SPXret19931996.mat';
    SPXret = SPXret19931996;
else
    load 'SPXret19972000.mat';
    SPXret = SPXret19972000;
end

beta = regress(percentReturn, SPXret);
warning('off');
adjustedAlpha = portalpha(percentReturn, SPXret);
result0 = [series2 res indicate s r0 r sumr percentReturnL portValue];

%if nargout == 0
    %% Plot results
    ax(1) = subplot(3,1,1);
    plot([series2(startPeriod-M:end,1),series2(startPeriod-M:end,2)]), grid on
    legend('Stock1','Stock2')
    title(['Pairs trading results, Sharpe Ratio = ',num2str(sh,3)])
    ylabel('Price (USD)')
    
    ax(2) = subplot(3,1,2);
    indicate_temp = indicate(startPeriod-M:end,1);
    plot([indicate_temp,spreadH*ones(size(indicate_temp)),-spreadH*ones(size(indicate_temp)),spreadL*ones(size(indicate_temp)),-spreadL*ones(size(indicate_temp))])
    grid on
    legend(['Indicator'],'Stock1: Over bought','Stock1: Over sold','Go Flat','Go Flat',...
        'Location','NorthWest')
    title(['Pairs indicator: rebalance every ' num2str(N)...
        ' days with previous ' num2str(M) ' days'' prices.'])
    ylabel('Indicator')
    
    ax(3) = subplot(3,1,3);
    plot([s(startPeriod-M:end,1),s(startPeriod-M:end,2),sumr(startPeriod-M:end)])
    grid on
    legend('Position for Stock1','Position for Stock2','Cumulative Return',...
        'Location', 'NorthWest')
    title(['Final Return = ',num2str(sum(r),3),' (',num2str(sum(r)/capital*100,3),'%)'])
    ylabel('Return (USD)')
    xlabel('Serial time number')
    linkaxes(ax,'x')
%else
    %% Return values
    for i = 1:nargout
        switch i
            case 1
                varargout{1} = s; % signal
            case 2
                varargout{2} = r; % return (pnl)
            case 3
                varargout{3} = sh; % sharpe ratio
            case 4
                varargout{4} = indicate; % indicator
            case 5
                varargout{5} = percentReturn;
            case 6
                varargout{6} = result0;                
            otherwise
                warning('PAIRS:OutputArg',...
                    'Too many output arguments requested, ignoring last ones');
        end 
    end
%end