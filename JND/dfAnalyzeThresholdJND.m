function Mean = dfAnalyzeThresholdJND(UD, varargin)
% modified from Palam
% Received from Ayoub Daliri
HighReversal = max(UD.reversal);
NumTrials = length(UD.response);

if ~isempty(varargin)
    NumOpts = length(varargin);
    for n = 1:2:NumOpts
        valid = 0;
        if strncmpi(varargin{n}, 'reversals',4)            
            criterion = 'reversals';
            number = varargin{n+1};
            LowReversal = HighReversal - number + 1;
            if LowReversal < 1
                warning('PALAMEDES:invalidOption','You asked for the last %s reversals. There are only %s reversals.',int2str(number),int2str(HighReversal));               
                LowReversal = 1;
            end
            valid = 1;
        end
        if strncmpi(varargin{n}, 'trials',4)            
            criterion = 'trials';
            number = varargin{n+1};
            LowTrial = NumTrials - number + 1;
            if LowTrial < 1
                warning('PALAMEDES:invalidOption','You asked for the last %s trials. There are only %s trials.',int2str(number),int2str(NumTrials));               
                LowTrial = 1;
            end
            valid = 1;
        end
        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})
        end        
    end            
else
    criterion = 'reversals';
    LowReversal = 3;
end

if strncmpi(criterion,'reversals',4)
    Mean = sum(UD.xStaircase(UD.reversal >= LowReversal))/(HighReversal-LowReversal+1);
else
    Mean = sum(UD.xStaircase(LowTrial:NumTrials))/(NumTrials-LowTrial+1);
end

figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',14,...
    'FontName','Arial');
box(axes1,'on');
hold(axes1,'all');

plot(UD.x,'LineWidth',3);
hold on;
plot(find(UD.response==1),UD.x(find(UD.response==1)),'o','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1],'MarkerSize',10);
plot(find(UD.response==0),UD.x(find(UD.response==0)),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerSize',10);
line([0 length(UD.response)], [Mean Mean],'LineStyle', '-.', 'LineWidth',3,'color',[1 0 1])
xlabel('Trials','FontSize',20,'FontName','Arial');
ylabel('Perturbations','FontSize',20,'FontName','Arial');

