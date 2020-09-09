function DCA_MOEAD(Global)
% <algorithm> <A>
% MOEA/D based on double hybrid
% nr    ---   2 --- Maximum number of solutions replaced by each offspring
% min   ---   0.05 --- Minimum probability for each operator
% alpha     ---  0.5 --- Learning Rate
% K         ---  5 --- Number of evolution step

    %% Parameter setting
    [nr,p_min,alpha,K] = Global.ParameterSet(2,0.05,0.5,5);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    for i = 1 : Global.N
        Population(i).adds(Population(i).dec);
    end
    Z = min(Population.objs,[],1);
    
    %% Define global parameters
    %对应的四个算子的顺序：SBX、DE、NCG、PCG
    probability = zeros(Global.maxgen,4);
    quality = zeros(Global.maxgen,4);
    moeadPerformance = zeros(Global.maxgen,4);
    nsgaPerformance = zeros(Global.maxgen,4);
    usedNumber = zeros(Global.maxgen,4);
    % External archive
    Archive = Population;            

    %% Optimization
    generation = 0;
    stepNumber = floor(Global.maxgen/K);
    while Global.NotTermination(Archive)
        tic;
        generation = generation + 1;% 第generation代进化开始       
        offspringPool = [];% 新产生的子代     
        strategyUsed = zeros(1,Global.N);% 个体在本代采用的进化策略               
        if(mod(generation-1,stepNumber)==0)
            probability(generation,:) = [0.25,0.25,0.25,0.25];
        else
            [probability(generation,:),quality(generation,:)] = CalcProbability...
            (quality(generation-1,:),moeadPerformance(generation-1,:),nsgaPerformance(generation-1,:),p_min,alpha);
        end
        
        for i = 1 : Global.N
            strategy = OperatorSelector(probability(generation,:));
            strategyUsed(i) = strategy;
            if strategy == 1 || strategy == 2
                if rand < 0.9
                    P = B(i,randperm(end));                  
                else
                    P = randperm(Global.N);
                end
                Offspring = GenerateOffspring(strategy,Population(i).dec,...
                    Population(P(1)).dec,Population(P(2)).dec);
            else                
                P = B(i,randperm(end));
                crossInd = Population(P(1));
                if strategy == 4
                    P = randperm(Global.N);
                end
                Offspring = GenerateOffspring(strategy,Population(i).dec,...
                    Population(P(2)).dec,Population(P(3)).dec,crossInd.dec,crossInd.add);
            end
           
            Z = min(Z,Offspring.obj);
            % Update the solutions in P by Tchebycheff approach
            g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
            g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);           

            Population(P(find(g_old > g_new,nr))) = Offspring;
            offspringPool = [offspringPool,Offspring];
            
            % Calculate the improve rate.
            IMPROVE_MAT = ((g_old - g_new)./g_old);
            impRate = mean(IMPROVE_MAT(g_old > g_new));
            if isnan(impRate)
                impRate = 0;
            end
            moeadPerformance(generation,strategy)=moeadPerformance(generation,strategy)+impRate;
        end
        [Archive,successful] = UpdateArchive(Archive,offspringPool);
        % 计算最终的平均奖励
        for s = 1 : 4
            if ~isnan(find(strategyUsed==s))
                usedNumber(generation,s)=length(find(strategyUsed==s));
                nsgaPerformance(generation,s) = ...
                    length(find(strategyUsed(successful)==s))/Global.N;
            end         
        end

    end
end