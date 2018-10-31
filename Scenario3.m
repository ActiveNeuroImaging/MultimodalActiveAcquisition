%A foolish consistency is the hobgoblin of little minds
%The best is the enemy of the good (apologies for mistakes in the code)


%A simple illustration of active acquisition for outlier detection using
%Bayesian optimization and summary measures of different modalities from
%the Camcan dataset

clear all
format compact;
rng(default);

%% Load data, remove missing data, find particpant ages, collect summary variables from the modalities
cd ~/Camcan/ 
camcan=readtable('CamcanSummaryFirstPassWithAgeEtc.csv');
camcanRedux= rmmissing(camcan);
sizecamcan=size(camcanRedux);
camcanTemp(:,1)=camcanRedux.T2_WM;
camcanTemp(:,2)=camcanRedux.Mean_FA;
camcanTemp(:,3)=camcanRedux.BoldAct_AudioVis;
camcanTemp(:,4)=camcanRedux.Rest6;
camcanTemp(:,5)=camcanRedux.Movie6;
camcanTemp(:,6)=camcanRedux.GM;


%Define a holdout dataset
camcanTemp=camcanTemp;
HoldoutP=0.2;
sizecamcan=size(camcanTemp2);
Indices = crossvalind('HoldOut', sizecamcan(1), HoldoutP);

trainIndices=find(Indices==1);
testIndices=find(Indices==0);


%Do a factor analysis on the z-transformed training data to define a
%one-dimensional ordering of modalities which can be used to define a (very
%simple) search space for the optimization based on covariance between
%modalities (i.e., nearby modalities are closer in the search space).


[Lambda,psi,T,stats]=factoran(zscore(camcanTemp(trainIndices,:)),1)
[SortedFactorScores,reorderTasks]=sort(Lambda); 
camcanTemp2=camcanTemp(:,reorderTasks);
imagesc(corr(camcanTemp2));  %Visualise the correlations between 

% Normalise the test data with mean and standard deviation of training data
A=camcanTemp2(trainIndices,:); 
C = bsxfun(@minus, A, mean(A));
camcanNormative = bsxfun(@rdivide, C, std(A));
B=camcanTemp3(testIndices,:);
C = bsxfun(@minus, B, mean(A));
camcanTest= bsxfun(@rdivide, C, std(A));

%Set up Bayesian optimization variable (modality) to be searched over.
ModalityParam= optimizableVariable('Modality',[1,6],'Type','integer');

%Run the optimization for each participant in the training set a 100 times,
%to compare active acquisition based on Bayesian optimization. The
%optimization is run several times because the burn in phase has a random element and the random
%ordering of the modalities in the comparison optimization is also random. 

for Iterations=1:100

    for i=1:size(testIndices,1)
    

        ModalityReorder=randperm(size(camcanTest,2)); %Create a random modality ordering to compare the true ordering against
   
        fun = @(x)camcanFuncForBOpt(i,x.Modality,camcanNormative, camcanTest); %Define Bayesian optimization with the true ordering (based on actual covariance between modalities)
        funWrong = @(x)camcanFuncForBOptWrong(i,x.Modality,camcanNormative, camcanTest,ModalityReorder); %Define the same optimization but with a random ordering not based on true covariance).

        results{i} = bayesopt(fun,[ModalityParam],'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',4,'NumSeedPoints',3,'PlotFcn' ,[],'IsObjectiveDeterministic',true); %Run the Bayesian optimization for the real ordering
        resultsWrong{i} = bayesopt(funWrong,[ModalityParam],'Verbose',0,'AcquisitionFunctionName','expected-improvement-plus','ExplorationRatio',0.5,'MaxObjectiveEvaluations',4,'NumSeedPoints',3,'PlotFcn' ,[],'IsObjectiveDeterministic',true); %Run the Bayesian optimization for the random ordering

   
        [BestScore(i) BestModality(i)]=min(camcanTest(i,:)); %Compare the true 
        [BestScoreWrong(i) BestModalityWrong(i)]=min(camcanTest(i,ModalityReorder));    

end

%Quantify average performance for each iteration in terms of proportion
%of most abnormal modality correctly estimated and the most abnormal
%measurement (in terms of z-score) calculated.
for k=1:size(results,2)
	MinIndex(k)=find(results{k}.ObjectiveMinimumTrace==min(results{k}.ObjectiveMinimumTrace),1);
   
    MinEstimated(k)=results{k}.EstimatedObjectiveMinimumTrace(end);
    MinEstimatedWrong(k)=resultsWrong{k}.EstimatedObjectiveMinimumTrace(end);
    MinSampled(k)=results{k}.MinObjective;
    MinSampledWrong(k)=resultsWrong{k}.MinObjective;
    tempTable=results{k}.bestPoint;
    MinEstimatedPoint(k)=tempTable.Modality;
    tempTable=resultsWrong{k}.bestPoint;
    MinEstimatedPointWrong(k)=tempTable.Modality;
    
    
    
    
end

KeepIndex=[];
KeepIndexWrong=[];
for i=1:size(results,2); if sum(unique(results{i}.XTrace.Modality)>0)==4; KeepIndex(i)=1; end; end
for i=1:size(resultsWrong,2); if sum(unique(resultsWrong{i}.XTrace.Modality)>0)==4; KeepIndexWrong(i)=1; end; end


MeanRightUnique(Iterations)=mean(MinEstimatedPoint(find(KeepIndex))==BestModality(find(KeepIndex)));
MeanWrongUnique(Iterations)=mean(MinEstimatedPointWrong(find(KeepIndexWrong))==BestModalityWrong(find(KeepIndexWrong)));
MeanRight(Iterations)=mean(MinEstimatedPoint==BestModality);
MeanWrong(Iterations)=mean(MinEstimatedPointWrong==BestModalityWrong);
KeepRight(Iterations)=sum(KeepIndex);
KeepWrong(Iterations)=sum(KeepIndexWrong);

MeanEstimate(Iterations)=mean(MinEstimated-MinEstimatedWrong);
CountEstimateRight(Iterations)=mean(MinEstimated<MinEstimatedWrong);
CountEstimateWrong(Iterations)=mean(MinEstimated>MinEstimatedWrong);
CountEstimateSame(Iterations)=mean(abs(MinEstimated-MinEstimatedWrong)<0.00001);



end


for k=1:size(find(KeepIndex),2); BestPoints(k)=results{k}.bestPoint.Modality; end



% Plot results (i.e., estimated most abnormal modality versus measured most
% abnormal modality) for the different approaches and different
% randomisations.

ii=16;disp(results{ii}.EstimatedObjectiveMinimumTrace');  disp(BestModality(ii)');disp(results{ii}.XTrace);
ObjectiveResults=camcanTest(ii,:)
plot(ObjectiveResults,'*k')
xlim([0.5 6.5]);

xticks([1.0 2.0 3.0 4.0 5.0 6.0]);

xticks('manual')

xticklabels({'\fontsize{12}Movie DMN','\fontsize{12}Rest DMN','\fontsize{14}Task BOLD','\fontsize{14}T2','\fontsize{14}FA','\fontsize{14}GM'})
ylabel('\fontsize{16}Z score')


objective = predictObjective(results{ii},results{ii}.XTrace)

T=results{ii}.XTrace;
T([4:end],:)=[];
objective = predictObjective(results{ii},T)


x1=MeanRightUnique;
y1=1:size(x1,2);


x2=MeanWrongUnique;
y2=[1:size(x1,2)]+50+size(x1,2);

x3=MeanRightUnique-MeanWrongUnique;
y3=[1:size(x1,2)]+50+y2(end);


x4=MeanEstimate;
y4=[1:size(x1,2)]+50+y2(end);




yyaxis left
ylabel('\fontsize{14}Proportion maximum modality')

yyaxis right
ylabel('\fontsize{14}Z-score improvement with oåptimization')
hold on

yyaxis left

scatter(y1,x1);
scatter(y2,x2);

yyaxis right

scatter(y4,x4);
xlim([-50 450]);
xticks([50 200 350]);
xticks('manual')

xticklabels({'\fontsize{14}True search space','\fontsize{14}Random search space','\fontsize{14}True-Random'})


hold off






