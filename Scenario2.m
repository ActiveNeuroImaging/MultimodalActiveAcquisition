%A foolish consistency is the hobgoblin of little minds
%The best is the enemy of the good (apologies for mistakes in the code)


%A simple illustration with age prediction, involving a decision tree
%approach to simulate active acquisition

clear all
format compact;
rng(2);

%% Load data, remove missing data, find particpant ages, collect summary variables from the modalities
cd ~/Camcan/ 
camcan=readtable('CamcanSummaryFirstPassWithAgeEtc.csv');
camcanRedux= rmmissing(camcan);
Ages=camcanRedux.age;
sizecamcan=size(camcanRedux);
camcanTemp(:,1)=camcanRedux.T2_WM;
camcanTemp(:,2)=camcanRedux.Mean_FA;
camcanTemp(:,3)=camcanRedux.BoldAct_AudioVis;
camcanTemp(:,4)=camcanRedux.Rest6;
camcanTemp(:,5)=camcanRedux.Movie6;
camcanTemp(:,6)=camcanRedux.GM;

%Define a holdout dataset
HoldoutP=0.2;
Indices = crossvalind('HoldOut', sizecamcan(1), HoldoutP);
trainIndices=find(Indices==1);
testIndices=find(Indices==0);

%Normalise the data, this could probably be done much better 
trainGM=zscore(camcanTemp(trainIndices,:));
testGM=zscore(camcanTemp(testIndices,:));


%Train the decision tree regression model with optimized hyperparameters
%and for comparison train a 
trainGM=[];
testGM=[];
DTModel=fitrtree(trainGM,Ages(trainIndices),'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));
SVRModel=fitrsvm(trainGM(:,:),Ages(trainIndices),'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));

view(DTModel,'Mode','graph'); %View the decision tree structure

%Calculate the predicted age differences for the decision tree and support
%vector regressions
MeanPAD_DT=mean(abs(predict(DTModel,testGM)-Ages(testIndices))); %Mean predicted age difference based on decision tree (i.e., simulating adaptive acquisition).
MedianAgesDiff=median(abs(predict(DTModel,testGM)-Ages(testIndices)));
  
MeanPAD_SVR=mean(abs(predict(SVRModel,testGM)-Ages(testIndices))); %Mean predicted age difference based on support vector regression (i.e., whole input dataset).
MedianPAD_SVR=median(abs(predict(SVRModel,testGM(:,:))-Ages(testIndices)));

