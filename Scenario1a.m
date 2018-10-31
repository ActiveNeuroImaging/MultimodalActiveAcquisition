%A foolish consistency is the hobgoblin of little minds
%The best is the enemy of the good (apologies for mistakes in the code)


%A very simple adaptive case where outlier is calculated and then a
%decision taken starting with lower resolution; can be used to decide
%whether a higher resolution scan is required. More 

clear all

format compact;

BigMatrix=[]; 
temp_img=[];
FileNames=[];
NScans=3; %Number of different resolution scans available
Scan1=[];
Scan2=[];
Scan3=[];

%Load nifti T1-weighted structural images, acquired at three resolutions,
%for controls and lesion patients and create a vector labelling
cd ~/T1s

count=0;
for i=1:13
    count=count+1;
    tempname=['S' num2str(i) '_1.nii'];
    temp_img=load_nii(tempname); 
    Scan1(:,:,:,count)=temp_img.img; 
    tempname=['S' num2str(i) '_2.nii'];
    temp_img=load_nii(tempname); 
    Scan2(:,:,:,count)=temp_img.img; 
    tempname=['S' num2str(i) '_3.nii'];
    temp_img=load_nii(tempname); 
    Scan3(:,:,:,count)=temp_img.img; 
    FileNames{count}=tempname;
    SubNum(count)=i;
    if i<8 
        Patient(count)=0;
    else
        Patient(count)=1;
    end
    
   
end
Scans{1}=Scan1;
Scans{2}=Scan2;
Scans{3}=Scan3;

temp_img=load_nii('MNI152_T1_2mm_brain.nii.gz'); 
brain_mask=temp_img.img; % imresize 



rng default;
%i=4;
N_SubBrains=3;


%Loop through each participant
for SubN=1:13;
 
Cord1=1;
Cord2=1;

KeepDataVect=[];
ToInclude=[1:count];
ToInclude(SubN)=[];

TempScans=Scans;
%Start at lowest resolution scan and then increase resolution
for i=NScans:-1:1
ShrinkingBrain=TempScans{i};
TempMask=brain_mask;
clear TempBigMatrixRedux TempBigMatrixHoldout TempBrainChunk;
Size_ShrinkingBrain=size(ShrinkingBrain);
%Split into N_SubBrains equally sized volumes  
temp_N_SubBrains=N_SubBrains;

        if Size_ShrinkingBrain(3)<N_SubBrains
             temp_N_SubBrains=Size_ShrinkingBrain(3);
        end

        if temp_N_SubBrains==0
            temp_N_SubBrains=1;
        end
        clear TempZLimits;
        for p=1:temp_N_SubBrains        
            tempStart=(p-1)*round(Size_ShrinkingBrain(3)./temp_N_SubBrains);
            if tempStart<1
                tempStart=1;
            end
            
            tempEnd=(p)*round(Size_ShrinkingBrain(3)./temp_N_SubBrains);
            if tempEnd>Size_ShrinkingBrain(3)
                tempEnd=Size_ShrinkingBrain(3);
            end
            if tempEnd==0
                tempEnd=Size_ShrinkingBrain(3);
            end
            for q=1:NScans
                TempBrainChunk{p,q}=TempScans{q}(:,:,tempStart:tempEnd,:,:); 
                TempMaskChunk{p}=TempMask(:,:,tempStart:tempEnd);
            end
            TempZLimits{p}=[tempStart tempEnd];
        end
        
        
        
        
    clear Outliers TempBigMatrixRedux; 
    
 
 %Assess outlierness within each subvolume for each individual.
 
    for p=1:temp_N_SubBrains    
        clear TempBigMatrixRedux;
        for j=1:size(ShrinkingBrain,4)
            TempBigMatrixRedux(:,:,:,j)=TempBrainChunk{p,i}(:,:,:,j);
        end
            
            for s=1:13
                NonStroke=find(Patient==0);
                NonStroke(find(NonStroke==s))=[];

                Outliers(s,p)=OutlierDetect3D(TempBigMatrixRedux(:,:,:,NonStroke),TempBigMatrixRedux(:,:,:,s),TempMaskChunk{p});
            end
            
    
    end
    
    %Convert outlierness into a z-score so can be compared even if
    %different average voxel intensities or different numbers of voxsels
      ControlIndices=[1:7];
    if find(ControlIndices==SubN)
       ControlIndices(3)=[]; 
    end
    if SubN>7
        tempRandIndices=randperm(7);
        ControlIndices=tempRandIndices(1:6);
    end
    MeanControl=squeeze(mean(Outliers(ControlIndices,:)));
    StdControl=squeeze(std(Outliers(ControlIndices,:)));
    ZScored(SubN,i,:)=(squeeze(Outliers(SubN,:))-MeanControl)./StdControl;
   
    [aa,bb]=max(ZScored(SubN,i,:));
   
   
   Cord3=bb;
  
  
   clear TempBigMatrixRedux; 
   
   
  %Keep subvolume which was quantified as maximum outlier 
   KeepDataCollected{SubN,i}=TempBrainChunk{Cord3,i}(:,:,:,SubN);
   for r=1:NScans
        
        TempScans{r}=TempScans{r};
   end


  KeepSlice(SubN,i)=Cord3;

  KeepZLimits{SubN,i}=TempZLimits{Cord3};
  
  
      
end

  
end

%Plot outlier values for maximum outlier subvolume for each participant
scatter(ones(1,7),max(squeeze(ZScored(1:7,3,:))'),'b','o')
hold on
scatter(ones(1,6)*1.1,max(squeeze(ZScored(8:13,3,:))'),'k','+')
scatter(ones(1,7)*2,max(squeeze(ZScored(1:7,2,:))'),'b','o')
scatter(ones(1,6)*2.1,max(squeeze(ZScored(8:13,2,:))'),'k','+')
scatter(ones(1,7)*3,max(squeeze(ZScored(1:7,1,:))'),'b','o')
scatter(ones(1,6)*3.1,max(squeeze(ZScored(8:13,1,:))'),'k','+')


scatter(ones(1,7)*4.50,(max(squeeze(ZScored(1:7,1,:))')+max(squeeze(ZScored(1:7,2,:))')+max(squeeze(ZScored(1:7,3,:))'))/3,'b','o')
scatter(ones(1,6)*4.6,(max(squeeze(ZScored(8:13,1,:))')+max(squeeze(ZScored(8:13,2,:))')+max(squeeze(ZScored(8:13,3,:))'))/3,'k','+')
hold off

xlim([0.5 5]);
xticks([1.05 2.05 3.05 4.55]);
xticks('manual')

xticklabels({'\fontsize{14}Scan 1','\fontsize{14}Scan 2','\fontsize{14}Scan 3','\fontsize{14}All scans'})
ylabel('\fontsize{16}Outlier strength')
legend('\fontsize{14}Control', '\fontsize{14}Patient');