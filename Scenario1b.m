%A foolish consistency is the hobgoblin of little minds
%The best is the enemy of the good (apologies for mistakes in the code)


%A simple approach similar to Scenario 1a, simulating actively "zooming" in to parts of an image which are
%more likely to be outliers, and "acquiring" higher resolution data.

clear all

format compact;

BigMatrix=[]; 
temp_img=[];
FileNames=[];
NScans=3; %Number of different scan resolutions acquired
Scan1=[];
Scan2=[];
Scan3=[];


%Load nifti T1-weighted structural images, acquired at three resolutions,
%for controls and lesion patients and create a vector labelling
%patients/controls
cd ~/T1s
count=0;

for i=1:13
    count=count+1;
    tempname=['S' num2str(i) '_1.nii']; %1 is highest resolution (1mmx1mmx1mm)
    temp_img=load_nii(tempname); 
    Scan1(:,:,:,count)=temp_img.img; 
    tempname=['S' num2str(i) '_2.nii']; %2 is middle resolution (2mmx2mmx2mm)
    temp_img=load_nii(tempname); 
    Scan2(:,:,:,count)=temp_img.img; 
    tempname=['S' num2str(i) '_3.nii']; %3 is middle resolution (2mmx2mmx2mm)
    temp_img=load_nii(tempname); 
    Scan3(:,:,:,count)=temp_img.img; 
    FileNames{count}=tempname;
    SubNum(count)=i;

    
    if i<8 
        Patient(count)=0; %Controls
    else
        Patient(count)=1; %Patients
    end
    
   
end
Scans{1}=Scan1;
Scans{2}=Scan2;
Scans{3}=Scan3;

temp_img=load_nii('MNI152_T1_2mm_brain.nii.gz'); 
brain_mask=temp_img.img; %Base brain mask on MNI template



rng default;

N_SubBrains=3; %Number of divisions to split brain into when "zooming


%Loop through subjects
for SubN=1:13;
 
Cord1=1;
Cord2=1;
%Start with whole brain, separate the volume (arbitraily) into thirds, then keep the reduced field which is quantified to be most likey at outlier and then repeat.

KeepDataVect=[];
ToInclude=[1:count];
ToInclude(SubN)=[];
TempScans=Scans;
%Iterate through scans, each iteration finding portion of the data that is greatest outlier and keeping for further analysis. 
for i=NScans:-1:1
    ShrinkingBrain=TempScans{i};
    TempMask=brain_mask;
    clear TempBigMatrixRedux TempBigMatrixHoldout TempBrainChunk;
    Size_ShrinkingBrain=size(ShrinkingBrain);
    temp_N_SubBrains=N_SubBrains;

        if Size_ShrinkingBrain(3)<N_SubBrains
             temp_N_SubBrains=Size_ShrinkingBrain(3);
        end

        if temp_N_SubBrains==0
            temp_N_SubBrains=1;
        end
        clear TempZLimits;
        %Split into 1/N_Subrain volumes
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
    
    %Loop through the three reduced sub-volumes
    for p=1:temp_N_SubBrains    
        clear TempBigMatrixRedux;
        for j=1:size(ShrinkingBrain,4)
            TempBigMatrixRedux(:,:,:,j)=TempBrainChunk{p,i}(:,:,:,j);
        end

        %Quantify outlierness for each volume compared to "normative" data
        %(i.e., 6 of the healthy controls) for each participant for the
        %specified volume. Do for each subject so that it can the outlier
        %distance can be adjusted (i.e., z-scored below) to correct for
        %differences in the number of voxels/ voxel intensity differences for different volumes. 
        
        for s=1:13
            NonStroke=find(Patient==0);
            NonStroke(find(NonStroke==s))=[];

            Outliers(s,p)=OutlierDetect3D(TempBigMatrixRedux(:,:,:,NonStroke),TempBigMatrixRedux(:,:,:,s),TempMaskChunk{p});
        end
    
    end
    
    
    %Calculate z-stat for the measure of outlierness specific to each
    %volume (this could be done in many other ways)
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

   
   %Keep the chosen (i.e., highest outlier) volume
   KeepDataCollected{SubN,i}=TempBrainChunk{Cord3,i}(:,:,:,SubN);
   for r=1:NScans
        TempScans{r}=TempBrainChunk{Cord3,r};
   end


  KeepSlice(SubN,i)=Cord3;
  KeepZLimits{SubN,i}=TempZLimits{Cord3};
  
  
      
end

%Save the chosen field of view using standard commands
temp_brain_Base=Scan3(:,:,:,SubN);
temp_nii_Base=make_nii(temp_brain_Base);
save_nii(temp_nii_Base,['SavedBrains/BrainBase_' num2str(SubN)])

temp_brain_3=Scan1(:,:,:,1).*0;
LowerZ=KeepZLimits{SubN,3}(1);
UpperZ=KeepZLimits{SubN,3}(2);
temp_brain_3(:,:,LowerZ:UpperZ)=KeepDataCollected{SubN,3};
temp_nii_3=make_nii(temp_brain_3);
save_nii(temp_nii_3,['SavedBrains/Brain_' num2str(SubN) '_' num2str(3)])

temp_brain_2=Scan1(:,:,:,1).*0;
LowerZ=KeepZLimits{SubN,3}(1)+KeepZLimits{SubN,2}(1);
UpperZ=KeepZLimits{SubN,3}(1)+KeepZLimits{SubN,2}(2);
temp_brain_2(:,:,LowerZ:UpperZ)=KeepDataCollected{SubN,2};
temp_nii_2=make_nii(temp_brain_2);
save_nii(temp_nii_2,['SavedBrains/Brain_' num2str(SubN) '_' num2str(2)])

temp_brain_1=Scan1(:,:,:,1).*0;
LowerZ=KeepZLimits{SubN,3}(1)+KeepZLimits{SubN,2}(1)+KeepZLimits{SubN,1}(1);
UpperZ=KeepZLimits{SubN,3}(1)+KeepZLimits{SubN,2}(1)+KeepZLimits{SubN,1}(2);
temp_brain_1(:,:,LowerZ:UpperZ)=KeepDataCollected{SubN,1};
temp_nii_1=make_nii(temp_brain_1);
save_nii(temp_nii_1,['SavedBrains/Brain_' num2str(SubN) '_' num2str(1)])

TotalBrain=temp_brain_Base+temp_brain_3+temp_brain_2+temp_brain_1;
temp_nii_total=make_nii(TotalBrain);
save_nii(temp_nii_total,['SavedBrains/TotalBrain_' num2str(SubN)])
 
end



%Display results of outlier analysis for part of the image that is
%quantified as most abnormal

scatter(ones(1,7),max(squeeze(ZScored(1:7,3,:))'),'b','o')
hold on
scatter(ones(1,6)*1.1,max(squeeze(ZScored(8:13,3,:))'),'k','+')
scatter(ones(1,7)*2,max(squeeze(ZScored(1:7,2,:))'),'b','o')
scatter(ones(1,6)*2.1,max(squeeze(ZScored(8:13,2,:))'),'k','+')
scatter(ones(1,7)*3,max(squeeze(ZScored(1:7,1,:))'),'b','o')
scatter(ones(1,6)*3.1,max(squeeze(ZScored(8:13,1,:))'),'k','+')

hold off

xlim([0.5 3.5]);
xticks([1.05 2.05 3.05]);
xticks('manual')
xticklabels({'\fontsize{14}Scan 1','\fontsize{14}Scan 2','\fontsize{14}Scan 3'})
ylabel('\fontsize{16}Outlier strength')
legend('\fontsize{14}Control', '\fontsize{14}Patient');
