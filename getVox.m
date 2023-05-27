% getVox 
% written by Yazan Kadkoy 6/2017
% Based off of FindPercentFill by Sangeeta Subramanian and Peter Michael
% This pogram prompts user for file directory and then extracts microCT grey scale photos.Photos should be arranged
% in folders denoting the sample number, all these sample folders should then be placed into a folder.These 
% photos are all analyzed and the highest grey scale(GSV) value it determined for all photos. .25*greatest GSV
% was used as image threshold. a binary 3-D matric was then created for each sample denoting voxels above the 
% specfied threshold valu. 3-D matix is then saved in as filename_vox, where filename is the name of the sample 
% folder


%% Clear all data
clear all
clc
close all


    pixel_noise=25;
    U_lim=.4;
    bone.vox=[];
   
for s_f=1:1

    selec_folder=s_f;
    
    if selec_folder == 2
    clearvars -except selec_folder bonepix zray_length percent
    end
    
% Load in directory and folder name
folder = uigetdir('C:\Users\POCLAB\');
dirListing = dir(folder);

for d = 3:length(dirListing)
    name{d-2}=dirListing(d).name;
end

% name=sort_nat(name);
end
%% Generate a vector of file names

for iii=1:numel(name)
 

 
fN1 = [];
cname=char(name(iii));   
subNames= fullfile(folder,cname);
dirsub=dir(subNames);
fN1= strvcat(fN1,subNames);

[n l]=size(fN1);

  fN2=[];
  fname=[];
  for d = 3:length(dirsub)
    fname{d-2}=dirsub(d).name;
  end 
for d=1:length(fname)
Names= fullfile(fN1,char(fname(d)));
fN2= strvcat(fN2,Names);
end


[num len]=size(fN2);
 

% %% Load in images and place them into structure "ROTATEDBONE"
for i = 1:num
 b=imread(fN2(i,:));
rotatedbone.slice(i).pic=b;

end

 
[R C]=size(rotatedbone.slice(1).pic);
thick=zeros(length(rotatedbone.slice),C);
for z=1:length(rotatedbone.slice)
    im=rotatedbone.slice(z).pic;
    %creates a histogam of frequencies of all GSV 0-255
    hist=imhist(im);
    %finds the largest GSV 
    locs=find(hist>10);
    %calculates 25% of max value
    %%THESH(z)=U_lim*(locs(end)+1);
    THESH(z)=111;
end
for z=1:length(rotatedbone.slice)
    im=rotatedbone.slice(z).pic;
    %applies Gaussian filter sigma=1
    im2=imgaussfilt(im,1);
    
    bw=im2>max(THESH);
    
    
    bw=bwareaopen(bw,pixel_noise);
    h=subplot(1,2,2);
    imshow(bw);
    pause(0.01)
    hold on
    subplot(1,2,1)
    imshow(im)
    
     %creates bone structure with matric representing slice
    bone(z).vox=bw;
    clear im im2 bw
 
end
close all
[pathstr,nam,ext] = fileparts(fN1);
bone(1).name=nam;
savename=[nam '_vox.mat']
save(savename,'bone')

clearvars -except fN FN FN2 iii fN1 name folder U_lim pixel_noise

 

end
 
