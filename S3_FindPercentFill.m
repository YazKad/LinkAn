%% Clear all data
clear all
clc
close all

%% Calculate full and partial data

    bonepix.full=[];
    zray_length.full=[];
    percent.full=[];
    bonepix.defect=[];
    zray_length.defect=[];
    percent.defect=[];


for s_f=1:2

    selec_folder=s_f;
    
    if selec_folder == 2
    clearvars -except selec_folder bonepix zray_length percent
    end
    
% Load in directory and folder name
folder = uigetdir('C:\Users\Peter Michael\Desktop\Bone Rotation');
dirListing = dir(folder);

for d = 3:length(dirListing)
    name{d-2}=dirListing(d).name;
end

name=sort_nat(name);
 
%% Generate a vector of file names
fN = [];
for d = 1:length(name)
fileNames = fullfile(folder,char(name(d)));
fN = strvcat(fN,fileNames);
end
 
[num len]=size(fN);
 
 
%% Load in images and place them into structure "ROTATEDBONE"
for i = 1:num
b=imread(fN(i,:));
rotatedbone.slice(i).pic=b;
end
 

%% Determine Threshold

imshow(rotatedbone.slice(round(num/3)).pic)  % display an image 1/3 of the way down
pause                 % hold code for user keystroke

L_lim=input('Define noise threshold (all below will be eliminated):  ');
U_lim=input('Define bone threshold (all above will be kept):  ');

 
%% Calculate percent fill
 
[R C]=size(rotatedbone.slice(1).pic);
 
for x=1:R
    for y=1:C
        for z=1:length(rotatedbone.slice)
            zray(z,:)=rotatedbone.slice(z).pic(x,y);
        end
        
        pix_bone=find(zray > U_lim); %find bone pixels
        pix_tissue=find(zray > L_lim); %find tissue pixels
        if length(pix_tissue)~=0
            zlen=(pix_tissue(end)- pix_tissue(1))+1; %first and last positive pixel (adjusted Z-ray length)
        elseif length(pix_tissue)==0
            zlen=1;  % to avoid dividing by zero
        end

        bone(x,y)=length(pix_bone);
        total_len(x,y)=zlen; %length of tissue in zray (ignoring )
        perc(x,y)=length(pix_bone)/zlen; %percent fill
        
    end
end
 
 
%% Color Map
figure
colormap('hot')
imagesc(perc*100)
colorbar;


%% Save into full or partial data
if selec_folder == 1
    bonepix(:).full=bone;
    zray_length(:).full=total_len;
    percent(:).full=perc;
end

if selec_folder == 2
    bonepix(:).defect=bone;
    zray_length(:).defect=total_len;
    percent(:).defect=perc;
end
end

uisave({'bonepix','zray_length','percent'},'percent.mat')

% profile off
% profile report
 
