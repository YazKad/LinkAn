function varargout = BoneRotate(varargin)
% BONEROTATE MATLAB code for BoneRotate.fig
%      BONEROTATE, by itself, creates a new BONEROTATE or raises the existing
%      singleton*.
%
%      H = BONEROTATE returns the handle to a new BONEROTATE or the handle to
%      the existing singleton*.
%
%      BONEROTATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BONEROTATE.M with the given input arguments.
%
%      BONEROTATE('Property','Value',...) creates a new BONEROTATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BoneRotate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BoneRotate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BoneRotate

% Last Modified by GUIDE v2.5 09-Sep-2014 15:39:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BoneRotate_OpeningFcn, ...
                   'gui_OutputFcn',  @BoneRotate_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BoneRotate is made visible.
function BoneRotate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BoneRotate (see VARARGIN)

% Choose default command line output for BoneRotate
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes BoneRotate wait for user response (see UIRESUME)
% uiwait(handles.figure1);

set(handles.imagefilt,'Visible','on')
set(handles.size,'string','0')
set(handles.thresh,'string','100')
set(handles.slider,'value',1)

evalin('base', 'exist frameskept;');
ans=evalin('base','ans');
if ans == 0
frameskept=[];
assignin('base','frameskept',frameskept)
end

evalin('base', 'exist FK;');
ans=evalin('base','ans');
if ans == 0
FK=[];
assignin('base','FK',FK)
end



% --- Outputs from this function are returned to the command line.
function varargout = BoneRotate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in loadfolder.
function loadfolder_Callback(hObject, eventdata, handles)
% hObject    handle to loadfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%profile on;

folder = uigetdir();
dirListing = dir(folder);

fN = [];
for d = 3:1:403
fileNames = fullfile(folder,dirListing(d).name);
fN = strvcat(fN,fileNames);
end

[num len]=size(fN);

assignin('base','fN',fN)

for i = 1:num
assignin('base','i',i)
set(handles.size,'string',num2str(i))
fNs = fN(i,end-16:end);
fNs = cellstr(fNs);
bone.name(i,1)=fNs;
evalin('base', 'b=imread(fN(i,:));');
b=evalin('base','b');
% bb=b(500:1500,500:1500);
bone.slice(i,1).pic=b;
end

set(handles.slider,'Min',1)
set(handles.slider,'Max',length(bone.slice))
set(handles.slider,'SliderStep',[1/(length(bone.slice)) 1/(length(bone.slice))])
f=get(handles.slider,'value');
f=round(f);
set(handles.frame,'string',num2str(f))
set(handles.slider,'value',f)

assignin('base','bone',bone)

%profile off;
%profile report



% --- Executes on button press in previewcrop.
function previewcrop_Callback(hObject, eventdata, handles)
% hObject    handle to previewcrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%profile on;

bone=evalin('base','bone');

Xmin=get(handles.Xmin,'string');
Xmin=str2num(Xmin);

Xmax=get(handles.Xmax,'string');
Xmax=str2num(Xmax);

Ymin=get(handles.Ymin,'string');
Ymin=str2num(Ymin);

Ymax=get(handles.Ymax,'string');
Ymax=str2num(Ymax);

axes(handles.imagefilt);

[a b]=size(bone.slice(1).pic);

for i=1:20:length(bone.slice)
    imshow(bone.slice(i).pic,'XData',[0 b],'YData',[0 a]);
    hold on
    line([Xmin Xmin],[Ymin Ymax],'LineWidth', 2)
    line([Xmax Xmax],[Ymin Ymax],'LineWidth', 2)
    line([Xmin Xmax],[Ymin Ymin],'LineWidth', 2)
    line([Xmin Xmax],[Ymax Ymax],'LineWidth', 2)
    hold off
    pause(0.1)
end

%profile off;
%profile report



% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bone=evalin('base','bone');

Xmin=get(handles.Xmin,'string');
Xmin=str2num(Xmin);

Xmax=get(handles.Xmax,'string');
Xmax=str2num(Xmax);

Ymin=get(handles.Ymin,'string');
Ymin=str2num(Ymin);

Ymax=get(handles.Ymax,'string');
Ymax=str2num(Ymax);

for i=1:length(bone.slice)
    b=bone.slice(i).pic;
    bb=b(Ymin:Ymax,Xmin:Xmax);
    bone.slice(i).pic=[];
    bone.slice(i).pic=bb;
end

assignin('base','bone',bone)



% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

bone=evalin('base','bone');

set(handles.slider,'Min',1)
set(handles.slider,'Max',length(bone.slice))
set(handles.slider,'SliderStep',[1/(length(bone.slice)) 1/(length(bone.slice))])
f=get(handles.slider,'value');
f=round(f);
set(handles.frame,'string',num2str(f))
set(handles.slider,'value',f)

thresh=get(handles.thresh,'string');
thresh=str2num(thresh);

axes(handles.imagefilt);

imgfilt=bone.slice(f).pic;
imgfilt(imgfilt<thresh)=0;
imshow(imgfilt)
assignin('base','frame',f)



% --- Executes on button press in keep.
function keep_Callback(hObject, eventdata, handles)
% hObject    handle to keep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frame=evalin('base','frame');
frameskept=evalin('base','frameskept');
FK=evalin('base','FK');

frameskept=[frameskept ' ' num2str(frame)];
FK=[FK frame];

assignin('base','frameskept',frameskept)
assignin('base','FK',FK)
set(handles.keptframes,'string',frameskept)



% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.keptframes,'string','')
frameskept=[];
assignin('base','frameskept',frameskept)
FK=[];
assignin('base','FK',FK)




% --- Executes on button press in edgedetect.
function edgedetect_Callback(hObject, eventdata, handles)
% hObject    handle to edgedetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bone=evalin('base','bone');
FK=evalin('base','FK');

for i=1:length(FK)
    
    set(handles.edgeframe,'string',num2str(FK(i)))
    ss=bone.slice(FK(i)).pic;
    thresh=str2num(get(handles.thresh,'string'));
    ss(ss<thresh)=0;
    
[B,L,N,A] = bwboundaries(ss,8);

 axes(handles.bounds);
 imshow(ss)
 hold on

for ii = 1:length(B)
    bound=B{ii};
    plot(bound(:,2), bound(:,1), 'w', 'LineWidth', 2)
end

for ii = 1:length(B)
    sB(ii,1:2)=size(B{ii});
end

[bsize bindex]=max(sB(:,1));
enclosed_boundaries = find(A(:,bindex));

for ii = 1:length(enclosed_boundaries)
    enc=B{enclosed_boundaries(ii)};
    size_enc(ii,1:2)=size(enc);
end


enc_max = max(size_enc(:,1));
enc_index = find(sB(:,1) == enc_max);
mmm=B{enc_index};
marrow{i}=mmm;

axes(handles.center);
imshow(ss)
hold on
m=marrow{i};
plot(m(:,2), m(:,1), 'r', 'LineWidth', 2)
means(i,:)=[mean(m(:,2)) mean(m(:,1)) FK(i)];
plot(means(i,1),means(i,2), 'g*')

end

assignin('base','centroids',means)



% --- Executes on button press in linfit.
function linfit_Callback(hObject, eventdata, handles)
% hObject    handle to linfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

centroids=evalin('base','centroids');
bone=evalin('base','bone');

x1=centroids(:,1);
x2=centroids(:,2);
x3=centroids(:,3);

mx1=round(mean(x1));
mx2=round(mean(x2));

[m n]=size(bone.slice);
for i=1:m
front(i,:)=(bone.slice(i).pic(mx2,:));
side(i,:)=(bone.slice(i).pic(:,mx1));
end

axes(handles.linearxz);
fitobjectxz = fit(x3,x1,'poly1');
plot(fitobjectxz,x3,x1,'.')
hold on
imshow(front','InitialMag', 'fit')
alpha(0.4)
hold off
legend('off')
xlabel('')
ylabel('')

axes(handles.linearyz)
fitobjectyz = fit(x3,x2,'poly1');
plot(fitobjectyz,x3,x2,'.')
hold on
imshow(side','InitialMag', 'fit')
alpha(0.4)
hold off
legend('off')
xlabel('')
ylabel('')


assignin('base','fitobjectxz',fitobjectxz)
assignin('base','fitobjectyz',fitobjectyz)



% --- Executes on button press in resampxz.
function resampxz_Callback(hObject, eventdata, handles)
% hObject    handle to resampxz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bone=evalin('base','bone');
fitobjectxz=evalin('base','fitobjectxz');

[s p]=size(bone.slice(1).pic);
[m n]=size(bone.slice);

for j=1:s
for i=1:m
front(i,:)=(bone.slice(i).pic(j,:));
end

% axes(handles.newxz)
% imshow(front')
% pause(1)
ang1=coeffvalues(fitobjectxz);
ang1=atan(ang1(1))*180/pi;
zzz=imrotate(front,-ang1);
axes(handles.newxz)
imshow(zzz')
% pause(1)
rotated.slice(j,1).pic=zzz;
end

assignin('base','rotated',rotated)



% --- Executes on button press in rsampyz.
function rsampyz_Callback(hObject, eventdata, handles)
% hObject    handle to rsampyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rotated=evalin('base','rotated');
fitobjectyz=evalin('base','fitobjectyz');

[s p]=size(rotated.slice(1).pic);
[m n]=size(rotated.slice);

for j=1:p

for i=1:m
side(:,i)=(rotated.slice(i).pic(:,j));
end

% axes(handles.newyz)
% imshow(side')
% pause(0.5)
ang1=coeffvalues(fitobjectyz);
ang1=atan(ang1(1))*180/pi;
zzz=imrotate(side,-ang1);
axes(handles.newyz)
imshow(zzz')
% pause(0.5)
rotated2.slice(j).pic=zzz;
end

assignin('base','rotated2',rotated2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function size_Callback(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of size as text
%        str2double(get(hObject,'String')) returns contents of size as a double


% --- Executes during object creation, after setting all properties.
function size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_Callback(hObject, eventdata, handles)
% hObject    handle to frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame as text
%        str2double(get(hObject,'String')) returns contents of frame as a double


% --- Executes during object creation, after setting all properties.
function frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function keptframes_Callback(hObject, eventdata, handles)
% hObject    handle to keptframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of keptframes as text
%        str2double(get(hObject,'String')) returns contents of keptframes as a double


% --- Executes during object creation, after setting all properties.
function keptframes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to keptframes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edgeframe_Callback(hObject, eventdata, handles)
% hObject    handle to edgeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edgeframe as text
%        str2double(get(hObject,'String')) returns contents of edgeframe as a double


% --- Executes during object creation, after setting all properties.
function edgeframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgeframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterate_Callback(hObject, eventdata, handles)
% hObject    handle to iterate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iterate as text
%        str2double(get(hObject,'String')) returns contents of iterate as a double


% --- Executes during object creation, after setting all properties.
function iterate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iterate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Xmin_Callback(hObject, eventdata, handles)
% hObject    handle to Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xmin as text
%        str2double(get(hObject,'String')) returns contents of Xmin as a double


% --- Executes during object creation, after setting all properties.
function Xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xmax_Callback(hObject, eventdata, handles)
% hObject    handle to Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xmax as text
%        str2double(get(hObject,'String')) returns contents of Xmax as a double


% --- Executes during object creation, after setting all properties.
function Xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ymin_Callback(hObject, eventdata, handles)
% hObject    handle to Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ymin as text
%        str2double(get(hObject,'String')) returns contents of Ymin as a double


% --- Executes during object creation, after setting all properties.
function Ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ymax_Callback(hObject, eventdata, handles)
% hObject    handle to Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ymax as text
%        str2double(get(hObject,'String')) returns contents of Ymax as a double


% --- Executes during object creation, after setting all properties.
function Ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
