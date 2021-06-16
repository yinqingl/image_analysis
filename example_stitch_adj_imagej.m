% adjust stitched tif image using imagej
% automatic adjust contrast
%
% Yinqing Li
% yinqingl@csail.mit.edu
% 2016

%%
%load functions
addpath('imread')
addpath('helper_function')
addpath('IAT_v0.9.2')
addpath('TSP')
addpath('bfmatlab')
addpath('bradley')
addpath('adaptivethreshold')

addpath('C:\\fiji-win64\\Fiji.app\\scripts')
ImageJ

%%
%list all tif files
datname_s = {};
%system path to the images
% path_sets = {'..\Slide1\','..\Slide2\'};
% path_sets = {'..\20170130\stitch\','..\20170207\stitch\','..\20170209\stitch\'};
path_sets = {'..\20170227\'};
for path_i = 1:length(path_sets),
    path = path_sets{path_i};
    path = strrep(path,'\',filesep);
    listing = dir(path);
    for i = 1:length(listing),
        filename = listing(i).name;
        ext = '.nd2.tif';
%         ext = '.nd2.tif.merge.tif';
        if length(filename)>length(ext) && strcmp(filename(end-length(ext)+1:end),ext),
            datname_s{end+1} = [path,filename];
        end
    end
end

i = cellfun(@(x) strcmp(x(end-7:end),'.adj.tif'), datname_s);
datname_s = datname_s(~i);

% i = cellfun(@(x) ~isempty(findstr(x,'PVT6')), datname_s);
% datname_s = datname_s(~i);

%print
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname))
end

%%
%imagej api configuration

%imagej process
imagej_macro_list = {
    'open("PATH2FILE");'
    'run("Channels Tool...");'
    'Stack.setDisplayMode("color");'
    'Stack.setChannel(CHANNEL);',
    'run("COLOR");' %Blue, Green, Red, Magenta
    'run("Brightness/Contrast...");'
    'run("Enhance Contrast", "saturated=0.35");'
    'run("Apply LUT");'
    'setMinAndMax(MIN_VAL, MAX_VAL);'
    'saveAs("Tiff", "PATH2FILE");'
    'run("Close All");'
    'run("Close");'
};
for i = 1:size(imagej_macro_list,1)
    disp(sprintf('%d, %s',i,imagej_macro_list{i}))
end

% imagej_ij.runMacro(str_macro);

%%
%iterate through all images
%write to composite images
imagej_w = ij.ImageJ;    %image window
imagej_ij = ij.IJ;    %imagej macro

color_str = {'Blue','Green','Red','Magenta','Grays'};

for i_img = [1],
% for i_img = 1:length(datname_s),
    datname = datname_s{i_img};
    datname = GetFullPath(datname);    
    display(sprintf('%d, %s',i_img, datname));
    
data_info = imfinfo(datname);
nchannels = length(data_info.BitsPerSample);

str_macro = imagej_macro_list{1};
str_macro = strrep(strrep(str_macro,'PATH2FILE',datname),'\','\\');
imagej_ij.runMacro(str_macro);

for c = 1:nchannels,
    str_macro = imagej_macro_list{2}; imagej_ij.runMacro(str_macro);
    str_macro = imagej_macro_list{3}; imagej_ij.runMacro(str_macro);
    
    str_macro = imagej_macro_list{4}; 
    str_macro = strrep(str_macro,'CHANNEL',num2str(c));
    imagej_ij.runMacro(str_macro);
    
	str_macro = imagej_macro_list{5};
    if c > length(color_str),
        c = length(color_str);
    end
    str_macro = strrep(str_macro,'COLOR',color_str{c});
    imagej_ij.runMacro(str_macro);   
    
    str_macro = imagej_macro_list{6}; imagej_ij.runMacro(str_macro);
    str_macro = imagej_macro_list{7}; imagej_ij.runMacro(str_macro);
%     str_macro = imagej_macro_list{8}; imagej_ij.runMacro(str_macro);
end


%write to tif
datname_output = [datname '.adj.tif'];
str_macro = imagej_macro_list{10};
str_macro = strrep(strrep(str_macro,'PATH2FILE',datname_output),'\','\\');

imagej_ij.runMacro(str_macro);

imagej_ij.runMacro('run("Close All");');
end

%quit imagej
str_macro = imagej_macro_list{2}; imagej_ij.runMacro(str_macro); pause(1)
str_macro = imagej_macro_list{12}; imagej_ij.runMacro(str_macro); pause(1)

str_macro = imagej_macro_list{6}; imagej_ij.runMacro(str_macro); pause(1)
str_macro = imagej_macro_list{12}; imagej_ij.runMacro(str_macro); pause(1)

imagej_ij.runMacro('run("Close All");');
imagej_w.quit();

%%
%compile thresholds
I_th_c_s = {};

% for i_img = [1:25],
for i_img = 1:length(datname_s),
    datname = datname_s{i_img};
    datname = [datname '.adj.tif'];
    datname = [datname '.mat'];
    display(sprintf('%d, %s',i_img, datname));
    try
        load(datname)
    catch
        display(sprintf('skip image: %s',datname))
        continue
    end
    I_th_c_s{i_img} = I_th_c;
end

I_th_c = cellfun(@(x) x(1,2), I_th_c_s);
figure;
plot(I_th_c(1,:))
hold on
plot(I_th_c(2,:),'r')

I_th_prev = I_th_c_s{20};

%%
%single channel images, using threshold determined for composite
% for i_img = [26:39],
for i_img = 1:length(datname_s),
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));

data = imread(datname);
data_mat = load([datname '.mat']);

% if sum(max(max(data))>0)<3
%     display(sprintf('skip image: %s',datname))
%     continue
% end

I_int_max = double(intmax('uint16'));

data_adj = data;
for c = 1:size(data,3),
    I = double(data(:,:,c));

    if is_bgsub_used
        background = imopen(I,strel('disk',25));
        I = I - background;
    end    
    
    if max(I(:))==0,
        continue
    end
    
    %use previous determined I_th
    I_th_c = I_th_c_s{i_img};
    if 1,
%         I_th_l = I_th_prev(1,c);
        I_th_l = 0;
        I_th_h = I_th_c(2,c);
    end
    
    %adjust contrast
    I_adj = I;
    I_adj(I_adj > I_th_h) = I_th_h;
    I_adj = I_adj - I_th_l;
    I_adj(I_adj<0) = 0;
    I_adj = I_adj/max(I_adj(:))*I_int_max;
    
    data_adj(:,:,c) = I_adj;
% end
% 
% %write to tif
% for c = 1:size(data_adj,3),
    datname_output = sprintf('%s.adj.c%d.tif',datname,c);
    data_w = data_adj;
    data_w = data_w(:,:,[c,c,c]);
    data_w = uint16(round(data_w));
    imwrite(data_w,datname_output);
end

end

%%
%sort files
datname_s = {};
%system path to the images
% path_sets = {'..\Slide1\','..\Slide2\'};
path_sets = {'..\20161118\sorted\'};
for path_i = 1:length(path_sets),
    path = path_sets{path_i};
    path = strrep(path,'\',filesep);
    listing = dir(path);
    for i = 1:length(listing),
        filename = listing(i).name;
        ext = '.nd2.tif.adj.tif';
%         ext = '.nd2.tif.merge.tif';
        if length(filename)>length(ext) && strcmp(filename(end-length(ext)+1:end),ext),
            datname_s{end+1} = [path,filename];
        end
    end
end

% i = cellfun(@(x) strcmp(x(end-7:end),'.adj.tif'), datname_s);
% datname_s = datname_s(~i);

%print
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname))
end

datname_s_idx = {};
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    i_slide = regexp(datname,'(?<=Slide)[0-9]+','match'); i_slide = str2num(i_slide{1});
    i_section = regexp(datname,'(?<=Section)[0-9]+','match'); i_section = str2num(i_section{1});
    datname_s_idx{i_slide,i_section} = datname;
end

s_idx = importdata('..\20161118\sorted\section_series.txt');
[~,s_tour] = sort(s_idx.data(:,3));

i = 1;
for i_img = hv(s_tour),
    datname = datname_s_idx{s_idx.data(i_img,1), s_idx.data(i_img,2)};
    f_name = strsplit(datname,filesep);
    f_name = f_name{end};
    display(sprintf('copy %s .%ssorted%s%03d_%s',f_name,filesep,filesep,i, f_name));
    i = i+1;
end