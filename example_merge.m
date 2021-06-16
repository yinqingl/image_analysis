% automatic merge of channels of FISH images
%
% use DAPI as stitching reference
% based on SIFT descriptor for matching
% also keep track of the dapi intensity as a mask for later processing
%
% Yinqing Li
% yinqingl@csail.mit.edu
% 2016

%%
%initialization

%analyze fish image
addpath('imread')
addpath('helper_function')
addpath('IAT_v0.9.2')
addpath('TSP')
addpath('bfmatlab')
addpath('adaptivethreshold')

%{
ATTENTION: run the following code once
iat_setup('forever')
%}
iat_setup('forever')

%%
%list tif files
path = '..\Slide1\'; %system path to the images
path = strrep(path,'\',filesep);
datname = [path,'20161021_Ecel1T4_Spp1T1_Slide1_Section1_TRNL.nd2'];

%list all tif files, already stitched
datname_s = {};
% path_sets = {'..\Slide1\','..\Slide2\'}; %system path to the images
path_sets = {'../20171215/',...
    '../20171220/',...
    '../20180102/',...
    '../20180111/',...
    '../20180112/',...
    '../20180115/'};
for path_i = 1:length(path_sets),
    path = path_sets{path_i};
    path = strrep(path,'\',filesep);
    listing = dir(path);
    for i = 1:length(listing),
        filename = listing(i).name;
        ext = '.nd2.tif';
        if length(filename)>length(ext) && strcmp(filename(end-length(ext)+1:end),ext),
            datname_s{end+1} = [path,filename];
        end
    end
end

%print
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname))
end

%group tif files for channel merge
%custom criteria apply

datname_s_grp = {};

%datname_s_i = cellfun(@(x) regexp(x,'(?<=_)\d+(?=_TRN)', 'match'), datname_s);
%datname_s_i = cellfun(@(x) str2num(x), datname_s_i);
datname_s_i = {};
for i_img = 1:length(datname_s),
   s = strsplit(datname_s{i_img},'_');
   s = join(s(2:end),'_');
   datname_s_i{i_img} = s{1};
end
[~,datname_s_i] = sort(datname_s_i);
datname_s(datname_s_i);

datname_s_i_excl = [45];
%for i_img = 1:length(datname_s),
%    if ismember(i_img,datname_s_i_excl)
%        continue
%    end
%    i = datname_s_i(i_img);
%    j_img = find(datname_s_i == i);
%    j_img = setdiff(j_img,[1:25]);
%    j_img = setdiff(j_img,datname_s_i_excl);
%    if ~isempty(j_img)
%        datname_s_grp{end+1} = [i_img,j_img];
%    end
%end
i_img_list = setdiff([1:length(datname_s)],datname_s_i_excl);

i_img_list = setdiff(datname_s_i,datname_s_i_excl,'stable');

for i = [1:2:length(i_img_list)],
    i_img = i_img_list(i);
    j_img = i_img_list(i+1);
    if length(regexp(datname_s{i_img},'cy7', 'match'))>0,
      datname_s_grp{end+1} = [j_img,i_img];
    else,
      datname_s_grp{end+1} = [i_img,j_img];
    end
end

%print
for i_grp = 1:length(datname_s_grp),
     i_img_s = datname_s_grp{i_grp};
     display(sprintf('%d, %s, %s, %s',i_grp,num2str(i_img_s),datname_s{i_img_s(1)},datname_s{i_img_s(2)}))
end


%%
%configuration

%set reference channel
ref_channel = 1;    %the reference channel is set to 1 by default

notes = {};
%%
%process all image files
for i_grp = 1:length(datname_s_grp),
%for i_grp = [12],    
    i_img_s = datname_s_grp{i_grp};
    
    %first image, use dapi from the first image
    i_img = i_img_s(1);
    
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));
    
    data = imread(datname);
    data_mat = load([datname '.mat']);
    z_ref_mask_ref = data_mat.z_ref_mask_ref;
    img_ref = double(data);
    
    %second image
    i_img = i_img_s(2);
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname)); 
    
    data = imread(datname);
    img = double(data);
    
    %process reference channel
    img_ref_c = img_ref(:,:,ref_channel);
    img_c = img(:,:,ref_channel);
    
    height = size(img(:,:,ref_channel),1);
    width = size(img(:,:,ref_channel),2);
    
    %remove background
    %may need to normalize itensity between images
%     I = img_ref_c;
%     background = imopen(I,strel('disk',25));
%     I = I - background;
%     img_ref_c = I;
% 
%     I = img_c;
%     background = imopen(I,strel('disk',25));
%     I = I - background;
%     img_c = I;

    is_find_sift_affine_t = 0;
    if_use_thresholded_c = 1;
    if if_use_thresholded_c,
        ref_channels = [1];
        c = 1;
        img_ref_mask = ones(size(img_ref_c));
        
            mask_adth = adaptivethreshold(img(:,:,ref_channels(c)),[25 25],0,0);
            mask_adth_ref = adaptivethreshold(img_ref_mask.*img_ref(:,:,ref_channels(c)),[25 25],0,0);
            try,
                [t,tform] = get_sift_affine_t(mask_adth.*img(:,:,ref_channels(c)), mask_adth_ref.*img_ref_mask.*img_ref(:,:,ref_channels(c)));
                note = sprintf('find matching between series %d and %d using channel %d, use thresholded reference channel for registration', s, s_pred, ref_channels(c));
                notes{end+1} = {i_img, datname, note};
                display(note);                
                is_find_sift_affine_t = 1;
                break
            catch
                try
                    [t,tform] = get_sift_affine_t(img(:,:,ref_channels(c)), img_ref_mask.*img_ref(:,:,ref_channels(c)));
                    note = sprintf('find matching between series %d and %d using channel %d', s, s_pred, ref_channels(c));
                    notes{end+1} = {i_img, datname, note};
                    display(note);
                    is_find_sift_affine_t = 1;                    
                    break;
                catch
                    note = sprintf('try matching between series %d and %d using channel %d, failed', s, s_pred, ref_channels(c));
                    notes{end+1} = {i_img, datname, note};
                    display(note);
                    is_find_sift_affine_t = 0;
                end
            end
    end
    
    if ~is_find_sift_affine_t,
        note = sprintf('cannot find matching between series %d and %d', s, s_pred);
        notes{end+1} = {i_img, datname, note};
        display(note);
        continue
    end
    
    %register two images
    [t,tform] = get_sift_affine_t(img_c,img_ref_c);
%     figure; imagesc(img_c)
%     figure; imagesc(img_ref_c)
    
    %construct a canvas for the reference and the image
    [~,xdata,ydata] = imtransform(img(:,:,ref_channel), tform,...
        'UData',[1,width], 'VData',[1,height]);
    xdata = ceil(abs(xdata)).*sign(xdata);
    ydata = ceil(abs(ydata)).*sign(ydata);
    t_ref = [1 0 0; 0 1 0; -min(0,xdata(1)), -min(0,ydata(1)), 1];
    t.T(3,1) = t.T(3,1) + t_ref(3,1);
    t.T(3,2) = t.T(3,2) + t_ref(3,2);
    xdata_max = max(size(img_ref,2),xdata(2)) + t_ref(3,1);
    ydata_max = max(size(img_ref,1),ydata(2)) + t_ref(3,2);

    Rcb = imref2d([ydata_max,xdata_max]);
    t_ref = affine2d(t_ref);
    
    %apply transformation to all channels and construct composite image
    img_w = imwarp(img,t,'outputview',Rcb);
    img_ref_w = imwarp(img_ref,t_ref,'outputview',Rcb);
    % figure, imshowpair(img_w(:,:,1), img_ref_w(:,:,1))
    
    %update reference mask
    z_ref_mask_ref = imwarp(z_ref_mask_ref,t_ref,'outputview',Rcb);
    
    %merge channels
    img_w_clist = [1:size(img_w,3)];
    img_w_clist = setdiff(img_w_clist,[ref_channel]);
    img_w_clist = setdiff(img_w_clist,find(squeeze(max(max(img_w))) == 0));
    img_ref = cat(3,img_ref_w,img_w(:,:,img_w_clist)); %the first channel is set to the reference channel by default
    
    %write to tiff
    i_img = i_img_s(1);
    datname = datname_s{i_img};    
    datname_output = [datname '.merge.tif'];
    % data = uint16(round(img_c));
    data = uint16(round(img_ref));
    while size(data,3)<3,   %imagej cannot open tif with less than 3 channels
        data(:,:,end+1) = zeros(size(data(:,:,1)));
    end
    t = Tiff(datname_output, 'w');
    tagstruct.ImageLength = size(data, 1);
    tagstruct.ImageWidth = size(data, 2);
    tagstruct.Compression = Tiff.Compression.None;
    % tagstruct.Compression = Tiff.Compression.LZW;        % compressed
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
    % tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;                        % float data
    tagstruct.SamplesPerPixel = size(data,3);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct);
    t.write(data);
    t.close();

    %write reference channel mask
    datname_output = [datname_output '.mat'];
    save(datname_output,'ref_channel','z_ref_mask_ref','tform','-v7.3');
    
end

