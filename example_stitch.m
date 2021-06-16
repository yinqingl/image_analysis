% automatic stitching of FISH images
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
%list nd2 files
path = '..\20171215\'; %system path to the images
path = strrep(path,'\',filesep);
datname = [path,'20161021_Ecel1T4_Spp1T1_Slide1_Section1_TRNL.nd2'];

%list all nd2 files
datname_s = {};
% path_sets = {'..\Slide1\','..\Slide2\'}; %system path to the images
% path_sets = {'..\20170130\','..\20170207\','..\20170209\'}; %system path to the images
% path_sets = {'..\20170227\'}; %system path to the images
% path_sets = {'../20170731/',...
%     '../20170809/',...
%     '../20170815/'};
path_sets = {'../20171215/',...
    '../20171220/',...
    '../20180102/',...
    '../20180111/',...
    '../20180112/',...
    '../20180115/'};
for path_i = 1:length(path_sets),
    path = path_sets{path_i};
    path = strrep(path,'\','*');
    path = strrep(path,'/','*');
    path = strrep(path,'*',filesep);
    listing = dir(path);
    for i = 1:length(listing),
        filename = listing(i).name;
        ext = '.nd2';
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

%38, 42, 48, 52
%%
%configuration

%set reference channel
ref_channel = 1;

%filter setup
h_stdfilt = ones(5); %std filter size

is_z_halfdepth_filter_used = 1; %use dapi to filter z stacks, set to 0 for weak dapi
z_halfdepth_around_ref_pl = 3; %the half depth of z stacks containing signal
med2filter_n = 5;   %number of x grids for z stack filter
med2filter_m = 5;   %number of y grids for z stack filter

%background adjustment
is_d_bk_sel = {'all'};
% is_d_bk_sel = {'on',[1]}; %only adjust background for the selected channels
% is_d_bk_sel = {'off',[1]}; %do not adjust background for the selected channels


%%
%process all image files
ch_name_seq = {'DAPI','FITC','TRITC','CY3','CY5','CY7'};

notes = {};
%for i_img = 39:length(datname_s),
for i_img = [38],
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));

% datname = [path,'2016-10-3_Slide1_Section2_TRN_3Color.nd2'];
% meta = imreadBFmeta_nd2(datname);
try
    [T,meta] = evalc('imreadBFmeta_nd2(datname);');
catch
    note = sprintf('cannot find open file');
    display(note);
    notes{end+1} = {i_img, datname, note};
    continue
end

meta

width = meta.width;
height = meta.height;
zsize = meta.zsize;
nframes = meta.nframes;
nseries = meta.nseries;
nchannels = meta.channels;
omemeta = meta.result{1,3};

%arrange channels by exitation wave length
ch_name = {};
for c = 1:nchannels,
    ch_name{c} = char(omemeta.getChannelName(0,c-1));
end
ch_arrange = [];

for c = 1:nchannels,
    for c2 = 1:length(ch_name_seq)
        if ~isempty(findstr(upper(ch_name{c}),ch_name_seq{c2}))
            ch_arrange(c) = c2;
            continue
        end
    end
end

for c = 1:nchannels,
    ch_name{c} = ch_name_seq{ch_arrange(c)};
end

[~,ch_arrange] = sort(ch_arrange);
ch_name = ch_name(ch_arrange);

disp(sprintf('arrange channels'))
ch_name
ch_arrange

if ~is_z_halfdepth_filter_used,
    z_halfdepth_around_ref_pl = zsize; %the half depth of z stacks containing signal
    med2filter_n = 1;   %number of x grids for z stack filter
    med2filter_m = 1;   %number of y grids for z stack filter    
end

if nseries <= 1,
    disp('nseries = 1')
    continue
end

%determine the order of series for alignment based on the x,y coordinates
pos_xy = double(cell2mat(meta.result(:,2)));
pos_xy_dist = dist(pos_xy', 'euclidean');
UG = sparse(pos_xy_dist);
[~,s_tour_pred] = graphminspantree(UG);
%BFS on series tour
s_tour = [];
queue = [find(s_tour_pred==0)];
while(~isempty(queue)),
    s_tour = [s_tour, queue(1)];
    queue(1) = [];
    queue = [queue, find(s_tour_pred == s_tour(end))];
end

%{
figure;
scatter3(pos_xy(:,1),pos_xy(:,2),[1:size(pos_xy,1)],21,[1:size(pos_xy,1)],'filled');view(2); grid off
hold on
for s = s_tour,
    s_pred = s_tour_pred(s);
    if s_pred == 0,
        continue
    end
    plot([pos_xy(s_pred,1),pos_xy(s,1)],[pos_xy(s_pred,2),pos_xy(s,2)])
end
daspect([1 1 1])
view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'))
view(biograph(s_tour,[],'ShowArrows','off','ShowWeights','on'))
s_tour
s_tour_pred(s_tour)
%}

%for each position
% s = 1;
% for s = [1:2],

img_mask_s = {};
    
for s = s_tour;
% for s = [9    14    16    15    17    18    19    20];
    disp(sprintf('loading frame: %d',s));
    
    %for each channel
    img = zeros(height, width, nchannels);
    for c = 1:nchannels,
        %read in all z stacks
%         [vol]=imreadBF_nd2(datname,[1:zsize],t,c);
        [T,vol] = evalc('imreadBF_nd2(datname,[1:zsize],s,ch_arrange(c));');
%         if meta.nd2,
%             [T,vol] = evalc('imreadBF_nd2(datname,[1:zsize],t,c);');
%         else
%             [T,vol] = evalc('imreadBF(datname,[1:zsize],t,c);');
%         end

        % figure;
        % imagesc(vol);
        % colormap(gray);

        %max projection
%         I_zmax = max(vol,[],3);
        %max projection
        if c == ref_channel,
            if zsize>1,
                [I_zmax,idx_zmax] = max(vol,[],3);
            else,
                I_zmax = vol;
                idx_zmax = ones(size(vol));
            end
            
            if min(min(I_zmax)) == max(max(I_zmax)),
                z_ref_mask = ones(size(I_zmax));
            else
                I = adaptivethreshold(I_zmax,[250 250],0,0);
                stdfilt_I = stdfilt(I,ones(3));
                z_ref_mask = I.*(stdfilt_I == 0);
                z_ref_mask = imdilate(z_ref_mask,ones(5));
            end
                
            if min(idx_zmax(:)) == max(idx_zmax(:)),
                z_ref_pl = ones(med2filter_n,med2filter_m);
                z_ref_pl(:,:) = min(idx_zmax(:));
            else,
%                 I = adaptivethreshold(I_zmax,[250 250],0,0);
%                 stdfilt_I = stdfilt(I,ones(3));
%                 z_ref_mask = I.*(stdfilt_I == 0);
%                 z_ref_mask = imdilate(z_ref_mask,ones(5));
%                 figure;imagesc(I_zmax)
%                 figure, imshowpair(z_ref_mask, I)
%                 stdfilt_J = stdfilt(idx_zmax,h_stdfilt); %find nuclei by looking at projected z, because the max z at no nuclei position should be random                
%                 [kmeans_idx_i,kmeans_c] = kmeans(stdfilt_J(i),2);
%                 kmeans_idx = zeros(size(i))-1;
%                 kmeans_idx(i) = kmeans_idx_i;
%                 kmeans_idx = reshape(kmeans_idx,height,width);
%         %         figure;imagesc(I_zmax);
%                 [~,kmeans_idx_min] = min(kmeans_c);
%                 z_ref_mask = (kmeans_idx == kmeans_idx_min);
                z_ref_pl = zeros(med2filter_n,med2filter_m);
                n_range = round(linspace(1,height,med2filter_n+1));
                m_range = round(linspace(1,width,med2filter_m+1));
                idx_zmax_msk = idx_zmax.*z_ref_mask;
                for n = 1:med2filter_n,
                    for m = 1:med2filter_m,
                        pl = idx_zmax_msk(n_range(n):n_range(n+1),m_range(m):m_range(m+1));
                        z_ref_pl(n,m) = median(pl(pl>0));
                    end
                end
                %clear  kmeans_idx idx_zmax_msk;
            end
        else
            %for other channels, take max projection within a range of z
            %the range is determined based on the reference channel
            if zsize>1,
                I_zmax = zeros(height,width);
                n_range = round(linspace(1,height,med2filter_n+1));
                m_range = round(linspace(1,width,med2filter_m+1));
                for n = 1:med2filter_n,
                    for m = 1:med2filter_m,
                        z_range = round([max(z_ref_pl(n,m)-z_halfdepth_around_ref_pl,1):1:...
                            min(z_ref_pl(n,m)+z_halfdepth_around_ref_pl,zsize)]);
                        vol_nm = vol(n_range(n):n_range(n+1),m_range(m):m_range(m+1),z_range);
                        [I_zmax_nm,idx_zmax] = max(vol_nm,[],3);
                        I_zmax(n_range(n):n_range(n+1),m_range(m):m_range(m+1)) = I_zmax_nm;
                    end
                end
            else
                I_zmax = vol;
            end
        end

        %substract background
        I = squeeze(I_zmax);
        background = imopen(I,strel('disk',25));
        I = I - background;
        I_bk_adj = I;
        img(:,:,c) = I_bk_adj;
    end
    %clear vol;
    %clear I_zmax I_bk_adj;
    
    

    %for the first position in the series, skip image alignment
    if s==s_tour(1),
        img_ref = img;
        img_ref_mask = ones(size(img_ref(:,:,ref_channel)));
        z_ref_mask_ref = z_ref_mask;
        
        img_mask_s{s} = img_ref_mask;
        continue
    end
    
    s_pred = s_tour_pred(s);
    img_ref_mask = img_mask_s{s_pred};
    
    % figure; imagesc(img(:,:,1)); daspect([1 1 1]);
    % figure; imagesc(img_ref_mask.*img_ref(:,:,1)); daspect([1 1 1]);

    
    %register two images
%     try
%         [t,tform] = get_sift_affine_t(img(:,:,ref_channel),img_ref_mask.*img_ref(:,:,ref_channel));
%         if_use_thresholded_c = 0;
%     catch
%         if_use_thresholded_c = 1;
%     end
    
    img_std_c = [];
    for c = 1:nchannels,
        img_std_c(c) = median([std(img(:,:,c),[],1), std(img(:,:,c),[],2)']);
    end
    ref_channels = [ref_channel];
    [~,img_std_c_i] = sort(img_std_c,'descend');
    for c = img_std_c_i,
        if c == ref_channel,
            continue
        end
        ref_channels(end+1) = c;
    end
    is_find_sift_affine_t = 0;
    if_use_thresholded_c = 1;
    if if_use_thresholded_c,
        for c = 1:length(ref_channels),
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
    end
    
    if ~is_find_sift_affine_t,
        note = sprintf('cannot find matching between series %d and %d', s, s_pred);
        notes{end+1} = {i_img, datname, note};
        display(note);
        %update s_pred
        s_tour_pred(s_tour_pred == s) = s_pred;
        disp('s_tour: ')
        s_tour
        disp('updated s_tour_pred: ')
        s_tour_pred(s_tour)        
        continue
    end
    
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
    
    %update series mask
    for ss = find(cellfun(@(x) ~isempty(x), img_mask_s)),
        img_mask_s{ss} = imwarp(img_mask_s{ss},t_ref,'outputview',Rcb);
    end
    
    for c = 1:length(ref_channels),
        if min(min(img_w(:,:,c))) ~= max(max(img_w(:,:,c))),
            img_mask_s{s} = img_w(:,:,c) ~= 0;
        end
    end
    
    % figure, imshowpair(img_w(:,:,1), img_ref_w(:,:,1))
    
    %adjust background based on the difference in the overlapping area
    %the background area is found by excluding the dapi high area
    img_mask_w = imwarp(~z_ref_mask,t,'outputview',Rcb);
    img_mask_ref_w = imwarp(~z_ref_mask_ref,t_ref,'outputview',Rcb);
    img_overlap_mask = (img_mask_w + img_mask_ref_w)==2;
%     figure; imagesc(img_w(:,:,2).*img_overlap_mask)
%     figure; imagesc(img_ref_w(:,:,2).*img_overlap_mask)
    
    for c = 1:nchannels,
        
        img_w_masked = img_w(:,:,c).*img_overlap_mask;
        img_ref_w_masked = img_ref_w(:,:,c).*img_overlap_mask;
        
        if min(min(img_w_masked)) == max(max(img_w_masked)),
            continue
        end
        if min(min(img_ref_w_masked)) == max(max(img_ref_w_masked)),
            continue
        end 
        
%         med_bk = median(img_w_masked(img_w_masked>0));
%         med_bk_ref = median(img_ref_w_masked(img_ref_w_masked>0));        
        
        img_w_bk = img_w_masked(img_w_masked>0);
        img_ref_w_bk = img_ref_w_masked(img_ref_w_masked>0);
        
        %norm using kde on linear scale hist local maximum
        I_bk = img_w_bk;
        I_th_5 = quantile(I_bk,0.95);
        [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_5)),256); n_bk = n/sum(n);
        [~,x_bk_i] = max(n_bk);
        med_bk = (x_bk(x_bk_i));
        
        I_bk = img_ref_w_bk;
        I_th_5 = quantile(I_bk,0.95);
        [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_5)),256); n_bk = n/sum(n);
        [~,x_bk_i] = max(n_bk);
        med_bk_ref = (x_bk(x_bk_i));
        
%         figure; plot(x_bk, n_bk); hold on; plot(x_bk,n_bk,'r'); 
        
%         %norm on log scale hist local maximum
%         [n,x] = hist(log(img_w_bk),256);
%         n_bk = histc(log(img_w_bk),x); n_bk = n_bk/sum(n_bk);
%         n_bk_ref = histc(log(img_ref_w_bk),x); n_bk_ref = n_bk_ref/sum(n_bk_ref);
% %         figure; plot(x,n_bk); hold on; plot(x,n_bk_ref,'r');
%         [~,x_bk_i] = max(n_bk);
%         med_bk = exp(x(x_bk_i));
%         [~,x_bk_i] = max(n_bk_ref);
%         med_bk_ref = exp(x(x_bk_i));
        
        d_bk = med_bk - med_bk_ref;
        
        if ~strcmp(is_d_bk_sel{1},'all'),
            if strcmp(is_d_bk_sel{1},'on'),
                if ~ismember(c,is_d_bk_sel{2})
                    d_bk = 0;
                end
            elseif strcmp(is_d_bk_sel{1},'off'),
                if ismember(c,is_d_bk_sel{2})
                    d_bk = 0;
                end
            end
        end
        
%         d_bk = (img_w_masked./img_ref_w_masked);
%         d_bk_mask = img_w_masked==0 | img_ref_w_masked==0;
%         d_bk = d_bk(~d_bk_mask);
% %         figure; hist(log(d_bk),256);
%         [n,x] = hist(log(d_bk),256);
%         [~,i_x] = max(n);
%         d_bk = x(i_x);
        
        %reduce background
        if d_bk > 0,
            img_d_bk = zeros(size(img_w_masked));
            img_d_bk = d_bk;
%             img_d_bk(img_mask_w) = d_bk;
            img_w(:,:,c) = max(img_w(:,:,c) - img_d_bk,0);
%             img_w(:,:,c) = max(img_w(:,:,c)/exp(img_d_bk),1);
%             img_ref_w(:,:,c) = max(img_ref_w(:,:,c)*exp(img_d_bk),1);
        else
            img_d_bk = zeros(size(img_ref_w_masked));
            img_d_bk = -d_bk;
%             img_d_bk(img_mask_ref_w) = -d_bk;            
            img_ref_w(:,:,c) = max(img_ref_w(:,:,c) - img_d_bk,0);
%             img_ref_w(:,:,c) = max(img_ref_w(:,:,c)/exp(img_d_bk),1);
%             img_w(:,:,c) = max(img_w(:,:,c)*exp(img_d_bk),1);
        end
    end
    
    %stitch two images
    img_c = max(img_w, img_ref_w);
    % figure, imagesc(img_c(:,:,2));
    
    %update reference mask
    img_mask_w = imwarp(z_ref_mask,t,'outputview',Rcb);
    img_mask_ref_w = imwarp(z_ref_mask_ref,t_ref,'outputview',Rcb);
    z_ref_mask_ref = max(img_mask_w, img_mask_ref_w);
    img_ref = img_c; 
    %clear img_c img_w img_ref_w;
end

% figure; imagesc(img_ref(:,:,1)); daspect([1 1 1]);
% figure; imagesc(z_ref_mask_ref(:,:)); daspect([1 1 1]);

% %adjust contrast
% I_int_max = double(intmax('uint16'));
% for c = 1:nchannels,
%     I = double(img_ref(:,:,c));
%     I = I/max(I(:))*I_int_max;
%     img_ref(:,:,c) = I;
% end

%write to tiff
datname_output = [datname '.tif'];

%redirect writing to folder with permission
% datname_output = strrep(datname_output,'slide_scanner','');
% datname_output = strrep(datname_output,'Confocal_retroFISH','');

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
save(datname_output,'ref_channel','z_ref_mask_ref','ch_name','-v7.3');

%memory
clear img_c img_w img_ref_w img_ref_mask img_mask_s img
clear img_mask_w img_mask_ref_w img_overlap_mask z_ref_mask_ref
clear img_w_masked img_ref_w_masked img_w_bk img_ref_w_bk

end

save([ 'image_stitching_notes_' date '_' num2str(round(rand(1)*1e6)) '.mat'], 'notes')

% fh = fopen('image_stitching_notes_17-Aug-2017_146729.tsv','w');
% for i = 1:length(notes)
%     fprintf(fh,'%s\t%s\n',notes{i}{2},notes{i}{3});
% end
% fclose(fh);

