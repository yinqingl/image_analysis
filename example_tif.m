% adjust tif image
% automatic adjust contrast
% identify nuclei area
% merge channels and write image file
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

%%
%list all tif files
datname_s = {};
%system path to the images
% path_sets = {'..\Slide1\','..\Slide2\'};
path_sets = {'..\20170130\stitch\','..\20170207\stitch\','..\20170209\stitch\'};
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

i = cellfun(@(x) ~isempty(findstr(x,'PVT6')), datname_s);
datname_s = datname_s(~i);

%print
for i_img = 1:length(datname_s)
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname))
end

%%
%configuration
% I_th_l_bk_lev = 1;  %1,2 for low or high background level
is_roi_used = 0; %use roi to select out defect in image
is_bgsub_used = 0; %use background substraction
is_I_th_prev_used = 0; %use previously determined contrast, 1 one for each individual image, 2 same for all images, set with path2I_th_prev
path2I_th_prev = '..\20161118_2\stitch\20161118_PVT6_Ecel1T4_Spp1T1_Slide18_Section7_7_TRNR002.nd2.tif.adj.tif.mat'; path2I_th_prev = strrep(path2I_th_prev,'\',filesep);
q_I_th_bk = 0.9; %percentile used to fit background distribution, 0.999
q_I_th_fg = 0.9; %percentile used to fit foreground distribution
i_composite = [3,4,0]; %red,green,blue, 0: zeros


%%
%iterate through all images
%write to composite images
% for i_img = [30],
for i_img = 1:length(datname_s),
    datname = datname_s{i_img};
    display(sprintf('%d, %s',i_img, datname));

data = imread(datname);
try,
    data_mat = load([datname '.mat']);
catch
    data_mat = struct();
    data_mat.ref_channel = 1;
end

% if sum(max(max(data))>0)<2
%     display(sprintf('skip image: %s',datname))
%     continue
% end

I_int_max = double(intmax('uint16'));

% I = data_mat.z_ref_mask_ref; %confocal max projection filter
% I_mask = bwareaopen(I, 25, 4);
%use reference channel as mask
I_zmax = data(:,:,data_mat.ref_channel);
I = adaptivethreshold(I_zmax,[250 250],0,0);
stdfilt_I = stdfilt(I,ones(3));
z_ref_mask = I.*(stdfilt_I == 0);
z_ref_mask = imdilate(z_ref_mask,ones(5));
I_mask = z_ref_mask;

%choose roi to mask the defect
if is_roi_used,
%     fig_i = [];
    for c = 1:size(data,3),
        I = double(data(:,:,c));
        try
            figure(fig_i{c});
        catch
            fig_i{c} = figure;
        end
        imagesc(I), daspect([1 1 1]);
        title(num2str(c))
    end
    x = input('which channel e.g. 1,2,3: ');
    figure(fig_i{x})
    I_roi = imrect();
    I_roi_mask=~(I_roi.createMask);
else
    I_roi_mask = ones(size(I_mask));
end

I_th_c = zeros(2,size(data,3));
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
    
    I_fg = I(I_mask & I_roi_mask); I_fg = (I_fg(I_fg>0));
    I_bk = I((~I_mask) & I_roi_mask); I_bk = (I_bk(I_bk>0));
    
%     I_th_q = quantile([I_fg;I_bk],0.99);
%     [~,n,x_fg,~] = kde((I_fg(I_fg<I_th_q)),256); n_fg = n/sum(n);
%     [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_q)),256); n_bk = n/sum(n);
%     
%     n_bk = interp1(x_bk,n_bk,x_fg)';
%     x = x_fg;
    
%     figure; plot(x, n_bk); hold on; plot(x,n_fg,'r'); 
    
%     [~,max_i_bk] = max(n_bk);
%     [~,max_i_fg] = max(n_fg);
%     i = min(max_i_bk, max_i_fg);
%     I_th_l_hi = (x(find(diff(n_fg(i:end) - n_bk(i:end)>0),1,'last')+i));  %likehood ratio test, produce high background    
    
%     [n,x] = hist((I_fg),256);
%     n_fg = histc((I_fg),x); n_fg = n_fg/sum(n_fg);
%     n_bk = histc((I_bk),x); n_bk = n_bk/sum(n_bk);

%     [n,x] = hist(log(I_fg),256);
%     n_fg = histc(log(I_fg),x); n_fg = n_fg/sum(n_fg);
%     n_bk = histc(log(I_bk),x); n_bk = n_bk/sum(n_bk);    

%     I_th_l_hi = (x(find(n_fg(i:ii) - n_bk(i:ii) <= 0,1,'last')+1));    %likehood ratio test, produce high background    
%     I_th_l_hi = exp(x(find(n_fg(1:i) - n_bk(1:i) <= 0,1,'last')+1));    %likehood ratio test, produce high background

%     I_th_h = quantile(I_fg(I_fg>I_th_l),0.95);
    
%     [n_bk,x] = hist((I_bk),256); n_bk = n_bk/sum(n_bk);
%     [~,n,x_bk,~] = kde((I_bk(I_bk<I_th_q)),256); n_bk = n/sum(n);
%     [~,x_bk_i] = max(n_bk);
%     med_bk = (x_bk(x_bk_i));
%     
% %     med_bk = median(I_bk);
% %     med_bk = exp(x(x_bk_i));
% 
%     sigma_bk = 1.4826*median(abs(I_bk-med_bk));
%     I_th_l_low = med_bk+2*sigma_bk; %keep background top 5% signal
%     I_th_l_low = med_bk+3*sigma_bk; %keep background top 1% signal
    
%     switch(I_th_l_bk_lev)
%         case 1, I_th_l = I_th_l_low;
%         case 2, I_th_l = I_th_l_hi;
%     end
    
%     %need to adjust foreground threshold based on its distribution
% %     I_fg = I_fg(I_fg>(med_bk+6*sigma_bk)); %use background dist to find foreground signal
%     [~,n,x_fg,~] = kde((I_fg(I_fg<I_th_q)),256); n_fg = n/sum(n);
%     [~,x_fg_i] = max(n_fg);
%     med_fg = (x_fg(x_fg_i));
%     sigma_fg = 1.4826*median(abs(I_fg-med_fg));
%     I_fg = I_fg(I_fg>(med_fg+3*sigma_fg)); 
   
    %use previous determined I_th
    if is_I_th_prev_used >0,
        if is_I_th_prev_used == 1,
            I_th_prev = load([datname '.adj.tif' '.mat']);
        elseif is_I_th_prev_used == 2,
            I_th_prev = load(path2I_th_prev);
        end
        I_th_prev = I_th_prev.I_th_c;
        I_th_l = I_th_prev(1,c);
        I_th_h = I_th_prev(2,c);
    else

        %background distribution, detect the tail using linear fit
        q_I_bk = [0.999:1e-4:0.9999];
        for i_f = [1:length(q_I_bk)],
            I_th_q = quantile(I_bk,q_I_bk(i_f));
            [~,n,x,c_n] = kde((I_bk(I_bk<I_th_q)),256); n = n/sum(n);
            n_polyfit = 100;
            xp = [ones(size(x')),x'];
            y = n';
            [~,i] = max(y);
        %     figure; plot(x,y)
            % p = polyfit(x(end-n_polyfit:end),y(end-n_polyfit:end),1);
            p = robustfit(xp(end-n_polyfit:end,:),y(end-n_polyfit:end),[],[],'off');
            [msglast, msgidlast] = lastwarn;
            if strcmp(msgidlast,'stats:statrobustfit:IterationLimit'),
                warning('FineId','fine') 
                display(sprintf('c = %d, quantile_th = %f',c,q_I_bk(i_f)));
    %             y = diff([0;n])';
    %             [~,i] = min(y);
    %             p = robustfit(xp(end-n_polyfit:end,:),y(end-n_polyfit:end),[],[],'off');
            end

            yfit = p(2)*x+p(1);
            yresid = y - yfit;
            v = sqrt(sum(yresid(end-n_polyfit).^2));
            z = yresid/v;
            pz = normcdf(-abs(z));
            I_th_l = x(find(pz(i+1:end)>0.05,1) + i);   %a 0.05
            I_th_l = quantile(I_bk(I_bk<I_th_l),q_I_th_bk);
            
            %approximiate the distribution using triangle
            %the tail contains less than 5% of the total distribution
            if y(i)/max(y(i)-yresid(i),1e-25) > 1/0.05    %good fit
                break
            end            
        end

    %     figure; plot(x,pz); figure; plot(x,yresid)
    %     figure; plot(x, y); hold on; plot(x, yfit, 'r')

        %saturate the top percentile of foreground distribution
        I_fg = I_fg(I_fg>I_th_l); %use background dist to find foreground signal
        if ~isempty(I_fg),
            I_th_h = quantile(I_fg,q_I_th_fg);   %allow a little saturation
        else
            I_th_h = max(I(:));
        end

        I_th_c(1,c) = I_th_l;
        I_th_c(2,c) = I_th_h;
    end
       
    %adjust contrast
    I_adj = I;
    I_adj(I_adj > I_th_h) = I_th_h;
    I_adj = I_adj - I_th_l;
    I_adj(I_adj<0) = 0;
    I_adj = I_adj/max(I_adj(:))*I_int_max;
    
    data_adj(:,:,c) = I_adj;
end

%write to tif
datname_output = [datname '.adj.tif'];
% data_adj = uint16(round(data_adj));
data_w = zeros(size(data_adj(:,:,[1:length(i_composite)])));
for i = [1:length(i_composite)],
    c = i_composite(i);
    if c>=1 && c<=size(data_adj,3),
        data_w(:,:,i) = data_adj(:,:,c);
    else
        data_w(:,:,i) = zeros(size(data_adj(:,:,1)));
    end
end
data_w = uint16(round(data_w));
imwrite(data_w,datname_output,'TIFF');

% t = Tiff(datname_output, 'w');
% tagstruct.ImageLength = size(data_w, 1);
% tagstruct.ImageWidth = size(data_w, 2);
% tagstruct.Compression = Tiff.Compression.None;
% % tagstruct.Compression = Tiff.Compression.LZW;        % compressed
% tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
% % tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
% tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
% tagstruct.BitsPerSample = 16;                        % float data
% tagstruct.SamplesPerPixel = size(data_w,3);
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
% t.setTag(tagstruct);
% t.write(data_w);
% t.close();


%write channel threshold intensity
if ~is_I_th_prev_used
    datname_output = [datname_output '.mat'];
    save(datname_output,'I_th_c','-v7.3');
end

end

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