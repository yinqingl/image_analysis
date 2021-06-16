function meta=imreadBFmeta_nd2(id)
%function meta=imreadBFmeta_nd2(id)
%
%returns metadata of nd2 image file using the BioFormats package
%
%
%OUT:
% meta.width : image width
% meta.height : image height
% meta.zsize : number of z slices
% meta.nframes : number of time frames
% meta.channels : number of channels
% meta.raw : all metadata as java hashtable
% 
% install bfmatlab
% http://downloads.openmicroscopy.org/bio-formats/
%
% Yinqing Li
% yinqingl@csail.mit.edu

% load the Bio-Formats library into the MATLAB environment
% status = bfCheckJavaPath(autoloadBioFormats);
% assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
%     'to the static Java path or add it to the Matlab path.']);

% initialize logging
% javaMethod('enableLogging', 'loci.common.DebugTools', 'INFO');


% Get the channel filler
r = bfGetReader(id, 0);

meta.width = r.getSizeX();
meta.height = r.getSizeY();
meta.zsize = r.getSizeZ();
meta.nseries = r.getSeriesCount();
meta.nframes = r.getSizeT();
meta.channels = r.getSizeC();

globalMetadata = r.getGlobalMetadata();
numSeries = r.getSeriesCount();
result = cell(numSeries, 2);

t_p = meta.nseries*meta.nframes*meta.zsize;

for s = 1:numSeries
    fprintf('Reading series #%d', s);
    r.setSeries(s - 1);
    
    numImages = r.getImageCount();
    
    %{
    pixelType = r.getPixelType();
    bpp = javaMethod('getBytesPerPixel', 'loci.formats.FormatTools', ...
                     pixelType);
    bppMax = power(2, bpp * 8);
    numImages = r.getImageCount();
    imageList = cell(numImages, 2);
    colorMaps = cell(numImages);
    
    for i = 1:numImages
        if mod(i, 72) == 1
            fprintf('\n    ');
        end
        fprintf('.');
        arr = bfGetPlane(r, i, varargin{:});

        % retrieve color map data
        if bpp == 1
            colorMaps{s, i} = r.get8BitLookupTable()';
        else
            colorMaps{s, i} = r.get16BitLookupTable()';
        end

        warning_state = warning ('off');
        if ~isempty(colorMaps{s, i})
            newMap = single(colorMaps{s, i});
            newMap(newMap < 0) = newMap(newMap < 0) + bppMax;
            colorMaps{s, i} = newMap / (bppMax - 1);
        end
        warning (warning_state);


        % build an informative title for our figure
        label = id;
        if numSeries > 1
            seriesName = char(r.getMetadataStore().getImageName(s - 1));
            if ~isempty(seriesName)
                label = [label, '; ', seriesName];
            else
                qs = int2str(s);
                label = [label, '; series ', qs, '/', int2str(numSeries)];
            end
        end
        if numImages > 1
            qi = int2str(i);
            label = [label, '; plane ', qi, '/', int2str(numImages)];
            if r.isOrderCertain()
                lz = 'Z';
                lc = 'C';
                lt = 'T';
            else
                lz = 'Z?';
                lc = 'C?';
                lt = 'T?';
            end
            zct = r.getZCTCoords(i - 1);
            sizeZ = r.getSizeZ();
            if sizeZ > 1
                qz = int2str(zct(1) + 1);
                label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
            end
            sizeC = r.getSizeC();
            if sizeC > 1
                qc = int2str(zct(2) + 1);
                label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
            end
            sizeT = r.getSizeT();
            if sizeT > 1
                qt = int2str(zct(3) + 1);
                label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
            end
        end

        save image plane and label into the list
        imageList{i, 1} = arr;
        imageList{i, 2} = label;
    end

    % save images and metadata into our master series list
%     result{s, 1} = imageList;
    %}
    % extract metadata table for this series
    seriesMetadata = r.getSeriesMetadata();
    javaMethod('merge', 'loci.formats.MetadataTools', ...
               globalMetadata, seriesMetadata, 'Global ');
%     result{s, 2} = seriesMetadata;
    p = r.getIndex(1-1, 1-1, 1-1)+1;
    p = numImages * (s-1) + p;
    
    str_xypos_ptn = {'#%d','#%02d','#%03d','#%04d','#%05d','endfor'};
    for i = 1:length(str_xypos_ptn),
        
        if strcmp(str_xypos_ptn{i},'endfor')
            display('cannot find x or y position')
            break;
        end
        
        posX_ptn = sprintf('Global X position for position %s', str_xypos_ptn{i});
        posY_ptn = sprintf('Global Y position for position %s', str_xypos_ptn{i});
        try
            posX = seriesMetadata.get(sprintf(posX_ptn, p)).value;
            posY = seriesMetadata.get(sprintf(posY_ptn, p)).value;
            break
        catch
            1;
        end
    end
    
%     try
%         posX = seriesMetadata.get(sprintf('Global X position for position #%d', p)).value;
%         posY = seriesMetadata.get(sprintf('Global Y position for position #%d', p)).value;
%     catch
%         try
%             posX = seriesMetadata.get(sprintf('Global X position for position #%02d', p)).value;
%             posY = seriesMetadata.get(sprintf('Global Y position for position #%02d', p)).value;
%         catch
%             try
%                 posX = seriesMetadata.get(sprintf('Global X position for position #%03d', p)).value;
%                 posY = seriesMetadata.get(sprintf('Global Y position for position #%03d', p)).value;
%             catch
%                 display('cannot find x or y position')
%             end
%         end
%     end
    

    
%     posZ = seriesMetadata.get(sprintf('Global Z position for position #%03d', p)).value;
    
    result{s, 1} = seriesMetadata;
    result{s, 2} = [posX,posY];
%     result{s, 3} = colorMaps;
%     result{s, 4} = r.getMetadataStore();
    result{s, 3} = r.getMetadataStore();
    fprintf('\n');
end
r.close();

meta.result = result;

    
end

            