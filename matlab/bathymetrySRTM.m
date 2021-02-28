function varargout = bathymetrySRTM(varargin)
% bathymetrySRTM Import/download NASA SRTM data files for bathymetry.
%
%	READHGT(LAT,LON) downloads the tiles and plots SRTM data corresponding 
%	to LAT and LON (in decimal degrees) coordinates (lower-left corner) 
%	from the USGS data server (needs an Internet connection and a companion  
%	file "readhgt_srtm_index.txt"). For better plot results, it is 
%	recommended to install DEM personal function available at author's 
%	Matlab page. 
%
%	LAT and/or LON can be vectors: in that case, tiles corresponding to all
%	possible combinations of LAT and LON values will be downloaded, and
%	optional output structure X will have as much elements as tiles.
%
%	READHGT(FILENAME) reads HGT data file FILENAME, must be in the form
%	"[e|w]xxx[n|s]yy.hgt[.zip]", as downloaded from SRTM data servers.
%
%	X=READHGT(...) returns a structure X containing: 
%		lat: coordinate vector of latitudes (decimal degree)
%		lon: coordinate vector of longitudes (decimal degree)
%		  z: matrix of elevations (meters, INT16 class)
%		hgt: downloaded filename(s)
%
%	X=READHGT(...,'plot') also plots the tile(s).
%
%	Additionnal options are available:
%
%	READHGT(LAT,LON,'merge'), in case of adjoining values of LAT and LON, 
%	will concatenate tiles to produce a single one.
%
%	READHGT(LAT,LON,'interp') linearly interpolates missing data.
%
%	READHGT(LAT,LON,...,'crop',[LAT1,lAT2,LON1,LON2]) crops the map using
%	latitude/longitude limits. READHGT(LAT,LON,...,'crop'), without limits
%	argument vector, crops the resulting map around existing land (reduces 
%	any sea or novalue areas at the borders).
%
%	READHGT(LAT,LON,OUTDIR) specifies output directory OUTDIR to write
%	downloaded files.
%
%	READHGT(LAT,LON,OUTDIR,URL) specifies the URL address to find HGT 
%	files (default is USGS).
%
%	Examples:
%	- to plot a map of the Paris region, France (single tile):
%		readhgt(48,2)
%
%	- to plot a map of Flores volcanic island, Indonesia (5 tiles):
%		readhgt(-9,119:123,'merge')
%
%	Informations:
%	- each file corresponds to a tile of ???? degree of a grid
%	  6000x4800 of elevation values (SRTM30 = 30 arc-seconds)
%	- elevations are of class INT16: sea level values are 0, unknown values
%	  equal -32768 (there is no NaN for INT class), use 'interp' option to
%	  fill the gaps.
%	- note that borders are included in each tile, so to concatenate tiles
%	  you must remove one row/column in the corresponding direction (this
%	  is made automatically with the 'merge' option).
%	- downloaded file is written in the current directory or optional  
%	  OUTDIR directory, and it remains there.
%	- NASA Shuttle Radar Topography Mission [February 11 to 22, 2000] 
%	  produced a near-global covering on Earth land, but still limited to 
%	  latitudes from 60S to 60N. Offshore tiles will be output as flat 0
%	  value grid.
%
%	Author: Fran�ois Beauducel <beauducel@ipgp.fr>
%		Institut de Physique du Globe de Paris
%   modified by Regis Cottereau <regis.cottereau@centralesupelec.fr>
%       CNRS, CentraleSupelec
%
%	References:
%		http://dds.cr.usgs.gov/srtm/version2_1
%
%	Acknowledgments: Yves Gaudemer, Jinkui Zhu, Greg
%
%	Created: 2012-04-22
%	Updated: 2014-05-17
%	Updated: 2015-12-04 (Regis Cottereau)

%	Copyright (c) 2014, Fran�ois Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

% Remarks (R�gis)
% 1) I did not check that the mergin was correctly done. In HGT files,
% apparently there is an overlap of 1 cell ... for these SRTM files, I do
% not know if there is overlap
% 2) I did not treat the southern part (below s60). The files are not
% labelled in the same way on topex.ucsd.edu as on 
% http://dds.cr.usgs.gov/srtm/version2_1 so a check should be made

fidx = 'readhgt_srtm30_index.txt';
% ATTENTION: this file must exist in the Matlab path
% since USGS delivers data continent-by-continent with nominative directories,
% this index file is needed to know the full path name of each tile.
urlftp = 'topex.ucsd.edu';
sz = [6000,4800]; % SRTM30 tile size
novalue = intmin('int16'); % -32768
n = 1;

makeplot = 0;
merge = 0;
decimflag = 0;
decim = 0;
inter = 0;
cropflag = 0;
crop = [];

if nargin > 0 
	makeplot = 1;%any(strcmp(varargin,'plot'));
	merge = any(strcmp(varargin,'merge'));
	kcrop = find(strcmp(varargin,'crop'));
    if ~isempty(kcrop)
        cropflag = 1;
        if (kcrop + 1) <= nargin && isnumeric(varargin{kcrop+1})
            crop = varargin{kcrop+1};
            if any(size(crop) ~= [1,4])
                error('CROP option arguments must be a 1x4 vector.')
            end
            cropflag = 2;
            crop = [minmax(crop(1:2)),minmax(crop(3:4))];
        end
    end
	inter = any(strcmp(varargin,'interp'));
	kdecim = find(strcmp(varargin,'decim'));
	if ~isempty(kdecim)
		decimflag = 1;
		if (kdecim + 1) <= nargin && isnumeric(varargin{kdecim+1})
			decim = round(varargin{kdecim+1});
			if ~isscalar(decim) || decim < 1
				error('DECIM option argument must be a positive integer.')
			end
			decimflag = 2;
		end
	end
end
nargs = makeplot + merge + cropflag + inter + decimflag;

if nargin < (2 + nargs)
	lat = str2double(filename(2:3));
	if filename(1) == 's'
		lat = -lat;
	end
	lon = str2double(filename(5:7));
	if filename(4) == 'w'
		lon = -lon;
	end
else
    % the tiles are given in 40(long) * 50(lat)
    % the min/max latitude of the tiles are -60/90
    lat = varargin{1}(:); 
    if length(lat)==1
        lat = [lat lat+1];
    end
    lon = varargin{2}(:); 
    if length(lon)==1
        lon = [lon lon+1];
    end
    lat = unique(floor((lat+10)/50));
    lon = unique(floor((lon-20)/40));
	if ~isnumeric(lon) || ~isnumeric(lat) || any(abs(lat) > 60) || any(lon < -180) || any(lat > 90) || isempty(lat) || isempty(lon)
		error('LAT and LON must be numeric and in valid SRTM interval (-60<LAT<90).');
	end
	if merge && (any(diff(lat) ~= 1) || any(diff(lon) ~= 1))
		error('With MERGE option, LAT and LON must be adjoining tiles.');
	end

	if nargin > (2 + nargs)
		out = varargin{3};
		if ~exist(out,'dir')
			error('OUTDIR must be a valid directory.')
		end
	else
		out = '.';
	end
	
	% if LAT/LON are vectors, NDGRID makes a grid of corresponding tiles
	[lat,lon] = ndgrid(lat,lon);
	f = cell(size(lat));
	for n = 1:numel(f)
		if lat(n) < 0
			slat = sprintf('s%02d',-(40+(lat(n)*50)));
		else
			slat = sprintf('n%02d',40+(lat(n)*50));
		end
		if lon(n) < 0
			slon = sprintf('w%03d',-(20+(lon(n)*40)));
		else
			slon = sprintf('e%03d',20+(lon(n)*40));
		end
		f{n} = sprintf('%s/%s%s.Bathymetry.srtm',out,slon,slat);
        [~,name] = fileparts(f{n});

		if ~exist(f{n},'file')
			if nargin > (3 + nargs)
				url = varargin{4};
				if ~ischar(url)
					error('URL must be a string.');
				end
			else
				%fsrtm = sprintf('%s/%s',fileparts(mfilename('fullpath')),fidx);
				fsrtm = fidx;
				if exist(fsrtm,'file')
					fid = fopen(fsrtm,'rt');
					idx = textscan(fid,'%s');
					fclose(fid);
					k = find(~cellfun('isempty',strfind(idx{1},sprintf('%s%s',slon,slat))));
					if isempty(k)
						%fprintf('READHGT: Warning! Cannot find %s tile in SRTM database. Consider it offshore...\n',ff);
						ff = '';
					else
                        ff = idx{1}{k(end)};
					end
				else
					error('Cannot find "%s" index file to parse SRTM database. Please download HGT file manually.',fsrtm);
				end
			end
			if isempty(ff)
				f{n} = '';
            else
                [pathname,filename,fileext] = fileparts(ff);
                filename = [filename fileext];
                nftp = ftp(urlftp);
                cd(nftp, pathname);
                mget(nftp, filename, out);
                close(nftp);
                fprintf('File "%s" downloaded from %s/%s\n',name,urlftp,ff)
            end
        else
            disp(['found file SRTM30:' name ' on local disk'])
		end
	end
end

% pre-allocates X structure (for each file/tile)
X = repmat(struct('hgt',[],'lat',[],'lon',[]),[n,1]);

if n == 1
	merge = 1;
end

for n = 1:numel(f)
	% unzips HGT file if needed
	if ~isempty(strfind(f{n},'.zip'));
		X(n).hgt = char(unzip(f{n}));
		funzip = 1;
	else
		X(n).hgt = f{n};
		funzip = 0;
	end

    if isempty(f{n})
		% offshore: empty tile...
		X(n).z = [];
    else
       
		% loads data from HGT file
        if strcmp(X(n).hgt(end-2:end),'hgt')
            fid = fopen(X(n).hgt,'rb','ieee-be');
            X(n).z = fread(fid,'*int16');
            fclose(fid);
            switch numel(X(n).z)
                case prod(sz1)
                    % srtm3 option: decimates the tile...
                    if srtm3
                        z = reshape(X(n).z,sz1);
                        X(n).z = z(1:3:end,1:3:end);
                        sz = sz3;
                    else
                        sz = sz1;
                    end
                case prod(sz3)
                    sz = sz3;
                otherwise
                    error('"%s" seems not a regular SRTM data file or is corrupted.',X(n).hgt);
            end
            X(n).z = rot90(reshape(X(n).z,sz));
            
		% loads data from DEM file
        elseif strcmp(X(n).hgt(end-2:end),'dem')
            fid = fopen(X(n).hgt,'rb','ieee-be');
            X(n).z = fread(fid,inf,'int16',0,'b');
            fclose(fid);
            X(n).z = reshape( X(n).z, sz(2), sz(1) )';
            X(n).z = X(n).z(end:-1:1,:);

		% loads data from SRTM file
        elseif strcmp(X(n).hgt(end-3:end),'srtm')
            fid = fopen(X(n).hgt,'rb','ieee-be');
            X(n).z = fread(fid,inf,'int16',0,'b');
            fclose(fid);
            X(n).z = reshape( X(n).z, sz(2), sz(1) )';
            X(n).z = X(n).z(end:-1:1,:);
        end

		% erases unzipped file if necessary
		if (funzip)
			delete(f{n});
		end
    end
    
    % builds latitude and longitude coordinates
	X(n).lon = linspace(20+(lon(n)*40),20+((lon(n)+1)*40),sz(2));
	X(n).lat = linspace(-10+(lat(n)*50),-10+((lat(n)+1)*50),sz(1))';
	% interpolates NaN (if not merged)
	if inter && ~merge
		X(n).z = fillgap(X(n).lon,X(n).lat,X(n).z,novalue);
	end
end

if merge
	% NOTE: cannot merge mixted SRTM1 / SRTM3 or discontiguous tiles
    minlat = Inf;
    maxlat = -Inf;
    for i1 = 1:numel(X)
        minlat = min(minlat,min(X(i1).lat));
        maxlat = max(maxlat,max(X(i1).lat));
    end
    minlon = Inf;
    maxlon = -Inf;
    for i1 = 1:numel(X)
        minlon = min(minlon,min(X(i1).lon));
        maxlon = max(maxlon,max(X(i1).lon));
    end
	Y.lat = linspace(minlat,maxlat,size(lat,1)*(sz(1)-1)+1)';
	Y.lon = linspace(minlon,maxlon,size(lon,2)*(sz(2)-1)+1);
	Y.z = zeros(length(Y.lat),length(Y.lon),'int16');
	for n = 1:numel(X)
		if ~isempty(X(n).z)
			Y.z((sz(1)-1)*(X(n).lat(1)-Y.lat(1)) + (1:sz(1)),(sz(2)-1)*(X(n).lon(1)-Y.lon(1)) + (1:sz(2))) = X(n).z;
		end
    end

	if cropflag
		if cropflag == 1 || isempty(crop)
			klat = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,2));
			klon = firstlast(any(Y.z ~= 0 & Y.z ~= novalue,1));
		else
			klat = find(Y.lat >= crop(1) & Y.lat <= crop(2));
			klon = find(Y.lon >= crop(3) & Y.lon <= crop(4));
		end			
		Y.lat = Y.lat(klat);
		Y.lon = Y.lon(klon);
		Y.z = Y.z(klat,klon);
	end

	if inter
		Y.z = fillgap(Y.lon,Y.lat,Y.z,novalue);
	end
end

if nargout == 0 || makeplot
	if merge
		fplot(Y.lon,Y.lat,Y.z,decim,urlftp,novalue)
	else
		for n = 1:numel(X)
			fplot(X(n).lon,X(n).lat,X(n).z,decim,urlftp,novalue)
		end
	end
end

if nargout == 3 % for backward compatibility...
	varargout{1} = X(1).lon;
	varargout{2} = X(1).lat;
	varargout{3} = X(1).z;
elseif nargout > 0
	if merge
		varargout{1} = Y;
	else
		varargout{1} = X;
	end
	if nargout == 2
		varargout{2} = f{1}; % for backward compatibility...
	end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fplot(x,y,z,decim,url,novalue)
%FPLOT plot the data using DEM function if exists, or IMAGESC

demoptions = {'latlon','legend','lake','nodecim'};

figure
if decim
	n = decim;
else
	n = 1;
end
if n > 1
	x = x(1:n:end);
	y = y(1:n:end);
	z = z(1:n:end,1:n:end);
	fprintf('READHGT: In the figure data has been decimated by a factor of %d...\n',n);
end

if exist('dem','file')
	dem(x,y,z,demoptions{:})
else
	warning('For better results you might install the function dem.m from http://www.ipgp.fr/~beaudu/matlab.html#DEM')
	z(z==novalue) = 0;
	imagesc(x,y,z);
	if exist('landcolor','file')
		colormap(landcolor(256).^1.3)
	else
		colormap(jet)
	end
	% aspect ratio (lat/lon) is adjusted with mean latitude
	xyr = cos(mean(y)*pi/180);
	set(gca,'DataAspectRatio',[1,xyr,1])

	orient tall
	axis xy, axis tight
end

title(sprintf('Data SRTM/NASA from %s',url),'FontSize',12,'Interpreter','none')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = firstlast(x)

k = find(x);
y = k(1):k(end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = minmax(x)

y = [min(x(:)),max(x(:))];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = fillgap(x,y,z,novalue)
% GRIDDATA is not efficient for large arrays, but has great advantage to be
% included in Matlab core functions! To optimize interpolation, we
% reduce the number of relevant data by building a mask of surrounding
% pixels of novalue areas... playing with linear index!

sz = size(z);
k = find(z == novalue);
k(k == 1 | k == numel(z)) = []; % removes first and last index (if exist)
if ~isempty(k)
	[xx,yy] = meshgrid(x,y);
	mask = zeros(sz,'int8');
	k2 = ind90(sz,k); % k2 is linear index in the row order
	% sets to 1 every previous and next index, both in column and row order
	mask([k-1;k+1;ind90(fliplr(sz),[k2-1;k2+1])]) = 1; 
	mask(k) = 0; % removes the novalue index
	kb = find(mask); % keeps only border values
	z(k) = int16(griddata(xx(kb),yy(kb),double(z(kb)),xx(k),yy(k)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k2 = ind90(sz,k)

[i,j] = ind2sub(sz,k);
k2 = sub2ind(fliplr(sz),j,i); % switched i and j: k2 is linear index in row order
