function d = swbd_shore(area,pfd,res)

% SWBD_SHORE  Create LandMask-Image(s) from SWBD-Tiles
%
% SWBD_SHORE( Area , Directory , Continent, Resolution )
%
% AREA = [ LonMin LonMax LatMin LatMax ]
%
% AREA = [ Lon Lat ]
%
% Directory  = Directory of Images to create
% Resolution = PixelResolution in ArcSeconds
% 
% Creates Images of LandCovering for [ 1° x 1° ] - SWBD-Tiles
%  with a Resolution of [ 1201 x 1201 ] (3arcsec). 
% 
% The Images are Gray-Scaled-PNG with the ColorValues: 
%   Land == 000, Ocean == 255, 100 < Lakes < 250, River == Blue
%
% The FileName of the Image refers to the SW-Corner
%  of the [ 1° x 1° ] - Tile: "Directory/*LON#LA_$$.png", 
%  where "*" is "e" or "w", "#" is "n" or "s", $$ is the Resolution
%
% The WaterCovering is read from SWBD-ShapeFiles using READ_SWBD.
%
%-----------------------------------------------------------------
%
% SWBD_SHORE( ZIP_Directory , Directory )
%
% Creates the Images for all SWBD-ZIP-Files ("*.zip"),
%  found in ZIP-Directory
%
%-----------------------------------------------------------------
%
% see also: READ_SWBD, SRTM_SHORE
%

Nin = nargin;

if Nin < 2
   pfd = '';
end

if Nin < 4
   res = [];
end

msg = cell(0,1);

if ~( ( isnumeric(area) & any(prod(size(area))==[2 4]) ) | chkstr(area,1) )
    msg = cat(1,msg,{'Area must be a DirectoryName or 4-Element Vector.'});
end

ok = chkstr(pfd,0);
if ~ok
    ok = ( isnumeric(pfd) & ( Nin <= 2 ) );
    if ok
       res = pfd;
       pfad = '';
    else
       msg = cat(1,msg,{'Invalid Input for Directory (or Resolution).'});
    end
end

if isempty(res)
   res = 1;
end

if ~( isnumeric(res) & ( prod(size(res)) == 1 ) )
    msg = cat(1,msg,{'Resolution must be a single or empty numeric'});
elseif ~( mod(60,res) == 0 )
    msg = cat(1,msg,{'Resolution must be a divisor of 60.'})
end

if ~isempty(pfd)
    if ~( exist(pfd,'dir') == 7 )
         msg = cat(1,msg,{sprintf('Directory "%s" doesn''t exist.',pfd)});
    end
end

if ~isempty(msg)
    msg = sprintf('%s\n',msg{:});
    error(sprintf('Invalid Inputs.\n%s',msg));
end

%*****************************************************************
% Check for Directory, look for ZIP-Files

if chkstr(area,1)

   w = dir(fullfile(area,'*.zip'));

   if ~isempty(w)
       ok = ~cat(1,w.isdir);
       if ~any(ok)
           w = {};
       else
           if ~all(ok);
               w = w(find(ok));
           end
           w = { w.name };
           w = w(:);
           w = w(:,[1 1]);
           s = size(w,1);
           ok = zeros(s,1);
           for ii = 1 : s
               [p,n,e] = fileparts(w{ii});
               [m,v] = str2vec(n);             % [ Lon Lat ]
               ok(ii) = ( isempty(m) & isequal(size(v),[1 2]) );
               if ok(ii) 
                  ok(ii) = all( isfinite(v) & ( mod(v,1) == 0 ) & ...
                               ( 0 <= v ) & ( v <= [ 180 90 ] ) );
               end
               if ok(ii)
                   w(ii,:) = { fullfile(area,[n e])  v };
               end
           end
           if ~any(ok)
               w = {};
           else
               if ~all(ok);
                   w = w(find(ok),:);
               end
           end
       end
   end

   if isempty(w)
       fprintf(1,'No valid ZIP-Files found in "%s".',area);
       return
   end

   nw = size(w,1);

else

   w = {};

   if prod(size(area)) == 2
      area = floor(area(:)');
      area = cat( 2 , area , area+1 );
      area = area([1 3 2 4]);
   end

   area([1 3]) = floor(area([1 3]));
   area([2 4]) =  ceil(area([2 4]));

   sz = area([2 4]) - area([1 3]);

   nw = prod(sz);

end

%*****************************************************************

off =  20;  % PixelOffset

siz = 3600/res;      % Size of Image

nn  = ceil(siz/600); % Number of Tiles

ext = siz/nn + 1;    % Size of Tile

pix = ext + 2 * off;

p01 = ( [0 pix] - 1/2 - off ) / siz;

ew = 'ew';
ns = 'ns';

form = '%s%3.3d%s%2.2d_%2.2d.%s';

typ = 'png';

%*****************************************************************
% Colors

ini = { 'Ocean' 'w'
        'Land'  'k'
        'River' 'b'
        'Lake'  ''
        'Isle'  'k' };

% Lakes

c = ( 100 : 250 )';
n = size(c,1);

cl = c(:,[1 1 1]);

for ii = 1 : 100
    jj = n - ii;
    cl = cat( 1 , cl , cl(1:jj,:)+ones(jj,1)*[0 0 ii] );
end

for ii = 1 : 100
    jj = n - ii;
    cl = cat( 1 , cl , cl(1:jj,:)+ones(jj,1)*[0 ii ii] );
end

cl = cl / 255;

%*****************************************************************
 
% fig = figure('position',[200 50 pix([1 1])], ...
%        'toolbar','none', ...
%        'menubar','none', ...
%        'color' , 'k' );
% 
% axe = axes('position',[0 0 1 1], ...
%      'visible','off', ...
%      'nextplot','add'  );

% Surrounding Tiles: E / NE / N / NW / W / SW / S / SE

dxy = cat( 1 , [ 1  1  0 -1 -1 -1  0  1 ] , ...
               [ 0  1  1  1  0 -1 -1 -1 ]       );

dxy = permute(dxy,[2 1]);

nc = size(dxy,1);

for ii = 1 : nw

    if isempty(w)
       iy = ceil(ii/sz(1));
       ix = ii - sz(1) * ( iy - 1 );
       f  = area([1 3]) + [ ix iy ] - 1;
       lm = f;
    else
       f  = w{ii,1};
       lm = w{ii,2};
    end

    try
       [d,f] = read_swbd(f,pfd);
    catch
        d = [];
    end
return
    if ~isempty(d)

        if 0 % Check for OceanBorder
           ok = zeros(1,nc);
           xy = cat( 2 , d.Ocean{:} );

           for jj = 1 : nc
            
           end
        end

        for jj = 1 : nc
            cc = lm + dxy(jj,:);
            try
               [c,f] = read_swbd(cc,pfd);
            catch
                warning(lasterr)
                c = [];
            end
            if ~isempty(c)
                fld = fieldnames(d);
                for ff = fld(:)'
                    v = cat(1,getfield(c,ff{1}),getfield(d,ff{1}));
                    d = setfield(d,ff{1},v);
                end
            else
                xy = cc(:);
                xy = xy(:,ones(1,4)) + [ 0 1 1 0 ; 0 0 1 1];
                d.Ocean = cat( 1 , d.Ocean , { xy } );
            end
        end

 %       delete(findobj(axe,'type','patch'));

        for jj = 1 : size(ini,1)
             c = getfield(d,ini{jj,1});
            if ~isempty(c)
                for kk = 1 : prod(size(c))
                    cc = ini{jj,2};
                    if isempty(cc)
                       cc = cl(kk,:);
                    end
                    patch(c{kk}(1,:),c{kk}(2,:),cc, ...
                      'edgecolor','none','parent',axe);
               end
            end
        end

% off =  20;  % PixelOffset
% siz = 3600/res;      % Size of Image
% nn  = ceil(siz/600); % Number of Tiles
% ext = siz/nn + 1;    % Size of Tile
% pix = ext + 2 * off;
% p01 = ( [0 pix] - 1/2 - off ) / siz;

        c = uint8(zeros(siz+1,siz+1,3));

        for ix = 1 : nn
            for iy = 1 : nn
                set(axe,'xlim',lm(1)+p01+(ix-1)/nn,'ylim',lm(2)+p01+(iy-1)/nn);
                figure(fig), drawnow, d = getframe(fig);
%%% try
                d = d.cdata(off+(1:ext),off+(1:ext),:);
%%% catch,keyboard,end
                d = d(ext:-1:1,:,:);
                jx = ( 1 : ext-(ix<nn) );
                jy = ( 1 : ext-(iy<nn) );
                c(jy+(ext-1)*(iy-1),jx+(ext-1)*(ix-1),:) = d(jy,jx,:);
            end
        end

        c = c((siz+1):-1:1,:,:);

        p = pfd;
        if isempty(p)
           p = fileparts(f);
        end

        f = sprintf(form,ew(1+(lm(1)<0)),abs(lm(1)),ns(1+(lm(2)<0)),abs(lm(2)),res,typ);

        f = fullfile(p,f);

        fprintf(1,'%s',f);

        imwrite(c,f,typ);

        fprintf(1,'\n');
        
    end

end
    
