function [d,fname] = read_swbd(file,outdir)

% READ_SWBD  Read SWBD-ShapeFiles (SRTM Water Body Data)
%
% [Struct,File] =  READ_SWBD( File )
%
% [Struct,File] =  READ_SWBD( [ Lon Lat ] )
%
% The FileName is build like "*LON#LA$.ext", 
%  where "*" is "e" or "w",
%        "#" is "n" or "s", 
%        "$" is the ContinentCode
%      "ext" is the Extension
%
% The Characters for Continent-Code are: 
%
%  "n"  NorthAmerica, "s" SouthAmerica, "e" Eurasia
%  "a" Australia, "f" Africa, "l" or "i" for Islands
%
% The Extension are "dbf" and "shp" (and "shx") or "zip"
%  for a ZIP-File, contains the single ShapeFiles.
%
% For a ZIP-File the UNIX-enviroment and UNZIP is required
%
%--------------------------------------------------------
% The FACC_CODE of SWBD-Files are:
%
% BA040 Ocean Struct.Ocean
% BH140 River Struct.River
% BH080 Lake  Struct.Lake
%
% Multiple Segments are Islands,
%
% For Oceans in Struct.Land and 
% for River and Lakes in Struct.Isle.
%
%--------------------------------------------------------
%
% see also: SHP_READ, READ_HGT, SRTM_SHORE
%

d = [];
url = 'http://dds.cr.usgs.gov/';
fsrtm = 'readhgt_swbd_index.txt';
%-------------------------------------------------------
% Check Input

if isnumeric(file) & ( prod(size(file)) == 2 )

   ew = 'ew';
   ns = 'ns';

   x = file(1);
   y = file(2);

   ns = ns(1+(y<0));
   ew = ew(1+(x<0));
   
   if exist(fsrtm,'file')
       fid = fopen(fsrtm,'rt');
       idx = textscan(fid,'%s');
       fclose(fid);
       k = find(~cellfun('isempty',strfind(idx{1},sprintf('%s%3.3d%s%2.2d',ew,abs(x),ns,abs(y)))),1);
       if isempty(k)
           %fprintf('READHGT: Warning! Cannot find %s tile in SRTM database. Consider it offshore...\n',ff);
           file = '';
       else
           file = idx{1}{k(end)};
       end
   else
       error('Cannot find "%s" index file to parse SRTM database. Please download HGT file manually.',fsrtm);
   end

elseif ~chkstr(file,1)

   error('File must be a String or [ Lat Lon ].');

end

%-------------------------------------------------------
% Check for ShapeFile

[~,fname] = fileparts(file);
[ok,name,fname,zip] = checkfile(fullfile(outdir,fname),'');

if ~ok
    try
        
%         ff = unzip([url,file],'./');

        tmp = tempname;
        mkdir(tmp)
        warning('no correct!!!')
        ftmp = sprintf('%s/%s%3.3d%s%2.2d.zip',tmp,ew,abs(x),ns,abs(y));
        warning('wget here')
        [~,aa]=system('which wget');
        [s,w] = system(sprintf([aa(1:end-1),' --no-check-certificate -O %s %s%s'],ftmp,url,file));
           if s
              disp(w)
           end
        ff = unzip(ftmp,'./');
        delete(ftmp)                
        if ~isempty(ff)
            fprintf('File "%s" downloaded from %s\n',file,url);
        end
    catch
        disp(['no file ' sprintf('%s%3.3d%s%2.2d',ew,abs(x),ns,abs(y)) '.zip on http server'])
    end
%     for c = 'einsalf'
%         [ok,name,fname,zip] = checkfile(file,c);
%         if ok, break, end
%     end
else
    disp(['found file SWBD:' name ' on local disk'])
end

% if ~ok
%     if exist(file,'file') == 2
%        m = 'Invalid';
%     else
%        m = 'Cann''t find';
%     end
%     error(sprintf('%s ShapeFile: "%s".',m,file));
% end
% 
%-------------------------------------------------------
% Check for ZIP-File, unzip

% pfd = '';  % TempDir
% 
% if ~isempty(zip)
% 
%     [pfd,tmp] = fileparts(tempname);
%     pfd = fullfile(pfd,sprintf('swbd_%s',tmp));
%     while exist(pfd,'dir') == 7
%          pfd = sprintf('%s0',pfd);
%     end
% 
%     ok = mkdir(pfd);
%     if ~ok
%         error(sprintf('Cann''t create TempDir "%s" for UNZIP.',pfd));
%     end
% 
%     cmd = sprintf('unzip -d "%s" "%s"',pfd,zip);
% 
%     [s,w] = unix(cmd);
% 
%     if ~( s == 0 )
% 
%         m = sprintf('Error call UNZIP for UNIX: %s',cmd);
%         if ~isempty(w)
%             m = sprintf('%s\n%s',m,w);
%         end
% 
%         [s,w] = unix(sprintf('rm -r "%s"',pfd));
%         if exist(pfd,'dir') == 7
%            m = sprintf('%s\nCann''t remove TempDir: "%s".',m,pfd);
%            if ~isempty(w)
%                m = sprintf('%s\n%s',m,w);
%            end
%         end     
% 
%         error(m)
% 
%     end
% 
%     fname = fullfile(pfd,name);
% 
% end

%-------------------------------------------------------
% Read ShapeFile, using SHP_READ

m = cell(0,1);

try
  shp = shp_read(fname);

catch

   m = cat(1,m,{sprintf('Error call SHP_READ( %s ).\n%s',fname,lasterr)});

end

%-------------------------------------------------------
% Remove TempDir if ZIP

% if ~isempty(zip) & ( exist(pfd,'dir') == 7 )
% 
%     [s,w] = unix(sprintf('rm -r "%s"',pfd));
%     if exist(pfd,'dir') == 7
%        m = cat(1,m,{sprintf('Cann''t remove TempDir: "%s".',pfd)});
%        if ~isempty(w)
%            m = cat(1,m,{sprintf('%s',w)});
%        end
%     end     
% 
% end
%  
% if ~isempty(m)
%     error(sprintf('%s\n',m{:}))
% end

%-------------------------------------------------------
% Check Data

ok = isstruct(shp);
if ok
   fld = fieldnames(shp);
   upp = upper(fld);
   for f = { 'X' 'Y' 'FACC_CODE' }; % 'BoundingBox' }
       k = strcmp(upp,f{1});
       kk = any(k);
       ok = ( ok & kk );
       if kk & ~any(strcmp(fld,f{1}))
          kk = find(k);
          shp = rnfield(shp,fld{kk},f{1});
       end
   end
end

if ~ok
    error('Invalid OutputStructure by SHP_READ.');
end

ini = { 'Ocean' 'BA040'
        'River' 'BH140'
        'Lake'  'BH080' };

d = cat( 1 , ini(:,1) , {'Land'} , {'Isle'} );

d = d(:,[1 1]);

d(:,2) = { {cell(0,1)} };

for ii = 1 : size(ini,1)

    [xy,isl] = get_coord(shp,ini{ii,2});

    d{ii,2} = {xy};

    if strcmp(ini{ii,1},'Ocean')
       d{end-1,2}{1} = cat(1,d{end-1,2}{1},isl);  % Land
    else
       d{end,2}{1} = cat(1,d{end,2}{1},isl);      % Isle
    end

end

d = permute(d,[2 1]);
d = struct(d{:});

if ~isempty(zip)
    fname = zip;
end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,name,file,zip] = checkfile(file,c);

zip = '';

[pfad,name,ext] = fileparts(file);

ext = { 'dbf' 'shp' };

name = [ name  c ];

file = name;
if ~isempty(pfad)
    file = fullfile(pfad,file);
end

ok = 1;
for e = ext
    ok = ( ok & ( exist(sprintf('%s.%s',file,e{1}),'file') == 2 ) );
end

if ok
   f = which(sprintf('%s.%s',file,ext{1}));
   if ~isempty(f)
       [p,f] = fileparts(f);
       file = fullfile(p,f);
   end
   return
end

zip = sprintf('%s.%s',file,'zip');

ok = ( exist(zip,'file') == 2 );
if ~ok
    return
end

f = which(zip);
if ~isempty(f)
    zip = f;
end

[s,w] = unix(sprintf('unzip -l "%s"',zip));

ok = ( s == 0 );
if ok
   for e = ext
       ok = ( ok & ~isempty(findstr(w,sprintf('%s.%s',name,e{1}))) );
   end
end

%************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [xy,isl] = get_coord(d,code)

xy = cell(0,1);
isl = cell(0,1);

ic = strcmp( { d.FACC_CODE } , code );

if ~any(ic)
    return
end

n = sum(ic);

ic = find(ic);

xy  = cell(n,1);
isl = cell(n,1);
isl(:) = { cell(0,1) };

ok  = zeros(n,1);

for ii = 1 : n

    zz = cat( 1 , d(ic(ii)).X , d(ic(ii)).Y );

    jj = isnan(zz(1,:));   % Seperator for Islands

    if any(jj)

       kk = cumsum(jj,2);  % Following Islands start with 1

       jj = find(jj);
       kk(jj) = NaN;       % Don't read NaN-Seperator

       m = max(kk);        % Number of Islands (following Segments)

       isl{ii} = cell(m,1);

       for il = 1 : m
           ll = find( kk == il );
           isl{ii}{il} = zz(:,ll);
           kk(ll) = NaN;
       end

       kk = ~isnan(kk);

       if ~any(kk)
           zz = [];
       else
           kk = find(kk);
           zz = zz(:,kk);
       end

    end

    ok(ii) = ~isempty(zz);
    if ok(ii)
       xy{ii} = zz;
    end

end

if ~all(ok)
    ok = find(ok);
    xy = xy(ok);
end

isl = cat(1,isl{:});

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function ok = chkstr(str,opt)


% CHKSTR  Checks Input for String
%
%  ok = chkstr(str,Option)
%
%  Option ~= 0 ==> only true for nonempty Strings
%
%   default:   Option == 0  (true for empty Strings)
%

 
if nargin < 2
   opt = 0;
end

ok = ( strcmp( class(str) , 'char' )      & ...
       ( prod(size(str)) == size(str,2) ) & ...
       ( isequal(opt,0) | ~isempty(str) )         );


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [ok,str] = chkcstr(str,opt)

% CHKCSTR  Checks Input for CellString, contains Strings !
%
%  [ok,str] = chkcstr(str,Option)
%
%  Option ~= 0 ==> CharacterArrays not allowed,
%
%   default: Option == 0   ==>  CharacterArrays --> CellString
%
 
if nargin < 2
   opt = 0;
end

if strcmp(class(str),'char') & isequal(opt,0)
   n = size(str,1);
   if n == 1
      str = strrep(str,char(32),char(1));
   end
   str = cellstr(str);
   if n == 1
      str = strrep(str,char(1),char(32));
   end
end

ok = iscellstr(str);
if ~ok
   return
end

try
  s = cat(2,str{:});
catch
  ok = 0;
  return
end
 
ok = ( strcmp(class(s),'char')  &  ( prod(size(s)) == size(s,2) ) );


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% SHP_READ
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%*********************************************************************

function varargout = shp_read(varargin)
%SHP_READ Read vector feature coordinates and attributes from a shapefile.
%
%   S = SHP_READ(FILENAME) returns an N-by-1 version 2 geographic data
%   structure (geostruct2) array, S, containing an element for each
%   non-null spatial feature in the shapefile. 
%
%   S combines feature coordinates/geometry and attribute values.
%
%   [S, A] = SHP_READ(FILENAME) returns an N-by-1 geostruct2 array,
%   S, and a parallel N-by-1 attribute structure array, A.  Each array
%   constains an element for each non-null spatial feature in the
%   shapefile.
%
%   S = shp_read(FILENAME,PARAM1,VAL1,PARAM2,VAL2,...) or
%   [S, A] = shp_read(FILENAME,PARAM1,VAL1,PARAM2,VAL2,...) returns a
%   subset of the shapefile contents in S or [S,A], as determined by the
%   parameters 'RecordNumbers','BoundingBox','Selector', or 'Attributes'.
%   In addition, the parameter 'UseGeoCoords' can be used in cases where
%   you know that the X- and Y-coordinates in the shapefile actually
%   represent longitude and latitude.
%
%   The shapefile format was defined by the Environmental Systems Research
%   Institute (ESRI) to store nontopological geometry and attribute
%   information for spatial features. A shapefile consists of a main file,
%   an index file, and an xBASE table. All three files have the same base
%   name and are distinguished by the extensions 'SHP', 'SHX', and 'DBF',
%   respectively (e.g., base name 'concord_roads' and filenames
%   'concord_roads.SHP', 'concord_roads.SHX', and 'concord_roads.DBF').
%
%   FILENAME can be the base name, or the full name of any one of the
%   component files.  SHP_READ reads all three files as long as they
%   exist in the same directory and have valid file extensions. If the
%   main file (with extension SHP) is missing, SHP_READ returns an
%   error. If either of the other files is missing, SHP_READ returns a
%   warning.
%
%   Supported shape types
%   ---------------------
%   SHP_READ supports the ordinary 2D shape types: 'Point', 'Multipoint',
%   'PolyLine', and 'Polygon'. ('Null Shape' features are may also be
%   present in a Point, Multipoint, PolyLine, or Polygon shapefile, but
%   are ignored.) SHP_READ does not support any 3D or "measured" shape
%   types: 'PointZ', 'PointM', 'MultipointZ', 'MultipointM', 'PolyLineZ',
%   'PolyLineM', 'PolygonZ', 'PolylineM', or 'Multipatch'.
%
%   Output structure
%   ----------------
%   The fields in the output structure array S and A depend on (1) the type
%   of shape contained in the file and (2) the names and types of the
%   attributes included in the file:
%
%     Field name            Field contents          Comment
%     ----------            -----------------       -------
%     'Geometry'            Shape type string
%
%     'BoundingBox'         [minX minY;             Omitted for shape type
%                            maxX maxY]             'Point'
%
%     'X' or 'Lon'          Coordinate vector       NaN-separators used
%                                                   in multi-part PolyLine
%     'Y' or 'Lat'          Coordinate vector       and Polygon shapes
%
%     Attr1                 Value of first          Included in output S
%                           attribute               if output A is omitted
%
%     Attr2                 Value of second         Included in output S
%                           attribute               if output A is omitted
%
%     ...                   ...                     ...
%
%   The names of the attribute fields (listed above as Attr1, Attr2, ...)
%   are determined at run-time from the xBASE table (with extension 'DBF')
%   and/or optional, user-specified parameters.  There may be many
%   attribute fields, or none at all.
%
%   'Geometry' field
%   -----------------
%   The 'Geometry' field will be one of the following values: 'Point',
%   MultiPoint', 'Line', or 'Polygon'.  (Note that these match the standard
%   shapefile types except for shapetype 'Polyline' the value of the
%   Geometry field is simply 'Line'.
%
%   'BoundingBox' field
%   -------------------
%   The 'BoundingBox' field contains a 2-by-2 numerical array specifying
%   the minimum and maximum feature coordinate values in each dimension
%   (min([x, y]); max([x, y] where x and y are N-by-1 and contain the
%   combined coordinates of all parts of the feature).
%
%   Coordinate vector fields ('X','Y' or 'Lon','Lat')
%   -------------------------------------------------
%   These are 1-by-N arrays of class double.  For 'Point' shapes, they
%   are 1-by-1.  In the case of multi-part 'Polyline' and 'Polygon' shapes,
%   NaN are added to separate the lines or polygon rings.  In addition,
%   terminating NaNs are added to support horizontal concatation of the
%   coordinate data from multiple shapes.
%
%   Attribute fields
%   ----------------
%   Attribute names, types, and values are defined within a given
%   shapefile. The following four types are supported: Numeric, Floating,
%   Character, and Date. SHP_READ skips over other attribute types.
%   The field names in the output shape structure are taken directly from
%   the shapefile if they contain no spaces or other illegal characters,
%   and there is no duplication of field names (e.g., an attribute named
%   'BoundingBox', 'PointData', etc. or two attributes with the names).
%   Otherwise the following 'name mangling' is applied: Illegal characters
%   are replaced by '_'. If the first character in the attribute name is
%   illegal, a leading 'Z' is added. Numerals are appended if needed to
%   avoid duplicate names. The attribute values for a feature are taken
%   from the shapefile and stored as doubles or character arrays:
%
%   Attribute type in shapefile     MATLAB storage
%   ---------------------------     --------------
%       Numeric                     double (scalar)
%       Float                       double (scalar)
%       Character                   char array
%       Date                        char array
%
%   Parameter-Value Options
%   -----------------------
%   By default, shp_read returns an entry for every non-null feature and
%   creates a field for every attribute.  Use the first three parameters
%   below (RecordNumbers, BoundingBox, and Selector) to be selective about
%   which features to read.  Use the 4th parameter (Attributes) to control
%   which attributes to keep.  Use the 5th (UseGeoCoords) to control the
%   output field names.
% 
%   Name            Description of Value        Purpose
%   ----            --------------------        -------
%   
%   RecordNumbers   Integer-valued vector,      Screen out features whose
%                   class double                record numbers are not
%                                               listed.
% 
%   BoundingBox     2-by-(2,3, or 4) array,     Screen out features whose
%                   class double                bounding boxes fail to
%                                               intersect the selected box.
%                                             
%   Selector        Cell array containing       Screen out features for
%                   a function handle and       which the function, when
%                   one or more attribute       applied to the the
%                   names.  Function must       corresponding attribute
%                   return a logical scalar.    values, returns false.
%                   
%   Attributes      Cell array of attribute     Omit attributes that are
%                   names                       not listed. Use {} to omit
%                                               all attributes. Also sets
%                                               the order of attributes in
%                                               the structure array.
%
%   UseGeoCoords    Scalar logical              If true, replace X and Y
%                                               field names with 'Lon' and
%                                               'Lat', respectively.
%                                               Defaults to false.
%
%   Examples
%   --------
%   % Read the entire concord_roads.shp shapefile, including the attributes
%   % in concord_roads.dbf.
%   S = shp_read('concord_roads.shp');
%
%   % Restrict output based on a bounding box and read only two
%   % of the feature attributes.
%   bbox = [2.08 9.11; 2.09 9.12] * 1e5;
%   S = shp_read('concord_roads','BoundingBox',bbox,...
%                 'Attributes',{'STREETNAME','CLASS'});
%
%   % Select the class 4 and higher road segments that are at least 200
%   % meters in length.  Note the use of an anonymous function in the
%   % selector.
%   S = shp_read('concord_roads.shp',...
%         'Selector',{@(v1,v2) (v1 >= 4) && (v2 >= 200),'CLASS','LENGTH'});
%   N = hist([S.CLASS],1:7)
%   hist([S.LENGTH])
%
%   % Read world-wide city names and locations in latitude and longitude.
%   % (Note presence of 'Lat' and 'Lon' fields.)
%   S = shp_read('worldcities.shp', 'UseGeoCoords', true)
%
%   See also SHP_INFO, UPDATEGEOSTRUCT.

%   Copyright 1996-2004 The MathWorks, Inc.  
%   $Revision: 1.1.10.7 $  $Date: 2004/12/18 07:46:34 $

%   Reference
%   ---------
%   ESRI Shapefile Technical Description, White Paper, Environmental
%   Systems Research Institute, July 1998.
%   (http://arconline.esri.com/arconline/whitepapers/ao_/shapefile.pdf)

nargoutchk(0,2,nargout);
switch(nargout)
    case {0,1}, separateAttributes = false;
    otherwise,  separateAttributes = true;
end

% Parse function inputs.
try
[filename, recordNumbers, boundingBox, selector, attributes, useGeoCoords] ...
    = parseInputs(varargin{:});
catch
  filename = varargin{1};
  recordNumbers = [];
  boundingBox   = [];
  selector      = [];
  attributes    = [];
  useGeoCoords  = []; 
end

% Try to open the SHP, SHX, and DBF files corresponding to  the filename
% provided. Selectively validate the header, including the shape type code.
[shpFileId, shxFileId, dbfFileId, headerTypeCode] ...
    = openShapeFiles(filename,'shp_read');

% Get the file offset for the content of each shape record.
if (shxFileId ~= -1)
    contentOffsets = readIndexFromSHX(shxFileId);
else
    contentOffsets = constructIndexFromSHP(shpFileId);
end

% Select which records to read.
records2read = selectRecords(shpFileId, dbfFileId, headerTypeCode,...
                   contentOffsets, recordNumbers, boundingBox, selector);
                   
% Read the shape coordinates from the SHP file into a cell array.
[shapeData, shapeDataFieldNames] ...
    = shpread(shpFileId, headerTypeCode, contentOffsets(records2read));
 
% Read the attribute data from the DBF file into a cell array.
[attributeData, attributeFieldNames] ...
    = dbfread(dbfFileId,records2read,attributes);

% Optionally rename coordinate field names.
if useGeoCoords
    shapeDataFieldNames{strmatch('X',shapeDataFieldNames,'exact')} = 'Lon';
    shapeDataFieldNames{strmatch('Y',shapeDataFieldNames,'exact')} = 'Lat';
end

% Concatenate the cell arrays, if necessary and convert to struct(s).
varargout = constructOutput(shapeData, attributeData,...
              shapeDataFieldNames, attributeFieldNames, separateAttributes);

% Clean up.
closeFiles([shpFileId, shxFileId, dbfFileId]);

%--------------------------------------------------------------------------
function outputs = constructOutput(shapeData, attributeData,...
              shapeDataFieldNames, attributeFieldNames, separateAttributes)

if separateAttributes
    if ~isempty(attributeData)
        A = cell2struct(attributeData,genvarname(attributeFieldNames),2);
    else
        A = [];
    end
    S = cell2struct(shapeData,shapeDataFieldNames,2);
    outputs = {S, A};
else
    if ~isempty(attributeData)
        % Concatenate the shape data field names for the current shape type
        % and the attribute field names from the DBF file (if available).
        % Ensure value, non-duplicate structure field names.
        %%% featureFieldNames = [shapeDataFieldNames,...
        %%%   genvarname(attributeFieldNames,feval(mapgate('geoReservedNames')))];

        featureFieldNames = chkname([shapeDataFieldNames,attributeFieldNames]);

        S = cell2struct([shapeData, attributeData],featureFieldNames,2);
    else
        S = cell2struct(shapeData,shapeDataFieldNames,2);
    end
    outputs = {S};
end

%--------------------------------------------------------------------------
function records2read = selectRecords(shpFileId, dbfFileId, headerTypeCode, ...
                          contentOffsets, recordNumbers, boundingBox, selector)
% Select record numbers to read as constrained by shapefile record types,
% user-specified record numbers, user-specified bounding box, and
% user-specified attribute-based selector function.

% Initialize selection to include all non-null shape records.
records2read = recordsMatchingHeaderType(shpFileId,contentOffsets,headerTypeCode);

% Narrow selection based on user-specified record numbers.
if ~isempty(recordNumbers)
    records2read = intersect(recordNumbers,records2read);
end

% Narrow selection based on bounding box.
if ~isempty(boundingBox)
    bbSubscripts = getshapetype(headerTypeCode,'BoundingBoxSubscripts');
    if hasBoundingBox(headerTypeCode)
        records2read = recordsIntersectingBox(shpFileId,...
            bbSubscripts,contentOffsets,boundingBox,records2read);
    else
        records2read = recordsWithPointsInbox(shpFileId,...
            bbSubscripts,contentOffsets,boundingBox,records2read);
    end
end

% Finalize selection based on selector function.
if (dbfFileId ~= -1) && ~isempty(selector)
    records2read = recordsMatchingSelector(dbfFileId,selector,records2read);
end

%---------------------------------------------------------------------------
function recs = recordsMatchingHeaderType(shpFileId,contentOffsets,headerTypeCode)
% Select the records that match the headerTypeCode.

totalNumRecords = length(contentOffsets);
recordMatchesHeaderType = false(1,totalNumRecords);
for n = 1:totalNumRecords
	fseek(shpFileId,contentOffsets(n),'bof');
	recordTypeCode = fread(shpFileId,1,'uint32','ieee-le');
	recordMatchesHeaderType(n) = (recordTypeCode == headerTypeCode);
end
recs = find(recordMatchesHeaderType);

%--------------------------------------------------------------------------
function answer = hasBoundingBox(shapeTypeCode)
fieldNames = getshapetype(shapeTypeCode,'ShapeDataFieldNames');
answer = ~isempty(strmatch('BoundingBox',fieldNames,'exact'));

%--------------------------------------------------------------------------
function recs = recordsIntersectingBox(...
    shpFileId, bbSubscripts, contentOffsets, box, recs)
% Select the records with bounding boxes intersecting the specified box.

currentNumberOfRecs = numel(recs);
intersectsBox = false(1,currentNumberOfRecs);
for k = 1:currentNumberOfRecs
    n = recs(k);
	fseek(shpFileId,contentOffsets(n) + 4,'bof');
    bbox = fread(shpFileId,8,'double','ieee-le');
    intersectsBox(k) = boxesIntersect(box,bbox(bbSubscripts));
end
recs(~intersectsBox) = [];

%--------------------------------------------------------------------------
function result = boxesIntersect(a,b)
result = ~(any(a(2,:) < b(1,:)) || any(b(2,:) < a(1,:)));

%--------------------------------------------------------------------------
function recs = recordsWithPointsInbox(...
    shpFileId, bbSubscripts, contentOffsets, box, recs)
% Select the point records for locations within the specified box.
% Note: This version assumes 2D-only.

currentNumberOfRecs = numel(recs);
insideBox = false(1,currentNumberOfRecs);
for k = 1:currentNumberOfRecs
    n = recs(k);
	fseek(shpFileId,contentOffsets(n) + 4,'bof');
    point = fread(shpFileId,[1 2],'double','ieee-le');
    insideBox(k) = all(box(1,:) <= point(1,:)) && all(point(1,:) <= box(2,:));
end
recs(~insideBox) = [];

%--------------------------------------------------------------------------
function recs = recordsMatchingSelector(fid,selector,recs)
% Apply selector to DBF file to refine list of records to read.

% The first byte in each record is a deletion indicator
lengthOfDeletionIndicator = 1;

% Initialize things...
info = dbfinfo(fid);
selectfcn  = selector{1};
fieldnames = selector(2:end);

% Determine the position, offset, and format string for each field 
% specified by the selector.  If any fieldnames fail to get a match,
% return without altering the list of records, and issue a warning.
allFieldNames = {info.FieldInfo.Name};
for l = 1:numel(fieldnames)
    m = strmatch(fieldnames{l}, allFieldNames ,'exact');
    if isempty(m)
        wid = sprintf('%s:%s:badSelectorFieldName',getcomp,mfilename);
        wrn = sprintf('Selector field name ''%s'' %s\n%s',...
                 fieldnames{l},'doesn''t match a shapefile attribute name.',...
                 '         Ignoring selector.');
        warning(wid,wrn)
        return;
    end
    position(l)  = m;
    offset(l)    = sum([info.FieldInfo(1:(m-1)).Length]) ...
                   + lengthOfDeletionIndicator;
   formatstr{l} = sprintf('%d*uint8=>char',info.FieldInfo(m).Length);
end

% Check each record in the current list to see if it satisifies the
% selector.
satisfiesSelector = false(1,numel(recs));
for k = 1:numel(recs)
    n = recs(k);
    for l = 1:numel(position)
        m = position(l);
        fseek(fid,info.HeaderLength + (n-1)*info.RecordLength + offset(l),'bof');
        data = fread(fid,info.FieldInfo(m).Length,formatstr{l});
        values(l) = feval(info.FieldInfo(m).ConvFunc,data');
    end
    satisfiesSelector(k) = feval(selectfcn,values{:});
end

% Remove records that don't satisfy the selector.
recs(~satisfiesSelector) = [];

%--------------------------------------------------------------------------
function contentOffsets = readIndexFromSHX(shxFileId)
% Get record content offsets (in bytes) from shx file.

fileHeaderLength    = 100;
recordHeaderLength  =   8;
bytesPerWord        =   2;
contentLengthLength =   4;

fseek(shxFileId,fileHeaderLength,'bof');
contentOffsets = recordHeaderLength + ...
    bytesPerWord * fread(shxFileId,inf,'uint32',contentLengthLength);

%--------------------------------------------------------------------------
function contentOffsets = constructIndexFromSHP(shpFileId)
% Get record content offsets (in bytes) from shp file.

bytesPerWord        =   2;
fileHeaderLength    = 100;
recordHeaderLength  =   8;
contentLengthLength =   4;

fseek(shpFileId,24,'bof');
fileLength = bytesPerWord * fread(shpFileId,1,'uint32','ieee-be');
lengthArray = [];
recordOffset = fileHeaderLength;
while recordOffset < fileLength,
    fseek(shpFileId,recordOffset + recordHeaderLength - contentLengthLength,'bof');
    contentLength = bytesPerWord * fread(shpFileId,1,'uint32','ieee-be');
    lengthArray(end + 1,1) = contentLength;
    recordOffset = recordOffset + recordHeaderLength + contentLength;
end
contentOffsets = fileHeaderLength ...
                 + cumsum([0;lengthArray(1:end-1)] + recordHeaderLength);

%--------------------------------------------------------------------------
function [shpdata, fieldnames] = ...
    shpread(shpFileId, headerTypeCode, contentOffsets)
% Read designated shape records.

shapeTypeLength = 4;
readfcn    = getshapetype(headerTypeCode,'ShapeRecordReadFcn');
fieldnames = getshapetype(headerTypeCode,'ShapeDataFieldNames');
shpdata = cell(numel(contentOffsets),numel(fieldnames));

for k = 1:numel(contentOffsets)
    fseek(shpFileId,contentOffsets(k) + shapeTypeLength,'bof');
    shpdata(k,:) = feval(readfcn,shpFileId);
end

%--------------------------------------------------------------------------
function [attributeData, attributeFieldNames] ...
    = dbfread(fid, records2read, requestedFieldNames)
% Read specified records and fields from a DBF file.  Fields will follow
% the order given in REQUESTEDFIELDNAMES.

% Return empties if there's no DBF file.
if (fid == -1)
    attributeData = [];
    attributeFieldNames = {};
    return;
end

info = dbfinfo(fid);
fields2read = matchFieldNames(info,requestedFieldNames);
attributeFieldNames = {info.FieldInfo(fields2read).Name};

% The first byte in each record is a deletion indicator
lengthOfDeletionIndicator = 1;

% Loop over the requested fields, reading the attribute data.
attributeData = cell(numel(records2read),numel(fields2read));
for k = 1:numel(fields2read),
    n = fields2read(k);
    fieldOffset = info.HeaderLength ...
                  + sum([info.FieldInfo(1:(n-1)).Length]) ...
                  + lengthOfDeletionIndicator;
    fseek(fid,fieldOffset,'bof');
    formatString = sprintf('%d*uint8=>char',info.FieldInfo(n).Length);
    skip = info.RecordLength - info.FieldInfo(n).Length;
    data = fread(fid,[info.FieldInfo(n).Length info.NumRecords],formatString,skip);
    attributeData(:,k) = feval(info.FieldInfo(n).ConvFunc,data(:,records2read)');
end

%--------------------------------------------------------------------------
function fields2read = matchFieldNames(info, requestedFieldNames)
% Determine which fields to read.

allFieldNames = {info.FieldInfo.Name};
if isempty(requestedFieldNames)
    if ~iscell(requestedFieldNames)
        % Default case: User omitted the parameter, return all fields.
        fields2read = 1:info.NumFields;
    else
        % User supplied '{}', skip all fields.
        fields2read = [];
    end
else
    % Match up field names to see which to return.
    fields2read = [];
    for k = 1:numel(requestedFieldNames)
        index = strmatch(requestedFieldNames{k},allFieldNames,'exact');
        if isempty(index)
            wid = sprintf('%s:%s:nonexistentAttibuteName',getcomp,mfilename);
            wrn = sprintf('Attribute name ''%s'' %s\n%s',requestedFieldNames{k},...
                     'doesn''t match a shapefile attribute name.',...
                     '         It will be ignored.');
            warning(wid,wrn)
        end
        for l = 1:numel(index)
            % Take them all in case of duplicate names.
            fields2read(end+1) = index(l);
        end
    end
end

%--------------------------------------------------------------------------
function [filename, recordNumbers, boundingBox, selector, attributes, useGeoCoords] ...
    = parseInputs(varargin)

validParameterNames = ...
    {'RecordNumbers','BoundingBox','Selector','Attributes','UseGeoCoords'};
checknargin(1, 1 + 2*numel(validParameterNames), nargin, mfilename);

% FILENAME is the only required input.
filename = varargin{1};
checkinput(filename, {'char'},{'vector'},mfilename,'FILENAME',1);

% Assign defaults for optional inputs.
recordNumbers = [];
boundingBox = [];
selector = [];
attributes = [];
useGeoCoords = false;

% Identify and validate the parameter name-value pairs.
for k = 2:2:nargin
    parName = checkstrs(varargin{k}, validParameterNames, mfilename, ...
        sprintf('PARAM%d',k/2), k);
    switch parName
      case 'RecordNumbers'
        checkExistence(k, nargin, mfilename, 'vector of record numbers', parName);
        recordNumbers = checkRecordNumbers(varargin{k+1},k+1);

      case 'BoundingBox'
        checkExistence(k, nargin, mfilename, 'bounding box', parName);
        checkboundingbox(varargin{k+1},mfilename,'BOUNDINGBOX',k+1);
        boundingBox = varargin{k+1};
        
      case 'Selector'
        checkExistence(k, nargin, mfilename, 'selector', parName);
        selector = checkSelector(varargin{k+1},k+1);
        
      case 'Attributes'
        checkExistence(k, nargin, mfilename,'attribute field names', parName);
        attributes = checkAttributes(varargin{k+1},k+1);
      
      case 'UseGeoCoords'
        checkExistence(k, nargin, mfilename, 'geo-coordinates flag (T/F)', parName);
        checkinput(varargin{k+1}, {'logical'}, {'scalar'}, mfilename, 'USEGEOCOORDS', k+1);
        useGeoCoords = varargin{k+1};
        
      otherwise
        eid = sprintf('%s:%s:internalProblem',getcomp,mfilename);
        error(eid,'Internal problem: unrecognized parameter name: %s',parName);
    end
end

%--------------------------------------------------------------------------

function checkExistence(position, nargs, fcnName, ...
                        propertyDescription, propertyName)
% Error if missing the property value following a property name.

if (position + 1 > nargs)
    eid = sprintf('%s:%s:missingParameterValue',getcomp,fcnName);
    error(eid,...
          'Expected %s to follow parameter name ''%s''.',...
          propertyDescription, propertyName);
end

%--------------------------------------------------------------------------

function recordNumbers = checkRecordNumbers(recordNumbers, position)

checkinput(recordNumbers, {'numeric'},...
              {'nonempty','real','nonnan','finite','positive','vector'},...
              mfilename, 'recordNumbers', position);

%--------------------------------------------------------------------------

function selector = checkSelector(selector, position)
% SELECTOR should be a cell array with a function handle followed by
% strings. Do NOT try to check the strings against actual attributes,
% because we won't know the attributes until we've read the DBF file.
% Instead, we issue a warning later on in matchFieldNames.

checkinput(selector, {'cell'}, {}, mfilename, 'SELECTOR', position);

if numel(selector) < 2
    eid = sprintf('%s:%s:selectorTooShort',getcomp,mfilename);
    error(eid,'Expected SELECTOR to have at least two elements.');
end

if ~isa(selector{1},'function_handle')
    eid = sprintf('%s:%s:selectorMissingFcnHandle',getcomp,mfilename);
    error(eid,'Expected the first element of SELECTOR to be a function handle.');
end

for k = 2:numel(selector)
    if ~ischar(selector{k})
        eid = sprintf('%s:%s:selectorHasNonStrAttrs',getcomp,mfilename);
        error(eid,'Expected the %s element of SELECTOR to be a string.',...
              num2ordinal(k));
    end
end

%--------------------------------------------------------------------------

function attributes = checkAttributes(attributes, position)

checkinput(attributes, {'cell'}, {}, mfilename, 'ATTRIBUTES', position);

for k = 1:numel(attributes)
    if ~ischar(attributes{k})
        eid = sprintf('%s:%s:nonCharAttribute',getcomp,mfilename);
        error(eid,'Expected ATTRIBUTES to contain only strings.');
    end
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function v = chkname(v,mode)

[ok,v] = chkcstr(v,0);

if ~ok
    return
end

n = prod(size(v));

if nargin < 2
   mode = 0;    % Both
end

%---------------------------------------------------------
if ~( mode == 2 )
%---------------------------------------------------------

% Replace invalid Characters by "_"

for ii = 1 : n
    w  = v{ii};
    jj = ~( ( ( '0' <= w ) & ( w <= '9' ) ) | ...
            ( ( 'A' <= w ) & ( w <= 'Z' ) ) | ...
            ( ( 'a' <= w ) & ( w <= 'z' ) )       );
    w(jj) = '_';
    v{ii} = w;
end

%---------------------------------------------------------
end
%---------------------------------------------------------

%---------------------------------------------------------
if ~( mode == 1 )
%---------------------------------------------------------

% Append duplicate Names with '_'

for ii = n : -1 : 1
    jj = cat(2,(1:(ii-1)),((ii+1):n));
    while any(strcmp(v{ii},v(jj)))
          v{ii} = cat(2,v{ii},'_');
    end
end

%---------------------------------------------------------
end
%---------------------------------------------------------


%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% PRIVATE Files for SHP_READ
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%*********************************************************************

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function checkboundingbox(...
               bbox, function_name, variable_name, argument_position)
%CHECKBOUNDINGBOX Check validity of bounding box array.
%   CHECKBOUNDINGBOX(...
%              BBOX, FUNCTION_NAME, VARIABLE_NAME, ARGUMENT_POSITION)
%   ensures that the bounding box array is a 2-by-2 array of double with
%   real, finite values, and that in each column the second value always
%   exceeds the first.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/12/13 02:52:25 $

try
    do_checkboundingbox(getcomp, bbox,...
                    function_name, variable_name, argument_position);
catch
    rethrow(lasterror);
end

%----------------------------------------------------------------------

function do_checkboundingbox(component, bbox,...
    function_name, variable_name, argument_position)

checkinput(bbox, {'double'},{'real','nonnan'},...
           function_name,variable_name,argument_position);
        
if size(bbox,1) ~= 2 || size(bbox,2) ~= 2 || ndims(bbox) ~= 2
    eid = sprintf('%s:%s:invalidBBoxSize',component,function_name);
    msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                   upper(function_name), num2ordinal(argument_position), ...
                   variable_name);
    msg2 = sprintf('to be a 2-by-2 array.');
    throwerr(eid, '%s\n%s', msg1, msg2);
end

if ~all(bbox(1,:) <= bbox(2,:))
    eid = sprintf('%s:%s:invalidBBoxOrder',component,function_name);
    msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
                   upper(function_name), num2ordinal(argument_position), ...
                   variable_name);
    msg2 = sprintf('to have element (2,k) greater than or equal to element (1,k).');
    throwerr(eid, '%s\n%s', msg1, msg2);
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function checkinput(a, classes, attributes, function_name, variable_name, ...
                    argument_position)
%CHECKINPUT Check validity of array.
%
%   CHECKINPUT(A,CLASSES,ATTRIBUTES,FUNCTION_NAME,VARIABLE_NAME,
%   ARGUMENT_POSITION) checks the validity of the array A and issues a
%   formatted error message if it is invalid.
%
%   CLASSES is a cell array of strings containing the set of classes that A
%   is expected to belong.  For example, CLASSES could be {'logical',
%   'cell'} if A is required to be either a logical array, cell array, or
%   cell.  The string 'numeric' is interpreted as an abbreviation for the
%   classes uint8, uint16, uint32, int8, int16, int32, single, double.
%
%   ATTRIBUTES is a cell array of strings containing the set of attributes
%   that A must satisfy.  For example, if ATTRIBUTES is {'real' 'nonempty'
%   'finite'}, then A must be real and nonempty, and it must contain only
%   finite values.  The supported list of attributes includes:
%   
%       real             vector              row            column
%       scalar           twod                2d             nonsparse
%       nonempty         integer             nonnegative    positive
%       nonnan           finite              nonzero        even
%       odd
%
%   FUNCTION_NAME is a string containing the function name to be used in
%   the formatted error message.
%
%   VARIABLE_NAME is a string containing the documented variable name to be
%   used in the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.

%   Copyright 1996-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/02/01 21:58:02 $

if ~iscell(classes)
  classes = {classes};
end
if ~iscell(attributes)
  attributes = {attributes};
end
checkinput_mex(a, classes, attributes, function_name, variable_name, ...
               argument_position);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function varargout = checknargin(varargin)
%CHECKNARGIN Check number of input arguments.
%   CHECKNARGIN(LOW,HIGH,NUM_INPUTS,FUNCTION_NAME) checks whether 
%   NUM_INPUTS is in the range indicated by LOW and HIGH.  If not, CHECKNARGIN 
%   issues a formatted error message using the string in FUNCTION_NAME.
%
%   LOW should be a scalar nonnegative integer.
%
%   HIGH should be a scalar nonnegative integer or Inf.
%
%   FUNCTION_NAME should be a string.
%
%   ERR = CHECKNARGIN(...) returns any error structure encounterd 
%   during the validation without rethrowing the error. 
%
%   [MSG, ID] = CHECKNARGIN(...) returns any error message in MSG
%   and error identifier in ID.  MSG and ID will be [] if no error
%   is encountered during the validation.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2003/08/23 05:55:02 $

% Input arguments are not checked for validity.

fcn = @do_checknargin;

[varargout, needRethrow] = ...
    checkfunction(fcn, nargout, varargin{:});

if needRethrow
   rethrow(lasterror);    
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_checknargin(component, low, high, numInputs, function_name)
% Main function for checknargin.
% COMPONENT is the name for the error ID's first component. 

if numInputs < low
  msgId = sprintf('%s:%s:tooFewInputs', component, function_name);
  if low == 1
    msg1 = sprintf('Function %s expected at least 1 input argument', ...
                   upper(function_name));
  else
    msg1 = sprintf('Function %s expected at least %d input arguments', ...
                   upper(function_name), low);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  throwerr(msgId, '%s\n%s', msg1, msg2);
  
elseif numInputs > high
  msgId = sprintf('%s:%s:tooManyInputs', component, function_name);

  if high == 1
    msg1 = sprintf('Function %s expected at most 1 input argument', ...
                   upper(function_name));
  else
    msg1 = sprintf('Function %s expected at most %d input arguments', ...
                   upper(function_name), high);
  end
  
  if numInputs == 1
    msg2 = 'but was called instead with 1 input argument.';
  else
    msg2 = sprintf('but was called instead with %d input arguments.', ...
                   numInputs);
  end
  
  throwerr(msgId, '%s\n%s', msg1, msg2);
end

  
%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function out = checkstrs(in, valid_strings, function_name, ...
                         variable_name, argument_position)
%CHECKSTRS Check validity of option string.
%   OUT = CHECKSTRS(IN,VALID_STRINGS,FUNCTION_NAME,VARIABLE_NAME, ...
%   ARGUMENT_POSITION) checks the validity of the option string IN.  It
%   returns the matching string in VALID_STRINGS in OUT.  CHECKSTRS looks
%   for a case-insensitive nonambiguous match between IN and the strings
%   in VALID_STRINGS.
%
%   VALID_STRINGS is a cell array containing strings.
%
%   FUNCTION_NAME is a string containing the function name to be used in the
%   formatted error message.
%
%   VARIABLE_NAME is a string containing the documented variable name to be
%   used in the formatted error message.
%
%   ARGUMENT_POSITION is a positive integer indicating which input argument
%   is being checked; it is also used in the formatted error message.

%   Copyright 1996-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2004/02/01 21:58:05 $

% Except for IN, input arguments are not checked for validity.

try
    if ~ischar(in) || ndims(in) > 2 || size(in,1) > 1
        id = sprintf('%s:%s:nonStrInput', getcomp, function_name);
        throwerr(id,...
          'Function %s expected its %s argument, %s,\nto be a character string.',...
          upper(function_name), num2ordinal(argument_position), variable_name);
    end

    matches = strncmpi(in,valid_strings,numel(in));
    if sum(matches) == 1
        out = valid_strings{matches};
    else
        out = substringMatch(valid_strings(matches));
        if isempty(out)
            failedToMatch(valid_strings, sum(matches), function_name,...
                          argument_position,variable_name,in);
        end
    end
catch
    rethrow(lasterror);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = substringMatch(strings)
%   STR = substringMatch(STRINGS) looks at STRINGS (a cell array of
%   strings) to see whether the shortest string is a proper substring of
%   all the other strings.  If it is, then substringMatch returns the
%   shortest string; otherwise, it returns the empty string.

if isempty(strings)
    str = '';
else
    len = cellfun('prodofsize',strings);
    [tmp,sortIdx] = sort(len);
    strings = strings(sortIdx);
    
    start = regexpi(strings(2:end), ['^' strings{1}]);
    if isempty(start) || (iscell(start) && any(cellfun('isempty',start)))
        str = '';
    else
        str = strings{1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function failedToMatch(valid_strings, num_matches, function_name,...
                       argument_position, variable_name, in)
% Convert valid_strings to a single string containing a space-separated list
% of valid strings.

list = '';
for k = 1:length(valid_strings)
    list = [list ', ' valid_strings{k}];
end
list(1:2) = [];

msg1 = sprintf('Function %s expected its %s input argument, %s,', ...
               upper(function_name), num2ordinal(argument_position), ...
               variable_name);
msg2 = 'to match one of these strings:';

if num_matches == 0
    msg3 = sprintf('The input, ''%s'', did not match any of the valid strings.', in);
    id = sprintf('%s:%s:unrecognizedStringChoice', getcomp, function_name);

else
    msg3 = sprintf('The input, ''%s'', matched more than one valid string.', in);
    id = sprintf('%s:%s:ambiguousStringChoice', getcomp, function_name);
end

throwerr(id,'%s\n%s\n\n  %s\n\n%s', msg1, msg2, list, msg3);

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function closeFiles(fileIds)
%CLOSEFILES Close multiple files.
%   Close any files in the list of file IDs that are actually open.

%   Copyright 1996-2003  The MathWorks, Inc.  
%   $Revision: 1.1.10.2 $ $Date: 2003/08/01 18:23:30 $

for k = 1:length(fileIds)
    if ~isempty(fopen(fileIds(k)))
        fclose(fileIds(k));
    end
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function fieldInfo = dbffieldinfo(S)
%DBFFIELDINFO Construct default DBF field info structure array from geostruct.
%
%   FIELDINFO = DBFFFIELDINFO(S) analyzes a geographic data structure S and
%   constructs a default fieldinfo structure suitable for use with SHAPEWRITE,
%   where it controls the layout of the DBF file. If a custom DBF layout is
%   needed or desired, you can use the output of DBFFIELDINFO as a starting
%   point and modify it before calling SHAPEWRITE.  If S has no attribute
%   fields, then FIELDINFO will be empty.

%   Copyright 2003-2004 The MathWorks, Inc.  
%   $Revision: 1.1.8.1 $  $Date: 2004/12/18 07:46:21 $

% Determine what types of fields to write for each attribute.
attributeNames = fields(S);
[t1,fIndex,t2] = setxor(attributeNames,{'Geometry','BoundingBox','X','Y','Lat','Lon'});
attributeNames = attributeNames(sort(fIndex));

% Return empty if there are no attributes
if isempty(attributeNames)
    fieldInfo = [];
end

% Keep it simple for now: Support only types 'N' and 'C'
% and assume blindly that all features are consistent with the
% first in terms of MATLAB storage classes.
for k = 1:numel(attributeNames)
    dataClass = class(S(1).(attributeNames{k}));
    switch(dataClass)
        
        case 'double'
            fieldType = 'N';
            v = [S.(attributeNames{k})];
            if all(v == 0)
                numRightOfDecimal = 0;
                fieldLength = 2;
            else
                numLeftOfDecimal = ceil(log10(max(abs(v))));
                if all(v == floor(v))
                    numRightOfDecimal = 0;
                    fieldLength = 1 + numLeftOfDecimal;
                else
                    numRightOfDecimal = 6;  % Hard-code for now
                    fieldLength = 1 + numLeftOfDecimal + 1 + numRightOfDecimal;
                end
            end
            if numRightOfDecimal == 0
                fmt = sprintf('%s%dd','%',fieldLength);
            else
                fmt = sprintf('%s%d.%df','%',fieldLength,numRightOfDecimal);
            end
                
        case 'char'
            fieldType = 'C';
            v = char({S.(attributeNames{k})});
            minCharFieldLength = 2;
            fieldLength = max(size(v,2), minCharFieldLength);
            numRightOfDecimal = 0;
            fmt = sprintf('%s%ds','%-',fieldLength);
            
        otherwise
            wid = sprintf('%s:%s:unsupportedDataClass',getcomp,mfilename);
            msg = sprintf('Omitting unsupported data class: %s',dataClass);
            warning(wid,msg);
            continue;
    end
    
    fieldInfo(k).AttributeName = attributeNames{k};
    fieldInfo(k).FieldName = attributeNames{k};
    fieldInfo(k).FieldType = fieldType;
    fieldInfo(k).FieldLength = fieldLength;
    fieldInfo(k).FieldDecimalCount = numRightOfDecimal;
    fieldInfo(k).Format = fmt;
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function info = dbfinfo(fid)
%DBFINFO Read header information from DBF file.
%   FID File identifier for an open DBF file.
%   INFO is a structure with the following fields:
%      Filename       Char array containing the name of the file that was read
%      DBFVersion     Number specifying the file format version
%      FileModDate    A string containing the modification date of the file
%      NumRecords     A number specifying the number of records in the table
%      NumFields      A number specifying the number of fields in the table
%      FieldInfo      A 1-by-numFields structure array with fields:
%         Name        A string containing the field name 
%         Type        A string containing the field type 
%         ConvFunc    A function handle to convert from DBF to MATLAB type
%         Length      A number of bytes in the field
%      HeaderLength   A number specifying length of the file header in bytes
%      RecordLength   A number specifying length of each record in bytes

%   Copyright 1996-2003  The MathWorks, Inc.  
%   $Revision: 1.1.10.2 $  $Date: 2003/08/01 18:23:31 $

[version, date, numRecords, headerLength, recordLength] = readFileInfo(fid);
fieldInfo = getFieldInfo(fid);

info.Filename     = fopen(fid);
info.DBFVersion   = version;
info.FileModDate  = date;
info.NumRecords   = numRecords;
info.NumFields    = length(fieldInfo);
info.FieldInfo    = fieldInfo;
info.HeaderLength = headerLength;
info.RecordLength = recordLength;

%----------------------------------------------------------------------------
function [version, date, numRecords, headerLength, recordLength] = readFileInfo(fid)
% Read from File Header.

fseek(fid,0,'bof');

version = fread(fid,1,'uint8');

year  = fread(fid,1,'uint8') + 1900;
month = fread(fid,1,'uint8');
day   = fread(fid,1,'uint8');

dateVector = datevec(sprintf('%d/%d/%d',month,day,year));
dateForm = 1;% dd-mmm-yyyy
date = datestr(dateVector,dateForm);

numRecords   = fread(fid,1,'uint32');
headerLength = fread(fid,1,'uint16');
recordLength = fread(fid,1,'uint16');

%----------------------------------------------------------------------------
function fieldInfo = getFieldInfo(fid)
% Form FieldInfo by reading Field Descriptor Array.
%
% FieldInfo is a 1-by-numFields structure array with the following fields:
%       Name      A string containing the field name 
%       Type      A string containing the field type 
%       ConvFunc  A function handle to convert from DBF to MATLAB type
%       Length    A number equal to the length of the field in bytes

lengthOfLeadingBlock    = 32;
lengthOfDescriptorBlock = 32;
lengthOfTerminator      =  1;
fieldNameOffset         = 16;  % Within table field descriptor
fieldNameLength         = 11;

% Get number of fields.
fseek(fid,8,'bof');
headerLength = fread(fid,1,'uint16');
numFields = (headerLength - lengthOfLeadingBlock - lengthOfTerminator)...
               / lengthOfDescriptorBlock;

% Read field lengths.
fseek(fid,lengthOfLeadingBlock + fieldNameOffset,'bof');
lengths = fread(fid,[1 numFields],'uint8',lengthOfDescriptorBlock - 1);

% Read the field names.
fseek(fid,lengthOfLeadingBlock,'bof');
data = fread(fid,[fieldNameLength numFields],...
             sprintf('%d*uint8=>char',fieldNameLength),...
             lengthOfDescriptorBlock - fieldNameLength);
data(data == 0) = ' '; % Replace nulls with blanks
names = cellstr(data')';

% Read field types.
fseek(fid,lengthOfLeadingBlock + fieldNameLength,'bof');
dbftypes = fread(fid,[numFields 1],'uint8=>char',lengthOfDescriptorBlock - 1);

% Convert DBF field types to MATLAB types.
typeConv = dbftype2matlab(upper(dbftypes));

% Return a struct array.
fieldInfo = cell2struct(...
    [names;  {typeConv.MATLABType}; {typeConv.ConvFunc}; num2cell(lengths)],...
    {'Name', 'Type',                'ConvFunc',          'Length'},1)';

%----------------------------------------------------------------------------
function typeConv = dbftype2matlab(dbftypes)
% Construct struct array with MATLAB types & conversion function handles.

typeLUT = ...
    {'N', 'double', @str2double2cell;...   % DBF numeric
     'F', 'double', @str2double2cell;...   % DBF float
     'C', 'char',   @cellstr;...           % DBF character
     'D', 'char',   @cellstr};             % DBF date

unsupported = struct('MATLABType', 'unsupported', ...
                     'ConvFunc',   @cellstr);
                     
% Unsupported types: Logical,Memo,N/ANameVariable,Binary,General,Picture

numFields = length(dbftypes);
if numFields ~= 0
  typeConv(numFields) = struct('MATLABType',[],'ConvFunc',[]);
end
for k = 1:numFields
    idx = strmatch(dbftypes(k),typeLUT(:,1));
    if ~isempty(idx)
        typeConv(k).MATLABType = typeLUT{idx,2};
        typeConv(k).ConvFunc   = typeLUT{idx,3};
    else
        typeConv(k) = unsupported;
    end
end

%----------------------------------------------------------------------------
function out = str2double2cell(in)
% Translate IN, an M-by-N array of class char, to an M-by-1 column vector
% OUT, of class double.  IN may be blank- or null-padded. If IN(k,:) does
% not represent a valid scalar value, then OUT(k) has value NaN.
out = num2cell(str2double(cellstr(in)));

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function comp=getcomp()
%GETCOMP Returns the error id component name. 
%
%   COMP  = GETCOMP Returns the component name for an error ID

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2003/08/01 18:23:40 $

comp = 'map';

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function result = getshapetype(shapeTypeCode,requestOrQuery)
%GETSHAPETYPE   Get information about a shape type.
%   Returns a single value, based on the second argument:
%     'TypeString'             -- Return a string
%     'IsValid'                -- Return a scalar logical
%     'IsSupported'            -- Return a scalar logical
%     'BoundingBoxSubscripts'  -- Return a 1-by-n double array
%     'ShapeRecordReadFcn'     -- Return a function handle
%     'ShapeDataFieldNames'    -- Return a cell array of string.

%   Copyright 1996-2003  The MathWorks, Inc.  
%   $Revision: 1.1.10.4 $  $Date: 2004/02/01 22:01:28 $

lutFields = { 'TypeCode',...
              'TypeString',...
              'IsValid',...
              'IsSupported',...
              'BoundingBoxSubscripts',...
              'ShapeRecordReadFcn',...
              'ShapeDataFieldNames' };

% Three kinds of bounding box subscripts
bbs2D = [1 2; 3 4];
bbsZ  = [1 2 5; 3 4 6]; % Ignore M for now, otherwise use [1 2 5 7; 3 4 6 8]
bbsM  = [1 2 7; 3 4 8];

fld1  = {'Geometry','X','Y'};                        % Single Point
fldn  = {'Geometry','BoundingBox','Length','X','Y'}; % Multiple Points

% Code  String      Valid  Supported  BBxSub Fcn  FieldNames
typeLUT = {...
   -1, 'Not Valid',   false, false, [],    [], {''};... 
    0, 'Null Shape',  true,  true,  [],    [], {''};... 
    1, 'Point',       true,  true,  bbs2D, @readPoint,      fld1 ; ...
    3, 'PolyLine',    true,  true,  bbs2D, @readPolyLine,   fldn ; ...
    5, 'Polygon',     true,  true,  bbs2D, @readPolygon,    fldn ; ...
    8, 'MultiPoint',  true,  true,  bbs2D, @readMultiPoint, fldn ; ...
   11, 'PointZ',      true,  false, bbsZ,  [], {''};...
   13, 'PolyLineZ',   true,  false, bbsZ,  [], {''};... 
   15, 'PolygonZ',    true,  true,  bbs2D, @readPolygon,    fldn ; ... %%% false, bbsZ,  [], {''} ; ... 
   18, 'MultiPointZ', true,  false, bbsZ,  [], {''};... 
   21, 'PointM',      true,  false, bbsM,  [], {''};... 
   23, 'PolyLineM',   true,  false, bbsM,  [], {''};... 
   25, 'PolygonM',    true,  false, bbsM,  [], {''};... 
   28, 'MultiPointM', true,  false, bbsM,  [], {''};... 
   31, 'MultiPatch',  true,  false, bbsZ,  [], {''};... 
  };
notValidRow = 1;
types = [typeLUT{:,1}];

% MAINTENANCE NOTE: To add support for additional types, add more rows
% to the type look up table (typeLUT), but be sure to keep 'Not Valid'
% in the first row.

row = find(shapeTypeCode == types);
if length(row) ~= 1
    row = notValidRow;
end

col = strmatch(lower(requestOrQuery),lower(lutFields));
if length(col) ~= 1;
    eid = sprintf('%s:%s:internalProblem',getcomp,mfilename);
    error(eid,'Internal error: Invalid second argument in private function.');
end

result = typeLUT{row,col};

%---------------------------------------------------------------------------
function shp = readPoint(fid)

point = fread(fid,[2 1],'double','ieee-le')';
shp = {'Point', point(1), point(2)};

%---------------------------------------------------------------------------
function shp = readMultiPoint(fid)

boundingBox    = fread(fid,4,'double','ieee-le')';
numPoints      = fread(fid,1,'uint32','ieee-le');
points         = fread(fid,[2 numPoints],'double','ieee-le')';
shp = {'MultiPoint', boundingBox([1 2; 3 4]),  numPoints, points(:,1)', points(:,2)'};

%---------------------------------------------------------------------------
function shp = readPolyLine(fid)

boundingBox    = fread(fid,4,'double','ieee-le')';
numParts       = fread(fid,1,'uint32','ieee-le');
numPoints      = fread(fid,1,'uint32','ieee-le');
partOffsets    = fread(fid,[1 numParts],'uint32','ieee-le');
points         = fread(fid,[2 numPoints],'double','ieee-le')';
[x,y] = organizeParts2D(partOffsets,points);
shp = {'Line', boundingBox([1 2; 3 4]), numPoints, x', y'};

%---------------------------------------------------------------------------
function shp = readPolygon(fid)

boundingBox    = fread(fid,4,'double','ieee-le')';
numParts       = fread(fid,1,'uint32','ieee-le');
numPoints      = fread(fid,1,'uint32','ieee-le');
partOffsets    = fread(fid,[1 numParts],'uint32','ieee-le');
points         = fread(fid,[2 numPoints],'double','ieee-le')';
[x,y] = organizeParts2D(partOffsets,points);
shp = {'Polygon', boundingBox([1 2; 3 4]), numPoints, x', y'};

%---------------------------------------------------------------------------
function [x,y] = organizeParts2D(partOffsets,points)

numParts  = size(partOffsets,2);
numPoints = size(points,1);
% Initialize x and y to be column vectors of NaN
% with length numPoints * numParts
x = NaN + zeros(numPoints + numParts, 1);
y = x;
partStart = 1 + [partOffsets numPoints];
partEnd   = [partOffsets(2:end) numPoints];
for k = 1:numParts
    rowsIn  = partStart(k):partEnd(k);
    rowsOut = rowsIn + k - 1; % Skip rows, leaving NaNs in the gaps
    x(rowsOut,1) = points(rowsIn,1);
    y(rowsOut,1) = points(rowsIn,2);
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function [shpFileId, shxFileId, dbfFileId, headerTypeCode] = ...
    openShapeFiles(filename,callingFcn)
%OPENSHAPEFILES Try to open .SHP, .SHX, and .DBF files.
%   Deconstruct a shapefile name that may include any of the standard
%   extensions (in either all lower of all upper case). Try to open the
%   corresponding SHP, SHX, and DBF files, returning a file ID for each one
%   if successful, -1 if not. Check the header file type, shapefile
%   version, and the shape type code found in the header.

%   Copyright 1996-2005  The MathWorks, Inc.  
%   $Revision: 1.1.10.3.14.1 $  $Date: 2005/01/14 21:54:16 $
        
% See if filename has an extension and extract basename.
[basename, shapeExtensionProvided] = deconstruct(filename);

% Open the SHP file, check shapefile code and version, and construct a
% qualified basename to ensure consistency in the event of duplicate
% filenames on the path.
[basename, ext] = checkSHP(basename,shapeExtensionProvided);

% Open the SHP file with the qualified basename and read/check the shape
% type code in the file header.
shpFileId = fopen([basename ext]);
headerTypeCode = readHeaderTypeCode(shpFileId,callingFcn);

% Open the SHX and DBF files with the qualified basename.
shxFileId = openSHX(basename,callingFcn);
dbfFileId = openDBF(basename);

%--------------------------------------------------------------------------
function [basename, shapeExtensionProvided] = deconstruct(filename)

shapefileExtensions = {'.shp','.shx','.dbf'};
[pathstr,name,ext] = fileparts(filename);

if isempty(ext)
    basename = filename;
    shapeExtensionProvided = false;
else
    if any(strcmp(lower(ext),shapefileExtensions))
        basename = fullfile(pathstr,name);
        shapeExtensionProvided = true;
    else
        % Make sure to allow filename  = 'test.jnk' where the full
        % shapefile name is actually 'test.jnk.shp'.
        basename = filename;
        shapeExtensionProvided = false;
    end
end

%--------------------------------------------------------------------------
function [basename, ext] = checkSHP(basename,shapeExtensionProvided)

shpFileId = fopen([basename '.shp']);
if (shpFileId == -1)
    shpFileId = fopen([basename '.SHP']);
end;

if (shpFileId == -1)
    if shapeExtensionProvided == false
        [pathstr,name,ext] = fileparts(basename);
        if ~isempty(ext)
            eid = sprintf('%s:%s:invalidExtension', getcomp, mfilename);
            msg = sprintf('Filename %s has an invalid extension. %s', basename,...
                  'If included must be .shp, .shx, or .dbf.');
            error(eid,'%s',msg)
        else
            eid = sprintf('%s:%s:failedToOpenSHP1', getcomp, mfilename);
            msg = sprintf('Failed to open both %s.shp and %s.SHP.',basename,basename);
            error(eid,'%s',msg)
        end
    else
        eid = sprintf('%s:%s:failedToOpenSHP2', getcomp, mfilename);
        msg = sprintf('Failed to open both %s.shp and %s.SHP.',basename,basename);
        error(eid,'%s',msg)
    end
end

standardShapefileCode = 9994;
fileCode = fread(shpFileId,1,'uint32','ieee-be');
if fileCode ~= standardShapefileCode
    fname = fopen(shpFileId);
    fclose(shpFileId);
    eid = sprintf('%s:%s:notAShapefile', getcomp, mfilename);
    msg = sprintf('Invalid shapefile %s.', fname);
    error(eid,'%s',msg)
end

versionSupported = 1000;
fseek(shpFileId,28,'bof');
version = fread(shpFileId,1,'uint32','ieee-le');
if version ~= versionSupported
    fclose(shpFileId);
    eid = sprintf('%s:%s:unsupportShapefileVersion', getcomp, mfilename);
    msg = sprintf('Unsupported shapefile version %d.', version);
    error(eid,msg)
end

% Construct fully qualified basename and get SHP extension.
[pathstr,name,ext] = fileparts(fopen(shpFileId));
if ~isempty(pathstr)
    basename = fullfile(pathstr,name);
else
    basename = fullfile('.',name);
end

fclose(shpFileId);

%--------------------------------------------------------------------------
function shxFileId = openSHX(basename,callingFcn)

if strcmp(callingFcn,'shaperead')
    msgForMissingSHX = 'Will build index from SHP file.';
else
    msgForMissingSHX = 'Depending on DBF file to get number of records.';
end

shxFileId = fopen([basename '.shx'],'r','ieee-be');
if (shxFileId == -1)
    shxFileId = fopen([basename '.SHX'],'r','ieee-be');
end

if (shxFileId == -1)
    wid = sprintf('%s:%s:missingSHX', getcomp, mfilename);
    msg = sprintf('Failed to open file %s.shx or file %s.SHX. %s',...
            basename, basename, msgForMissingSHX);
    warning(wid,'%s',msg)
end

%--------------------------------------------------------------------------
function dbfFileId = openDBF(basename)

dbfFileId = fopen([basename '.dbf'],'r','ieee-le');
if (dbfFileId == -1)
    dbfFileId = fopen([basename '.DBF'],'r','ieee-le');
end

if (dbfFileId == -1)
    wid = sprintf('%s:%s:missingDBF', getcomp, mfilename);
    msg = sprintf('Failed to open file %s.dbf or file %s.DBF. %s',...
            basename, basename,...
            'Shape output structure will have no attribute fields.');
    warning(wid,'%s',msg)
end

%--------------------------------------------------------------------------
function headerTypeCode = readHeaderTypeCode(shpFileId,callingFcn)

% Read the type code from a shapefile header and check to see if it's (1)
% valid and (2) supported by shaperead.  If it's not supported by
% shaperead, generate an error if called from shaperead but only a warning
% if called by shapeinfo.

fseek(shpFileId,32,'bof');
headerTypeCode = fread(shpFileId,1,'uint32','ieee-le');

if ~getshapetype(headerTypeCode,'IsValid')
    eid = sprintf('%s:%s:invalidShapeTypeCode', getcomp, mfilename);
    msg = sprintf('Invalid shape type (type code = %g).',headerTypeCode);
    fclose(shpFileId);
    error(eid, msg);
end

if ~getshapetype(headerTypeCode,'IsSupported')
    typeString = getshapetype(headerTypeCode,'TypeString');
    eid = sprintf('%s:%s:unsupportedShapeType', getcomp, mfilename);
    msg = sprintf('Unsupported shape type %s (type code = %g).',...
                  typeString,headerTypeCode);
    if strcmp(callingFcn,'shaperead')
        fclose(shpFileId);
        error(eid, msg);
    else
        warning(eid, msg);
    end
end

%*********************************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function throwerr(id,varargin)
%THROWERR constructs and rethrows an error.
%   THROWERR(ID,VARGARGIN) constructs an error structure
%   with the ID and inputs from varargin and rethrows this structure.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2003/08/01 18:19:46 $

err.message = sprintf(varargin{:});
err.identifier = id;
rethrow(err);

