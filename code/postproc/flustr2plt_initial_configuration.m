%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: flustr2plt_v2(flu,str,filename,opt)
% This routine exports pairs of flu and str to an ASCII plt file
% Usage: flustr2plt_v2(flu,str,filename,opt)
% flu and str are the structures that are going to be exported to the file
% given by "filename". opt selects which variables are written to the file
% Possible options
% For exporting STR
% 'strx': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, ltens]
%
% 'strx0': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the UNdeformed coordinates. 
%       Variables exported: [x0,y0, vx,vy, stress, F, ltens]
%
% 'stradvx': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, FX]
%       where FX=F from eulerian material position
%
% 'strx_xmat': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, ltens, Xmat, Ymat]
%
% For exporting FLU
% 'flu_avgvx_p_avgf': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Variables exported: [xmed,ymed, avgvx, avgvy, p, avgfx, avgfy]
%
% 'flu_avgvx_p_avgf_vs': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Spatial positions are shared with first zone (concatenate ONLY)
%       Variables exported: [avgvx, avgvy, p, avgfx, avgfy]
%
% 'flu_mz_vx_p_fx': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell EDGES
%       Variables are exported in multiple zones without averaging.
%       Variables exported: {[x,y,p] [x,y,vx,fx],[x,y,vy,fy]]
%
% 'flu_mz_advx': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Variables exported: {[xmed,ymed, avgvx, avgvy, p, avgfx, avgfy]
%                             [x,y,xmat],[x,y,ymat]}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [str,flup,strp]                                                   = flustr2plt_initial_configuration(flu,str,filename,opt,append,interpolate)


if nargin < 5
    append                                                                 = 0;
end
if nargin < 6
    interpolate                                                            = 0; % by default, we don't interpolate to nodes.
end
header                                                                     = 'Example: mixed cell-centred and point data on 2D/3D MAC/FE mesh in BLOCK format';
str_export                                                                 = 0;
flu_export                                                                 = 0;
fem_export                                                                 = 0;

if isempty(flu)
    switch str.data.dim
        case 2
            flu.control.z                                                  = 0;
        case {3,13,23}
            flu.control.z                                                  = 1;
    end
    flu.control.outputfilename                                             = str.control.outputfilename;
    flu.control.t                                                          = str.time;
end

if ~isfield(flu.control,'z')
    flu.control.z                                                          = 0;
end

if isnumeric(opt)
    switch opt
        case {1,11,12,7,10,14}            
            switch opt                
                case {1,7,11}
                    options                                                = 'strx';
                case 12
                    options                                                = 'strx0';                
                case 10
                    options                                                = 'stradvx';
                case 14
                    options                                                = 'strx_xmat';
            end
        case {2,6,13,9}            
            switch opt
                case {2,6}
                    options                                                = 'flu_avgvx_p_avgf';                
                case 13
                    options                                                = 'flu_mz_vx_p_fx';
                case 9
                    options                                                = 'flu_mz_advx';
            end
    end
else
    options = opt;
end

name_p = {'p'};
name_residual = {'normR'};

if flu.control.z
    name_pos                                                               = {'X','Y','Z'};
    name_xmat                                                              = {'Xmat','Ymat','Zmat'};
    name_pos0                                                              =  {'X0','Y0','Z0'};
    eltype                                                                 = 'TETRAHEDRON';
else
    name_pos                                                               = {'X','Y'};
    name_xmat                                                              = {'Xmat','Ymat'};
    name_pos0                                                              =  {'X0','Y0'};
    eltype                                                                 = 'TETRAHEDRON';    eltype                                                                 = 'TRIANGLE';
end

passive                                                                    = {[]};
var_sharing                                                                = {[]};
con_sharing                                                                = {[]};

numzones                                                                   = 1;
size_var                                                                   = {@extract_pressure};

 
nodes                                                              = @(s) s.postproc.Lagrangian_X;
nodes0                                                             = @(s) s.postproc.Lagrangian_X;
eldata                                                             = @(s) s.postproc.connectivity;
[str]  =  plot_variables_treatment(str);
%----------------------------------
% Data. 
%----------------------------------
data                                                               = {{nodes,nodes0,eldata}}; 


%----------------------------------
% Names of the variables associated to data.
%----------------------------------
type                                                               = ones(numel(data{1}),1);
type(end)                                                          = 0;
variables                                                          = [name_pos name_pos0];
fem_export                                                          = 1;
switch str.data.dim
    case 2
        switch str.data.shape
            case 1 % quadrilateral
                eltype                                             = 'QUADRILATERAL';
            case 0 % triangular
                eltype                                             = 'TRIANGLE';
        end
    case {3,13,23}
        switch str.data.shape
            case 1 % hexahedral
                eltype                                             = 'BRICK';
            case 0 % tetrahedral
                eltype                                             = 'TETRAHEDRON';
        end
        
end
datatype = 'FEBLOCK';

        


loc                                                                        = 'NODAL';
loc                                                                        = 'CELLCENTERED';

if str_export 
switch str.kindelem.eltype
    case {'3dtetra4_nodal','2dtria3nodal_opt','2dtria3nodal_gq','3dtetra4_nodal_gq'}
        ifem                                                               = 0; % I need to finish detecting the type
    otherwise
        ifem                                                               = 1;
end
else
    ifem                                                                   = 0;
end

% Variable location depends on the method. For the ISPM all str-exported
% variables are nodal, whereas for the ifem node positions and velocities
% are nodal, but the rest of the variables are cell-centered.
% For mac flus variables are exported as "nodal", even though strictly
% speaking they are not.
if (ifem && ~flu_export)
    location                                                               = cellfun(@(t) 'CELLCENTERED', variables,'UniformOutput',false);
    mult                                                                   = 2;
    switch options
        case {'strx_x0_eigs','strx_x0_eigs_cs','strx_x0_eigs_csvs'}
            mult                                                           = 3;
    end
    for i = 1:(mult*size(str.x,1))
        location{i}                                                        = 'NODAL';
    end
    switch options
        case 'strx_x0_eigs_vinc'
            location{end}                                                  = 'NODAL';
    end
else
    location                                                               = cellfun(@(t) 'NODAL', variables,'UniformOutput',false);
%     switch options
%         case 'dgfem_u_residual'
%             switch str.d
%                 case 2
%                     location{5}                                          = 'CELLCENTERED';    
%                 case 3
%                     location{7}                                          = 'CELLCENTERED';    
%             end
%     end
end

if isunix()
    dirsep                                                                 = '/';
else
    dirsep                                                                 = '\';
end
a=find(flu.control.outputfilename==dirsep);
if isunix()
    zonenumber                                                             = flu.control.outputfilename((a(end-2)+1):(a(end)-1));
    zonenumber                                                             = strrep(zonenumber,dirsep,'\\');
else
    if ~isempty(a)
        zonenumber                                                         = flu.control.outputfilename(1:(a(end)-1));
        zonenumber                                                         = strrep(zonenumber,dirsep,'\\');
    else
        zonenumber                                                         = 'not automatically assigned';
    end
end



switch str.data.example_comparison
       case 0
            jobfolder                                                      = [str.jobfolder '\results'];
            cd(jobfolder)
       case 1            
            cd(str.destination_jobfolder)
end

%datatype = 'FEPOINT';
% Open file for writing

if append
    f                                                                      = fopen(filename,'a+');
else
    f                                                                      = fopen(filename,'w+');
end

% Write header
writeheader(f,header);

% Write variable header
writevariableheader(f,variables);

if fem_export
     for zon = 1:numzones
        % Write zone header
        writezoneheader_t(f,zonenumber, size(str.postproc.nodes,2),size(str.postproc.t,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 

        % Write variable location
        writevarlocation(f,location);
        
        % Now get the data ready to output
    % STR ISPM
        % In the case of a structure
%         if ispm_str
%             [flup,strp]                                                  = peek_structure(flu,str,15);
%         end
        % write node positions and velocities; always in block format and
        % nodal
        
        % Now get the data ready to output
        for i = 1:numel(data{zon})
            switch type(i)
                case 1
                    formato                                                = '%.16g';
                case 0
                    formato                                                = '%d';
            end
%             switch use_peeked(i)
%                 case 1
%                     s                                                    = strp;
%                 case 0
                    s                                                      = str;
%             end
            writepointdata_faster(f,data{zon}{i}(s),formato);
        end 
        

     end
end

if str_export
     for zon = 1:numzones
             % STR ISPM
        % In the case of a structure
%         if ispm_str
            [flup,strp]                                                    = peek_structure(flu,str,15);
%         end
        % write node positions and velocities; always in block format and
        % nodal
        if interpolate
            if sum(use_peeked) > 0
                str_interp                                                 = resample_str_to_nodes(strp);
            else
                error('Cannot interpolate not peeking the structure!!')
            end
            writezoneheader_t(f,zonenumber, size(str_interp.kindelem.x,2),size(str_interp.kindelem.lnods,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 
        else
            writezoneheader_t(f,zonenumber, size(str.kindelem.x,2),size(str.kindelem.lnods,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 
        end
        
        % Write zone header
        
        % Write variable location
        writevarlocation(f,location);
        
        % Now get the data ready to output

        
        % Now get the data ready to output
        for i = 1:numel(data{zon})
            switch type(i)
                case 1
                    formato                                                = '%.16g';
                case 0
                    formato                                                = '%d';
            end
            switch use_peeked(i)
                case 1
                    s                                                      = strp;
                case 0
                    s                                                      = str;
            end
            if interpolate
                s                                                          = str_interp;
            end
            writepointdata_faster(f,data{zon}{i}(s),formato);
        end
        
% %         blockdata = [nodes(str); vel(str)];
% %         % write data
% %         writepointdata_faster(f,blockdata);
% %         % write the other variables; always in block format. nodal or
% %         % cell-centered depending on ispm/ifem
% %         blockdata = [S(strp) ; ...
% %                     F(strp);...
% %                     L(strp)];                
% %         writepointdata_faster(f,blockdata);
% %         % get element data ready
% %         elementdata = (str.kindelem.lnods');
% %         % write element data
% %         writeelementdata_faster(f,elementdata);
% %     % STR IFEM
     end
end

if flu_export
    for zon = 1:numzones
        P                                                                  = size_var{zon}(flu);
        if flu.control.z
            % mac flu 3d
            writezoneheader_ijkt(f,zonenumber, size(P,1),size(P,2), size(P,3),flu.control.t,datatype, passive{zon},var_sharing{zon});
        else
            % mac flu 2d
            writezoneheader_ijt(f,zonenumber, size(P,1),size(P,2),flu.control.t,datatype, passive{zon},var_sharing{zon});
        end 
        % Write variable location
        writevarlocation(f,location);

        % Now get the data ready to output
        for i = 1:numel(data{zon})
            writepointdata_faster(f,data{zon}{i}(flu));
        end
        flup                                                               = flu;
        strp                                                               = str;
    end
end
            
% close file
fclose(f);

end


function writeheader(f,s)
headline                                                                   = ['TITLE = "' s '"'];
fprintf(f,'%s\n',headline);
end
function writevariableheader(f,s)

a                                                                          = cell2mat(cellfun(@(t) ['"' t '", '],s,'UniformOutput',false));
variables                                                                  = ['VARIABLES = ' a(1:(end-2))];
fprintf(f,'%s\n',variables);
end
function writezoneheader(f,zonename,numpoints,numel,type,eltype,varargin)
if nargin > 6
    %passive = 1;  
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive = 0;
    pass                                                                   = '';
end

zones                                                                      = ['ZONE T ="' zonename '", N=' num2str(numpoints) ', E=' num2str(numel) ', F=' type ', ET=' eltype pass ];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ij(f,zonename,numpoints,numel,type,varargin)
if nargin > 6
    %passive                                                               = 1;
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(numpoints) ', J=' num2str(numel) ', F=' type pass];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijt(f,zonename,numpoints,numel,t,type,varargin)
if nargin > 6
    %passive                                                               = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 7
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(numpoints) ', J=' num2str(numel) ', F=' type ', SOLUTIONTIME=' num2str(t,20) pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijk(f,zonename,i,j,k,type, varargin)
if nargin > 6
    %passive = 1;
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 7
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end

zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(i) ', J=' num2str(j) ', K=' num2str(k) ', F=' type pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijkt(f,zonename,i,j,k,t,type, varargin)
if nargin > 7
    %passive = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 8
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(i) ', J=' num2str(j) ', K=' num2str(k) ', F=' type ', SOLUTIONTIME=' num2str(t,20) pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_t(f,zonename,numpoints,numel,t,type,eltype,varargin)
if nargin > 7
    %passive                                                               = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 8
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
if nargin > 9
    if ~isempty(varargin{3})
        cs                                                                 = [', CONNECTIVITYSHAREZONE=' varargin{3}];
    else
        cs                                                                 = '';
    end
else    
    cs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", N=' num2str(numpoints) ', E=' num2str(numel) ', F=' type ', ET=' eltype ', SOLUTIONTIME=' num2str(t,20) pass vs cs];
fprintf(f,'%s\n',zones);
end
function writepointdata(f,elementdata)
% save(filename,'-ascii','-append','-double','pointdata');
c=1;
for i=1:size(elementdata,1)
    for j=1:size(elementdata,2)
        fprintf(f,'%g ',elementdata(i,j));
        c = c+1;
        if mod(c,50)==1 % This is just to avoid extremely long lines and tecplot crashing with them
            fprintf(f,'\n');
        end
    end
    fprintf(f,'\n');
end
end
function writepointdata_fast(filename,elementdata)
save(filename,'-ascii','-append','-double','elementdata');
end
function writeelementdata(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% fprintf(f,'\n');
for i=1:size(elementdata,1)
    for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i,j));
    end
    fprintf(f,'\n');
end
end
function writepointdata_faster(f,elementdata,varargin)
% save(filename,'-ascii','-append','elementdata');
% % fprintf(f,'\n');
if nargin > 2
    formato                                                                = varargin{1};
else
    formato                                                                = '%.16g';
end
stride                                                                     = 100;
elementdata                                                                = elementdata';
 j                                                                         = 0:(stride-1);
 if numel(elementdata) > stride
for i=1:stride:(numel(elementdata)-stride)
%     for j=1:size(elementdata,2)
        fprintf(f,[formato ' '],elementdata(i+j));
%     end
    fprintf(f,'\n');
end
fprintf(f,[formato ' '],elementdata(i+stride:end));
 else
     fprintf(f,[formato ' '],elementdata);
 end
% fprintf(f,'%d ',elementdata'); % this is really fast, but gives long lines that fail with tecplot
fprintf(f,'\n\n');
end
function writeelementdata_fast(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% fprintf(f,'\n');
for i=1:size(elementdata,1)
%     for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i,:));
%     end
    fprintf(f,'\n');
end
end
function writeelementdata_faster(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% % fprintf(f,'\n');
stride                                                                     = 100;
elementdata                                                                = elementdata';
 j = 0:(stride-1);
 if numel(elementdata) > stride
for i=1:stride:(numel(elementdata)-stride)
%     for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i+j));
%     end
    fprintf(f,'\n');
end
fprintf(f,'%d ',elementdata(i+stride:end));
 else
     fprintf(f,'%d ',elementdata);
 end
% fprintf(f,'%d ',elementdata'); % this is really fast, but gives long lines that fail with tecplot
fprintf(f,'\n\n');
end
function writedata(filename,elementdata)
save(filename,'-ascii','-append','elementdata');
end
function writevarlocation(f,s)
a                                                                          = cell2mat(cellfun(@(t) [t ', '],s,'UniformOutput',false));
variables                                                                  = ['VARLOCATION = (' a(1:(end-2)) ')'];
fprintf(f,'%s\n',variables);
end
function pointdata                                                         = pointformat(c)
pointdata                                                                  = cell2mat(c);
end
function T                                                                 = eT(tensor)
    if size(tensor,1) == 2
        T                                                                  = [squeeze(tensor(1,1,:))'; ...
                                                                              squeeze(tensor(1,2,:))'; ...
                                                                              squeeze(tensor(2,1,:))'; ...
                                                                              squeeze(tensor(2,2,:))'];
    elseif size(tensor,1) == 3
        T                                                                  = [squeeze(tensor(1,1,:))'; ...
                                                                              squeeze(tensor(1,2,:))'; ...
                                                                              squeeze(tensor(1,3,:))'; ...
                                                                              squeeze(tensor(2,1,:))'; ...
                                                                              squeeze(tensor(2,2,:))'; ...
                                                                              squeeze(tensor(2,3,:))'; ...
                                                                              squeeze(tensor(3,1,:))'; ...
                                                                              squeeze(tensor(3,2,:))'; ...
                                                                              squeeze(tensor(3,3,:))'];
    end
end
function v                                                                 = ext_fluid_vel(f)
    if f.control.z
        vx                                                                 = f.vx(f.index1,f.pindex2,f.pindex3);
        vy                                                                 = f.vy(f.pindex1,f.index2,f.pindex3);
        vz                                                                 = f.vz(f.pindex1,f.pindex2,f.index3);
        v                                                                  = {vx;vy;vz};
    else
        vx                                                                 = f.vx(f.index1,f.pindex2);
        vy                                                                 = f.vy(f.pindex1,f.index2);
        v                                                                  = {vx;vy};
    end
end
function v= ext_forces(f)
    if f.control.z
        vx                                                                 = f.fx(f.index1,f.pindex2,f.pindex3);
        vy                                                                 = f.fy(f.pindex1,f.index2,f.pindex3);
        vz                                                                 = f.fz(f.pindex1,f.pindex2,f.index3);
        v                                                                  = {vx;vy;vz};
    else                                                                   
        vx                                                                 = f.fx(f.index1,f.pindex2);
        vy                                                                 = f.fy(f.pindex1,f.index2);
        v                                                                  = {vx;vy};
    end   
end
function p                                                                 = extract_pressure(f)
     if f.control.z
         p                                                                 = f.p(f.pindex1,f.pindex2,f.pindex3);
     else
         p                                                                 = f.p(f.pindex1,f.pindex2);
     end
     %p={p};
end
function po = mctmat(p)
    switch numel(p)
        case 3
            avx                                                            = p{1}; avy = p{2}; avz = p{3};
            po                                                             = [avx(:) avy(:) avz(:)]';
        case 2
            avx                                                            = p{1}; avy = p{2};
            po                                                             = [avx(:) avy(:)]';
        case 1
            avx                                                            = p{1};
            po                                                             = [avx(:)]';
    end
end
function vo                                                                = avg(v)
    if numel(v)==3
        vx  = v{1}; vy = v{2}; vz = v{3};
        avx                                                                = 0.5*(vx(1:end-1,:,:)+vx(2:end,:,:));
        avy                                                                = 0.5*(vy(:,1:end-1,:)+vy(:,2:end,:));
        avz                                                                = 0.5*(vz(:,:,1:end-1)+vz(:,:,2:end));
        vo                                                                 = {avx;avy;avz};        
    else
        vx                                                                 = v{1}; vy = v{2};% vz = v{3};
        avx                                                                = 0.5*(vx(1:end-1,:)+vx(2:end,:));
        avy                                                                = 0.5*(vy(:,1:end-1)+vy(:,2:end)); 
        vo                                                                 = {avx;avy};
    end    
end
function xo                                                                = ext_fluid_midgrid(f)
    if f.control.z
        % CORRECT THIS IF NECESSARY FOR 3D
        % uncomment if there is any problem
%         xmed=((f.x(1:end-1,1:end-1,1:end-1)+f.x(2:end,1:end-1,1:end-1))/2);
%         ymed=((f.y(1:end-1,1:end-1,1:end-1)+f.y(1:end-1,2:end,1:end-1))/2);
%         zmed=((f.z(1:end-1,1:end-1,1:end-1)+f.z(1:end-1,1:end-1,2:end))/2);
%         xo = [xmed(:) ymed(:) zmed(:)]';
        xo                                                                 = [f.xmed(:) f.ymed(:) f.zmed(:)]';
    else        
        xmed                                                               = ((f.x(1:end-1,1:end-1)+f.x(2:end,1:end-1))/2);
        ymed                                                               = ((f.y(1:end-1,1:end-1)+f.y(1:end-1,2:end))/2);
        xo                                                                 = [xmed(:) ymed(:)]';        
    end
end
function xo                                                                = ext_fluid_xgrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                  (f.xmin(1)-f.deltax*f.nxghostcells):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells),...
                                                                                  (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay),...
                                                                                  (f.xmin(3)-f.deltaz*f.nzghostcells+0.5*f.deltaz):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells-0.5*f.deltaz));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        [xr,yr]                                                            = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay));
        xo                                                                 = [xr(:) yr(:)]';                
    end
end
function xo                                                                = ext_fluid_ygrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells),...
                                                                                    (f.xmin(3)-f.deltaz*f.nzghostcells+0.5*f.deltaz):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells-0.5*f.deltaz));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        [xr,yr]                                                            = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells));
        xo = [xr(:) yr(:)]';                
    end
end
function xo                                                                = ext_fluid_zgrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay),...
                                                                                    (f.xmin(3)-f.deltaz*f.nzghostcells):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        xo                                                                 = [];                
    end
end

function e                                                                 = comp_eigenvals_and_eigenvects(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        e                                                                  = zeros(6,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            [V,D]                                                          = eig(A);
            d                                                              = diag(D);
            [e(1:2,i),p]                                                   = sort(d);
            V                                                              = V(:,p);
            e(3:6,i)                                                       = V(:);             
        end
    else
        e                                                                  = zeros(12,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i);... 
                                                                              S(4,i) S(2,i) S(6,i);...
                                                                              S(5,i) S(6,i) S(3,i)];
            [V,D]                                                          = eig(A);
            d                                                              = diag(D);
            [e(1:3,i),p]                                                   = sort(d);
            V                                                              = V(:,p);
            e(4:12,i)                                                      = V(:);
        end
    end
end

function e                                                                 = comp_eigenvals(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        e                                                                  = zeros(2,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    else
        e                                                                  = zeros(3,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i); S(4,i) S(2,i) S(6,i); S(5,i) S(6,i) S(3,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    end
end
function v                                                                 = comp_eigenvects(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        v                                                                  = zeros(4,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            [e,vv]                                                         = sort(eig(A));
        end
    else
        v                                                                  = zeros(3,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i); S(4,i) S(2,i) S(6,i); S(5,i) S(6,i) S(3,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    end
end
