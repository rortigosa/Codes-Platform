%--------------------------------------------------------------------------
% Division of the boundary for contact problems
%--------------------------------------------------------------------------
function str                              =  SpecificBoundaryDivision(str)

zmin                                      =  min(str.mesh.volume.x.nodes(3,:));
zmax                                      =  max(str.mesh.volume.x.nodes(3,:));
xmin                                      =  min(str.mesh.volume.x.nodes(1,:));
xmax                                      =  max(str.mesh.volume.x.nodes(1,:));
xmed                                      =  0.5*(xmin + xmax);
n_elements_boundary                       =  size(str.mesh.surface.x.boundary_edges,2);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Identify number of edges within each region
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
areaA_entries                             =  0;
areaB_entries                             =  0;
areaC_entries                             =  0;
areaD_entries                             =  0;
areaE_entries                             =  0;
areaF_entries                             =  0;
for iedge=1:n_elements_boundary
    xelem                                 =  (str.mesh.volume.x.nodes(:,str.mesh.surface.x.boundary_edges(:,iedge)))';
    xcenter                               =  (sum(xelem)/size(xelem,1))';
    %----------------------------------------------------------------------
    % Region A
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmax)<1e-5
       if xcenter(1)<xmed
          areaA_entries                   =  areaA_entries + 1;
       end
    end
    %----------------------------------------------------------------------
    % Region B
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmax)<1e-5
       if xcenter(1)>xmed
          areaB_entries                   =  areaB_entries + 1;          
       end
    end    
    %----------------------------------------------------------------------
    % Region C
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmin)<1e-5
       if xcenter(1)<xmed
          areaC_entries                   =  areaC_entries + 1;
       end
    end    
    %----------------------------------------------------------------------
    % Regiobn D
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmin)<1e-5
       if xcenter(1)>xmed
          areaD_entries                   =  areaD_entries + 1;          
       end
    end    
    %----------------------------------------------------------------------
    % Region E
    %----------------------------------------------------------------------
    if abs(xcenter(1) - xmin)<1e-5
          areaE_entries                   =  areaE_entries + 1;
    end
    %----------------------------------------------------------------------
    % Region F
    %----------------------------------------------------------------------
    if abs(xcenter(1) - xmax)<1e-5
          areaF_entries                   =  areaF_entries + 1;
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Specify the different regions
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
areaA                                     =  zeros(areaA_entries,1);
areaB                                     =  zeros(areaB_entries,1);
areaC                                     =  zeros(areaC_entries,1);
areaD                                     =  zeros(areaD_entries,1);
areaE                                     =  zeros(areaE_entries,1);
areaF                                     =  zeros(areaF_entries,1);
areaA_entrie                              =  0;
areaB_entrie                              =  0;
areaC_entrie                              =  0;
areaD_entrie                              =  0;
areaE_entrie                              =  0;
areaF_entrie                              =  0;
for iedge=1:n_elements_boundary
    xelem                                 =  (str.mesh.volume.x.nodes(:,str.mesh.surface.x.boundary_edges(:,iedge)))';
    xcenter                               =  (sum(xelem)/size(xelem,1))';
    %----------------------------------------------------------------------
    % Region A
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmax)<1e-5
       if xcenter(1)<xmed
          areaA_entrie                    =  areaA_entrie + 1;
          areaA(areaA_entrie)             =  iedge;
       end
    end
    %----------------------------------------------------------------------
    % Region B
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmax)<1e-5
       if xcenter(1)>xmed
          areaB_entrie                    =  areaB_entrie + 1;
          areaB(areaB_entrie)             =  iedge;
       end
    end    
    %----------------------------------------------------------------------
    % Region C
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmin)<1e-5
       if xcenter(1)<xmed
          areaC_entrie                    =  areaC_entrie + 1;
          areaC(areaC_entrie)             =  iedge;
       end
    end    
    %----------------------------------------------------------------------
    % Regiobn D
    %----------------------------------------------------------------------
    if abs(xcenter(3) - zmin)<1e-5
       if xcenter(1)>xmed
          areaD_entrie                    =  areaD_entrie + 1;
          areaD(areaD_entrie)             =  iedge;
       end
    end    
    %----------------------------------------------------------------------
    % Region E
    %----------------------------------------------------------------------
    if abs(xcenter(1) - xmin)<1e-5
          areaE_entrie                    =  areaE_entrie + 1;
          areaE(areaE_entrie)             =  iedge;
    end
    %----------------------------------------------------------------------
    % Region F
    %----------------------------------------------------------------------
    if abs(xcenter(1) - xmax)<1e-5
          areaF_entrie                    =  areaF_entrie + 1;
          areaF(areaF_entrie)             =  iedge;
    end
end

str.mesh.surface.x.boundary_division.edges{1}  =  areaA;
str.mesh.surface.x.boundary_division.edges{2}  =  areaB;
str.mesh.surface.x.boundary_division.edges{3}  =  areaC;
str.mesh.surface.x.boundary_division.edges{4}  =  areaD;
str.mesh.surface.x.boundary_division.edges{5}  =  areaE;
str.mesh.surface.x.boundary_division.edges{6}  =  areaF;





