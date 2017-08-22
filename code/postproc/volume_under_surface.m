clear u volume pressure factor VA PA
n_loads                    =  53;

for ifile=1:n_loads
    ifile
filename  =  ['Load_increment_' num2str(ifile) '.mat'];
load(filename)
str=new_str;
u(ifile)  =  (max(abs(str.Eulerian_x(3,:)-str.Lagrangian_X(3,:))));
factor(ifile)  =  str.data.acumulated_factor;
end 
% 
% 

plot(u,factor,'-o')

figure(1)
% % 
% n_loads                    =  34;
% volume                     =  zeros(n_loads,1);
% pressure                     =  zeros(n_loads,1);
% for iload=1:n_loads    
%     filename               =  ['Load_increment_' num2str(iload) '.mat'];
%     load(filename)
%     str                    =  new_str;
%     pressure(iload)        =  str.data.acumulated_factor;
%     volume(iload)          =  0;
%     for ielem=1:str.n_elem
%         nodes              =  str.connectivity(ielem,:);
%         xelem              =  str.Eulerian_x(:,nodes);
%         Xelem              =  str.Lagrangian_X(:,nodes);
%         phielem            =  str.phi(nodes);
%         str                =  gradients(xelem,Xelem,phielem,str);
%         for igauss=1:size(str.quadrature.Chi,1)
%             DX_chi         =  str.grad.DX_chi(:,:,igauss);
%             J_t            =  abs(det(DX_chi));
%             W              =  str.quadrature.W_v(igauss);
%             volume(iload)  =  volume(iload) + str.grad.J(igauss)*W*J_t;
%         end
%     end
% end
% 
% 
% plot(volume,pressure)
% 
% 
% 


volume                     =  zeros(n_loads,1);
area                     =  zeros(n_loads,1);
VA                       =  zeros(n_loads,1);
PA                       =  zeros(n_loads,1);
ppressure                     =  zeros(n_loads,1);
for iload=1:n_loads    
    iload
    filename               =  ['Load_increment_' num2str(iload) '.mat'];
    load(filename)
    str                    =  new_str;
    pressure(iload)        =  str.data.acumulated_factor;
    volume(iload)          =  0;
    area(iload)             =  0;
    for ielem=1:str.n_elem
        if str.solid.bc.pressure_load.active_element(ielem)
        local_nodes                                              =  str.solid.bc.pressure_load.local_nodes(ielem,:);
        global_nodes                                             =  str.connectivity(ielem,local_nodes);
        Xelem                                                    =  str.Lagrangian_X(:,global_nodes);
        xelem                                                    =  str.Eulerian_x(:,global_nodes);
        Total_nodes                                              =  str.connectivity(str.ielem,:);
        X0                                                       =  sum(str.nodes(Total_nodes,:))'/size(Total_nodes,2);

        for igauss=1:size(str.quadrature.surface.Chi,1)
            %----------------------------------------------------------------------
            % Shape functions in the surface
            %----------------------------------------------------------------------
            N                                                    =  str.f_e.surface.N(:,igauss);
            DN_chi                                               =  str.f_e.surface.DN_chi(:,:,igauss);
            x                                                    =  xelem*N;
            X                                                    =  Xelem*N;
            R                                                    =  X - X0;
            %----------------------------------------------------------------------
            % DX_chi, Dx_chi and DN_X
            %----------------------------------------------------------------------
            DX_chi                                               =  Xelem*DN_chi';
            Dx_chi                                               =  xelem*DN_chi';
            DN_X                                                 =  DX_chi'\DN_chi;
            %----------------------------------------------------------------------
            % Normal vector at the Gauss point in the surface
            %----------------------------------------------------------------------
            N_vector                                             =  cross(DX_chi(:,1),DX_chi(:,2));
            N_vector                                             =  N_vector/norm(N_vector);
            N_vector                                             =  N_vector*sign(N_vector'*R);
            n_vector                                             =  cross(Dx_chi(:,1),Dx_chi(:,2));
            n_vector                                             =  n_vector/norm(n_vector);
            DX_chi                                               =  [DX_chi N_vector];
            J_t                                                  =  abs(det(DX_chi));
            W                                                    =  str.quadrature.surface.W_v(igauss);
            %----------------------------------------------------------------------
            % Compute gradients of the deformation in the surface
            %----------------------------------------------------------------------
            F                                                    =  xelem*DN_X' + n_vector*N_vector';
            H                                                    =  0.5*Javier_double_cross_product(F,F,1,1,str.data.dim);
            HN                                                   =  H*N_vector;
            %----------------------------------------------------------------------
            % Compute volume
            %----------------------------------------------------------------------
            volume(iload)                                        =  volume(iload) - 1/3*(x'*(H*N_vector))*(W*J_t);
            area(iload)                                          =  area(iload) + norm(HN)*W*J_t;
        end
        end
    end
    VA(iload)                                                    =  volume(iload)/area(iload);
    PA(iload)                                                     =  pressure(iload)*area(iload);
end
 
volume =  [0;volume];
VA =  [0;VA];
PA =  [0;PA];
factor=[0;factor'];
pressure  =  [0;pressure'];
figure(2)
plot(volume,pressure,'-o')        
figure(3)
plot(VA,factor,'-o')        
figure(4)
plot(VA,PA,'-o')        
figure(5)
plot(volume,factor,'-o')        
        


asdf=98
