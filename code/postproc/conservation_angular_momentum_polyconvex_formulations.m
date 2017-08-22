function   [integrated_Kirchhoff,integrated_Kirchhoff_element,momentum_vanishing,check1,check2,check3,check4]    =  conservation_angular_momentum_polyconvex_formulations(str)

integrated_Kirchhoff               =  zeros(3,3);
integrated_Kirchhoff_element       =  zeros(3,3,str.n_elem);
momentum_vanishing              =  zeros(3,1);
check1                          =  zeros(3,3);
check2                          =  zeros(3,3);
check3                          =  0;
check4                          =  zeros(3,3,str.n_nodes);
for ielem=1:str.n_elem
    str.ielem                   =  ielem;
    %------------------------------------------------------------
    % Compute information at gauss level. 
    %------------------------------------------------------------
    nodes                       =  str.connectivity(ielem,:);
    xelem                       =  str.Eulerian_x(:,nodes);
    Xelem                       =  str.Lagrangian_X(:,nodes);
    phielem                     =  str.phi(nodes,1);
    str                         =  gradients(xelem,Xelem,phielem,str);
    gauss_level_information     =  gauss_level_information_mixed_formulations(str);
    Piola                       =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);    
    for igauss=1:size(str.grad.F,3)
        F                       =  str.grad.F(:,:,igauss);
        J                       =  str.grad.J(igauss);
        %----------------------------------------------------------------------
        % Isoparametric information 
        %----------------------------------------------------------------------
        J_t                     =  abs(det(str.grad.DX_chi(:,:,igauss)));
        W                       =  str.quadrature.W_v(igauss);        
        integrated_Kirchhoff       =  integrated_Kirchhoff + Piola(:,:,igauss)*F'*(W*J_t);
        integrated_Kirchhoff_element(:,:,...
            ielem)              =  integrated_Kirchhoff_element(:,:,ielem) + Piola(:,:,igauss)*F'*(W*J_t);
        momentum_vanishing      =  momentum_vanishing + levi_civita(F*Piola(:,:,igauss)',3)*(W*J_t);
        check1                  =  check1 + str.properties.stabilisation_parameter_TJ(1)*(F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*F'*(W*J_t);
        %check1                  =  check1 + str.properties.stabilisation_parameter_TJ(1)*(F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*eye(3)'*(W*J_t);
        check2                  =  check2 + (F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*gauss_level_information.SigmaF(:,:,igauss)'*(W*J_t);

        
        
        check3                  =  check3 + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*F')*(W*J_t);
        for inode=1:4
            NF                  =  str.f_e.N_F(inode,igauss);
            SigmaF              =  str.SigmaF(:,:,nodes(inode));
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*NF)*(W*J_t);
            check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + ((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*NF*SigmaF)*(W*J_t);
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*NF*eye(3))*(W*J_t);
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*eye(3))*(W*J_t);
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*(F'-reshape(gauss_level_information.DGDSF(:,igauss),3,[])))*(W*J_t);
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*F')*(W*J_t);
            %check4(:,:,nodes(inode))   =  check4(:,:,nodes(inode)) + str.properties.stabilisation_parameter_TJ(1)*((F - reshape(gauss_level_information.DGDSF(:,igauss),3,[])')*(-reshape(gauss_level_information.DGDSF(:,igauss),3,[])))*(W*J_t);
        end
    end
end
