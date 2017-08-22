function    str               =  additional_constraints(str)




% n_node1                            =  [];
% for i=1:size(str.nodes,1)
%     if norm(str.nodes(i,3) - 0.25)<1e-4
%        n_node1                     =  [n_node1;i];
%     end
% end
% n_node2                            =  [];
% for i=1:size(str.nodes,1)
%     if norm(str.nodes(i,3) - 0.0)<1e-4
%        n_node2                     =  [n_node2;i];
%     end
% end
% prescribed_potential                =  1e7;
% 
% 
% str.elec.freedof_prescribed         =  n_node2;
% str.elec.freedof_prescribed_value   =  prescribed_potential*ones(size(str.elec.freedof_prescribed,1),1);
% str.mec.freedof_prescribed          =  [];
% str.mec.freedof_prescribed_value    =  [];
 
%--------------------------------------------------------------------------
% Dirichlet type bc's
%--------------------------------------------------------------------------
% str.fixdof_arclength                =  [str.mec.fixdof; 3*str.n_nodes + str.elec.fixdof;3*str.n_nodes+n_node1];
% switch str.data.formulation.mixed_type
%     case {'u_p_phi_incompressible','full_mixed_formulation_electroelasticity_F_D0_V_incompressible'}
%          str.freedof_arclength               =  (1:str.n_nodes*4+size(str.nodes_p,1))';
%     case {'u_phi_D0_c_formulation'}
%          str.freedof_arclength               =  (1:str.n_nodes*4+3*size(str.nodes_E0,1))';
%     case {'full_mixed_formulation_electroelasticity_F_D0_V_c'}
%          str.freedof_arclength               =  (1:str.n_nodes*4 + 2*9*size(str.nodes_F,1) + 2*9*size(str.nodes_H,1) + ...
%                                                                    2*size(str.nodes_J,1) + 3*size(str.nodes_E0,1) + 2*3*size(str.nodes_V,1))';
%     case {'full_mixed_formulation_electroelasticity_F_D0_V_p_c'}
%          str.freedof_arclength               =  (1:str.n_nodes*4 + 2*9*size(str.nodes_F,1) + 2*9*size(str.nodes_H,1) + ...
%                                                                    2*size(str.nodes_J,1) + 3*size(str.nodes_E0,1) + 2*3*size(str.nodes_V,1) + size(str.nodes_p,1))';
%     otherwise
%          str.freedof_arclength               =  (1:str.n_nodes*4)';
% end
% str.freedof_arclength(str.fixdof_arclength,...
%     :)                              =  []; 
%str.freedof_prescribed              =  3*str.n_nodes + n_node2;

 
% str.elec_freedof_arclength          =  (1:str.n_nodes)';
% freedof_arclength                   =  str.freedof_arclength(str.freedof_arclength<=str.n_nodes*4);
% str.elec_freedof_arclength          =  freedof_arclength(freedof_arclength>str.n_nodes*3) - 3*str.n_nodes;

% elec_fixdof                         =  [n_node1;n_node2];
% str.elec.fixdof                     =  [str.elec.fixdof;elec_fixdof];
% str.elec.cons_val                   =  [str.elec.cons_val;zeros(size(elec_fixdof,1)/2,1);  str.elec.freedof_prescribed_value];
% 
% str.elec.freedof_prescribed         =  n_node2;
% str.elec.freedof_prescribed_value   =  prescribed_potential*ones(size(str.elec.freedof_prescribed,1),1);
% str.mec.freedof_prescribed          =  [];
% str.mec.freedof_prescribed_value    =  [];
 
%--------------------------------------------------------------------------
% Newmann type bc's
%--------------------------------------------------------------------------
mec_node_1 =  [];
for i=1:size(str.nodes,1)
    if norm(str.nodes(i,1) - 0)<1e-4
       if norm(str.nodes(i,3) - 0)<1e-4
          mec_node_1               =  [mec_node_1;i];
       end
    end
%     if norm(str.nodes(i,1) - 10)<1e-4
%        if norm(str.nodes(i,3) - 0.025)<1e-4
%           mec_node_1               =  [mec_node_1;i];
%        end
%     end
%     if norm(str.nodes(i,2) - 0)<1e-4
%        if norm(str.nodes(i,3) - 0.025)<1e-4
%           mec_node_1               =  [mec_node_1;i];
%        end
%     end
%     if norm(str.nodes(i,2) - 1)<1e-4
%        if norm(str.nodes(i,3) - 0.025)<1e-4
%           mec_node_1               =  [mec_node_1;i];
%        end
%     end
end
str.mec.fixdof                     =  sort(unique([str.mec.fixdof;3*mec_node_1]));
str.mec.cons_val                   =  [str.mec.cons_val;zeros(size(unique(mec_node_1),1),1);];
% % 


% 
n_node1                            =  [];
for i=1:size(str.nodes,1)
    if norm(str.nodes(i,3) - 0.025)<1e-4
       n_node1                     =  [n_node1;i];
    end
end
% n_node2                            =  [];
% for i=1:size(str.nodes,1)
%     if norm(str.nodes(i,3) - 0.0)<1e-4
%        n_node2                     =  [n_node2;i];
%     end
% end
 prescribed_potential                =  0;


% str.elec.freedof_prescribed         =  n_node2;
% str.elec.freedof_prescribed_value   =  prescribed_potential*ones(size(str.elec.freedof_prescribed,1),1);
% str.mec.freedof_prescribed          =  [];
% str.mec.freedof_prescribed_value    =  [];
 

elec_fixdof                         =  n_node1;
str.elec.fixdof                     =  [str.elec.fixdof;elec_fixdof];
str.elec.cons_val                   =  [str.elec.cons_val;prescribed_potential*ones(size(elec_fixdof,1),1);];
% 
