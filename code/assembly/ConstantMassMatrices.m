
function assembly                      =  ConstantMassMatrices(geom,mesh,fem,quadrature)

%--------------------------------------------------------------------------
%  Mx0x0
%--------------------------------------------------------------------------
dim                                    =  geom.dim;
Mass_xx                                =  zeros(dim*mesh.volume.x.n_node_elem,...
                                                dim*mesh.volume.x.n_node_elem,...
                                                size(quadrature.volume.mass.Chi,1));
for igauss=1:size(quadrature.volume.mass.Chi,1)
    for anode=1:mesh.volume.x.n_node_elem
        for bnode=1:mesh.volume.x.n_node_elem
            adofs                      =  dim*(anode-1)+1:dim*anode;
            bdofs                      =  dim*(bnode-1)+1:dim*bnode;
            Mass_xx(adofs,...
      bdofs,igauss)                    =  Mass_xx(adofs,bdofs,igauss)  +  (fem.volume.mass.x.N(anode,igauss)*...
                                                                             fem.volume.mass.x.N(bnode,igauss))*eye(dim);
        end
    end
end 
assembly.mass_matrices.constant_Mxx    =  Mass_xx;

