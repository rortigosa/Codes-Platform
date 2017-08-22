function str                              =  nodal_average_paraview(str)

nodes_counter                             =  str.nodes_counter;
%--------------------------------------------------------------------------
% Stress conjugates.
%--------------------------------------------------------------------------
for i = 1:size(str.postproc.SigmaF,1)
    for j = 1:size(str.postproc.SigmaF,2)
        str.postproc.SigmaF(i,j,:)        =  reshape(str.postproc.SigmaF(i,j,:),size(str.postproc.SigmaF(i,j,:),3),1)'./nodes_counter;
    end
end

for i = 1:size(str.postproc.SigmaH,1)
    for j = 1:size(str.postproc.SigmaH,2)
        str.postproc.SigmaH(i,j,:)        =  reshape(str.postproc.SigmaH(i,j,:),size(str.postproc.SigmaH(i,j,:),3),1)'./nodes_counter;
    end
end

for i = 1:str.postproc.n_nodes
    str.postproc.SigmaJ(i,1)              =  str.postproc.SigmaJ(i,1)/nodes_counter(i);
end

for i = 1:size(str.postproc.SigmaF,1)
    for j = 1:size(str.postproc.SigmaF,2)
        str.postproc.First_Piola(i,j,:)   =  reshape(str.postproc.First_Piola(i,j,:),size(str.postproc.First_Piola(i,j,:),3),1)'./nodes_counter;
    end
end
%--------------------------------------------------------------------------
% Cauchy and pressure.
%--------------------------------------------------------------------------
for i = 1:size(str.postproc.sigma,1)
    str.postproc.sigma(i,:)               =  str.postproc.sigma(i,:)./nodes_counter;
end
str.postproc.sigma_pressure               =  str.postproc.sigma_pressure'./nodes_counter;
%--------------------------------------------------------------------------
% Lagrangian Electric displacement.
%--------------------------------------------------------------------------
for i = 1:size(str.postproc.D0,1)
    str.postproc.D0(i,:)                  =  str.postproc.D0(i,:)./nodes_counter;
end
%--------------------------------------------------------------------------
% Eulerian Electric displacement.
%--------------------------------------------------------------------------
for i = 1:size(str.postproc.D,1)
    str.postproc.D(i,:)                   =  str.postproc.D(i,:)./nodes_counter;
end
%--------------------------------------------------------------------------
% Electric field.
%--------------------------------------------------------------------------
for i = 1:size(str.postproc.E0,1)
    str.postproc.E0(i,:)                  =  str.postproc.E0(i,:)./nodes_counter;
end

for i = 1:size(str.postproc.E,1)
    str.postproc.E(i,:)                   =  str.postproc.E(i,:)./nodes_counter;
end
%--------------------------------------------------------------------------
% Deformation variables.
%--------------------------------------------------------------------------
for i = 1:str.postproc.n_nodes
    str.postproc.F(:,:,i)                 =  str.postproc.F(:,:,i)/nodes_counter(i);
end
for i = 1:str.postproc.n_nodes
    str.postproc.H(:,:,i)                 =  str.postproc.H(:,:,i)/nodes_counter(i);
end
for i = 1:str.postproc.n_nodes
    str.postproc.J(i,1)                   =  str.postproc.J(i,1)/nodes_counter(i);
end

end