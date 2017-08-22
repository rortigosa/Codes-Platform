% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %
% % %  Provisional postprocessing in Matlab
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % figure(1)
% for ielem=1:str.n_elem
%     %----------------------------------------------------------------------
%     %  first face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [1 2 3];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  second face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [1 2 4];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  third face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [2 3 4];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  fourth face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [3 1 4];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
% end 
% hv         =  [h1;h2;h3;h4];
% colorbar('vert')
% set(hv,'edgecolor','none'); 
% set(hv,'MarkerEdgeColor','none');
% set(hv); 
% drawnow;
% % % 
% % % 
 figure(4)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 3];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    %h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 3 4];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    %h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 1 4];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    %h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
end 
hv         =  [h1;h2;h3;h4];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;
% % % 
% % % 
% % % figure(1)
% % % for ielem=1:str.n_elem
% % % %for ielem=80:80
% % %     %----------------------------------------------------------------------
% % %     %  first face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [1 2 6 5];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % %     %----------------------------------------------------------------------
% % %     %  second face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [2 4 8 6];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % %     %----------------------------------------------------------------------
% % %     %  third face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [4 3 7 8];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % %     %----------------------------------------------------------------------
% % %     %  fourth face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [3 1 5 7];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % %     %----------------------------------------------------------------------
% % %     %  fifth face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [1 2 4 3];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h5                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % %     %----------------------------------------------------------------------
% % %     %  sixth face
% % %     %----------------------------------------------------------------------
% % %     face_connectivity         =  [5 6 8 7];
% % %     face_nodes                =  str.connectivity(ielem,face_connectivity);
% % %     phi_face                  =  str.phi(face_nodes);
% % %     x_face                    =  str.Eulerian_x(:,face_nodes);
% % %     h6                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
% % %     hold on
% % % end 
% % % hv         =  [h1;h2;h3;h4;h5;h6];
% % % colorbar('vert')
% % % %set(hv,'edgecolor','none'); 
% % % set(hv,'MarkerEdgeColor','none');
% % % set(hv); 
% % % drawnow;
% % % 
% % % 
% % % N_loads         =  105;
% % % u               =  zeros(N_loads,1);
% % % phi             =  zeros(N_loads,1);
% % % for iload=1:N_loads
% % %     filename    =  ['Load_increment_' num2str(iload) '.mat'];
% % %     str         =  new_str;
% % %     load(filename)
% % %     u(iload)    =  norm(str.Eulerian_x(:,165) - str.Lagrangian_X(:,165));
% % %     phi(iload)  =  max(str.phi);
% % % end
% % % plot(phi,u)
% %     
% %      
figure(5)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face 
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 3];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
%    h1                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face,'EdgeColor','none');    
    h1                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
%    h2                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face,'EdgeColor','none');    
    h2                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 3 4];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
%    h3                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face,'EdgeColor','none');    
    h3                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 1 4];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
%    h4                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face,'EdgeColor','none');    
    h4                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
end 
hv         =  [h1;h2;h3;h4];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;
% % 
% % 
% figure(3)
% for ielem=1:str.n_elem
%     %----------------------------------------------------------------------
%     %  first face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [1 3 6];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  second face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [1 3 10];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  third face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [3 6 10];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
%     %----------------------------------------------------------------------
%     %  fourth face
%     %----------------------------------------------------------------------
%     face_connectivity         =  [6 1 10];
%     face_nodes                =  str.connectivity(ielem,face_connectivity);
%     phi_face                  =  str.phi(face_nodes);
%     x_face                    =  str.Eulerian_x(:,face_nodes);
%     h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),phi_face);    
%     hold on
% end 
% hv         =  [h1;h2;h3;h4];
% colorbar('vert')
% set(hv,'edgecolor','none'); 
% set(hv,'MarkerEdgeColor','none');
% set(hv); 
% drawnow;
% % 
% % 
% % 
% % 

figure(1)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4 3];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 4 8 6];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 3 7 5];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 4 8 7];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fifth face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 6 5];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h5                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  sixth face
    %----------------------------------------------------------------------
    face_connectivity         =  [5 6 8 7];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.p(ielem)*ones(4,1);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h6                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
end  
hv         =  [h1;h2;h3;h4;h5;h6];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;



figure(1)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4 3];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 4 8 6];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h2                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 3 7 5];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h3                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 4 8 7];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h4                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fifth face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 6 5];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h5                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  sixth face
    %----------------------------------------------------------------------
    face_connectivity         =  [5 6 8 7];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face                    =  str.phi(face_nodes);
    x_face                    =  str.Eulerian_x(:,face_nodes);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h6                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face);    
    hold on
end 
hv         =  [h1;h2;h3;h4;h5;h6];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;




figure(1)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4 3];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h1                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 4 8 6];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h2                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 3 7 5];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h3                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 4 8 7];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h4                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fifth face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 6 5];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h5                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  sixth face
    %----------------------------------------------------------------------
    face_connectivity         =  [5 6 8 7];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h6                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
end 
hv         =  [h1;h2;h3;h4;h5;h6];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;



figure(1)
for ielem=1:str.n_elem
    %----------------------------------------------------------------------
    %  first face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 4 3];
    p_face_connectivity         =  [1 2 4 3];
    face_nodes                =  str.connectivity(ielem,face_connectivity);
    p_face_nodes              =  str.connectivity_p(ielem,p_face_connectivity);
    p_face                    =  str.p(p_face_nodes);
    x_face                    =  str.nodes(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h1                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  second face
    %----------------------------------------------------------------------
    face_connectivity         =  [2 4 8 6];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h2                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  third face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 3 7 5];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h3                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fourth face
    %----------------------------------------------------------------------
    face_connectivity         =  [3 4 8 7];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h4                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  fifth face
    %----------------------------------------------------------------------
    face_connectivity         =  [1 2 6 5];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h5                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
    %----------------------------------------------------------------------
    %  sixth face
    %----------------------------------------------------------------------
    face_connectivity         =  [5 6 8 7];
    face_nodes                =  str.connectivity_p(ielem,face_connectivity);
    p_face                    =  str.p(face_nodes);
    x_face                    =  str.nodes_p(face_nodes,:);
    %h1                        =  fill3(x_face(1,:),x_face(2,:),x_face(3,:),p_face,'EdgeColor','none');    
    h6                        =  fill3(x_face(:,1),x_face(:,2),x_face(:,3),p_face);    
    hold on
end 
hv         =  [h1;h2;h3;h4;h5;h6];
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;







