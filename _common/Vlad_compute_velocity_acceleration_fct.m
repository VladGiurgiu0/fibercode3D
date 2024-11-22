function [x_raw,y_raw,z_raw,...
            x_fitted,y_fitted,z_fitted,...
            x_dot,y_dot,z_dot,...
            x_dot_dot,y_dot_dot,z_dot_dot]=Vlad_compute_velocity_acceleration_fct(AllFibers,p)

        for ij=1:size(AllFibers.Centroid(:,1),1)
            index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';
        
        
            positions=cell2mat(AllFibers.Centroid_Refined(ij,index_time)') * p.dx; %%% stored as y,x,z
            timesteps=index_time-index_time(1)+1;
            time=(timesteps-1)*p.dt;
        
            x=positions(:,1);
            y=positions(:,2);
            z=positions(:,3);
       
        
            %%%%%%% ----positions and translational velocity---- %%%%%%%%
            [x_fitted_temp,~]=fit_data(time,x,p.fitting_type,p.kernel_trajectory_positions);
            [y_fitted_temp,~]=fit_data(time,y,p.fitting_type,p.kernel_trajectory_positions);
            [z_fitted_temp,~]=fit_data(time,z,p.fitting_type,p.kernel_trajectory_positions);
        
            x_dot_temp=Vlad_compute_derivative(time,x_fitted_temp,p.derivative_stencil,p.disable_edge_points);
            y_dot_temp=Vlad_compute_derivative(time,y_fitted_temp,p.derivative_stencil,p.disable_edge_points);
            z_dot_temp=Vlad_compute_derivative(time,z_fitted_temp,p.derivative_stencil,p.disable_edge_points);

            x_dot_dot_temp=Vlad_compute_derivative(time,x_dot_temp,p.derivative_stencil,p.disable_edge_points);
            y_dot_dot_temp=Vlad_compute_derivative(time,y_dot_temp,p.derivative_stencil,p.disable_edge_points);
            z_dot_dot_temp=Vlad_compute_derivative(time,z_dot_temp,p.derivative_stencil,p.disable_edge_points);
  

            %% check
%             figure(1); clf;
%             subplot(1,3,1); grid on; box on; hold on;
%             plot(time,x,'r.')
%             plot(time,y,'g.')
%             plot(time,z,'b.')
%             plot(time,x_fitted_temp,'r-')
%             plot(time,y_fitted_temp,'g-')
%             plot(time,z_fitted_temp,'b-.')
% 
%             legend('x','y','z','x fit','y fit', 'z fit')
%             xlabel('time'); ylabel('x,y,z')
% 
%             subplot(1,3,2); grid on; box on; hold on;
%             plot(time,x_dot_temp,'r-.')
%             plot(time,y_dot_temp,'g-.')
%             plot(time,z_dot_temp,'b-.')
% 
%             legend('x_dot','y_dot','z_dot')
%             xlabel('time'); ylabel('velocity')
% 
%             subplot(1,3,3); grid on; box on; hold on;
%             plot(time,x_dot_dot_temp,'r-.')
%             plot(time,y_dot_dot_temp,'g-.')
%             plot(time,z_dot_dot_temp,'b-.')
%             legend('x_dot_dot','y_dot_dot','z_dot_dot')
%             xlabel('time'); zlabel('accelerations')
            

            %%
            for it=index_time'
                if it==1
                    tt=1;
                else
                    tt=it-index_time(1)+1;
                end

            x_raw(ij,it)=x(tt);
            y_raw(ij,it)=y(tt);
            z_raw(ij,it)=z(tt);

            x_fitted(ij,it)=x_fitted_temp(tt);
            y_fitted(ij,it)=y_fitted_temp(tt);
            z_fitted(ij,it)=z_fitted_temp(tt);

            x_dot(ij,it)=x_dot_temp(tt);
            y_dot(ij,it)=y_dot_temp(tt);
            z_dot(ij,it)=z_dot_temp(tt);

            x_dot_dot(ij,it)=x_dot_dot_temp(tt);
            y_dot_dot(ij,it)=y_dot_dot_temp(tt);
            z_dot_dot(ij,it)=z_dot_dot_temp(tt);

            end
        end

end