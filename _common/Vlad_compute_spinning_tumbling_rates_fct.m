function [Spinning_rate, Tumbling_rate]=Vlad_compute_spinning_tumbling_rates_fct(AllFibers,p)

        for ij=1:size(AllFibers.Centroid(:,1),1)
            index_time=find(~cellfun('isempty',AllFibers.Centroid(ij,:)))';
        
            if p.correct_orientations==1
                for it=index_time'
                    switch p.use_which_vectors_for_orientation
                        case 'poly'
                            AllFibers.red_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(1,1),AllFibers.EigenVectors_Tensor{ij,it}(1,2),AllFibers.EigenVectors_Tensor{ij,it}(1,3)];
                            AllFibers.green_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(2,1),AllFibers.EigenVectors_Tensor{ij,it}(2,2),AllFibers.EigenVectors_Tensor{ij,it}(2,3)];
                            AllFibers.blue_tensor{ij,it}=[AllFibers.EigenVectors_Tensor{ij,it}(3,1),AllFibers.EigenVectors_Tensor{ij,it}(3,2),AllFibers.EigenVectors_Tensor{ij,it}(3,3)];
                        case 'region'
                            AllFibers.red_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,1),AllFibers.EigenVectors{ij,it}{1,1}(2,1),AllFibers.EigenVectors{ij,it}{1,1}(3,1)];
                            AllFibers.green_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,2),AllFibers.EigenVectors{ij,it}{1,1}(2,2),AllFibers.EigenVectors{ij,it}{1,1}(3,2)];
                            AllFibers.blue_tensor{ij,it}=[AllFibers.EigenVectors{ij,it}{1,1}(1,2),AllFibers.EigenVectors{ij,it}{1,1}(2,3),AllFibers.EigenVectors{ij,it}{1,1}(3,3)];
                        case 'svd'
                            PrincAxis = Beppe_PricipalAxis(full(AllFibers.Object{ij,it}));
                            AllFibers.red_tensor{ij,it}=PrincAxis(:,1);
                            AllFibers.green_tensor{ij,it}=PrincAxis(:,2);
                            AllFibers.blue_tensor{ij,it}=PrincAxis(:,3);
                    end            
                end
        
                %%% set the red vector in the x direction for the first time-step
                if abs(Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,index_time(1)},[1 0 0]))>90
                    AllFibers.red_tensor{ij,index_time(1)}=-AllFibers.red_tensor{ij,index_time(1)};
                end
        
                for it=index_time(1:end-1)'
        
                    AllFibers.blue_tensor{ij,it}=cross(AllFibers.red_tensor{ij,it},AllFibers.green_tensor{ij,it});
        
        
                    if Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,it},AllFibers.red_tensor{ij,it+1})>p.angle
                        AllFibers.red_tensor{ij,it+1}=-AllFibers.red_tensor{ij,it+1};
                        if Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                            AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                            AllFibers.blue_tensor{ij,it+1}=cross(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                        end
                    elseif Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                        AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                        if Vlad_angle_between_two_vectors(AllFibers.red_tensor{ij,it},AllFibers.red_tensor{ij,it+1})>p.angle
                            AllFibers.red_tensor{ij,it+1}=-AllFibers.red_tensor{ij,it+1};
                            AllFibers.blue_tensor{ij,it+1}=cross(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                        end
                    elseif Vlad_angle_between_two_vectors(AllFibers.blue_tensor{ij,it},AllFibers.blue_tensor{ij,it+1})>p.angle
                        AllFibers.blue_tensor{ij,it+1}=-AllFibers.blue_tensor{ij,it+1};
                        if Vlad_angle_between_two_vectors(AllFibers.green_tensor{ij,it},AllFibers.green_tensor{ij,it+1})>p.angle
                            AllFibers.green_tensor{ij,it+1}=-AllFibers.green_tensor{ij,it+1};
                            AllFibers.red_tensor{ij,it+1}=-cross(AllFibers.blue_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1});
                        end
        
                    end
        
                    %AllFibers.Orientation_Tensor{ij,it+1}=Vlad_find_euler_angles(AllFibers.red_tensor{ij,it+1},AllFibers.green_tensor{ij,it+1},AllFibers.blue_tensor{ij,it+1});
                end
                AllFibers.blue_tensor{ij,index_time(end)}=cross(AllFibers.red_tensor{ij,index_time(end)},AllFibers.green_tensor{ij,index_time(end)});
        
            end
        
            timesteps=index_time-index_time(1)+1;
            time=(timesteps-1)*p.dt;
   
        
            red=[]; green=[]; blue=[]; R=[];
            for it=index_time'
                if it==1
                    tt=1;
                else
                    tt=it-index_time(1)+1;
                end
                % extract the vectors of the fibre
                red(tt,:)=AllFibers.red_tensor{ij,it};
                green(tt,:)=AllFibers.green_tensor{ij,it};
                blue(tt,:)=AllFibers.blue_tensor{ij,it};
        
            end
        
        
            if p.use_fitted_vectors_for_rotation_matrix == 1
                % fit the red vector
                [red_fitted_1,~]=fit_data(time,red(:,1),p.fitting_type,p.kernel_trajectory_vectors);
                [red_fitted_2,~]=fit_data(time,red(:,2),p.fitting_type,p.kernel_trajectory_vectors);
                [red_fitted_3,~]=fit_data(time,red(:,3),p.fitting_type,p.kernel_trajectory_vectors);
                % fit the green vector
                [green_fitted_1,~]=fit_data(time,green(:,1),p.fitting_type,p.kernel_trajectory_vectors);
                [green_fitted_2,~]=fit_data(time,green(:,2),p.fitting_type,p.kernel_trajectory_vectors);
                [green_fitted_3,~]=fit_data(time,green(:,3),p.fitting_type,p.kernel_trajectory_vectors);
        
                for it=index_time'
                        if it==1
                            tt=1;
                        else
                            tt=it-index_time(1)+1;
                        end
                
                        % compute a temporary vector perpendicular to the red and green fitted
                        % vectors
                        green_temp=cross([red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)],[green_fitted_1(tt); green_fitted_2(tt); green_fitted_3(tt)]);
                        green_temp=green_temp/norm(green_temp);
                
                        % compute a vector perpendicular to the red and the temporary
                        % vector (this is in the same plane as the fitted red and green vectors, but it's perpendicular to the fitted red one)
                        green_fitted=cross(green_temp,[red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)]); green_fitted=green_fitted/norm(green_fitted);
                        green_fitted_1(tt)=green_fitted(1);
                        green_fitted_2(tt)=green_fitted(2);
                        green_fitted_3(tt)=green_fitted(3);
                    
                        blue_fitted=cross([red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)],green_fitted); blue_fitted=blue_fitted/norm(blue_fitted);
                        blue_fitted_1(tt)=blue_fitted(1);
                        blue_fitted_2(tt)=blue_fitted(2);
                        blue_fitted_3(tt)=blue_fitted(3);
        
                        % generate the rotation matrix - each column corresponds to one vector of the fibre
                        R(:,1,tt)=[red_fitted_1(tt); red_fitted_2(tt); red_fitted_3(tt)];
                        R(:,2,tt)=[green_fitted_1(tt); green_fitted_2(tt); green_fitted_3(tt)];
                        R(:,3,tt)=[blue_fitted_1(tt); blue_fitted_2(tt); blue_fitted_3(tt)];
                end
            else             
                for it=index_time'
                    if it==1
                        tt=1;
                    else
                        tt=it-index_time(1)+1;
                    end
                    % generate the rotation matrix - each column corresponds to one vector of the fibre
                    R(:,1,tt)=red(tt,:)';
                    R(:,2,tt)=green(tt,:)';
                    R(:,3,tt)=blue(tt,:)';
                end
            end
        
            red_dot_1=rad2deg(Vlad_compute_derivative(time,red_fitted_1,p.derivative_stencil,p.disable_edge_points));
            red_dot_2=rad2deg(Vlad_compute_derivative(time,red_fitted_2,p.derivative_stencil,p.disable_edge_points));
            red_dot_3=rad2deg(Vlad_compute_derivative(time,red_fitted_3,p.derivative_stencil,p.disable_edge_points));
        
        
            % Computation with rotation matrix
            % skip some time-steps to lower the effect of noise on the result
            % compute the derivative of the rotation matrix
            R_dot=[]; R_q = []; av3=[];     OMEGA_s_R_dot=[]; omega_s_R_dot=[]; OMEGA_b_R_dot=[]; omega_b_R_dot=[];
            for tt=1:size(R,3)
                av2(tt,:)=[NaN, NaN, NaN];
                omega_b_R_dot(:,tt)=[NaN; NaN; NaN];
                omega_s_R_dot(:,tt)=[NaN; NaN; NaN];
            end
        
            for it=1:size(R,3)-p.skip_timesteps_rotation_matrix
        
                R_dot(:,:,it) = (R(:,:,it+p.skip_timesteps_rotation_matrix) - R(:,:,it))/p.skip_timesteps_rotation_matrix;
        
                R_q(:,:,1) = R(:,:,it);
                R_q(:,:,2) = R(:,:,it+p.skip_timesteps_rotation_matrix);
                q3 = quaternion(R_q,'rotmat','frame');
                av2_temp = rad2deg(angvel(q3,p.skip_timesteps_rotation_matrix,'frame')); % units in deg/frame
                av2(it,:) = av2_temp(2,:);  % disregard the first data point, because it is always wrong -> outlier
        
            end
        
            for tt=1:size(R,3)-p.skip_timesteps_rotation_matrix
                OMEGA_s_R_dot(:,:,tt)=rad2deg(R_dot(:,:,tt)*R(:,:,tt)');
            
                omega_s_R_dot(1,tt)=OMEGA_s_R_dot(3,2,tt);
                omega_s_R_dot(2,tt)=OMEGA_s_R_dot(1,3,tt);
                omega_s_R_dot(3,tt)=OMEGA_s_R_dot(2,1,tt);
           
            
                OMEGA_b_R_dot(:,:,tt)=rad2deg(R(:,:,tt)'*R_dot(:,:,tt));
            
                omega_b_R_dot(1,tt)=OMEGA_b_R_dot(3,2,tt);
                omega_b_R_dot(2,tt)=OMEGA_b_R_dot(1,3,tt);
                omega_b_R_dot(3,tt)=OMEGA_b_R_dot(2,1,tt);
                
            end
        
            omega_q2_x=av2(:,1);
            omega_q2_y=av2(:,3);
            omega_q2_z=-av2(:,2);

            for it=index_time'
                if it==1
                    tt=1;
                else
                    tt=it-index_time(1)+1;
                end
    
            Spinning_rate(ij,it)=omega_b_R_dot(1,tt);
            Tumbling_rate(ij,it)=sqrt(omega_b_R_dot(2,tt)^2 + omega_b_R_dot(3,tt)^2);
       
            end
        end

end