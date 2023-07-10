%noise =[0 0.05 1];
%noise_1 = [0.00001]
noise = [0.01];

for vnoiseloop=1:length(noise)
    
    clearvars -except noise vnoiseloop order_vs_time_final* order_vs_time_meanstd* Cellnumber_final* num_radg_* Tot_number_* radG_av_* theta* coordsx_* coordsy_* orderparam_center_final*;
    close all;
     
    rng('shuffle');
    
    nSteps = 5000 %65000;
    Polar=zeros(nSteps,1);
    
    scat_fun=zeros(nSteps,1);  
     

    cyclenum =1;
     
    ds_msd=zeros(nSteps,cyclenum); 
     
    radG=zeros(nSteps,cyclenum); 
     
    den=zeros(nSteps,cyclenum);
     
    df=zeros(nSteps,cyclenum);
     
    ds_msdav=zeros(nSteps,1); 
     
    radG_av=zeros(nSteps,1); 
     
    den_av=zeros(nSteps,1);
     
    df_av=zeros(nSteps,1);
    
    liferows=1;
    visdata=zeros(liferows,15);
    visdata_row=zeros(liferows,15);
    speed_sum_temp=0;
    theta_sum_temp=0;
    dt = 5;         % Integration ticme
    dt2 = dt*dt;        % Integration time, squared
    
    vnoise =(noise(vnoiseloop));

     eval(['theta' strrep(num2str(vnoise), '.', 'p') '(:,1)= zeros(1,1);' ]);
    %eval(['theta' strrep(num2str(vnoise), '.', 'p') '_' num2str(cyclemsd) '=zeros(2000,1);' ]);
    
    %cyclemsd=2;
    for cyclemsd=1:cyclenum
    
%         L=200;
%         ab=[L L]; % rectangle dimensions; [width height]
%         R_min=3;      % minimum circle radius
%         R_max=4.5;     % maximum circle radius
%         cnst=false;
%         vis=false;
%         [C,R]=random_circle_packing_rectangle(ab,R_min,R_max,cnst,vis)
%         coords = C'; 
%         rad=R;
%         nPart = length(C);
%         numin = length(C);
        
        theta_sum_temp=0;  
        cs_theta=0;
        sn_theta=0;
        
        numin=200;
       
        % ===================
        %     Initialize
        % ===================
        
        % Set configuration parameters
        nPart = numin;        % Number of particles
        density = 0.001;     % Density of particles
        mass = 1;           % Particles' mass
        nDim = 2;           % The dimensionality of the system (3D in our case)
        
        % Set simulation parameters
%         dt = 10.0;         % Integration time
%         dt2 = dt*dt;        % Integration time, squared
        
        %nSteps = 3800;       % Total simulation time (in integration steps)
        sampleFreq = 100;    % Sampling frequency
        sampleCounter = 1;  % Sampling counter
        printFreq = 1000;     % Printing frequency
        plotFreq = 10;  %plot frequency
        %volrate = 1.4; %rate of vol growth dv/dt
        radmitosis = 5.0; %mitotic radius aka critical size
        samplecounter = 1;
        eta=0.005; %ECM Viscosity
        gfac=1; %08/09/22 ***********Changing to smaller values makes cell divison occur faster (0.5, 0.2, 1)
        taumin = gfac*1000; % second - minimum cell cycle time
        
        track(nDim,nSteps+1,numin)=zeros;
          
        % Set initial configuration
        %[coords L] = initCubicGrid(nPart,density,nDim);
       [coords L] = initnoverlap(nPart, density,nDim);



        
        % Set initial velocities with random numbers
        uo=0.1; % particle speed
       % vnoise =noise(vnoiseloop);
        ra=10;
        vels = zeros(nDim,nPart);
        theta=2*pi*(rand(1,nPart)-0.5);
        velsx=uo*cos(theta);
        velsy=uo*sin(theta);
        vels_tot = vels+[velsx; velsy];
        %vels(1, :) = velsx(1, :);
        %vels(2, :) = velsy(1, :);
    
        figure
        axis([-60*L 60*L -60*L 60*L])
        axis('square')
        hold on
            
        %randomly assign a radius to each particle
        rad = zeros(nPart,1);
        modulus = zeros(nPart,1);
        poisson = zeros(nPart,1);
        receptor = zeros(nPart,1);
        ligand = zeros(nPart,1);
        countdr=0; 
        
        eval(['theta_timecourse_temp' strrep(num2str(vnoise), '.', 'p') '{1}= theta;']);
        
        for part = 1:nPart
    
            rad(part,1) = randgaussrad(4.5,0.5); 
            modulus(part,1) = randgaussrad(10^-3,10^-4);
            poisson(part,1) = randgaussrad(0.5,0.02);
            receptor(part,1) = randgaussrad(0.9,0.02);
            ligand(part,1) = randgaussrad(0.9,0.02);
            lifetime(part,1) = 0;
            label(part,1)=part;
    
        end

       
        
        % ===================
        % Molecular Dynamics
        % ===================
        
        time = 0; % Following simulation time
        
        %volrate = (1*2*pi*(radmitosis)^3)/(3*taumin);
        volrate = (pi*(radmitosis)^2)/(2*taumin);
        for avini=1:numin
        initial(:,avini)=coords(:,avini);
        end
        
        for avini=1:numin
        track(:,1,avini)=coords(:,avini);
        end
          
    
        no_of_cells_ever_involved =  nPart;  
        
        for step = 1:nSteps

            oldcoords= coords;
            oldv =[velsx; velsy];
    
            theta_sum_temp=0;  
            cs_theta=0;
            sn_theta=0;
          
            centerM = zeros(nDim,1);
            % === Calculate new forces ===
            [forces,gamma3,pressure, theta] = Forcepara2(coords,rad,poisson,modulus,...
                nPart,receptor,ligand, vels_tot,theta, vnoise, ra, dt);
    
    %         coords(1,:) = coords(1,:) + velsx.*dt;
    %         coords(2,:) = coords(2,:) + velsy.*dt;
            % === Second integration step ===
            
            gammat = 0;
            % Implement the Andersen thermostat
            
	          for part =1:nPart
                %theta= (2.*pi).*0.05*(randn(1,1)-0);
    
                %vel_mag = sqrt(sum(bsxfun(@minus, vels(:,part), 0).^2, 1));
                %theta(1, part)= vels(1, part)/(vel_mag);
                %theta(2, part)= vels(2, part)/(vel_mag);
                %theta_noise(:, part) = theta(:,part)  (2.*pi).*0.05*(rand(1,1)-0.5); 
                gammat = 6*pi*eta*rad(part,1) + gamma3(part,1);
                %dropping the inertial term update the coords
                %coords(:,part) = coords(:,part) + (dt*forces(:,part))/(gammat);
                % Update velocities - All velocities are updated at once
                %vels(:,part) = vels(:,part)+ forces(:,part)/(gammat);
                %vels(1,part) = vels(1,part)+cos(theta_noise(1, part))*0.0 +0.*forces(1,part)/(gammat);
                %vels(2,part) = vels(2,part)+cos(theta_noise(2, part))*0.0 +0.*forces(2,part)/(gammat);
                %vels(3,part) = vels(3,part)+forces(3,part)/(gammat).*0.1;
                %dropping the inertial term update the coords
    %             coords(1,part) = coords(1,part) + u0*cos(theta(1,part));
    %             %(velsx(1, part).*1+vels(1, part)).*dt;
    %             coords(2,part) = coords(2,part) + u0*sin(theta(1,part));
                %(velsy(1, part).*1+vels(2, part)).*dt; 

                 velsx(1,part)=(uo*(cos(theta(1,part))));
                 velsy(1,part)=(uo*(sin(theta(1,part)))); 

                 coords(1,part) = (coords(1,part) + (velsx(1,part)*dt)+(dt*(forces(1,part)))/(gammat));
                 coords(2,part) = (coords(2,part) + (velsy(1,part)*dt)+(dt*(forces(2,part)))/(gammat));     
                
                 vels(:,part) = forces(:,part)/(gammat); 
                 
                 vels_tot(:,part) = (vels(:,part)+ [velsx(1,part) velsy(1,part)]');

                 theta(:,part) = atan2(vels_tot(2, part),vels_tot(1, part) );
           
                 lifetime(part,1) = lifetime(part,1)+ dt;
              end
             
            % === First integration step ===
            
            %update the size of particles - all radii updated at one time
            death = 10^(-20);
            deadpart=0;
	          deadnumin=0;
            
            for part=1:nPart
                
                if rand <= death*dt
                    
                    deadpart=deadpart+1;
                    deadindex(deadpart,1) = part;
                    
                elseif (rad(part,1) < radmitosis) && (pressure(part,1) < 0.0001)
                    
                    %grate = (1*volrate/(4*pi*rad(part,1)*rad(part,1)));
                    grate = (volrate/(2*pi*rad(part,1)));
                    rad(part,1) = rad(part,1) + dt*randgaussrad(grate,10^-7);
                    
                    
                elseif (rad(part,1) >= radmitosis) && (pressure(part,1) < 0.0001)
                    
                    %rad(end+1,1)=(2^(-1/3))*rad(part,1); %%****** FIX THIS *********
                    rad(end+1,1)=(2^(-1/2))*rad(part,1); 
                    rad(part,1) = (2^(-1/2))*rad(part,1); %new radius after division
                    modulus(end+1,1) = randgaussrad(10^-3,10^-4);
                    poisson(end+1,1) = randgaussrad(0.5,0.02);
   
                    receptor(end+1,1) = randgaussrad(0.9,0.02);
                    ligand(end+1,1) = randgaussrad(0.9,0.02);

                    lifetime(end+1,1) = 0; 
                    no_of_cells_ever_involved= no_of_cells_ever_involved +1;
                    label(end+1,1)= no_of_cells_ever_involved; 
                    
                    %must now add elements to vels, forces as well
                    %due to coords=coords+dt*vels + 0.5*dt2*forces equation
                    %vels(:,end+1)= ([0,0,0]'); %%TURNED OFF 12 Sept. 2022
                    vels(:,end+1)= ([0,0]'); 
                    theta(1, end+1) = pi.*2.*rand;
                    velsx(1, end+1)=uo*cos(theta(1, end));
                    velsy(1, end+1)=uo*sin(theta(1, end));
                    vels_tot(:,end+1) = vels(:,end)+ [velsx(1,end) velsy(1,end)]';
   
                    %generating random numbers between 0,1
                    a=0;
                    b=1;
                    r3=(b-a).*rand(1) + a;
                    r4=pi*(b-a).*rand(1) + pi*a;
                    r5=2*pi*(b-a).*rand(1) + 2*pi*a;
                    psi=size(coords,2)+1;
                    
                    %coords(1,psi) = coords(1,part)+radmitosis*(1-2^(-1/3))*sin(r4)*cos(r5);
                    %coords(2,psi) = coords(2,part)+radmitosis*(1-2^(-1/3))*sin(r4)*sin(r5);
                    %coords(3,psi) = coords(3,part)+radmitosis*(1-2^(-1/3))*cos(r4);c
                    coords(1,psi) = coords(1,part)+(radmitosis)*(1-1*2^(-1/2))*cos(r5);
                    coords(2,psi) = coords(2,part)+(radmitosis)*(1-1*2^(-1/2))*sin(r5);

                    coords(1,part)=coords(1,part)-(radmitosis)*(1-1*2^(-1/2))*cos(r5);
                    coords(2,part)=coords(2,part)-(radmitosis)*(1-1*2^(-1/2))*sin(r5);
                    %coords(3,part) = coords(3,part)-radmitosis*(1-2^(-1/3))*cos(r4);
                    
                end    
            end
    
            if deadpart > 0 
                
                    coords(:,deadindex(1:deadpart))=[];
                    rad(deadindex(1:deadpart))=[];
                    modulus(deadindex(1:deadpart)) = [];
                    poisson(deadindex(1:deadpart)) = [];
                    receptor(deadindex(1:deadpart)) =[];
                    ligand(deadindex(1:deadpart)) = [];
                    vels(:,deadindex(1:deadpart))= [];
		            lifetime(deadindex(1:deadpart))=[];
                    label(deadindex(1:deadpart))= [];  
                    theta(:, deadindex(1:deadpart)) =[];
                    velsx(:,deadindex(1:deadpart))= [];
                    velsy(:,deadindex(1:deadpart))= [];
	    
	            for part=1:deadpart
		       
		                    if deadindex(part,1)<= numin
                                deadnumin=deadnumin+1;
                                deadindexn(deadnumin,1)=deadindex(part,1);
                            end 
                end
            
            end
            
            if deadnumin > 0
               
                       track(:,:,deadindexn(1:deadnumin,1))= [];
                       initial(:,deadindexn(1:deadnumin,1))=[];            
            end
    
            numin=size(initial,2);
            nPart=size(coords,2);
            %this is the radius of gyration squared
	        centerM(1,1) = mean(coords(1,:));
            centerM(2,1) = mean(coords(2,:));
            %centerM(3,1) = mean(coords(3,:));
    
	          for part=1:nPart
                  radG(step,cyclemsd) = radG(step,cyclemsd) + ...
                    ((norm(coords(:,part) - centerM(:,1)))^2)/nPart; 
              end
            
            den(step,cyclemsd) =nPart/((1/3)*4*pi*radG(step,cyclemsd)^(3/2));  
            
            radG_av(step) = radG_av(step) + (radG(step,cyclemsd)/cyclenum);
            
            den_av(step) = den_av(step) + (den(step,cyclemsd)/cyclenum);
            
            % === Move time forward ===
            time = time + dt;
            %velsx=uo*(rand(1, nPart)-0.5)+vnoise*(rand(1,nPart)-0.5);
            %velsy=uo*(rand(1, nPart)-0.5)+vnoise*(rand(1,nPart)-0.5);
            if mod(step,printFreq) == 0
                step % Print the step
                cyclemsd  
            end
            
            if mod(step,plotFreq) == 0
                numcells(samplecounter,cyclemsd) = nPart;
                samplecounter = samplecounter + 1;
            end
            
            
            for avini=1:numin
                ds_msd(step,cyclemsd)=(norm(coords(:,avini)-initial(:,avini)))^2;  
                ds_msdav(step) = ds_msdav(step) + (ds_msd(step,cyclemsd))/(cyclenum*numin);
                
                df(step,cyclemsd)=(norm(coords(:,avini)-initial(:,avini)))^4;
                df_av(step)=df_av(step)+(df(step,cyclemsd))/(cyclenum*numin);
            end
            
            for avini=1:numin
                track(:,step+1,avini)=coords(:,avini);  
            end
            
           
            save('fnumin.txt','numin','-ascii','-append');
            Tot_number(step)=nPart;

	         if mod(step,plotFreq) == 0
        		    %countdr = countdr+1;
        		    %deltart(countdr,cyclemsd) = rtumor(coords,centerM);
		                for part=1:nPart
                             visdata(1,1:nDim)=coords(1:nDim,part);
                             visdata(1,nDim+1)= label(part,1);
                             visdata(1,nDim+2)= lifetime(part,1);
                             visdata(1,nDim+3)= step*dt;
                             visdata(1,nDim+4)= receptor(part,1);
                             visdata(1,nDim+5)= ligand(part,1);
                             visdata(1,nDim+6)= modulus(part,1);
                             visdata(1,nDim+7)= poisson(part,1);
                             visdata(1,nDim+8)= rad(part,1);
                             visdata(1,nDim+9)=velsx(1,part);
                             visdata(1,nDim+10)=velsy(1,part);
                             visdata(1,nDim+11)=theta(1,part);
                             visdata(1,nDim+12)=sqrt(velsx(1,part).^2+velsy(1,part).^2);
                             visdata(1,nDim+13)=cyclemsd;
                             visdata_row= visdata(1,:);
                             %timecourse_order(nPart,:) = visdata_row(1,:);
                             %timecourse_order(step,2) = sum(visdata(:,13))/(nPart*(sum(visdata(:,14)))) ;
                             save('lifetime1.txt','visdata_row', '-ascii','-append');
                             speed_sum_temp = speed_sum_temp+visdata_row(1, 14);
                             cs_theta = (cs_theta +cos((visdata_row(1,13)))/1)/(nPart*1);
                             sn_theta = (sn_theta+sin((visdata_row(1,13)))/1)/(nPart*1);
                             %plotcell(part) = nsidedpoly(600, 'Center', [coords(1,part) coords(2,part)], 'Radius', 5);                          
                        end
    
                       % order_vs_time(step,1) =step;
                       % order_vs_time(step,2) =theta_sum_temp/(nPart*1);
                        eval(['order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(step,1)= step;']);
                        eval(['order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(step,cyclemsd+1)= abs(cs_theta+sn_theta);']);
                        %eval(['theta' strrep(num2str(vnoise), '.', 'p') '_' num2str(cyclemsd) '(:,step)= theta ;']);
                        %Tot_number(step)=nPart; %***USED 
                      %  eval(['Tot_number' strrep(num2str(vnoise), '.', 'p') '(step,cyclemsd)=nPart/radG_av ;' ]);
                      %  eval(['theta' strrep(num2str(vnoise), '.', 'p') '(:,1)= part;']);
                        %eval(['theta' strrep(num2str(vnoise), '.', 'p') '(:,cyclemsd+1)= theta ;']);
                      %  eval(['theta' strrep(num2str(vnoise), '.', 'p') '=' 'theta' strrep(num2str(vnoise), '.', 'p') '( find(' 'theta' strrep(num2str(vnoise), '.', 'p') '(:,1)) ,:)'  ]);
                      eval(['theta_timecourse_temp' strrep(num2str(vnoise), '.', 'p') '{step+1}= theta;']);
                      eval(['velxy_deltat_' strrep(num2str(cyclemsd), '.', 'p') ' = oldcoords-coords(:,length(oldcoords));' ]);
%                       velxy_deltat = (oldcoords-coords(:,length(oldcoords)))./norm(oldcoords-coords(:,length(oldcoords)));
%                       velxy_deltat =(oldcoords-coords(:,length(oldcoords)))./norm(oldcoords-coords(:,length(oldcoords)));

                      eval(['orderparam_center_' strrep(num2str(vnoise), '.', 'p') '{cyclemsd, step}=mean(dot((coords-centerM.*(ones(2,length(coords)))), [vels_tot])./(vecnorm(coords-centerM.*(ones(2,length(coords)))).*vecnorm([vels_tot])));']);

%                                                       norm_velsx = (velsx)./sqrt(velsx.^2+velsy.^2);
%                       norm_velsy = (velsy)./sqrt(velsx.^2+velsy.^2);

%                       angle_delta = atan2(velxy_deltat(2,:), velxy_deltat(1,:));
                      %diff_rdeltai_vi= abs(velxy_deltat-[norm_velsx(1:length(velxy_deltat )); norm_velsy(1:length(velxy_deltat))]); 
                      %diff_rdeltai_vi= abs([velsx; velsy]./sqrt(velsx.^2+velsy.^2)-[oldv(1,:); oldv(2,:)]./sqrt(oldv(1,:).^2+oldv(2,:).^2)); 
                      %diff_angle = atan2(velsy-oldv(2,:), velsx-oldv(1,:));
                      %total_cos_theta = sum(cos(angle_delta), 2);
                      %total_sin_theta = sum(sin(angle_delta), 2);
                      %total_add_cos_sin = sqrt(total_sin_theta.^2+total_cos_theta.^2);
                      %average_order_param = total_add_cos_sin/(length(velxy_deltat));
                      %eval(['orderparameter_timecourse' strrep(num2str(cyclemsd), '.', 'p') '(step,cyclemsd) = sum(velxy_deltat_' strrep(num2str(cyclemsd), '.', 'p') ',2);' ]);
                      %eval(['orderparameter_angle_timecourse' strrep(num2str(cyclemsd), '.', 'p') '(step,cyclemsd) = average_order_param;' ]);
             end
%              cla
%              set(gcf,'doublebuffer','on')         
%              plot(coords(1,:),coords(2,:), 'bo' );
%              hold on
%              quiver(coords(1,:),coords(2,:),10*velsx(1,:)./sqrt(velsx(1,:).^2+velsy(1,:).^2),10*velsy(1,:)./sqrt(velsy(1,:).^2+velsy(1,:).^2),'Color','red', 'AutoScale','off', 'LineWidth',1.5, 'AutoScaleFactor',2, 'MaxHeadSize',10);
%              drawnow
            
             %hold off
    
    %          figure
    %          hold on
    %          %cla
    %          %set(gcf,'doublebuffer','on')  
    %          plot(step, sum(visdata(:,13))/(nPart*(sum(visdata(:,14)))));
    %          drawnow

        end

    eval(['order_vs_time_final_noise' strrep(num2str(vnoise), '.', 'p') '(:,1)=' 'order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(find(' 'order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(:,1)),1);']);
    eval(['order_vs_time_final_noise' strrep(num2str(vnoise), '.', 'p') '(:,cyclemsd+1)=' 'order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(find(' 'order_vs_time_noise' strrep(num2str(vnoise), '.', 'p') '(:,1)),cyclemsd+1);']);
    
    eval(['order_vs_time_meanstd_noise' strrep(num2str(vnoise), '.', 'p') '(:,1)=' '(order_vs_time_final_noise' strrep(num2str(vnoise), '.', 'p') '(:,1));']);
    eval(['order_vs_time_meanstd_noise' strrep(num2str(vnoise), '.', 'p') '(:,2)=' 'mean(order_vs_time_final_noise' strrep(num2str(vnoise), '.', 'p') '(:,2:cyclemsd+1),2);']);
    eval(['order_vs_time_meanstd_noise' strrep(num2str(vnoise), '.', 'p') '(:,3)=' 'std(order_vs_time_final_noise' strrep(num2str(vnoise), '.', 'p') '(:,2:cyclemsd+1),0,2);']);
    

     eval(['Cellnumber_final(vnoiseloop,1)=' 'vnoise;']);
     eval(['Cellnumber_final(vnoiseloop,cyclemsd+1)=' 'nPart;']);
     theta_transpose = theta';
     eval(['theta' strrep(num2str(vnoise), '.', 'p') '= vertcat(theta' strrep(num2str(vnoise), '.', 'p') ', theta_transpose) ;']);   

     eval([ 'theta_timecourse' strrep(num2str(vnoise), '.', 'p') '{cyclemsd,:}= theta_timecourse_temp' strrep(num2str(vnoise), '.', 'p') '(~cellfun(''isempty'',' 'theta_timecourse_temp' strrep(num2str(vnoise), '.', 'p') '));']);

       
        % Simulation results
        % ===================
        radG_inst(1,1) = sqrt(radG(step,cyclemsd));
        rho = nPart/((4/3)*pi*radG_inst^3);
        
        L=2*radG_inst(1,1);
        dL = 1.0;

        save('initial.txt','initial','-ascii','-append');    

     Tot_number= Tot_number(find(Tot_number));
     eval([ 'Tot_number' strrep(num2str(vnoise), '.', 'p') '(cyclemsd,:)=Tot_number ;']);
     eval([ 'Tot_number_by_RadGav' strrep(num2str(vnoise), '.', 'p') '(:, cyclemsd)=Tot_number''./radG_av ;' ]);

      eval(['coordsxy_combined_' strrep(num2str(cyclemsd), '.', 'p') ' = coords;' ]);
      %eval(['coordsy_combined_' strrep(num2str(cyclemsd), '.', 'p') ' = coords(2,:);'  ]);

    end
    eval(['idx=not(cellfun(@isempty,orderparam_center_' strrep(num2str(vnoise), '.', 'p')  '));']);
    eval(['orderparam_center_final_' strrep(num2str(vnoise), '.', 'p') '=arrayfun(@(x) cell2mat(orderparam_center_' strrep(num2str(vnoise), '.', 'p')  '(x,idx(x,:))),[1:size( orderparam_center_' strrep(num2str(vnoise), '.', 'p') ' ,1)]'',''un'',0);']);

    
      save('numcellf.txt','numcells','-ascii');
      save('fds_msd.txt','ds_msd','-ascii');
      save('fds_msdav.txt','ds_msdav','-ascii');
      save('fradG.txt','radG','-ascii');
      save('fradG_av.txt','radG_av','-ascii');
      save('fden.txt','den','-ascii');
      save('fden_av.txt','den_av','-ascii');
      save('radius.txt','rad','-ascii');
      %save('deltart.txt','deltart','-ascii'); 
end
