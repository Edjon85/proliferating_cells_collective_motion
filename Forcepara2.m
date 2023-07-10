function [forces,gamma3,pressure, theta] = Forcepara2(coords,rad,poisson,modulus,...
                                             nPart,receptor,ligand, vels_tot, theta, vnoise, ra,dt)


etab=0.0001;
% Initialize all forces to 0close 
forces = zeros(size(coords));
gamma3 = zeros(size(rad));
pressure = zeros(size(rad));
velsmean = zeros(size(coords));
%lengthint=0;
%sumtheta = zeros(1,length(coords));
%theta_vtot=(atan2((vels_tot(2,:)),(vels_tot(1,:))));

% Get the number of particles
%nPart = size(coords,2);

f = 0.0001; %adhesion f^{ad} parameter

for part=1:nPart %change to parfor to call parallel looping in Matlab
%instead of looping over all pairs in the force calculation, find the distance between all
%pairs

dlist= sqrt(sum(bsxfun(@minus, coords(:,part), coords).^2, 1));
%dlist= sqrt(sum(bsxfun(@minus, coords(1:2,part), coords(1:2,:)).^2, 1));

[d_dist, ind_dist] = sort(dlist);
%index each set of nearest distances
[d, ind] = sort(rad(part,1)+rad'-dlist);
                
                %begin_index=find(d==2*rad(part,1));
                begin_index=find(d>0.0);
                %ind_closest = ind(2:(min(nPart,26))); %find the n nearest neighbors
                ind_closest=ind((begin_index));
                coords_closest = coords(:,ind_closest);
                %rad_closest=rad(ind_closest,1);
                %poission_closest=poisson(partA,1)
                for partA=1:(size(ind_closest,2))
                    dr =  coords(:,part) - coords_closest(:,partA);
                    
                    if norm(dr) > 0.0 
                        Rij = norm(dr);
                        RijHat = dr/Rij;     
                        hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
                        Eij = ((1 - poisson(ind_closest(partA),1)^2)/(modulus(ind_closest(partA),1)) ...
                               + (1 - poisson(part,1)^2)/(modulus(part,1)));  
                        Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));   
                        areaint = pi*(1/Rijf)*hij;
                         lengthint = (2*sqrt((rad(part,1)^2-((rad(part,1)^2-rad(ind_closest(partA),1)^2+Rij^2)^2)/(2*Rij)^2)))  ;
                        invDr2 = (hij)^(3/2); % 1/r^2
                        %forceFact = (invDr2/(0.75*(Eij)*sqrt(Rijf)))-areaint*f...
                        %*0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                         %     receptor(part,1)*ligand(ind_closest(partA),1));
                        %pressure(part,1) = pressure(part,1)+ ...
                        %abs(forceFact/areaint);

                        forceFact = (invDr2/(0.75*(Eij)*sqrt(Rijf)))-lengthint*f...
                        *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                              receptor(part,1)*ligand(ind_closest(partA),1));

                        pressure(part,1) = pressure(part,1)+ ...
                                        abs(forceFact/lengthint);
                             %abs(forceFact/areaint);
                        forces(:,part) = forces(:,part) + (forceFact*RijHat);
                        %velsmean(:, part) = vels(:, part) + vels(:, ind_closest(partA),1);  
                    end 
                end
               %velsmean(:, part) = velsmean(:,part)./(size(ind_closest,2)+1); 

               begin_index_theta=find(d_dist<ra); 
               begin_index_nearneighbor_order=find(d_dist<2*ra);
                %find the n nearest neighbors
                %clear ind_closest theta_closest;
                ind_closest_theta=ind_dist((begin_index_theta));
                %vels_closest = vels(:, ind_closest_theta);
                theta_closest = theta(1,ind_closest_theta);
               % theta_closest = theta_vtot(1,ind_closest_theta);
                %sumtheta(1,part) = mod(atan2d(sum((sin(theta_closest))),sum((cos(theta_closest)))), 360);
                sumtheta(1,part) = (atan2(mean((sin(theta_closest))),mean((cos(theta_closest)))));
                %sumtheta(1,part) = (atan2(mean((sin(theta_vtot(part)))),mean((cos(theta_vtot(part))))));
              

               for partA=1:(size(ind_closest,2))
                
                    dr =  coords(:,part) - coords_closest(:,partA);
                    
                    if norm(dr)>0.0
                        Rij = norm(dr);
                        RijHat = dr/Rij;
                        hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
                        Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
                        mag=norm(forces(:,part));
                        areaint = pi*(1/Rijf)*hij;    
                        lengthint = (2*sqrt((rad(part,1)^2-((rad(part,1)^2-rad(ind_closest(partA),1)^2+Rij^2)^2)/(2*Rij)^2))) ;
                        %gamma2 = etab*areaint*0.5*(1+(sum(bsxfun(@times, forces(:,part),RijHat)))/mag)...
                        %*0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                        %receptor(part,1)*ligand(ind_closest(partA),1));
                         gamma2 = etab*lengthint*0.5*(1+(sum(bsxfun(@times, forces(:,part),RijHat)))/mag)...
                        *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                        receptor(part,1)*ligand(ind_closest(partA),1));
         
                        gamma3(part,1) = gamma3(part,1) + gamma2;
                    
                    end
              
                end
                
 end
%theta(1, :) = (sumtheta(1,:)+ sqrt(vnoise*dt)*1*pi*(rand(1,nPart)-0.5)); 
theta(1, :) = (sumtheta(1,:)+ vnoise*sqrt(dt)*pi*(rand(1,nPart)-0.5));
%theta(1, :) = (sumtheta(1,:)+ vnoise*2*pi*(rand(1,nPart)-0.5)) ; 
                
end
