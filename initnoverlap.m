function [coords, L] = initoverlap(nPart, density, nDim)
ncell=nPart;
radius=5;
tot_time =10;
L=250;
lbox=L;
%r=15;
r_rep=10;
k=5;
vrep_init=10;
vrep=vrep_init;
onebd=ones(1,ncell)';
ang=pi.*2.*(rand(ncell,1));
init_R = lbox;
init_theta = 2*pi*rand(ncell,1);
init_r = init_R*(rand(ncell,1));
xb = init_r.*cos(init_theta);
yb = init_r.*sin(init_theta);
vxb=cos(ang);
vyb=sin(ang);
ang=pi.*2.*(rand(ncell,1));
ang_force = zeros(ncell,1);
dt=3;
for nsteps=1:tot_time;
    nsteps;
    %  if(nsteps<t_relax)
    %      vs = 0; ve=0; vrep = vrep_init;
    for cell1=1:ncell
        %find mean angle of neigbours (include cell1)
        sep(1:ncell)=sqrt((onebd.*xb(cell1)-xb(1:ncell)).^2+...
        (onebd.*yb(cell1)-yb(1:ncell)).^2);
        dist_x(1:ncell)=(-onebd.*xb(cell1)+ xb(1:ncell));
        dist_y(1:ncell)=(-onebd.*yb(cell1)+yb(1:ncell));
        ang_force=atan2(dist_y(find(sep<r_rep & sep>0)), dist_x(find(sep<r_rep & sep>0) ));
        ang_force_mean(cell1)=mean(ang_force);
        forcex(cell1)=mean(-k.*((2*radius-(sep(sep<r_rep & sep>0)'))).*cos(ang_force'));
        forcey(cell1)=mean(-k.*((2*radius-(sep(sep<r_rep & sep>0)'))).*sin(ang_force'));
        forcex(isnan(forcex))=0;
        forcey(isnan(forcey))=0;
        mag_force(cell1)=sqrt(forcex(cell1).^2+forcey(cell1).^2);
        ang_force2(cell1)=atan2(forcey(cell1),forcex(cell1) );
        %nearang=ang(sep<r);  %'r' this is alignment radius
        %mang(cell1)=mean(nearang);
        %plotcell(cell1) =nsidedpoly(1000, 'Center', [xb(cell1) yb(cell1)], 'Radius', 5);
    end
    ang_force2(isnan(ang_force2))=0;
    ang_force_mean(isnan(ang_force_mean))=0;
    ang_force2_trans=ang_force2';
    %ang=mang' ;
    %mag_force(find(sep>r_rep))=0;
    vxb= (mag_force.*(cos(ang_force2)))' ;
    vyb=(mag_force.*(sin(ang_force2)))' ;
    xb=xb+vxb.*dt;
    yb=yb+vyb.*dt;
    coords(1,:) = xb';
    coords(2,:) = yb';
end         
  end


       