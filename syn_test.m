% This is the benchmark test appeared in paper "Mapping the near-field
% scattered energy in sedimentary basins from ambient noise correlations"
% by Li, Clayton and Villa(2024)
% One can directly run this file to get all the results.
% alpha: the damping variable
% beta: the smoothing variable


load('scatter_syn_interarr_c3_point.mat')
% Replace the file for different cases:
% scatter_syn_interarr_c3_point.mat: inter array with point source
% scatter_syn_intraarr_c3_point.mat: intra array with point source
% scatter_syn_interarr_c3_horiz.mat: inter array with horiontal sources
% scatter_syn_interarr_c3_vert.mat:  inter array with vertical sources
% Variables in the *.mat files:
% c:       velocity model, in ny*nx
% mesh:    the (regular) mesh coordinates ny*nx*2, 3rd dimension: 1 for x, 2 for y
% SampleT: Sampling rate in s.
% sgn:     the correlation functions, recorded in ns*nt. ns is the number of
%          station pairs; nt is the number of time steps in eaca correlagram.
% x1, x2:  virtual source and virtual receiver in ns*2. First column for x
%          and second column for y. Each row represent a station-pair.
alpha = 1;
beta  = 1;
stat_x1=unique(x1,'rows');
stat_x2=unique(x2,'rows');
ind=ones([size(x1,1),1]);
mx=mesh(:,:,1);
my=mesh(:,:,2);
mx=mx(:);
my=my(:);
for i=1:length(x1)
    if all((x1(i,:)-x2(i,:))==0)||any(((x1(i,1)-mx).^2+(x1(i,2)-my).^2)==0)||any(((x2(i,1)-mx).^2+(x2(i,2)-my).^2)==0)
        ind(i)=0;
    end
end
x1=x1(ind==1,:);
x2=x2(ind==1,:);
sgn=sgn(ind==1,:);
sgnmax=max(sgn(:));
dsmp=10;
sourceT=5;
mesh=mesh(1:dsmp:end,1:dsmp:end,:);
c=c(1:dsmp:end,1:dsmp:end);
sgnamp=zeros(size(sgn));
for i=1:size(sgn,1)
    zxcv=fft(sgn(i,:));
    tmp=2*abs(ifft(zxcv.*(-sign([0:length(zxcv)-1]-floor(length(zxcv)/2))+1)/2));
    sgnamp(i,:)=tmp;
end 
nt=round((size(sgnamp,2)+1)/2);
nx=size(mesh,2);
ny=size(mesh,1);
m1=zeros(ny,nx);
dx=mesh(1,2,1)-mesh(1,1,1);
dy=mesh(2,1,2)-mesh(1,1,2);
pre_m=ones(size(m1));
source=exp(-([(-1*sourceT:SampleT:1*sourceT)]*3/sourceT).^2);
[COR1,submesh,ind,wt,~]=calc_f(mesh,c,m1,x1,x2,SampleT,nt);
[dL_dm0,RHS]=calc_dL_dm(c,m1,x1,x2,SampleT,nt,submesh,ind,wt,sgnamp,source,pre_m(:)',mesh);
K=(0.25*diag(ones(nx*ny-1,1),1)+0.25*diag(ones(nx*ny-1,1),-1)+...
   0.25*diag(ones(nx*ny-nx,1),nx)+0.25*diag(ones(nx*ny-nx,1),-nx)-diag(ones(nx*ny,1),0))/dx/dy;
lambda=alpha*ones(size(m1));
dL_dm=dL_dm0+diag(lambda(:)).^2+beta.^2*K*K';
[m,wt,~]=nnls(dL_dm,RHS,['Accy',0]);
m=reshape(m,size(m1));
[COR,submesh,ind,wt,~]=calc_f(mesh,c,m,x1,x2,SampleT,nt);
for i=1:size(COR,1)
    cor(i,:)=conv(COR(i,:),source,'same');
end
sensitivity=reshape(2*dL_dm*m(:),size(m1));


% make figure
figure
source_i=1; %change this index for virtual source index
S_max=max(max(sgnamp));
cormax=max(cor(:));
cmap=bluewhitered;
iii=1;
src_ind=ones([length(x1),1]);
for i=2:length(x1)
    if x1(i,1)~=x1(i-1,1) || x1(i,2)~=x1(i-1,2)
        iii=iii+1;
    end
    src_ind(i)=iii;
end
tt=((-nt+1):(nt-1))*SampleT;
ax1=subplot('Position',[0.08,0.1,0.35,0.385]);
imagesc(mesh(1,:,1),mesh(:,1,2),m*(tt(2)-tt(1)))
title(ax1,'b. output')
colormap(ax1,bluewhitered);
cb1=colorbar();
cb1.Label.Position = [cb1.Label.Position(1) - 0.5, cb1.Label.Position(2), cb1.Label.Position(3)];
ylabel(cb1,'Source intensity (mm/s)')
axis xy
hold on
scatter(stat_x1(:,1),stat_x1(:,2),30,'k')
scatter(stat_x2(:,1),stat_x2(:,2),30,'k')
scatter(stat_x1(source_i,1),stat_x1(source_i,2),300,'k','p','filled')
xlabel('x (km)','FontWeight','bold')
ylabel('y (km)','FontWeight','bold')
ax2=subplot('Position',[0.08,0.565,0.35,0.385]);
imagesc(mesh(1,:,1),mesh(:,1,2),c)
title(ax2,'a. input')
axis xy
colormap(ax2,'parula');
cb2=colorbar();
ylabel(cb2, 'velocity (km/s)')
hold on
scatter(stat_x1(:,1),stat_x1(:,2),30,'k')
scatter(stat_x2(:,1),stat_x2(:,2),30,'k')
scatter(stat_x1(source_i,1),stat_x1(source_i,2),300,'k','p','filled')
hold off
ylabel('y (km)','FontWeight','bold')
ax3=subplot('Position',[0.52,0.1,0.45,0.85]);
rec_sta_coord1=x2(src_ind==source_i,1);
rec_sta_coord2=x2(src_ind==source_i,2);
plot(repmat(tt,[sum(src_ind==source_i),1])',(3*sgn(src_ind==source_i,:)/S_max+repmat(sqrt((rec_sta_coord1-rec_sta_coord1(1)).^2+(rec_sta_coord2-rec_sta_coord2(1)).^2),[1,size(sgn,2)]))','k')
ylim([-2,max(sqrt((rec_sta_coord1-rec_sta_coord1(1)).^2+(rec_sta_coord2-rec_sta_coord2(1)).^2))+1])
hold on
plot(repmat(tt,[sum(src_ind==source_i),1])',(3*(cor(src_ind==source_i,:)/S_max)+repmat(sqrt((rec_sta_coord1-rec_sta_coord1(1)).^2+(rec_sta_coord2-rec_sta_coord2(1)).^2),[1,size(cor,2)]))','r')
ylim([-2,max(sqrt((rec_sta_coord1-rec_sta_coord1(1)).^2+(rec_sta_coord2-rec_sta_coord2(1)).^2))+1])
xlim([-150,150])
hold off
title(ax3,'c. waveform')
hold off
xlabel('t (s)','FontWeight','bold')
yl3=ylabel('distance (km)','FontWeight','bold');
yl3.Position(1)=yl3.Position(1)+10;

function [t,submesh,ind,wt,theta]=calc_f(mesh,c,m,x1,x2,dt,nt)
    %mesh: ny by nx by 2
    %x: nsta by 2
    %no auto-correlation should be included
    %avoid station being exactly on the grid
    [field,submesh,ind,wt]=interp_reg2D(mesh,m,30);
    freq=2;%Hz
    t=zeros(length(x1),2*nt-1);
    t1=zeros(length(x1),2*nt-1);
    t2=zeros(length(x1),2*nt-1);
    theta=zeros(length(x1),2*nt-1);
%     tic
    for i=1:length(x1)
        if size(c,1)==1&&size(c,2)==1
            t_submesh=(sqrt((submesh(:,:,1)-x2(i,1)).^2+(submesh(:,:,2)-x2(i,2)).^2))/c-(sqrt((submesh(:,:,1)-x1(i,1)).^2+(submesh(:,:,2)-x1(i,2)).^2))/c;
        else
            sx1=mesh(:,:,1);sx1=sx1(:);sx2=repmat(x2(i,1),size(sx1));
            sx=[sx1,sx2];
            sy1=mesh(:,:,2);sy1=sy1(:);sy2=repmat(x2(i,2),size(sy1));
            sy=[sy1,sy2];
            t_mesh1=calc_trave_t(mesh(1,:,1),mesh(:,1,2)',sx,sy,c);
            sx1=mesh(:,:,1);sx1=sx1(:);sx2=repmat(x1(i,1),size(sx));
            sx=[sx1,sx2];
            sy1=mesh(:,:,2);sy1=sy1(:);sy2=repmat(x1(i,2),size(sy));
            sy=[sy1,sy2];
            t_mesh2=calc_trave_t(mesh(1,:,1),mesh(:,1,2)',sx,sy,c);
            t_mesh=t_mesh1-t_mesh2;
            t_mesh=reshape(t_mesh,size(mesh(:,:,1)));
            [t_submesh,~,~,~]=interp_reg2D(mesh,t_mesh,30);
            [sub_c,submesh,ind,wt]=interp_reg2D(mesh,c,30);
        end
        tmpc=(x1(i,1)-mesh(:,:,1)).^2+(x1(i,2)-mesh(:,:,2)).^2;
        c1=c(tmpc==min(min(tmpc)));c1=c1(1);
        tmpc=(x2(i,1)-mesh(:,:,1)).^2+(x2(i,2)-mesh(:,:,2)).^2;
        c2=c(tmpc==min(min(tmpc)));c2=c2(1);
        amp_submesh=sub_c./sqrt(c1*c2).*1./sqrt(sqrt((submesh(:,:,1)-x2(i,1)).^2+(submesh(:,:,2)-x2(i,2)).^2).*sqrt((submesh(:,:,1)-x1(i,1)).^2+(submesh(:,:,2)-x1(i,2)).^2));
        cos1=((submesh(:,:,1)-x1(i,1))*(x2(i,1)-x1(i,1))+(submesh(:,:,2)-x1(i,2))*(x2(i,2)-x1(i,2)))./sqrt(((submesh(:,:,1)-x1(i,1)).^2+(submesh(:,:,2)-x1(i,2)).^2).*((x2(i,1)-x1(i,1)).^2+(x2(i,2)-x1(i,2)).^2));
        cos2=((submesh(:,:,1)-x2(i,1))*(x2(i,1)-x1(i,1))+(submesh(:,:,2)-x2(i,2))*(x2(i,2)-x1(i,2)))./sqrt(((submesh(:,:,1)-x2(i,1)).^2+(submesh(:,:,2)-x2(i,2)).^2).*((x2(i,1)-x1(i,1)).^2+(x2(i,2)-x1(i,2)).^2));
%        amp_submesh=amp_submesh.*abs(cos1.*cos2); %for Love wave only
        amp_submesh=amp_submesh(:);
        t_submesh=t_submesh(:);
        t(i,:)=accumarray(floor(t_submesh(t_submesh>-(nt-1)*dt&t_submesh<nt*dt)/dt)+nt,field(t_submesh>-(nt-1)*dt&t_submesh<nt*dt).*amp_submesh(t_submesh>-(nt-1)*dt&t_submesh<nt*dt),[2*nt-1,1]);
        theta(i,:)=acos(t1(i,:)./t(i,:)).*sign(t2(i,:));
    end
    theta(isnan(theta))=0;
end

function [dL_dm,RHS]=calc_dL_dm(c,m,x1,x2,dt,nt,submesh,ind,wt,Signal,source,pre_m,mesh)
    %dF_dm: 2nt-1 by nx*ny 

    dL_dm=zeros(length(m(:)));
    RHS=zeros(length(m(:)),1);
    tic
    for i=1:length(x1)
        if size(c,1)==1&&size(c,2)==1
            t_mesh=sqrt((x1(i,1)-x2(i,1)).^2+(x1(i,2)-x2(i,2)).^2)/c;
            t_submesh=(sqrt((submesh(:,:,1)-x2(i,1)).^2+(submesh(:,:,2)-x2(i,2)).^2))/c-(sqrt((submesh(:,:,1)-x1(i,1)).^2+(submesh(:,:,2)-x1(i,2)).^2))/c;
        else
            sx1=mesh(:,:,1);sx1=sx1(:);sx2=repmat(x2(i,1),size(sx1));
            sx=[sx1,sx2];
            sy1=mesh(:,:,2);sy1=sy1(:);sy2=repmat(x2(i,2),size(sy1));
            sy=[sy1,sy2];
            t_mesh1=calc_trave_t(mesh(1,:,1),mesh(:,1,2)',sx,sy,c);
            sx1=mesh(:,:,1);sx1=sx1(:);sx2=repmat(x1(i,1),size(sx));
            sx=[sx1,sx2];
            sy1=mesh(:,:,2);sy1=sy1(:);sy2=repmat(x1(i,2),size(sy));
            sy=[sy1,sy2];
            t_mesh2=calc_trave_t(mesh(1,:,1),mesh(:,1,2)',sx,sy,c);
            t_mesh=t_mesh1-t_mesh2;
            t_mesh=reshape(t_mesh,size(mesh(:,:,1)));
            [t_submesh,~,~,~]=interp_reg2D(mesh,t_mesh,30);
        end
        t_submesh=t_submesh(:);
        tmpc=(x1(i,1)-mesh(:,:,1)).^2+(x1(i,2)-mesh(:,:,2)).^2;
        c1=c(tmpc==min(min(tmpc)));c1=c1(1);
        tmpc=(x2(i,1)-mesh(:,:,1)).^2+(x2(i,2)-mesh(:,:,2)).^2;
        c2=c(tmpc==min(min(tmpc)));c2=c2(1);
        amp_mesh=c./sqrt(c1*c2).*1./sqrt(sqrt((mesh(:,:,1)-x2(i,1)).^2+(mesh(:,:,2)-x2(i,2)).^2).*sqrt((mesh(:,:,1)-x1(i,1)).^2+(mesh(:,:,2)-x1(i,2)).^2));
        %amp_mesh=amp_mesh.*abs(cos1.*cos2); %for Love wave only
        amp_mesh=amp_mesh(:);
        t_ind=t_submesh>-(nt-1)*dt&t_submesh<nt*dt;
        submat1=floor(t_submesh(t_ind)/dt)+nt+(2*nt-1)*(ind(t_ind,1)-1);
        submat2=floor(t_submesh(t_ind)/dt)+nt+(2*nt-1)*(ind(t_ind,2)-1);
        submat3=floor(t_submesh(t_ind)/dt)+nt+(2*nt-1)*(ind(t_ind,3)-1);
        submat4=floor(t_submesh(t_ind)/dt)+nt+(2*nt-1)*(ind(t_ind,4)-1);
        pre_OF=0*ones(size(Signal(i,:)));
        if size(c,1)==1&&size(c,2)==1
        	pre_OF(abs([1-nt:nt-1]*dt)<sqrt((x1(i,1)-x2(i,1)).^2+(x1(i,2)-x2(i,2)).^2)/c-30)=1;
        else
            pre_OF(abs([1-nt:nt-1]*dt)<calc_trave_t(mesh(1,:,1),mesh(:,1,2)',[x1(i,1),x2(i,1)],[x1(i,2),x2(i,2)],c)-30)=1;
        end
        dF_dm=reshape(accumarray(submat1,wt(t_ind,1),[length(m(:))*(2*nt-1),1])+accumarray(submat2,wt(t_ind,2),[length(m(:))*(2*nt-1),1])+accumarray(submat3,wt(t_ind,3),[length(m(:))*(2*nt-1),1])+accumarray(submat4,wt(t_ind,4),[length(m(:))*(2*nt-1),1]),[(2*nt-1),length(m(:))]);
        for j=1:size(dF_dm,2)
            dF_dm(:,j)=conv(dF_dm(:,j),source,'same');
            dF_dm(:,j)=dF_dm(:,j).*pre_OF';
        end
        for k=1:size(dF_dm,1)
            dF_dm(k,:)=dF_dm(k,:).*amp_mesh';
        end
        dL_dm=dL_dm+dF_dm'*dF_dm;
        RHS=RHS+dF_dm'*(pre_OF'.*Signal(i,:)');
        if mod(i,100)==1
            disp(['i:',num2str(i),', time=',num2str(toc)]);
        end
    end
end

function [out,submesh,ind,wt]=interp_reg2D(mesh,f,N)
%out:field on submesh, nsy by nsx
%ind:nsx*nsy by 4, ll,lu,rl,ru index for each submesh
%wt:nsx*nsy by 4, ll,lu,rl,ru weight
%%%%%%%%%%%test%%%%%%%%%%%%%%%%%%
% % % mesh=zeros(201,201,2);
% % % mesh(:,:,1)=repmat(-1:0.01:1,[201,1]);
% % % mesh(:,:,2)=repmat([-1:0.01:1]',[1,201]);
% % % f=exp(-mesh(:,:,1).^2-mesh(:,:,2).^2);
% % % [out,submesh,ind,wt]=interp_reg2D(mesh,f,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ny=size(mesh,1);
nx=size(mesh,2);
nsy=(ny-1)*N+1;
nsx=(nx-1)*N+1;
xmin=mesh(1,1,1);
xmax=mesh(1,end,1);
ymin=mesh(1,1,2);
ymax=mesh(end,1,2);
dx=(xmax-xmin)/(nx-1);
dy=(ymax-ymin)/(ny-1);
xsmax=xmax-dx/N/2;
xsmin=xmin+dx/N/2;
ysmax=ymax-dy/N/2;
ysmin=ymin+dy/N/2;
dsx=(xsmax-xsmin)/(nsx-1);
dsy=(ysmax-ysmin)/(nsy-1);
submesh=zeros(nsy,nsx,2);
submesh(:,:,1)=repmat(xsmin:dsx:xsmax,[nsy,1]);
submesh(:,:,2)=repmat([ysmin:dsy:ysmax]',[1,nsx]);
llx=floor((submesh(:,:,1)-xmin)/dx);%lower left corner x index(0,1,...)
lly=floor((submesh(:,:,2)-ymin)/dy);%lower left corner y index
ll=lly(:)+llx(:)*ny+1;
lu=ll+1;
rl=ll+ny;
ru=rl+1;
ind=[ll,lu,rl,ru];
meshx=mesh(:,:,1);
meshy=mesh(:,:,2);
submeshx=submesh(:,:,1);
submeshy=submesh(:,:,2);
meshx=meshx(:);
meshy=meshy(:);
submeshx=submeshx(:);
submeshy=submeshy(:);
wtll=(meshx(ru)-submeshx).*(meshy(ru)-submeshy)/dx/dy;
wtlu=(meshx(rl)-submeshx).*(submeshy-meshy(rl))/dx/dy;
wtrl=(submeshx-meshx(lu)).*(meshy(lu)-submeshy)/dx/dy;
wtru=(submeshx-meshx(ll)).*(submeshy-meshy(ll))/dx/dy;
wt=[wtll,wtlu,wtrl,wtru];
out=zeros(nsy,nsx,size(f,3));
for i=1:size(f,3)
ff=f(:,:,i);
ff=ff(:);
tmp=ff(ll).*wtll+ff(lu).*wtlu+ff(rl).*wtrl+ff(ru).*wtru;
out(:,:,i)=reshape(tmp,[nsy,nsx]);
end
end

function t=calc_trave_t(gx,gy,sx,sy,c)
    dx=gx(2)-gx(1);
    dy=gy(2)-gy(1);
    nx=length(gx);
    ny=length(gy);
    ns=length(sx(:,1));
    dist=sqrt((sx(:,1)-sx(:,2)).^2+(sy(:,1)-sy(:,2)).^2);
    nsegment = 1000;
    G=zeros(ns,nx*ny);
    for i=1:ns
        for j=1:nsegment
            index=ceil((sx(i,1)+((j-0.5)/nsegment)*(sx(i,2)-sx(i,1))-gx(1)+0.5*dx)/dx)+ceil((sy(i,1)+((j-0.5)/nsegment)*(sy(i,2)-sy(i,1))-gy(1)-0.5*dy)/dy)*nx;
            G(i,index)=G(i,index)+dist(i)/nsegment;
        end
    end
    m=reshape(1./c',[nx*ny,1]);
    t=G*m;
end