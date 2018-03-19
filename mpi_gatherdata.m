load mpi_param.data;

SIZEX = mpi_param(1);
SIZEY = mpi_param(2);

dims(1) = mpi_param(3);
dims(2) = mpi_param(4);

DSTEP = mpi_param(5);

x = SIZEX/dims(1);
y = SIZEY/dims(2);

SIZEX = SIZEX/DSTEP;
SIZEY = SIZEY/DSTEP;
x = x/DSTEP;
y = y/DSTEP;

data_hole = zeros(SIZEX, SIZEY);
data_hole_etendu = zeros(SIZEX, SIZEY);

for i=1:dims(1)
    for j=1:dims(2)
        source = (i-1)*dims(2)+j;
        
        file_name=sprintf('hole%d.txt', source);
        data=load(file_name);
        
        data_hole(1+(i-1)*x:i*x,1+(j-1)*y:j*y) = data;
        
        file_name=sprintf('hole_etendu%d.txt', source);
        data=load(file_name);
        
        data_hole_etendu(1+(i-1)*x:i*x,1+(j-1)*y:j*y) = data;
    end
end

%save data_hole.txt data_hole
%save data_hole_etendu.txt data_hole_etendu

%system("grep -vE '#' data_hole.txt > tmp.txt");
%system("mv tmp.txt data_hole.txt");
%system("grep -vE '#' data_hole_etendu.txt > tmp.txt");
%system("mv tmp.txt data_hole_etendu.txt");

%hole = data_hole(1:10:4950,1:10:4950);
%hole_etendu = data_hole_etendu(1:10:4950,1:10:4950);
hole = data_hole;
hole_etendu = data_hole_etendu;
clear data_hole data_hole_etendu;

% DEBUT CA_extraire_contour
Nx=4958;
Ny=Nx;
Nwound=300;%round(Nx/2);          % # of points defining the wound  

Lx=1; % dcm => 10 cm => 100 mm
Ly=1; % dcm => 10 cm => 100 mm

Ansys_scaling=100;% 1 USI <=> 100mm

hx=Lx/(Nx-1); hy=Ly/(Ny-1);

N_wound=Nwound-1;

Area_wound=sum(sum(hole));

% extraction of the wound level set:   
if Area_wound>0,
	clear C HC
	[C HC]=contour(hole,0.5);
	dim_C=size(C);Nc=dim_C(2)-1;
	clear XC YC xy
	
	XC=C(2,1:Nc+1); YC=C(1,1:Nc+1);
	XC(1,Nc+2)=XC(1,1); YC(1,Nc+2)=YC(1,1);
	for j=1:Nc+2,
		xy(j)=j-1;
	end
	for j=1:N_wound,
		newxy(j)=1+(Nc-1)*(j-1)/(N_wound-1);
	end
	%save truc % for debug
	%clear C;
	XC=(XC-1)*hx-Lx/2; YC=(YC-1)*hy-Ly/2;
	xC=interp1(xy,XC,newxy);
	yC=interp1(xy,YC,newxy);
	
	XC0=mean(xC); YC0=mean(yC);
	
	clear C HC
	[C HC]=contour(hole_etendu,0.5);
	dim_C=size(C);Nc=dim_C(2)-1;
	clear XC YC xy
	
	XC=C(2,1:Nc+1); YC=C(1,1:Nc+1);
	XC(1,Nc+2)=XC(1,1); YC(1,Nc+2)=YC(1,1);
	for j=1:Nc+2,
		xy(j)=j-1;
	end
	for j=1:N_wound,
		newxy(j)=1+(Nc-1)*(j-1)/(N_wound-1);
	end
	
	% clear C;
	XC=(XC-1)*hx-Lx/2; YC=(YC-1)*hy-Ly/2;
	xC_E=interp1(xy,XC,newxy);
	yC_E=interp1(xy,YC,newxy);
	
end

clear hole_etendu

if Area_wound>0,
	for j=1:N_wound-1
		X_ansys(j)=Ansys_scaling*(xC(j)); % *Ansys_scaling sCALING FROM DCM TO MM
		Y_ansys(j)=Ansys_scaling*(yC(j)); % and recentering completly the wound to a fix coord system
	end
	for j=1:N_wound-1
		X_ansys_E(j)=Ansys_scaling*(xC_E(j)); % *Ansys_scaling sCALING FROM DCM TO MM
		Y_ansys_E(j)=Ansys_scaling*(yC_E(j)); % and recentering completly the wound to a fix coord system
	end
else
	X_ansys=zeros(N_wound-1,1);  Y_ansys=X_ansys; X_ansys_E=X_ansys; Y_ansys_E=X_ansys;
end

f=fopen('/home/remi/Desktop/CA/code/Area_wound.txt', 'w');

fprintf(f, '%d', Area_wound);

fclose(f);

f=fopen('/home/remi/Desktop/CA/code/trace_E.out', 'w');

for i=1:N_wound-1
	fprintf(f, '%f %f\n', X_ansys_E(i), Y_ansys_E(i));
end

fclose(f);

f=fopen('/home/remi/Desktop/CA/code/X_ansys.txt', 'w');

for i=1:N_wound-1
	fprintf(f, '%f ', X_ansys(i));
end

fclose(f);

f=fopen('/home/remi/Desktop/CA/code/Y_ansys.txt', 'w');

for i=1:N_wound-1
	fprintf(f, '%f ', Y_ansys(i));
end

fclose(f);

f=fopen('/home/remi/Desktop/CA/code/X_ansys_E.txt', 'w');

for i=1:N_wound-1
	fprintf(f, '%f ', X_ansys_E(i));
end

fclose(f);

f=fopen('/home/remi/Desktop/CA/code/Y_ansys_E.txt', 'w');

for i=1:N_wound-1
	fprintf(f, '%f ', Y_ansys_E(i));
end

fclose(f);


%% FIN CA_extraire_contour.m

clear all;

system('rm hole*.txt 2>/dev/null');
