function cercle5(rayon1,rayon2,nbpoints1,nbpoints2)
%angle=pi/12;
angle1=2*pi/nbpoints1;
angle2=2*pi/nbpoints2;
nbpoints=2*nbpoints1+nbpoints2;
%nbpoints=24;

%coordonnées des noeuds

for i=0:2:2*nbpoints1-2
  X(:,i+1)=[i+1;rayon1*cos(i*angle1/2);rayon1*sin(i*angle1/2);1];
end
for i=1:2:2*nbpoints1-1
  X(:,i+1)=[i+1;rayon1*(cos((i-1)*angle1/2)+cos((i+1)*angle1/2))/2;rayon1*(sin((i-1)*angle1/2)+sin((i+1)*angle1/2))/2;1];
end
  for i=0:nbpoints2-1
  X(:,2*nbpoints1+i+1)=[2*nbpoints1+i+1;rayon2*cos(i*angle2);rayon2*sin(i*angle2);3];
end


%aretes
Y(:,1)=[1;1;2;1];
Y(:,2*nbpoints1)=[2*nbpoints1;2*nbpoints1;1;1];
Y(:,2*nbpoints1+1)=[2*nbpoints1+1;2*nbpoints1+1;2*nbpoints1+2;3];
Y(:,nbpoints)=[nbpoints;nbpoints;2*nbpoints1+1;3];

for i=2:2*nbpoints1-1
  Y(:,i)=[i;i;i+1;1];
end
  for i=2*nbpoints1+2:nbpoints-1
  Y(:,i)=[i;i;i+1;3];
end

debut=[nbpoints,2,0,2];
milieu=[nbpoints,1];
fin2=[1,0,(rayon1+rayon2)/2,1];
fin3=[1,0,0,2];

%ecriture du fichier cercle.poly
fid=fopen(['cercle3_julien.poly'],'wt');
fprintf(fid,'%.16g %.16g %.16g %.16g\n',debut);
fprintf(fid, '%.16g %22.15E %22.15E %.16g\n',X);
fprintf(fid,'\n');
fprintf(fid,'%.16g %.16g\n',milieu);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',Y);
fprintf(fid, '\n');
fprintf(fid,'%.16g\n',0);
fprintf(fid,'\n');
fprintf(fid,'%.16g\n',2);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin2);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin3);
fclose(fid);


