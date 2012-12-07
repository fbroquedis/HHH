function cercle(rayon1,rayon2,angle1)
%angle=pi/12;
angle=pi/angle1
nbpoints=2*angle1;
%nbpoints=24;

%coordonnées des noeuds

for i=0:nbpoints-1
  X(:,i+1)=[i+1;rayon1*cos(i*angle);rayon1*sin(i*angle);1];
  X(:,nbpoints+i+1)=[nbpoints+i+1;rayon2*cos(i*angle);rayon2*sin(i*angle);3];
end


%aretes
Y(:,1)=[1;1;2;3];
Y(:,nbpoints)=[nbpoints;nbpoints;1;3];
Y(:,nbpoints+1)=[nbpoints+1;nbpoints+1;nbpoints+2;1];
Y(:,2*nbpoints)=[2*nbpoints;2*nbpoints;nbpoints+1;1];
for i=2:nbpoints-1
  Y(:,i)=[i;i;i+1;3];
  Y(:,nbpoints+i)=[nbpoints+i;nbpoints+i;nbpoints+i+1;1];
end

debut=[2*nbpoints,2,0,2];
milieu=[2*nbpoints,1];
fin2=[1,0,(rayon1+rayon2)/2,1];
fin3=[1,0,0,2];

%ecriture du fichier cercle.poly
fid=fopen(['cercle2.poly'],'wt');
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


