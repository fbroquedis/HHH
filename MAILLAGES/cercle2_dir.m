function cercle2_dir(rayon1,rayon2,nbpoints1,nbpoints2,a1,a2)
%angle=pi/12;
angle1=2*pi/nbpoints1;
angle2=2*pi/nbpoints2;
nbpoints=nbpoints1+nbpoints2;
%nbpoints=24;

%coordonnées des noeuds

for i=0:nbpoints1-1
  X(:,i+1)=[i+1;rayon1*cos(i*angle1);rayon1*sin(i*angle1);1];
end
  for i=0:nbpoints2-1
  X(:,nbpoints1+i+1)=[nbpoints1+i+1;rayon2*cos(i*angle2);rayon2*sin(i*angle2);2];
end


%aretes
Y(:,1)=[1;1;2;1];
Y(:,nbpoints1)=[nbpoints1;nbpoints1;1;1];
Y(:,nbpoints1+1)=[nbpoints1+1;nbpoints1+1;nbpoints1+2;2];
Y(:,nbpoints)=[nbpoints;nbpoints;nbpoints1+1;2];

for i=2:nbpoints1-1
  Y(:,i)=[i;i;i+1;1];
end
  for i=nbpoints1+2:nbpoints-1
  Y(:,i)=[i;i;i+1;2];
end

debut=[nbpoints,2,0,2];
milieu=[nbpoints,1];
%fin2=[1,0,(rayon1+rayon2)/2,1,a1];
%fin3=[1,0,0,2,a2];
fin2=[1,0,(rayon1+rayon2)/2,1];
fin3=[1,0,0,2];

%ecriture du fichier cercle.poly
fid=fopen(['cercle_dir.poly'],'wt');
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


