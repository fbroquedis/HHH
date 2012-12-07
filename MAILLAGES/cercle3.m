function cercle3(rayon1,rayon2,nbpoints1,nbpoints2)
%angle=pi/12;
angle1=2*pi/nbpoints1;
angle2=2*pi/nbpoints2;
nbpoints=3*nbpoints1+nbpoints2;
%nbpoints=24;
delta=2*angle1*rayon1;
%coordonnées des noeuds
rayon3=rayon1-delta;
rayon4=rayon1+delta;
for i=0:nbpoints1-1
    toto=rand;
  X(:,i+1)=[i+1;rayon1*cos((i+0.01*toto)*angle1);rayon1*sin((i+0.01*toto)*angle1);1];
end
for i=0:nbpoints1-1
    toto=rand;
  X(:,nbpoints1+i+1)=[nbpoints1+i+1;rayon3*cos((i+0.01*toto)*angle1);rayon3*sin((i+0.01*toto)*angle1);1];
end
for i=0:nbpoints1-1
    toto=rand;
  X(:,2*nbpoints1+i+1)=[2*nbpoints1+i+1;rayon4*cos((i+0.01*toto)*angle1);rayon4*sin((i+0.01*toto)*angle1);1];
end
  for i=0:nbpoints2-1
  X(:,3*nbpoints1+i+1)=[3*nbpoints1+i+1;rayon2*cos(i*angle2);rayon2*sin(i*angle2);3];
end


%aretes
Y(:,1)=[1;1;2;1];
Y(:,nbpoints1)=[nbpoints1;nbpoints1;1;1];
Y(:,nbpoints1+1)=[nbpoints1+1;nbpoints1+1;nbpoints1+2;1];
Y(:,2*nbpoints1)=[2*nbpoints1;2*nbpoints1;nbpoints1+1;1];
Y(:,2*nbpoints1+1)=[2*nbpoints1+1;2*nbpoints1+1;2*nbpoints1+2;1];
Y(:,3*nbpoints1)=[3*nbpoints1;3*nbpoints1;2*nbpoints1+1;1];
Y(:,3*nbpoints1+1)=[3*nbpoints1+1;3*nbpoints1+1;3*nbpoints1+2;3];
Y(:,nbpoints)=[nbpoints;nbpoints;3*nbpoints1+1;3];

for i=2:nbpoints1-1
  Y(:,i)=[i;i;i+1;1];
end
for i=nbpoints1+2:2*nbpoints1-1
  Y(:,i)=[i;i;i+1;1];
end
for i=2*nbpoints1+2:3*nbpoints1-1
  Y(:,i)=[i;i;i+1;1];
end
  for i=3*nbpoints1+2:nbpoints-1
  Y(:,i)=[i;i;i+1;3];
end

debut=[nbpoints,2,0,2];
milieu=[nbpoints,1];
fin2=[1,0,(rayon1+rayon2)/2,1];
fin3=[1,0,0,2];
fin1=[1,0,rayon1+delta/2,1];
fin4=[1,0,rayon1-delta/2,2];

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
fprintf(fid,'%.16g\n',4);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin2);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin3);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin1);
fprintf(fid,'%.16g %.16g %.16g %.16g\n',fin4);
fclose(fid);


