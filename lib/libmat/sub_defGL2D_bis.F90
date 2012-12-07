SUBROUTINE sub_defGL2D_bis
  Use m_mat
  implicit none
  integer :: i,indice
  real*16 :: aa,bb
  NGL2D=91
  
  allocate(wGL2D(NGL2D))
  allocate(PtGL2D(NGL2D,2))
  indice=1
  WGL2D(indice)=  0.5Q0*0.0275622569528764809669070448245143Q0
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0220602154134885011913507340331164Q0
  END DO
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0234600159386714884930134449523000Q0
  END DO
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0003268895950471905462145575015465Q0
  END DO
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0032653194629399682343353040958667Q0
  END DO
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0117564629154127977043079692133821Q0
  END DO
  DO I=1,3
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0117807684199115168455575790986761Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0022688108188011408053357043343043Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0025960109644363200606737836654882Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0046345297858718602123478905615969Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0047943360545488579348574487199119Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0057124788367236115672506383429634Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0058658276043221216369557987000023Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0094137630590915875898182685203471Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0134149437966564249100220266108931Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0157169180920832459435000011378462Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0168636830144369045916509638861999Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0213900270853200983778322980803590Q0
  END DO

  DO I=1,6
     indice=indice+1
     WGL2D(indice)=0.5Q0*0.0230767921894926813678808755218915Q0
  END DO



  indice=1
  PtGL2D(indice,1)=1/3Q0
  PtGL2D(indice,2)=1/3Q0

  aa=.2009352770650852798729618515641637Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa


  aa=0.4376591659619271797318338441880541Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa

  aa=.0034339564905961768509599122096049Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa

  aa=0.0466434847753067534951762404321419Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa

  aa=0.3864222517630714909403520241677264Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa

  aa=0.0954354711085309101085716810414760Q0
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-2Q0*aa
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-2Q0*aa
  PtGL2D(indice,2) = aa


  aa=0.9555138033504563605013147251467712Q0
  bb=0.0357186278731633582380416089754387Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa


  aa=0.8866388134288682261249005746914376Q0
  bb=0.1081432249156462115273886110463127Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa


  aa=0.7842628458804341542966439903981954Q0
  bb=0.2074644495998764568243804295157274Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.8829239550502000327113489873168897Q0
  bb=0.0856847087203169400000000000000100Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.6689919644410772404913224832098946Q0
  bb=0.3214940030142888168816832126834860Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.5520721210355609641571609652527788Q0
  bb=0.4379422187933413835523680769629170Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.7975929655965685676293142232957258Q0
  bb=0.1619164530635778567510067702038591Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.6775147151197714846349911663441326Q0
  bb=0.2745047674019949038590029729073332Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.5429974155890916053311361168391934Q0
  bb=0.4053359980750069279498908953763256Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.7054599055699685616588563415406017Q0
  bb=0.1877376806564353427728167439451200Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa

  aa=0.5748005730665084622159824505498500Q0
  bb=0.3056968347660551665127925566498432Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa


  aa=0.4717788085046148166039770401349242Q0
  bb=0.3121444668708908816708046058155764Q0
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = aa
  indice=indice+1
  PtGL2D(indice,1) = bb
  PtGL2D(indice,2) = 1Q0-aa-bb
  indice=indice+1
  PtGL2D(indice,1) = aa
  PtGL2D(indice,2) = bb
  indice=indice+1
  PtGL2D(indice,1) = 1Q0-aa-bb
  PtGL2D(indice,2) = aa



END SUBROUTINE sub_defGL2D_bis
