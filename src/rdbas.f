      subroutine rdbad(idebug,dfree,istats,ityp)
      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxonh=100)
      parameter (mxexpo=100)
      parameter (mxshl=20)
      logical fndshl,eledon
      character*3 cnam
      character*2 ele
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp,dfree,rdbuf
      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      common /coord / xyz(3,numatm)
      common /oniomh/ ioni,nion,ionih(mxonh),natonh(mxonh),
     &                iomap(numatm),icntat(mxonh),fct(mxonh),
     &                xyzi(3,mxonh)
      common /slagau/ ihasd,isgau,ido5d,ido7f,ido9g,ihasg
      common /rdwr/  iun1,iun2,iun3,iun4,iun5
      integer*2 ityp
      equivalence (nshell,i)
      character*137 line,str
      common /curlin/ line
      integer getlin,getele
      logical gnreal
      dimension r(4),ityp(*),expo(mxexpo),coefo(mxexpo)
      dimension ishll(numatm,mxshl),ieler(numatm)

      do i=1,natoms
         do j=1,mxshl
            ishll(i,j)  = -1
         end do
         ieler(i) = -1
      end do

      istats = 1
      ig03 = 0
      ieleol = -1

      if (idebug.eq.1)
     &    call inferr('looking for gaussian basis-set',0)

1030  format(4d18.10)
1040  format(4x,4d18.10)
      if (dfree.ne.0) then
         if (dfree.eq.1) then
            call searchu(line,'[GTO]',istat)
         elseif (dfree.eq.2) then
            call srchmf(line,'$BASIS',istat)
         elseif (dfree.eq.3) then
            rewind(iun2)
            call search(line,'bas_input',istat)
            call search(line,'basis',istat)
            if (index(line,"cartesian").ne.0) then
               ido5d = 0
               ido7f = 0
               ido9g = 0
            endif
            if (index(line,"spherical").ne.0) then
               ido5d = 1
               ido7f = 1
               ido9g = 1
               print*,"mulliken charges of spherical F and G ",
     &                "dont add up !!"
              
            endif
         endif
      else
         call seardu(line,'Basis set in the form',
     &                    'AO basis set:',istat)
         if (istat.eq.0) then
             call inferr(
     &       'Use keyword gfinput in gaussian inputfile',1)
             goto 1000
         endif
         if (index(line,'AO basis set').ne.0) ig03 = 1
         if (index(line,'AO basis set:').ne.0) ig03 = 2
      endif

      i = 0
      k = 0
      if (ioni.eq.1) then
         do m=1,numatm
            iomap(m) = 0
            ityp(m) = 20000
         end do
      endif
      mm = 1
      mmdf = 1
      jann = 0
      jmode = 1
      ns = 0

10    if (dfree.ne.0) then
         if (getlin(0).eq.1) idum = 0
      else
         call nxtlin(line,jstat)
         if (jstat.eq.1.or.jstat.eq.2) goto 1000
      endif

      if (dfree.ne.0) then

         ktype = nxtwrd(str,nstr,itype,rtype)

         if (ktype.eq.1) then
            if (nstr.ge.1) then
               if (index(str,'[').ne.0) return
               if (dfree.eq.2) then
                  if (icdex(str,'$END').ne.0) return
               endif
               if (str(1:1).eq.'#'.or.str(1:1).eq.'*') then
                  if (dfree.eq.2) jmode = 1
                  goto 10
               endif
               if (jmode.eq.1.and.dfree.eq.2) then
                  jann = jann + 1
                  jmode = 0
                  goto 10
               endif
            endif
         elseif (ktype.eq.2) then
            jann = itype
            goto 10
         elseif (ktype.eq.0) then
            goto 10
         endif

      else

         if (ig03.eq.2.and.line(2:5).ne.'Atom') return

         if (i.eq.0.or.line(2:5).eq.'****'.or.line(2:5).eq.'Atom') then

            if (line(2:5).eq.'Atom') then
               ig03 = 2

c example: Atom B3       Shell     8 SP   2    bf   20 -    23

               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)

c break up string and extract jann

               do kk=1,nstr
                  ik = ichar(str(kk:kk))
                  if (ik.lt.48.or.ik.gt.57) str(kk:kk) = ' '
               end do

               do kk=1,nstr
                  if (str(kk:kk).ne.' ') istrt = kk
               end do

               jann = reada(str,1,nstr)
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)
               ktype = nxtwrd(str,nstr,itype,rtype)

               if (ktype.eq.1) then
                  cnam = '   '
                  if (nstr.gt.3) nstr = 3
                  cnam(1:nstr) = str(1:nstr)
                  call tocap(cnam,3)
               else
                  goto 1000
               endif

               ktype = nxtwrd(str,nstr,itype,rtype)
               if (ktype.eq.2) then
                  nsh = itype
               else
                  goto 1000
               endif

            else

               if (line(2:5).eq.'****') then
                   call nxtlin(line,jstat)
                   if (jstat.eq.1.or.jstat.eq.2) goto 1000
                   if (linlen(line).eq.0) return
               endif

               if (line(7:7).ne.' ') then
                  read(line,'(i7)',err=1000) jann
               elseif (line(4:4).ne.' ') then
                  read(line,'(i4)',err=1000) jann
               else
                  read(line,'(i3)',err=1000) jann
               endif

               if (ioni.eq.1) then
                  k = k + 1
                  iomap(k) = jann
                  nat(k) = nat(jann)
                  ityp(jann) = ityp(jann) - 20000

                  do ll=1,3
                     xyz(ll,k) = xyz(ll,jann)
                  end do

                  do jj=1,nion
                      if (jann.eq.ionih(jj)) then
                         nat(k) = natonh(jj)
                         do ll=1,3
                            xyz(ll,k) = xyzi(ll,jj)
                         end do
                      endif
                  end do

               endif

               call nxtlin(line,jstat)
               if (jstat.eq.1.or.jstat.eq.2) goto 1000
            endif
         endif
      endif

      if (dfree.eq.0.and.ig03.eq.1) then
         ktype = nxtwrd(str,nstr,itype,rtype)
      endif

      if (dfree.ne.0.or.ig03.eq.1) then
      
c         if (dfree.ne.3.or.ig03.eq.1) then
         if (ig03.eq.1) then
            cnam = '   '
            if (nstr.gt.3) nstr = 3
            cnam(1:nstr) = str(1:nstr)
            call tocap(cnam,3)
            ktype = nxtwrd(str,nstr,itype,rtype)
            if (ktype.eq.2) then
               nsh = itype
            endif
         else if (dfree.eq.3) then
            call bckfil
         endif

      else
         if (ig03.eq.0) then
            read(line,'(1x,a3,1x,i2)',err=1000) cnam, nsh
         endif
      endif

      i = i + 1

      if (dfree.eq.3) then

c ish is the overall counter for shells
c ns  is the element counter for shells

         ish = 0

c read in element basis set

         do j=1,natoms

            ns = 0

            if (ieler(j).eq.-1) then

               call rewmf
               call search(line,'bas_input',istat)
               call search(line,'basis',istat)

               nprim = 0

               do while (nprim.ne.-2)
                  nprim = rdbuf(nat(j),expo,coefo,cnam)
                  if (nprim.eq.-1) goto 1000
                  if (j.eq.natoms) then
                     if (nprim.eq.-2) then
                        nshell = ish
c                        call basprt(6,.false.,.false.)
                        return
                     endif
                  else
                     if (nprim.eq.-2) goto 30
                  endif

                  ish = ish + 1
                  ns  = ns + 1

c store shell in common block /b/:

                  call pttshl(ish,nprim,j,ns,ishll,mm,mmdf,
     &                        expo,coefo,cnam)

               end do

               ieler(j) = ish

            else

c copy already read in basis for element into the 

               eledon = fndshl(j,ish,ishll,mm,mmdf)

            endif

30          continue
            
         end do

100      continue

      else

         if (i.gt.numprm) print*,"exceded parameter numprm=",numprm
         call settc(cnam,shellt(i),shellc(i))

         shelln(i)  = nsh
         shella(i)  = mm
         shladf(i)  = mmdf

         call rdshl(idebug,dfree,ioni,jann,ig03,i,k,nsh,mm,
     &                    mmdf,istats,cnam)
         if (istats.eq.0) goto 1000

      endif

      goto 10

1000  call inferr('error while reading Basis set!',1)
      istats = 0
      return
      end

      subroutine pttshl(ish,nprim,jann,ns,ishll,mm,mmdf,expo,coefo,cnam)
      implicit double precision (a-h,o-z)

      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxexpo=100)
      parameter (mxshl=20)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp
      character*3 cnam

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp

      common /coord / xyz(3,numatm)
      dimension expo(mxexpo),coefo(mxexpo)
      dimension ishll(numatm,mxshl)

      call settc(cnam,shellt(ish),shellc(ish))

      ishll(jann,ns)  = ish

      shelln(ish)  = nprim
      shella(ish)  = mm
      shladf(ish)  = mmdf

      jan(ish)     = jann
      gx(ish)      = xyz(1,jann)
      gy(ish)      = xyz(2,jann)
      gz(ish)      = xyz(3,jann)

      call putshl(ish,nprim,mm,mmdf,expo,coefo)

      return
      end

      integer function iseqa(isq,ns,ishll)

c find the first occurance (sequence number) of element number isq 
c in coordinate set, that has not been initialised

      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxshl=20)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension ishll(numatm,mxshl)

      iseqa = -1

      do j=1,natoms
          if (isq.eq.nat(j).and.ishll(j,ns).eq.-1) then
             iseqa = j
             return
          endif
      end do

      return
      end

      logical function fndshl(jann,ishell,ishll,mm,mmdf)

c Finds element thats has previously been read and put into
c exponent and contraction coefficient arrays by putshl
c and calls copy shell (cpshl) for other instances of this element

      implicit double precision (a-h,o-z)
      parameter (numatm=2000)
      parameter (mxshl=20)
      common /moldat/ natoms, norbs, nelecs,nat(numatm)
      dimension ishll(numatm,mxshl)

      fndshl = .false.

      jold = -1

      do j=1,numatm
         if (nat(jann).eq.nat(j).and.ishll(j,1).ne.-1) then
            jold = j
         endif
      end do

      if (jold.ne.-1) then

          do k=1,mxshl

             if (ishll(jold,k).ne.-1) then
                 fndshl = .true.
                 call cpshl(jann,ishell,ishll(jold,k),mm,mmdf)
                 ishll(jann,k) = ishell
                 ishell = ishell + 1
             endif

          end do

      else

          print*,"presumed read in shell does not exist!"

      endif

      return
      end

      integer function getele(ele)
      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      character*2 ele,tocapf,elemnt
      common /elem/elemnt(mxel)

      getele = 0
      do i=1,mxel
         if (tocapf(ele).eq.tocapf(elemnt(i))) then
             getele = i
         endif
      end do 

      return
      end

      integer function rdbuf(iele,expo,coefo,cnam)

c Reads a NWCHEM element (iele) shell and returns the name of the shell in cnam,
c the coefficients and exponents
c in coefo and expo arrays and returns the number of primitives as rdbuf

      implicit double precision (a-h,o-z)
      parameter (mxel=100)
      parameter (mxexpo=100)
      character*137 line,tlin,str
      character*3 cnam
      character*2 ele,tocapf
      common /curlin/ line
      character*2 elemnt
      common /elem/elemnt(mxel)
      integer getlin,getele
      logical gnreal
      dimension r(2),expo(mxexpo),coefo(mxexpo)

      ifndel = 0
      mm = 0

      do while (.true.)

20      if (getlin(0).eq.1) then

            tlin = line
            ktype = nxtwrd(str,nstr,itype,rtype)

c line starts with a string: read in element and shell name

            if (ktype.eq.1) then

               if (nstr.eq.3) then
                  if (str(1:3).eq.'end') then
                     if (ifndel.eq.1) goto 50
                     goto 200
                  endif
               endif

               if (ifndel.eq.1) goto 50

c convert element string to integer

               ele(1:2) = '  '
               if (nstr.eq.1) ele = ' '//str(1:1)
               if (nstr.eq.2) ele = str(1:2)

               iel = getele(ele)

c unknown element:

               if (iel.eq.0) goto 100


               if (ifndel.ne.1) then
c element is not the desired element:
                  if (iel.ne.iele) then
c                     mm = 0
                     goto 20
                  endif
               else
c done
                  call bckfil
                  goto 200
               endif

c found the desired element

               ifndel = 1
               
               ktype = nxtwrd(str,nstr,itype,rtype)

               if (ktype.eq.1) then
                  if (nstr.eq.3) then
                     cnam = str(1:nstr)
                  else if (nstr.eq.2) then
                     cnam = str(1:2)//' '
                  else if (nstr.eq.1) then
                     cnam = str(1:1)//'  '
                  else
                     goto 100
                  endif
                  call tocap(cnam,nstr)
                  mm = 0
               else
                  goto 100
               endif

            else

c line starts not with a string: read in exponent and coefficient

               line = tlin
               if (gnreal(r,2,.false.)) then
                  mm = mm + 1
                  expo(mm)   = r(1)
                  coefo(mm)  = r(2)
               else
                  goto 100
               endif
            endif

         else
            goto 100
         endif

      end do

50    rdbuf = mm
      call bckfil
      return

c     error
100   rdbuf = -1
      return

c     end
200   rdbuf = -2
      return
      end

      subroutine rdshl(idebug,dfree,ioni,jann,ig03,ishell,ksh,nshl,mm,
     &                 mmdf,istats,cnam)
      implicit double precision (a-h,o-z)

      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      logical gnreal
      character*3 cnam
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp,dfree

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /coord / xyz(3,numatm)
      character*137 line,str
      common /curlin/ line
      dimension r(4)

1030  format(4d18.10)
1040  format(4x,4d18.10)

      istats = 1

      if (ioni.eq.1) then
         jan(ishell)     = ksh
         gx(ishell)      = xyz(1,ksh)
         gy(ishell)      = xyz(2,ksh)
         gz(ishell)      = xyz(3,ksh)
      else
         jan(ishell)     = jann
         gx(ishell)      = xyz(1,jann)
         gy(ishell)      = xyz(2,jann)
         gz(ishell)      = xyz(3,jann)
      endif

      do j = 1, nshl

           if (shellt(ishell).eq.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c1(mm)
                  else
                     read(line,1030,err=1000) exx(mm), c1(mm)
                  endif
               endif
           endif

           if (shellt(ishell).eq.1.and.shellc(ishell).eq.1) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c2(mm)  = r(2)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c2(mm)
                  else
                     read(line,1030,err=1000) exx(mm), c2(mm)
                  endif
               endif
           endif

           if (shellt(ishell).eq.1.and.shellc(ishell).ne.1) then
               if (dfree.ne.0) then
                  if (gnreal(r,3,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                     c2(mm)  = r(3)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c1(mm), c2(mm)
                  else
                     read(line,1030,err=1000) exx(mm), c1(mm), c2(mm)
                  endif
               endif
           endif

           if (shellt(ishell).eq.2.and.shellc(ishell).eq.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,4,.true.)) then
                     exx(mm) = r(1)
                     c1(mm)  = r(2)
                     c2(mm)  = r(3)
                     c3(mmdf)= r(4)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) 
     &                        exx(mm), c1(mm), c2(mm), c3(mmdf)
                  else
                     read(line,1030,err=1000) 
     &                        exx(mm), c1(mm), c2(mm), c3(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif

           if (shellt(ishell).eq.2.and.shellc(ishell).ne.0) then
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c3(mmdf)= r(2)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c3(mmdf)
                  else
                     read(line,1030,err=1000) exx(mm), c3(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif

           if (shellt(ishell).eq.3) then 
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c4(mmdf)= r(2)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c4(mmdf)
                  else
                     read(line,1030,err=1000) exx(mm), c4(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif

           if (shellt(ishell).eq.4) then 
               if (dfree.ne.0) then
                  if (gnreal(r,2,.true.)) then
                     exx(mm) = r(1)
                     c5(mmdf)= r(2)
                  endif
               else
                  call nxtlin(line,jstat)
                  if (jstat.eq.1.or.jstat.eq.2) goto 1000

                  if (ig03.gt.0) then
                     read(line,1040,err=1000) exx(mm), c5(mmdf)
                  else
                     read(line,1030,err=1000) exx(mm), c5(mmdf)
                  endif
               endif
               mmdf = mmdf + 1
           endif

           mm = mm + 1
      end do

      return

1000  istats = 0
      return
      end

      subroutine cpshl(jann,ishell,jshell,mm,mmdf)

c copies previously read shell (rdbuf) and parsed by putshl
c to same shell and element combination for another center

      implicit double precision (a-h,o-z)

      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      logical gnreal
      character*3 cnam
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      common /coord / xyz(3,numatm)

      jan(ishell)    = jann
      gx(ishell)     = xyz(1,jann)
      gy(ishell)     = xyz(2,jann)
      gz(ishell)     = xyz(3,jann)

      shelln(ishell) = shelln(jshell)
      shellt(ishell) = shellt(jshell)
      shellc(ishell) = shellc(jshell)
      shella(ishell) = mm
      shladf(ishell) = mmdf
      mmo            = shella(jshell) 
      mmdfo          = shladf(jshell) 

      isht = shellt(jshell)
      ishc = shellc(jshell)

      do j = 1, shelln(jshell)

           if (isht.eq.0) then
               exx(mm) = exx(mmo)
               c1(mm)  = c1(mmo)
           endif

           if (isht.eq.1.and.ishc.eq.1) then
               exx(mm) = exx(mmo)
               c2(mm)  = c2(mmo)
           endif

           if (isht.eq.1.and.ishc.ne.1) then
               exx(mm) = exx(mmo)
               c1(mm)  = c1(mmo)
               c2(mm)  = c2(mmo)
           endif

           if (isht.eq.2.and.ishc.eq.0) then
               exx(mm) = exx(mmo)
               c1(mm)  = c1(mmo)
               c2(mm)  = c2(mmo)
               c3(mmdf) = c3(mmdfo)

               mmdf  = mmdf  + 1
               mmdfo = mmdfo + 1
           endif

           if (isht.eq.2.and.ishc.ne.0) then
               exx(mm) = exx(mmo)
               c3(mmdf) = c3(mmdfo)
               mmdf  = mmdf  + 1
               mmdfo = mmdfo + 1
           endif

           if (isht.eq.3) then 
               exx(mm) = exx(mmo)
               c4(mmdf) = c4(mmdfo)
               mmdf  = mmdf  + 1
               mmdfo = mmdfo + 1
           endif

           if (isht.eq.4) then 
               exx(mm) = exx(mmo)
               c5(mmdf) = c5(mmdfo)
               mmdf  = mmdf  + 1
               mmdfo = mmdfo + 1
           endif

           mm  = mm  + 1
           mmo = mmo + 1
      end do

      return

      end

      subroutine putshl(ishell,nshl,mm,mmdf,expo,coefo)

c putshl: takes temporary exponents and contraction coefficients
c         read by rdbuf, and puts them in the permanent arrays
c         exx and c1 to c5 in common block /b/

      implicit double precision (a-h,o-z)

      parameter (numatm=2000)
      parameter (numprm=1600)
      parameter (numcex=numprm*3)
      parameter (mxexpo=100)
      integer jan,shella,shelln,shellt,shellc,shladf,aos,
     &        nshell,maxtyp

      common/b/exx(numcex),
     &         c1(numcex),c2(numcex),c3(numprm),c4(numprm),c5(numprm),
     &         shladf(numprm),gx(numprm),gy(numprm),gz(numprm),
     &         jan(numprm),shella(numprm),shelln(numprm),shellt(numprm),
     &         shellc(numprm),aos(numprm),nshell,maxtyp
      dimension expo(mxexpo),coefo(mxexpo)

      mmo = 1

      do j = 1, nshl

           if (shellt(ishell).eq.0) then

               exx(mm) = expo(mmo)
               c1(mm)  = coefo(mmo)
               mm  = mm  + 1
           endif

           if (shellt(ishell).eq.1) then
               exx(mm) = expo(mmo)
               c2(mm)  = coefo(mmo)
               mm  = mm  + 1
           endif

           if (shellt(ishell).eq.2) then
               exx(mm)  = expo(mmo)
               c3(mmdf) = coefo(mmo)
               mm  = mm  + 1
               mmdf  = mmdf  + 1
           endif

           if (shellt(ishell).eq.3) then 
               exx(mm)  = expo(mmo)
               c4(mmdf) = coefo(mmo)
               mm  = mm  + 1
               mmdf  = mmdf  + 1
           endif

           if (shellt(ishell).eq.4) then 
               exx(mm)  = expo(mmo)
               c5(mmdf) = coefo(mmo)
               mm  = mm  + 1
               mmdf  = mmdf  + 1
           endif

           mmo = mmo + 1
      end do

      return

      end

