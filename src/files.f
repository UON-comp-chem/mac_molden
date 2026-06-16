      subroutine files(oempty,idocub,idogrd)
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      common /rdwr/   iun1,iun2,iun3,iun4,iun5
      character curdir*1024
      common /cwd/ curdir
      common /cmdfil/ icmdf
      logical oempty
      character fniun*256
c      character temps*80
      common /fnunit/ fniun
      common /keywrd/ keywrd,keyori

c
c----- process FILE or CUBE directive --------------------
c
      idocub = 0
      idogrd = 0
      i = index(keywrd,'FILE=')
      if (i.eq.0) then
         i = index(keywrd,'GAUCUB')
         if (i.ne.0) then
            idocub = 1
            icmdf = 1
         endif
      endif
      if (i.eq.0) then
         i = index(keywrd,'RDBAS')
         if (i.ne.0) then
            idogrd = 1
            icmdf = 1
         endif
         return
      endif
      if (i.ne.0) then
         ind = i + 4
         i   = index(keywrd(ind:320),'=')
         do j=i+ind,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j
         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1

         if (index(keyori(ind1:ind2),'/').eq.1) then
            open(unit=iun2,form='formatted',file=keyori(ind1:ind2),
     &        status='old',err=2200)
            fniun = keyori(ind1:ind2)
         else
            ld = linlen(curdir)
            open(unit=iun2,form='formatted',
     &        file=curdir(1:ld)//keyori(ind1:ind2),
     &        status='old',err=2200)
            fniun = curdir(1:ld)//keyori(ind1:ind2)
         endif
         do j=ind1,ind2
            keyori(j:j) = ' '
            keywrd(j:j) = ' '
         end do
         icmdf = 1
      else
         call inferr('Keyword FILE not supplied !!',1)
         oempty = .true.
         icmdf = 0
      endif

      return

2200  write(iun3,*) 
     &   'File '//keyori(ind1:ind2)//' does not exist !!!!!!!'
c      call inferr(temps,1)
      oempty = .true.
      icmdf = 0
      return
      end

      subroutine filkey(fname,iun,iform)
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      character fname*(*)
      character fniun*256
      common /fnunit/ fniun
      common /keywrd/ keywrd,keyori
      character curdir*1024
      common /cwd/ curdir
      common /cmdfil/ icmdf

      ln = len(fname)
      i = index(keywrd,fname(1:ln))

      if (i.ne.0) then

         ind = i + ln
         i   = index(keywrd(ind:320),'=')

         do j=i+ind,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j

         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1

         if (iform.eq.1) then
            if (index(keyori(ind1:ind2),'/').eq.1) then
              open(unit=iun,form='formatted',
     &           file=keyori(ind1:ind2),status='new',err=2200)
              fniun = keyori(ind1:ind2)
            else
              if (icmdf.eq.1) then
                ld = linlen(curdir)
                open(unit=iun,form='formatted',
     &           file=curdir(1:ld)//keyori(ind1:ind2),
     &           status='new',err=2200)
                 fniun = curdir(1:ld)//keyori(ind1:ind2)
              else
                open(unit=iun,form='formatted',
     &           file=keyori(ind1:ind2),
     &           status='new',err=2200)
                 fniun = keyori(ind1:ind2)
              endif
            endif
         else
            if (index(keyori(ind1:ind2),'/').eq.1) then
                open(unit=iun,form='unformatted',
     &           file=keyori(ind1:ind2),status='new',err=2200)
                fniun = keyori(ind1:ind2)
            else
              if (icmdf.eq.1) then
                ld = linlen(curdir)
                open(unit=iun,form='unformatted',
     &           file=curdir(1:ld)//keyori(ind1:ind2),
     &           status='new',err=2200)
                fniun = curdir(1:ld)//keyori(ind1:ind2)
              else
                open(unit=iun,form='unformatted',
     &           file=keyori(ind1:ind2),
     &           status='new',err=2200)
                fniun = keyori(ind1:ind2)
              endif
            endif
         endif

         do j=ind1,ind2
            keyori(j:j) = ' '
            keywrd(j:j) = ' '
         end do
      endif

2200  continue
      return
      end

      subroutine filmap(mapit)
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil
      character curdir*1024
      common /cwd/ curdir

c
c----- process MAPFIL directive --------------------
c
      mapit = 0
      i = index(keywrd,'MAPFIL')
      if (i.ne.0) then
         mapit = 1
         ind = i + 4
         i   = index(keywrd(ind:320),'=')
         do j=i+ind,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j
         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1
         if (index(keyori(ind1:ind2),'/').eq.1) then
            mapfil = keyori(ind1:ind2)
         else
            ld = linlen(curdir)
            mapfil = curdir(1:ld)//keyori(ind1:ind2)
         endif
      endif

      call filgrd()

      return
      end

      subroutine filgrd()
      implicit double precision (a-h,o-z)
      character keywrd*320, keyori*320
      common /keywrd/ keywrd,keyori
      character*80 mapfil
      character*80 grdfil
      common /maphlp/ mapfil,grdfil
      character curdir*1024
      common /cwd/ curdir

c
c----- process RDBAS directive --------------------
c
      i = index(keywrd,'RDBAS=')
      if (i.ne.0) then
         ind = i + 5
         i   = index(keywrd(ind:320),'=')
         do j=ind+i,320
            if (ichar(keywrd(j:j)).ne.32) goto 2170
         end do
2170     ind1 = j
         do j=ind1,320
            if (ichar(keywrd(j:j)).eq.32) goto 2190
         end do
2190     ind2 = j-1

         if (index(keyori(ind1:ind2),'/').eq.1) then
            grdfil = keyori(ind1:ind2)
         else
            ld = linlen(curdir)
            grdfil = curdir(1:ld)//keyori(ind1:ind2)
         endif
      endif

      return
      end

