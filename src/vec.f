      subroutine vec(small,ohoh,u,c,j,k)
      implicit double precision (a-h,o-z)
      logical ohoh
 
      dimension c(*),r(3),u(3)
 
      data dzero/0.0d0/
 
      jtemp = (j-1)*3
      ktemp = (k-1)*3
      r2 = dzero

      do i=1,3
         r(i) = c(i+jtemp) - c(i+ktemp)
         r2   = r2 + r(i)*r(i)
      end do

      r2 = dsqrt(r2)
      ohoh = r2.lt.small
      if (ohoh) return

      do i=1,3
         u(i) = r(i)/r2
      end do

      return
      end
