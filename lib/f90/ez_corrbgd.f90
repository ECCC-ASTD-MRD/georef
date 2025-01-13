
      subroutine ez_corrbgd(zout, ni,nj, hem)
      implicit none
      integer ni,nj,hem
      real zout(ni,nj)

      integer i
      real moyenne, somme


      if (hem.eq.0.or.hem.eq.2) then
         somme = 0.
         do i=1,ni
            somme = somme + zout(i,1)
         enddo
         moyenne = somme / (ni *1.0)
         
         do i=1,ni
            zout(i,1) = moyenne
         enddo
      endif
         
      if (hem.eq.0.or.hem.eq.1) then
         somme = 0.
         do i=1,ni
            somme = somme + zout(i,nj)
         enddo
         moyenne = somme / (ni *1.0)
         
         do i=1,ni
            zout(i,nj) = moyenne
         enddo

      endif
      return
      end
      
      
