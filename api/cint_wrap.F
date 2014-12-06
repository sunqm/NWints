      module cint_wrap
c store libcint environments in this module
      integer CINTnatm
      integer CINTnbas
      integer,allocatable CINTatm(:,:)
      integer,allocatable CINTbas(:,:)
      double precision,allocatable CINTenv(:,:)
c initialize and delete opt by
c call f2e_optimizer(CINTopt, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
c call CINTdel_optimizer(CINTopt)
      integer(8) CINTopt
      end module


      subroutine init_cint_env()
      use cint_wrap
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "geom.fh"
      integer offset = 0
      integer geom, basis       ! [input] handles
      integer natoms, nshell
      integer i, j, ptr, ish
      double precision,allocatable envbuf(:)

      status = geom_ncent(geom, natoms)
      status = status .and. bas_numcont(basis, nshell)
      if (.not. status) call errquit('schwarz_init: info failed', 0,
     &       INPUT_ERR)
      CINTnatm = natoms
      CINTnbas = nshell

      allocate (CINTatm(6,natoms))
      allocate (CINTbas(8,nshell))
c a buf which assumes 50->20 contraction for each shell should be enough
      allocate(envbuf(natoms*4+nshell*50*20))
      do i = 1, natoms
        CINTatm(1,i) = 0 !TODO: nuclear charge
        CINTatm(2,i) = offset
        envbuf(offset+1) = coords(1)
        envbuf(offset+2) = coords(2)
        envbuf(offset+3) = coords(3)
        offset = offset + 3
        CINTatm(3,i) = 1 ! point nuclear model
        offset = offset + 1
c        CINTatm(4,i) = offset ! pointer to nuclear mass
c        envbuf(offset+1) = mass
c        offset = offset + 1
      end do

      do ish = 1, nshell
        CINTbas(1,ish) = ibs_cn2ce (ncont,ish) - 1 ! atom_id, 0-based
        CINTbas(2,ish) = infbs_cont(CONT_TYPE,uncont,ish)  ! angular momentum
        CINTbas(3,ish) = infbs_cont(CONT_NPRIM,ucont,ish) ! num of pGTO
        CINTbas(4,ish) = infbs_cont(CONT_NGEN ,ucont,ish) ! num of contracted GTO
        CINTbas(5,ish) = 0 ! kappa for spinor basis
        CINTbas(6,ish) = offset ! pointer to exp of primitive GTO
        ptr = mb_exndcf(infbs_cont(CONT_IEXP,ucont,ish), ish)
        do j = 1, CINTbas(3,ish)
          envbuf(offset+j) = dbl_mb(ptr+j)
        end do
        offset = offset + CINTbas(3,ish)
        CINTbas(7,ish) = offset ! pointer to contraction coefficients
        ptr = mb_exndcf(infbs_cont(CONT_ICFP,ucont,ish), ish)
        do j = 1, CINTbas(3,ish)*CINTbas(4,ish)
          envbuf(offset+j) = dbl_mb(ptr+j)
        end do
        offset = offset + CINTbas(3,ish)*CINTbas(4,ish)
      end do

      allocate (CINTenv(offset))
      do i = 1, offset
        CINTenv(i) = envbuf(i)
      end do
      deallocate (envbuf)
      end subroutine init_cint_env


      subroutine del_cint_env()
      use cint_wrap
      deallocate (CINTatm)
      deallocate (CINTbas)
      deallocate (CINTenv)
      end subroutine del_cint_env


      subroutine cint_2e4c(brain, ish, jsh, ketin, ksh, lsh,
     &                     lscr, scr, leri, eri)
c same API as int_2e4c
      use cint_wrap
      implicit none
c:: passed
c:tex-\begin{verbatim}
      integer brain              !< [Input] no use
      integer ish                !< [Input] the shell index
      integer jsh                !< [Input] the shell index
      integer ketin              !< [Input] no use
      integer ksh                !< [Input] the shell index
      integer lsh                !< [Input] the shell index
      integer lscr               !< [Input] no use
      double precision scr(lscr) !< [Scratch] no use
      integer leri               !< [Input] the length of integral array
      double precision eri(leri) !< [Output] the 2-electron integrals
c:tex-\end{verbatim}
c:: local
      integer shls(4)
      shls(1) = ish - 1
      shls(2) = jsh - 1
      shls(3) = ksh - 1
      shls(4) = lsh - 1
      cint2e_sph(eri, shls, CINTatm, CINTnatm, CINTbas, CINTnbas,
     &           CINTenv, CINTopt) 
      end subroutine cint_2e4c
