      module cint_wrap
c store libcint environments in this module
      integer CINTnatm
      integer CINTnbas
      integer,allocatable :: CINTatm(:,:)
      integer,allocatable :: CINTbas(:,:)
      double precision,allocatable :: CINTenv(:) ! Use only as 1D array
c initialize and delete opt by
c call f2e_optimizer(CINTopt, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
c call CINTdel_optimizer(CINTopt)
      integer(8) CINTopt
      end module


      subroutine init_cint_env(rtdb,basis)
      use cint_wrap
      implicit none
cHACK
c     Ideally these two include files should not be here.
c     The <name>P.fh indicates that this is a Private include file of
c     the module indicated by <name>. Instead we should be calling
c     the regular API functions of the BAS module. Nevertheless, let's
c     get this stuff working first, then we can make it nice.
#include "nwc_const.fh"
#include "mafdecls.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "bas_exndcf.fh"
#include "bas_ibs_dec.fh"
#include "geobasmapP.fh"
#include "bas_ibs_sfn.fh"
cHACK
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
      integer offset
      integer rtdb        !< [Input] The RTDB handle
      integer geom, basis       ! [input] handles
      integer natoms, nshell
      integer i, j, k, ptr, iat
      integer ishlo, ishhi, ish, myprim, mycont, ucont
      logical status
      character*16 gtag         ! geometry tag 
      double precision chg
      double precision cxyz(3)
      double precision,allocatable :: envbuf(:)
      data offset/0/
      integer ibas, expoffset
      double precision exptmp
      double precision,external :: CINTgto_norm

      status = bas_geom(basis,geom)
      status = status .and. geom_ncent(geom, natoms)
      status = status .and. bas_numcont(basis, nshell)
      if (.not. status) call errquit('cint_wrap: info failed', 0,
     &       INPUT_ERR)
      CINTnatm = natoms
      CINTnbas = nshell

      allocate (CINTatm(6,natoms))
      allocate (CINTbas(8,nshell))
c a buf which assumes 50->20 contraction for each shell should be enough
      allocate(envbuf(natoms*4+nshell*50*20))
      ibas = basis + BASIS_HANDLE_OFFSET
      do iat = 1, natoms
        if (.not.geom_cent_get(geom, iat, gtag, cxyz, chg)) call errquit
     &      ('cint_wrap: geom_cent_get failed',911, GEOM_ERR)
        CINTatm(1,iat) = chg
        CINTatm(2,iat) = offset
        envbuf(offset+1) = cxyz(1)
        envbuf(offset+2) = cxyz(2)
        envbuf(offset+3) = cxyz(3)
        offset = offset + 3
        CINTatm(3,iat) = 1 ! point nuclear model, 2 = gaussian nuc mod
c        CINTatm(4,iat) = offset ! pointer to nuclear mass
c        envbuf(offset+1) = mass
c        offset = offset + 1

        if (.not. bas_ce2cnr(basis, iat, ishlo, ishhi))
     $        call errquit('cint_wrap: bas_ce2cnr failed', ish,
     &       BASIS_ERR)
        do ish = ishlo, ishhi
          ucont = sf_ibs_cn2ucn(ish, ibas)
          CINTbas(1,ish) = sf_ibs_cn2ce(ish,ibas) - 1 ! atom_id, 0-based
          CINTbas(2,ish) = infbs_cont(CONT_TYPE ,ucont,ibas) ! angular momentum
          myprim = infbs_cont(CONT_NPRIM,ucont,ibas) ! num of primitive GTO
          mycont = infbs_cont(CONT_NGEN ,ucont,ibas) ! num of contracted GTO
          CINTbas(3,ish) = myprim
          CINTbas(4,ish) = mycont
          CINTbas(5,ish) = 0 ! kappa for spinor basis
          CINTbas(6,ish) = offset ! pointer to exp of primitive GTO

          ptr = infbs_cont(CONT_IEXP,ucont,ibas) - 1
          do i = 1, myprim
            envbuf(offset+i) = sf_exndcf(ptr+i,ibas)
          end do
          expoffset = offset
          offset = offset + CINTbas(3,ish)

          CINTbas(7,ish) = offset ! pointer to contraction coefficients
          ptr = infbs_cont(CONT_ICFP,ucont,ibas) - 1
          do j = 1, mycont
            do i = 1, myprim
              exptmp = envbuf(expoffset+i)
              envbuf(offset+i) = sf_exndcf(ptr+i,ibas)
     $                         * CINTgto_norm(CINTbas(2,ish), exptmp)
            end do
            ptr = ptr + myprim
            offset = offset + myprim
          end do
        end do
      end do

      allocate (CINTenv(offset))

      do i = 1, offset
        CINTenv(i) = envbuf(i)
      end do
      deallocate (envbuf)

      CINTopt = 0
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
      call cint2e_sph(eri, shls, CINTatm, CINTnatm, CINTbas, CINTnbas,
     &                CINTenv, CINTopt) 
      end subroutine cint_2e4c
