	subroutine nei_find(ndir,ic,jc,isc,jsc,in,jn,isn,jsn,sldbound)
	
	use globals
	implicit none
	
	integer,intent(in)::ic,jc,isc,jsc,ndir
	integer,intent(out)::in,jn,isn,jsn,sldbound
	integer::lev, ms, lev_nei, msn
	
	sldbound=0; lev=c_lev(ic,jc); ms=2**lev
	
	select case(ndir)
! 1 = upper cell, 3: down cell
		case(1,3)
			if(ndir==1)then
				if(jsc<ms)then
				  in=ic; jn=jc; isn=isc; jsn=jsc+1
				  return
				endif
				
				in=ic; jn=jc+1
				if(jn>jm)then
				  sldbound=1; isn=0; jsn=0
				  return
				else
				  lev_nei=c_lev(in,jn); jsn=1
				endif
			else
				if(jsc>1)then
				  in=ic; jn=jc; isn=isc; jsn=jsc-1
				  return
				endif
				
				in=ic; jn=jc-1
				if(jn<1)then
				  sldbound=1; isn=0; jsn=0
				  return
				else
				  lev_nei=c_lev(in,jn); jsn=2**lev_nei
				endif
			endif
			
			if(lev_nei==lev)then
			  isn=isc
			  return
			elseif(lev_nei==lev+1)then
			  isn=2*isc-1
			  return
			elseif(lev_nei==lev-1)then
				if(isc/2*2==isc)then
				  isn=isc/2
				else
				  isn=(isc+1)/2
				endif
				return
			else
			  sldbound=1; isn=0; jsn=0
			  return
			endif
! 2= right cell, 4 = left cell			
		case(2,4)
			if(ndir==2)then
				if(isc<ms)then
				  in=ic; jn=jc; isn=isc+1; jsn=jsc
				  return
				endif
				
				in=ic+1; jn=jc
				if(in>im)then
				  sldbound=1; isn=0; jsn=0
				  return
				else
				  lev_nei=c_lev(in,jn); isn=1
				endif
			else
				if(isc>1)then
				  in=ic; jn=jc; isn=isc-1; jsn=jsc
				  return
				endif
				in=ic-1; jn=jc
				if(in<1)then
				  sldbound=1; isn=0; jsn=0
				  return
				else
				  lev_nei=c_lev(in,jn); isn=2**lev_nei
				endif
			endif
			
			if(lev_nei==lev)then
			  jsn=jsc
			  return
			elseif(lev_nei==lev+1)then
			  jsn=2*jsc-1
			  return
			elseif(lev_nei==lev-1)then
				if(jsc/2*2==jsc)then
				  jsn=jsc/2
				else
				  jsn=(jsc+1)/2
				endif
				return
			else
			  sldbound=1; isn=0; jsn=0
			  return
			endif
 ! 5= right upper cell           
        case(5)
          if(isc<Ms.and.jsc<Ms) then
            in=ic; jn=jc; isn=isc+1; jsn=jsc+1
            return
          end if
          
          if(isc==Ms.and.jsc==Ms) then
              in=ic+1; jn=jc+1
            if(in>im.or.jn>jm) then !boundary
              isn=0; jsn=0
              return
            else
              isn=1; jsn=1
              return        
            end if
          
          elseif(isc==Ms.and.jsc/=Ms) then !isc defining local boundary
            in=ic+1; jn=jc
            if(in>jm) then !boundary
            isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn)
              if(lev_nei==lev+1) then
                isn=1; jsn=2*jsc+1
                return          
              elseif(lev_nei==lev) then
                isn=1; jsn=jsc+1
                return
              elseif(lev_nei==lev-1) then
                if(jsc/2*2==jsc) then !jsc is even
                  isn=1; jsn=jsc/2+1
                  return
                else
                  isn=0; jsn=0
                  return
                end if
              end if
            end if
        
          elseif(isc/=Ms.and.jsc==Ms) then !jsc defining local boundary
              in=ic; jn=jc+1
            if(jn>jm) then !boundary
			  isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn)
              if(lev_nei==lev+1) then
                isn=2*isc+1; jsn=1
                return
              elseif(lev_nei==lev) then
                isn=isc+1; jsn=1
                return
              elseif(lev_nei==lev-1) then
                if(isc/2*2==isc) then !jsc is even
                  isn=isc/2+1; jsn=1
                  return
                else
                  isn=0; jsn=0
                  return
                end if
              end if
            end if
          end if
! 6 = left down cell
        case(6)
          if(isc<Ms.and.jsc>1) then
            in=ic; jn=jc; isn=isc+1; jsn=jsc-1
            return
          end if
          
          if(isc==Ms.and.jsc==1) then
              in=ic+1; jn=jc-1
            if(in>im.or.jn<1) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              isn=1; jsn=Msn
              return        
            end if
          
          elseif(isc==Ms.and.jsc/=1) then !isc defining local boundary
              in=ic+1; jn=jc
            if(in>im) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn)
              if(lev_nei==lev+1) then
                isn=1; jsn=2*jsc-2
                return          
              elseif(lev_nei==lev) then
                isn=1; jsn=jsc-1
                return
              elseif(lev_nei==lev-1) then
                if(jsc/2*2.EQ.jsc) then !jsc is even
                  isn=0; jsn=0
                  return
                else
                  isn=1; jsn=(jsc-1)/2
                  return
                end if
              end if
            end if
        
          elseif(isc/=Ms.and.jsc==1) then !jsc defining local boundary
            in=ic; jn=jc-1
            if(jn<1) then !boundary
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              if(lev_nei==lev+1) then
                isn=2*isc+1; jsn=Msn
                return
              elseif(lev_nei==lev) then
                isn=isc+1; jsn=Msn
                return
              elseif(lev_nei==lev-1) then
                if(isc/2*2.EQ.isc) then !jsc is even
                  isn=isc/2+1; jsn=Msn
                  return
                else
                  isn=0; jsn=0
                  return
                end if
              end if
            end if
          end if
! 7= left down cell
        case(7)
          if(isc>1.and.jsc>1) then
            in=ic; jn=jc; isn=isc-1; jsn=jsc-1
            return
          end if
          
          if(isc==1.and.jsc==1) then
              in=ic-1; jn=jc-1
            if(in<1.or.jn<1) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              isn=Msn; jsn=Msn
              return        
            end if
          
          elseif(isc==1.and.jsc/=1) then !isc defining local boundary
              in=ic-1; jn=jc
            if(in<1) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              if(lev_nei==lev+1) then
                isn=Msn; jsn=2*jsc-2
                return          
              elseif(lev_nei==lev) then
                isn=Msn; jsn=jsc-1
                return
              elseif(lev_nei==lev-1) then
                if(jsc/2*2==jsc) then !jsc is even
                  isn=0; jsn=0
                  return
                else
                  isn=Msn; jsn=(jsc-1)/2
                  return
                end if
              end if
            end if
        
          elseif(isc/=1.and.jsc==1) then !jsc defining local boundary
              in=ic; jn=jc-1
            if(jn<1) then !boundary
			  isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              if(lev_nei==lev+1) then
                isn=2*isc-2; jsn=Msn
                return
              elseif(lev_nei==lev) then
                isn=isc-1; jsn=Msn
                return
              elseif(lev_nei==lev-1) then
                if(isc/2*2==isc) then !jsc is even
                  isn=0; jsn=0
                  return
                else
                  isn=(isc-1)/2; jsn=Msn
                  return
                end if
              end if
            end if
          end if
! 8 = left up cell
        case(8)
          if(isc>1.and.jsc<Ms) then
            in=ic; jn=jc; isn=isc-1; jsn=jsc+1
            return
          end if
          
          if(isc==1.and.jsc==Ms) then ! A real corner cell
              in=ic-1; jn=jc+1
            if(in<1.or.jn>jm) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              isn=Msn; jsn=1
              return        
            end if
          
          elseif(isc==1.and.jsc/=Ms) then !isc defining local boundary
              in=ic-1; jn=jc
            if(in<1) then !boundary
			  isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn); Msn=2**lev_nei
              if(lev_nei==lev+1) then
                isn=Msn; jsn=2*jsc+1
                return          
              elseif(lev_nei==lev) then
                isn=Msn; jsn=jsc+1
                return
              elseif(lev_nei==lev-1) then
                if(jsc/2*2==jsc) then !jsc is even
                  isn=Msn; jsn=jsc/2+1
                  return
                else
                  isn=0; jsn=0
                  return
                end if
              end if
            end if
        
          elseif(isc/=1.and.jsc==Ms) then !jsc defining local boundary
              in=ic; jn=jc+1
            if(jn>jm) then !boundary
              isn=0; jsn=0
              return
            else
              lev_nei=c_lev(in,jn)
              if(lev_nei==lev+1) then
                isn=2*isc-2; jsn=1
                return
              elseif(lev_nei==lev) then
                isn=isc-1; jsn=1
                return
              elseif(lev_nei==lev-1) then
                if(isc/2*2==isc) then !jsc is even
                  isn=0; jsn=0
                  return
                else
                  isn=(isc-1)/2; jsn=1
                  return
                end if
              end if
            end if
          end if
	end select
	
	end subroutine nei_find
