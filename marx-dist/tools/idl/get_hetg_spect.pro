pro get_hetg_spect,path,order,rm,xm,ym,rh,xh,yh,p,xp,yp
;-----------------------------------------------------------------------
; Name: GET_HETG_SPECT
;
; Purpose: Uses MARX internal simualtion variables to extract the
;          MEG and HEG dispersed spectra from a given simulation.
;          
; Inputs: 
;         
;         
; Comments: 
;           
;           
; Revision history:
;       written by Michael Wise, 06-06-97
;-----------------------------------------------------------------------
;
; Check for correct number of parameters
;
np=n_params(0)
if ((np lt 3) or (np gt 20)) then begin
   print,string(7B),'CALLING SEQUENCE: ', $
   'get_hetg_spect,path,order,rm,[xm,ym,rh,xh,yh,p,xp,yp]'
   return   
endif

;
; Read in necessary files
;
y=read_marx_file(path+'/ypos.dat')
z=read_marx_file(path+'/zpos.dat')
o=read_marx_file(path+'/order.dat')
m=read_marx_file(path+'/mirror.dat')
p=read_marx_file(path+'/pha.dat')

r=sqrt(y*y+z*z)


io=where(abs(o) eq 0) 
im=where( (abs(o) eq order) and ((m eq 0) or (m eq 1)) )
ih=where( (abs(o) eq order) and ((m eq 2) or (m eq 3)) )

rp=p(io)
if (order eq 0) then begin
   rm=r(io)
endif else begin
   rm=r(im)
   rh=r(ih)
endelse


;
; Convert to energy grids
;
ym=histogram(rm,binsize=0.024,min=0)
nm=n_elements(ym)
xm=0.024*(findgen(nm)+1.0)
xm=(12.3985/xm)*(8634./4001.)

yh=histogram(rh,binsize=0.024,min=0)
nh=n_elements(yh)
xh=0.024*(findgen(nh)+1.0)
xh=(12.3985/xh)*(8634./2002.)

yp=histogram(rp,binsize=1.0,min=0)
np=n_elements(yp)
xp=1.0*(findgen(np)+1.0)
xp=xp*0.0045

;
; Return to IDL
;
return
end
