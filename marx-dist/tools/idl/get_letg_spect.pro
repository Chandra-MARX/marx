pro get_letg_spect,path,order,rl,xl,yl,p,xp,yp
;-----------------------------------------------------------------------
; Name: GET_LETG_SPECT
;
; Purpose: Uses MARX internal simualtion variables to extract the
;          LETG dispersed spectra from a given simulation.
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
   'get_letg_spect,path,order,rl,[xl,yl,p,xp,yp]'
   return   
endif

;
; Read in necessary files
;
y=read_marx_file(path+'/ypos.dat')
z=read_marx_file(path+'/zpos.dat')
o=read_marx_file(path+'/order.dat')
p=read_marx_file(path+'/pha.dat')

r=sqrt(y*y)

io=where(abs(o) eq 0) 
il=where( (abs(o) eq order) )

rp=p(io)
if (order eq 0) then begin
   rl=r(io)
endif else begin
   rl=r(il)
endelse


;
; Convert to energy grids
;
yl=histogram(rl,binsize=0.0064,min=0)
nl=n_elements(yl)
xl=0.0064*(findgen(nl)+1.0)
xl=(12.3985/xl)*(8634./9912.)

yp=histogram(rp,binsize=1.0,min=0)
np=n_elements(yp)
xp=1.0*(findgen(np)+1.0)
xp=xp/107.0

;
; Return to IDL
;
return
end
