PRO tg_extract, hml, maxdist, y, z, y0, z0, tg_r, tg_d, list_selected
;+
; NAME:
;	tg_extract
;
; PURPOSE:
;	Routine to "extract" a spectrum - distance from zero-order within
;	some distance perpendicular from dispersion line
;
; CATEGORY:
;	Image processing, data analysis
;
; CALLING SEQUENCE:
;	tg_extract, hml, maxdist, y, z, y0, z0, tg_r, tg_d, list_selected
;
; INPUTS:
;	hml = 'h' or 'm' or 'l' (or uppercase) = HEG or MEG or LEG;
;		specifies which arm of the 'x' to use for distance from
;		dispersion line calculation.
;
;	maxdist = mm max distance from dispersion line (h or m) to allow.
;	       (could be other units, as long as y,z,y0,z0 are same units)
;
;	y,z = input coordinates in focal plane (typically mm, but
;	could be angle or pixel, as long as y0,z0,maxdist are the same)
;
;	y0,z0 - the Origin - hopefully the zero-order coord any point
;		on the line will do, but w/ arbitrary offset in
;		result.  (give in same units as y,z
;
;
; KEYWORD PARAMETERS:
;	NONE
;
; OUTPUTS:
;	tg_r = [mm], output vector of photon distances from zero order parallel to
;		dispersion.  [or same units as y,z]
;
;	tg_d = mm, output vector of photon distances from dispersion line in
;		perpendicular direction. [or same units as y,z]
;
;	list_selected = optional output vector of indices within input
;		vectors of selected photons
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;
;      Mainly a quick-and-dirty extraction, for examination of
;      simulations, though it also will work on real data.
;
;      Doesn't do 3D geometric transformations - uses 2D
;      approximation, which is very close for HETG, a bit less close
;      for LETG, due to larger departures from Rowland geometry.
;      Correct treatment will be done in ASC pipeline s/w, from pixel
;      3D position to diffraction angles.
;
;      Hard-coded clocking angles; doesn't read from simulation
;      parameter file (or from data header; in fact, doesn't even know
;      about a header).
;
;      Doesn't assign wavelengths.  Left as exercise for the user.
;      (well, see examples)
;
;	NOTE: To get a list of all photons distances when treated as either
;             all MEG or HEG, use a maxdist larger than a CCD, say
;             50mm.  Otherwise, the tg_r,tg_d pairs no longer
;             correspond to the input vectors or associated vectors
;             (such as nord).  (added optional output of
;             list_selected)
;
; PROCEDURE:
;	A (crude) diagram of the geometry:
;
;                                                                 o P4
;                         P1  o (y,z)= photon position            ^
;                              \                                  |
;                               \                                 | +Z
;                                \                                |
;                                 o P2= intersection of perp from |
;                                       P1 to dispersion line     |
;                                                                 |
;                                                                 |
; o P0=(0,0) --------------------- +Y --------------------------> o P3
;
;         = origin == zero-order position;
;                      dispersion direction is from P0 to P2
;
;	X = P0_P1
;	tg_d = P1_P2
;	tg_r = P0_P2
;	P0_P3 = y-axis of detector array
;	P3_P4 = z-direction
;
;	X^2 = y^2 + z^2
;	X^2 + tg_d^2 = tg_r^2
;
; EXAMPLE:
;
; - Run MARX to produce and output directory (on-axis source).
; - Use read_marx_file to get coordinate vectors:
;   yp = read_marx_file('./marx/ypos.dat')
;   zp = read_marx_file('./marx/zpos.dat')
; - Look at "image":
;   plot, yp, zp, psy=3
; - Use tg_extract to extract spectra:
;   tg_extract,'H', 0.1, yp, zp, 0.0, 0.0, tg_rh, tg_dh, lh
;   tg_extract,'M', 0.1, yp, zp, 0.0, 0.0, tg_rm, tg_dm, lm
; - See if you got the right ones:
;   oplot, yp[lm], zp[lm], psy=3, color=99
;   oplot, yp[lh], zp[lh], psy=3, color=44
; - Convert MEG to 1D:
;   spm = make_image(tg_rm, tg_dm, xra=[-80,80], yra=[-.1,.1], $
;         xbin=0.024, ybin=0.2, xaxis=axxm, yaxis=axym)
;   spm = reform(spm) ; remove superficial index of 1
; - Look at it:
;   plot, axxm, spm, psy=10, /yty, yra=[1,max(spm)]
;
; - etc. (like, use PHA to pick orders, or cheat and use .marx/order.dat...)
;
; - Want wavelengths?
;   tg_lambda = 4001.d0 * abs(tg_rm) / 8634.d0  ; m*lambda = P*sin(theta)
;
; MODIFICATION HISTORY:
; 	Written by:
;	d. huenemoerder 6/96
;
;       971014 d.huenemoerder: sanitized - removed vestigial peices of
;                              code; edited prologue and provided
;                              examples.

;; check params

   IF n_params() EQ 0 THEN BEGIN
      print,'tg_extract, hml, maxdist, y, z, y0, z0, tg_r, tg_d, [list_selected]'
      print, 'WARNING: unsupported code.  Use and/or modify at your own discretion!'
      return
   ENDIF

;;;;;;;;; establish variables and commons ;;;;;;;;;;;;;;;;
;
;;; these need to be tied to MARX parameters
                                ; params have been set to default
                                ; values. They will be applicable to
                                ; marx simulations w/ default param
                                ; file.  But if params are changed
                                ; (not recommended), they will be
                                ; inconsistent.

alpha_inner = -5.0d*!dtor       ; HEG clocking angle
alpha_outer = 5.d0*!dtor        ; MEG clocking angle, radians
alpha_middle = 0.0d0            ; LEG clocking angle
                                ; ("middle", though all shells,
                                ; because retro-fit into old bad code)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; select photons within specified distance perp. from dispersion line:
; (we don't have to select by HRMA shell, since HEG or MEG dispersion
; does this)  alpha_dispersion is the angle of the dispersion line
; relative to the -S array's long axis (the axaf y coordinate)
;
IF strupcase(hml) EQ 'M' THEN  BEGIN
   alpha_dispersion = alpha_outer
ENDIF ELSE IF strupcase(hml) EQ 'H' THEN BEGIN
   alpha_dispersion = alpha_inner
ENDIF ELSE IF strupcase(hml) EQ 'L' THEN BEGIN
   alpha_dispersion = alpha_middle
ENDIF


; the dispersion line by definition goes through (x0,y0)
;

   tg_d=-(y-y0)*sin(alpha_dispersion)+(z-z0)*cos(alpha_dispersion)
   list_selected = where(abs(tg_d) LT maxdist)

; again, since we have defined the origin at (x0,y0) as the position of
; the zero-order (not true for flight in general when there is aspect
; motion or an offset in pointing, but let's not consider it now),
; then the distance of the point from the zero order parallel to the
; dispersion is:

   tg_r = sqrt((float(y)-y0)^2+(float(z)-z0)^2 - float(tg_d)^2)

; adjust sign:

   neg_ords =  where((y-y0) LT 0)
   IF (size(neg_ords))(0) NE 0 THEN tg_r(neg_ords) = -tg_r(neg_ords)

;  select photons only within the specified distance:

   tg_d = tg_d(list_selected)
   tg_r = tg_r(list_selected)

END


