FUNCTION  make_image, vx, vy, $
            XRANGE=xminmax, XBINSIZE = xbin, XAXIS = xax,$
                      YRANGE=yminmax, YBINSIZE = ybin, YAXIS = yax,$
                      INDEX_LIST = idx, REVERSE_INDICES = revidx
;+
; NAME:
;	make_image
;
; PURPOSE:
;	Return the density function (histogram) of two variables.
;
; CATEGORY:
;	Image processing, statistics, probability.
;
; CALLING SEQUENCE:
;	Result = make_image(vx, vy, XRANGE=[xmin,xmax],XAXIS=xax,
;	XBINSIZE=xbin, YRANGE=[ymin,ymax], YAXIS=yax, YBINSIZE=ybin,
;	INDEX=idx, REVERSE_INDICES=revidx)
; INPUTS:
;	Vx and Vy = arrays containing the variables. 
;
; KEYWORD PARAMETERS:
;    XRANGE = [xmin,xmax]  --- input range in input array to make into image
;                              xmin refers to the left edge of the
;                              first bin; xmax is the RIGHT edge of
;                              the last bin.
;     YRANGE = [ymin,ymax]  --- input range in input array to make into image
;                              ymin refers to the LEFT EDGE of the
;                              first bin; ymax is the RIGHT EDGE of
;                              the last bin.
;     XBINSIZE = xbin       --- binsize for image x
;     YBINSIZE = ybin       --- binsize for image y
;     XAXIS = xax           --- output x-coordinate vector;
;                               xax(i) is the LEFT EDGE of bin i.
;     YAXIS = yax           --- output y-coordinate vector;
;                               yax(i) is the LEFT EDGE of bin i.
;     INDEX_LIST = idx      --- output vector of array indices of
;                               array elements used from vx, vy; if
;                               clipped with XRANGE or YRANGE, not all
;                               are used in the image.
;     REVERSE_INDICES = revidx --- output reverse-index vector, but
;                                  only of the sub-array actually
;                                  specified by XRANGE, YRANGE.
;                                  (use map_idx function --- see below).
;
; OUTPUTS:
;	The two dimensional density function of the two variables,
;	a longword array of dimensions (MAX(v1)+1, MAX(v2)+1).  Result(i,j)
;	is equal to the number of sumultaneous occurences of V1 = i,
;	and V2 = j, at the same element.
;	
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.

; RESTRICTIONS:
;	Reverse indexing is not complete - it refers to indices in the
;	sub-array, not the original array. --- this needs to be fixed.
;       As a work-around, the index list is returned which will allow
;       selection of the same points in the input vx, vy or any other
;       corresponding vector.   The function map_idx included here
;       will allow calculation of the proper index, given the limits
;       and binning used to make the image.
;
;       Note that the internal IDL histogram fuction treats bins
;       slightly differently - it's maximum parameter refers to the
;       LEFT EDGE of the last bin, so the real maximum data value
;       included is really maximum+binsize.
;
; PROCEDURE:
;	Creates a combines array from the two variables, equal to the
;	linear subscript in the resulting 2D histogram, then applies
;	the standard histogram function.
;
;	The following pseudo-code shows what the result contains,
;	not how it is computed:
;		r = LONARR(MAX(v1)+1, MAX(v2)+1)  ;Result
;		FOR i=0, N_ELEMENTS(v1)-1 DO $
;		  r(v1(i), v2(i)) = r(v1(i), v2(i)) +1
; EXAMPLE:
;	Return the 2D histogram made from two floating point images
;	with range of -1 to +1, and with 100 bins:
;	img =
;	    make_image(x,y,xra=[-1,1],yra=[-1,1],xbin=0.02,ybin=0.02,
;	               xaxis=xax, yaxis=yax)
;	Display the image with coordinates:
;	contour,img,xax,yax
;
; MODIFICATION HISTORY:
; 	Written by:
;	d. huenemoerder 6/96, based on hist_2d in the RSI/IDL lib.
;-


;; check params
   
   IF n_params() NE 2 THEN BEGIN
      message,$
       'USAGE: Result = make_image(vx, vy, XRANGE=[xmin,xmax], XAXIS=xax,'+$
       'YRANGE=[ymin,ymax], XBINSIZE=xbin, YBINSIZE=ybin, ' + $
       'YAXIS=yax,'+$
       'INDEX_LIST=idx, REVERSE_INDICES=revidx)'
   ENDIF
   
;;; determine limits

   xmin = min(vx,max = xmax)
   ymin = min(vy,max = ymax)
   
   IF n_elements(xminmax) NE 0 THEN BEGIN
      IF n_elements(xminmax) NE 2 THEN BEGIN
         message,'bad XRANGE.'
      ENDIF ELSE BEGIN
         xmin = xminmax(0)
         xmax = xminmax(1)
      ENDELSE
   ENDIF

   IF n_elements(yminmax) NE 0 THEN BEGIN
      IF n_elements(yminmax) NE 2 THEN BEGIN
         message,'bad YRANGE.'
      ENDIF ELSE BEGIN
         ymin = yminmax(0)
         ymax = yminmax(1)
      ENDELSE
   ENDIF
   
   IF n_elements(xbin) EQ  0 THEN xbin = 1
   IF n_elements(ybin) EQ  0 THEN ybin = 1
   
;;;;;
; test bin sizes
   IF xbin LE 0 THEN message,'XBINSIZE must be >0'
   IF ybin LE 0 THEN message,'YBINSIZE must be >0'
   
;;; determine array element list in clipping region:

   ll = (vx GE xmin) AND (vx LT xmax) AND (vy GE ymin) AND (vy LT ymax)
   idx = where(ll)

;;; truncate and scale vectors to longs:

   sx = long( ( vx(where(ll))-xmin) / xbin)
   sy = long( ( vy(where(ll))-ymin) / ybin)
   
   pxmax = long((xmax-xmin)/xbin) ; number of x pixels
   pymax = long((ymax-ymin)/ybin) ; number of y pixels
   
;;; transform coords to 1d and make 1D histogram
;;
   sum = pxmax * sy  +  sx

;   h = histogram(sum, min=0, max=m1 * m2 -1)  ; orig idl
   h = histogram(sum, min=0, max=pxmax*pymax-1, reverse=revidx)
   h =  reform(h, pxmax, pymax, /overwrite) ;and make it 2D
   
;;; compute coordinate axes - "left" edge of bin
;;
   xax = findgen(pxmax) * xbin + xmin
   yax = findgen(pymax) * ybin + ymin

   return, h
   
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;
;
FUNCTION map_idx, x, y, xmin, xmax, ymin, xbin, ybin
   
;; compute indices of x,y in the reverse index array.
   
;; x,y input coordinates
;; xmin,xmax,ymin = input, ranges used in make_image
;; xbin,ybin = bin sizes used in make_image   

   
   IF n_params() NE 7 THEN BEGIN
      message,'USAGE:  index=map_idx(x,y,xmin,xmax,ymin,xbin,ybin)'
   ENDIF
   

   idx = long ((x-xmin) / xbin)  + $
    long ( (y-ymin) / ybin)  * $
    (long ((xmax-xmin)/xbin))
   
   return,idx
   
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;
;
FUNCTION unmap_idx, pix, array
   
;; given the serial index, "pix", of a bin 2D array, "array", unmap
;; the index to the (x,y) bin pair.
   
   IF n_params() NE 2 THEN message, 'USAGE:  xy=unmap_idx(pix,array)'
   
   sz = size(array)
   
   IF sz(0) NE 2 THEN message, 'Array must be 2D.'
   
   xy = [pix MOD sz(1)  , pix/sz(1)]
   
   return, xy
   
END

   
      
      
