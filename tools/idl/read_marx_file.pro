function read_marx_file,inpfile,type,name,rlen,clen,reserv,magic
;-----------------------------------------------------------------------
; Name: READ_MARX_FILE
;
; Purpose: Returns the contents of a binary output file 
;          from the MARX simulator
;          
; Inputs:  inpfile -- string containing name of input file
;         
; Outputs:  header -- character array containing embedded 
;                     parameter file
;
; Comments: 
;           
; Revision history:
;       written					M. Wise	02-23-96
;       added BYTEORDER patch for LINUX		M. Wise	01-03-97
;       modified to read new PHAS format	M. Wise 05-13-97
;	added MAGIC number checking 
;
;-----------------------------------------------------------------------
;
; Check for correct number of parameters
;
np=n_params(0)
if ((np lt 1) or (np gt 7)) then begin
   print,string(7B),'CALLING SEQUENCE: ', $
   'array=read_marx_file(inpfile[,type,name,rlen,clen,reserv,magic]'
   return,-1
endif
;
; Open binary file
;
get_lun,unit
openr,unit,inpfile
;
; Read in data header info
;
;
; Offset  length   Interpretation
; ----------------------------------------------------------------------
; 0       4      Magic number: 0x83 0x13 0x89 0x8D
; 4       1      Data type: 
;                    'A'    8 bit signed integer (character)
;		     'I'   16 bit signed integer
;		     'J'   32 bit signed integer
;		     'E'   32 bit float
;		     'D'   64 bit float
;
; 5      15      Data column name.  If the length of the name is less than 15
;		 characters, it will be padded with 0.  If the name is 15
;		 characters, there will be no padding.
;
; 20      4       Number of Rows
; 24      4       Number of Columns, if 0 it is a vector
; 28      4       Reserved
; 32      ??      Data
; 
;
;
;magic=0L
magic=bytarr(4)
type=' '
name=bytarr(15)
rlen=0L
clen=0L
reser=bytarr(4)
readu,unit,magic,type,name,rlen,clen,reser
name=string(name)
;
; Check to see if this is a valid MARX file
;
marx_default=[131,19,137,141]
if (total(magic-marx_default) ne 0.0) then begin
   print,string(7B),'ERROR: Input file is not a valid MARX binary'
   return,-1
endif
;
; If on a LINUX machine do byte-swapping
;
if (!version.os eq 'linux') then byteorder,rlen,/NTOHL 
if (!version.os eq 'linux') then byteorder,clen,/NTOHL 
;
; Create array to hold data
;
if (clen eq 0) then begin
	case type of
	     'A': data=bytarr(rlen)
	     'I': data=intarr(rlen)
	     'J': data=lonarr(rlen)
	     'E': data=fltarr(rlen)
	     'D': data=dblarr(rlen)
	    else: begin
	          print,string(7B),'ERROR: Invalid data type input'
	          return,-1
	          end
	endcase
endif else begin
	case type of
	     'A': data=bytarr(clen,rlen)
	     'I': data=intarr(clen,rlen)
	     'J': data=lonarr(clen,rlen)
	     'E': data=fltarr(clen,rlen)
	     'D': data=dblarr(clen,rlen)
	    else: begin
	          print,string(7B),'ERROR: Invalid data type input'
		  return,-1
	          end
	endcase
endelse
;
; Read data into memory
;
readu,unit,data
;
; Close binary file
;
close,unit
free_lun,unit
;
; If on a LINUX machine do byte-swapping
;
if (!version.os eq 'linux') then begin
   case type of
       'A': 
       'I': byteorder,data,/NTOHS
       'J': byteorder,data,/NTOHL
       'E': byteorder,data,/XDRTOF
       'D': byteorder,data,/XDRTOD
      else: begin
            print,string(7B),'ERROR: Invalid data type input'
            return,-1
            end
   endcase
endif
;
; Convert byte data into signed integer
;
if (type eq 'A') then begin
   sign=-fix(data/128)
   data=data+sign*256
endif
;
; Return to IDL
;
return,data
end
