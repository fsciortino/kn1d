; $Id: reverse.pro,v 1.10 2002/02/06 21:45:53 scottm Exp $
;
; Copyright (c) 1991-2002, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.

;+
; NAME:
;	REVERSE
;
; PURPOSE:
;	Reverse the order of rows or columns in an array or vector.
;
; CATEGORY:
;	Array manipulation.
;
; CALLING SEQUENCE:
;	Result = REVERSE(Array [, Subscript_Index])
;
; INPUTS:
;	Array:	The array or vector containing the original data.
;
; OPTIONAL INPUT PARAMETERS:
; Subscript_Index:  If this parameter is omitted or 1, the first subscript is
;		reversed (i.e., rows are reversed).  Set this parameter to
;		2 to reverse columns.
;
; KEYWORD PARAMETERS:
;	OVERWRITE = Set this keyword to do the transformation "in-place".
;             The result overwrites the previous contents of the variable.
;
; OUTPUTS:
;	REVERSE returns the input array, but reversed about
;	one of its dimensions.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	Uses the REFORM function.
;
; MODIFICATION HISTORY:
;	Old.
;	Apr, 1991, DMS,	Added 3D reversing.
;       Sept, 1992 Mark L. Rivers, added simple return for scaler argument
;	Sept, 1994. Added default for 3D case.
;   May 2000, CT, Rewrote to handle any dimensions, added OVERWRITE keyword.
;-

function reverse, a, subscript, $
	OVERWRITE=overwrite

	on_error,2                             ;Return to caller if an error occurs

	ndims = SIZE(a,/N_DIMENSIONS)
	if ndims eq 0 then return, a

	if n_elements(subscript) le 0 then subscript = 1  ;Default case
	IF (subscript GT ndims) THEN MESSAGE, $
		"Subscript_index must be less than or equal to number of dimensions."


; handle 1 or 2 dimensions using ROTATE for efficiency
	IF (ndims EQ 1) THEN RETURN, ROTATE(a,5)
	IF (ndims EQ 2) THEN BEGIN
		CASE (subscript) OF
		1: RETURN, ROTATE(a, 5)
		2: RETURN, ROTATE(a, 7)
		ENDCASE
	ENDIF


; for 3 or more dimensions, collapse down to 3 dimensions & loop over index
	b = KEYWORD_SET(overwrite) ? TEMPORARY(a) : a
	; compress the smaller (inner) & larger (outer) dimensions
	; so we only have to deal with a 3-dimensional array
	dimensions = SIZE(b,/DIMENSIONS)
	nDo = dimensions[subscript-1]   ; dimension size for Subscript_index
; array size for dimensions smaller than Subscript_index
	nLess = 1L
	FOR i=0,subscript-2 DO nLess = nLess*dimensions[i]
; array size for dimensions greater than Subscript_index
	nMore = 1L
	FOR i=subscript,ndims-1 DO nMore = nMore*dimensions[i]
	b = REFORM(b, nLess, nDo, nMore, /OVERWRITE)


	; manually loop over the middle dimension (could do this using an
	; index array, but it might use too much memory)
	FOR i=0ull,(nDo-1)/2 DO BEGIN
		temp = b[*,nDo-1-i,*]
		b[0,nDo-i-1,0] = b[*,i,*]
		b[0,i,0] = temp
	ENDFOR

	return,REFORM(b,dimensions,/OVERWRITE)  ; restore original dimensions

end
