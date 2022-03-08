
; function dis_surf, a0=a0, b0=b0, a1=a1, b1=b1
; 	dis = ((a0-a1)^2 + (b0-b1)^2)^0.5
;     return, dis
; end
; function inv_sq, x=x
; 	xinv = 1/x^2
;     return, xinv
; end
; function get_nst, an0=an0, bn0=bn0, array_x=array_x, array_y=array_y
; 	print,'success import get_nearst'
;     dis_all = dis_surf(an0, bn0, array_x, array_y)
;     dex_array = sort(dis_all)[0:(4-1)]
;     dis_array = dis_all[dex_array]
; 	ret_array = [transpose(dis_array),transpose(dex_array)]
;     return, ret_array
; end
