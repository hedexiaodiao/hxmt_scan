function modulate_curve_le_psf_modify, RA=RA, DEC=DEC, Flux=Flux, Qt1=Qt1, Qt2=Qt2, Qt3=Qt3, BoxID=BoxID, DL=DL, euler=euler, par=par
	
	L0_le = 6.d  *!dpi/180
	L1_le = 1.6d *!dpi/180
	L2_le = 4.d  *!dpi/180
	
	L0_le = L0_le + DL[0] *!dpi/180
	L1_le = L1_le + DL[1] *!dpi/180
		
	if BoxID eq 0 then R_b = [[0.5, -0.866, 0.d], [-0.866, -0.5, 0], [0, 0, -1]]
	if BoxID eq 1 then R_b = [[1.d,      0, 0.d], [     0, -1.d, 0], [0, 0, -1]]
	if BoxID eq 2 then R_b = [[0.5,  0.866, 0.d], [ 0.866, -0.5, 0], [0, 0, -1]]
	
	euler_rad = euler * !dpi/180
	
	rota_00 = cos(euler_rad[0])*cos(euler_rad[1])
	rota_01 = cos(euler_rad[0])*sin(euler_rad[1])*sin(euler_rad[2])-sin(euler_rad[0])*cos(euler_rad[2])
	rota_02 = cos(euler_rad[0])*sin(euler_rad[1])*cos(euler_rad[2])+sin(euler_rad[0])*sin(euler_rad[2])
	rota_10 = sin(euler_rad[0])*cos(euler_rad[1])
	rota_11 = sin(euler_rad[0])*sin(euler_rad[1])*sin(euler_rad[2])+cos(euler_rad[0])*cos(euler_rad[2])
	rota_12 = sin(euler_rad[0])*sin(euler_rad[1])*cos(euler_rad[2])-cos(euler_rad[0])*sin(euler_rad[2])
	rota_20 =-sin(euler_rad[1])
	rota_21 = cos(euler_rad[1])*sin(euler_rad[2])
	rota_22 = cos(euler_rad[1])*cos(euler_rad[2])
	
	rota = [[rota_00, rota_10, rota_20], $
			[rota_01, rota_11, rota_21], $
			[rota_02, rota_12, rota_22]]
	
	;------------
	Qt0 = (1-Qt1^2-Qt2^2-Qt3^2)^0.5
	dex = where(abs(Qt1^2+Qt2^2+Qt3^2-1) lt 1d-06)
	if dex[0] ne -1 then Qt0[dex] = 0.d
	
	x_target = cos(RA*!dpi/180) * cos(DEC*!dpi/180)
	y_target = sin(RA*!dpi/180) * cos(DEC*!dpi/180)
	z_target = sin(DEC*!dpi/180)
	
	r_target = transpose([x_target, y_target, z_target])
	
	n = n_elements(Qt1)
	alpha_sf = dblarr(n)  &  beta_sf = dblarr(n)
	for i=0L, n-1 do begin
		R = [[Qt0[i]^2+Qt1[i]^2-Qt2[i]^2-Qt3[i]^2, 2*(Qt1[i]*Qt2[i]+Qt0[i]*Qt3[i]), 2*(Qt1[i]*Qt3[i]-Qt0[i]*Qt2[i])], $
			 [2*(Qt1[i]*Qt2[i]-Qt0[i]*Qt3[i]), Qt0[i]^2+Qt2[i]^2-Qt1[i]^2-Qt3[i]^2, 2*(Qt2[i]*Qt3[i]+Qt0[i]*Qt1[i])], $
			 [2*(Qt1[i]*Qt3[i]+Qt0[i]*Qt2[i]), 2*(Qt2[i]*Qt3[i]-Qt0[i]*Qt1[i]), Qt0[i]^2+Qt3[i]^2-Qt1[i]^2-Qt2[i]^2]]
		
		rm_b = [[0,0,1],[0,1,0],[-1,0,0]] ## rota ## R ## r_target
		;print, R_b ## [[0,0,1],[0,1,0],[-1,0,0]]
		;------------------------
		r_d = R_b ## rm_b
		
		alpha_sf[i] = atan(r_d[0]/r_d[2])  &  beta_sf[i] = atan(r_d[1]/r_d[2])
		
	endfor
	
	model_sf = Flux*(1-abs(tan(alpha_sf))/tan(L0_le)) * (1-abs(tan(beta_sf))/tan(L1_le)) / (tan(alpha_sf)^2 + tan(beta_sf)^2 + 1)^0.5
	
	alpha_sf_deg = alpha_sf * 180/!dpi
	beta_sf_deg = beta_sf * 180/!dpi
	
	;print,'LE:',alpha_sf_deg,beta_sf_deg
	
	model_sf = model_sf * (par[0]*alpha_sf_deg^2*beta_sf_deg^2 + par[1]*alpha_sf_deg^2 + par[2]*beta_sf_deg^2 + par[3])
		
	dex_sf = where(abs(alpha_sf) gt L0_le or abs(beta_sf) gt L1_le)
	if dex_sf[0] ne -1 then model_sf[dex_sf] = 0
	
	return, model_sf
	
end

function get_z_zerr, a0=a0, b0=b0, array_x=array_x, array_y=array_y, array_psf=array_psf, array_zerr=array_zerr
	;print,'successs import get_z_zerr'
    ;ret_array = get_nst(an0=a0, bn0=b0, array_x=data_x, array_y=data_y)
	dis_all = ((a0-array_x)^2 + (b0-array_y)^2)^0.5
	sort_dis = sort(dis_all)
    dex_array = sort_dis[0:(4-1)]
    dis_array = dis_all[dex_array]
	; dis_array = transpose(ret_array[0,*])
	; dex_array = transpose(ret_array[1,*])
    nrst4z= array_psf[dex_array]
    nrst4zerr = array_zerr[dex_array]
    if min(dis_array) eq 0 then begin
        calc_z = array_psf[dex_array[0]]
        ;calc_weight = data_zerr[dex_array[0]]
	endif
    if min(dis_array) ne 0 then begin
        calc_z = total((1/nrst4zerr^2)*(1/dis_array^2)*nrst4z)/total((1/nrst4zerr^2)*(1/dis_array^2))
        ;calc_weight = total(inv_sq(nrst4zerr)*inv_sq(dis_array))/total(inv_sq(dis_array))
	endif
    return, calc_z
end

function get_psf, x=x, y=y, instru=instru, BoxID=BoxID
	;print, 'success import get_psf'
    data_x = transpose(read_data(STRCOMPRESS(instru+'_psf_alpha_box'+strtrim(BoxID)+'.txt', /REMOVE_ALL)))
    data_y = transpose(read_data(STRCOMPRESS(instru+'_psf_beta_box'+strtrim(BoxID)+'.txt', /REMOVE_ALL)))
    data_psf = transpose(read_data(STRCOMPRESS(instru+'_psf_value_box'+strtrim(BoxID)+'.txt', /REMOVE_ALL)))
    data_psferr = make_array(n_elements(data_psf),value=1)
    
	n = n_elements(x)
	calc_z = dblarr(n)
	for i=0L, n-1 do begin
		;print, x[i], y[i],data_x,data_y,data_psf,data_psferr
    	calc_z[i] = get_z_zerr(a0=x[i], b0=y[i], array_x=data_x, array_y=data_y, array_psf=data_psf, array_zerr=data_psferr)
	endfor
	psf_value = calc_z

	;print, 'success excu get_psf'
    return, psf_value
end

function modulate_curve_me_psf_modify, RA=RA, DEC=DEC, Flux=Flux, Qt1=Qt1, Qt2=Qt2, Qt3=Qt3, BoxID=BoxID
	;print, 'success import psf modify'
	instru = 'ME'

	if BoxID eq 0 then R_b = [[0.5, -0.866, 0.d], [-0.866, -0.5, 0], [0, 0, -1]]
	if BoxID eq 1 then R_b = [[1.d,      0, 0.d], [     0, -1.d, 0], [0, 0, -1]]
	if BoxID eq 2 then R_b = [[0.5,  0.866, 0.d], [ 0.866, -0.5, 0], [0, 0, -1]]
	
	
	;------------
	Qt0 = (1-Qt1^2-Qt2^2-Qt3^2)^0.5
	dex = where(abs(Qt1^2+Qt2^2+Qt3^2-1) lt 1d-06)
	if dex[0] ne -1 then Qt0[dex] = 0.d
	
	x_target = cos(RA*!dpi/180) * cos(DEC*!dpi/180)
	y_target = sin(RA*!dpi/180) * cos(DEC*!dpi/180)
	z_target = sin(DEC*!dpi/180)
	
	r_target = transpose([x_target, y_target, z_target])
	
	n = n_elements(Qt1)
	alpha_sf = dblarr(n)  &  beta_sf = dblarr(n)
	for i=0L, n-1 do begin
		R = [[Qt0[i]^2+Qt1[i]^2-Qt2[i]^2-Qt3[i]^2, 2*(Qt1[i]*Qt2[i]+Qt0[i]*Qt3[i]), 2*(Qt1[i]*Qt3[i]-Qt0[i]*Qt2[i])], $
			 [2*(Qt1[i]*Qt2[i]-Qt0[i]*Qt3[i]), Qt0[i]^2+Qt2[i]^2-Qt1[i]^2-Qt3[i]^2, 2*(Qt2[i]*Qt3[i]+Qt0[i]*Qt1[i])], $
			 [2*(Qt1[i]*Qt3[i]+Qt0[i]*Qt2[i]), 2*(Qt2[i]*Qt3[i]-Qt0[i]*Qt1[i]), Qt0[i]^2+Qt3[i]^2-Qt1[i]^2-Qt2[i]^2]]
		
		rm_b = [[0,0,-1],[0,1,0],[1,0,0]] ## R ## r_target
		;print, R_b ## [[0,0,1],[0,1,0],[-1,0,0]]
		;------------------------
		r_d = R_b ## rm_b
		
		alpha_sf[i] = atan(r_d[0]/r_d[2])  &  beta_sf[i] = atan(r_d[1]/r_d[2])
		
	endfor
	
	psf_value = get_psf(x=alpha_sf,y=beta_sf,instru=instru,BoxID=BoxID)
	model_sf = Flux*psf_value
	
	; alpha_sf_deg = alpha_sf * 180/!dpi
	; beta_sf_deg = beta_sf * 180/!dpi
	; print,'ME:',alpha_sf_deg,beta_sf_deg
	;print, 'success excu psf modify'
	return, model_sf
	
end


pro position_test_v5, energy_LE=energy_LE, energy_ME=energy_ME, burst_dex=burst_dex, print_result=print_result, use_LE=use_LE, use_ME=use_ME,n_step=n_step,n_burn=n_burn
	
	if energy_ME eq '7-20' then begin
		mepsf_para = transpose(read_data('Parabolic_ME_2020_2201_7_20.txt'))
	endif
	if energy_ME eq '7-12' then begin
		mepsf_para = transpose(read_data('Parabolic_ME_2020_2201_7_12.txt'))
	endif
	DL_ME_7_20_box0 = [mepsf_para[9*0+3,0],mepsf_para[9*0+4,0]]
	DL_ME_7_20_box1 = [mepsf_para[9*1+3,0],mepsf_para[9*1+4,0]]
	DL_ME_7_20_box2 = [mepsf_para[9*2+3,0],mepsf_para[9*2+4,0]]
	
	euler_ME_7_20_box0 = [mepsf_para[9*0+0,0],mepsf_para[9*0+1,0],mepsf_para[9*0+2,0]]
	euler_ME_7_20_box1 = [mepsf_para[9*1+0,0],mepsf_para[9*1+1,0],mepsf_para[9*1+2,0]]
	euler_ME_7_20_box2 = [mepsf_para[9*2+0,0],mepsf_para[9*2+1,0],mepsf_para[9*2+2,0]]
	
	par_ME_7_20_box0 = [mepsf_para[9*0+5,0],mepsf_para[9*0+6,0],mepsf_para[9*0+7,0],mepsf_para[9*0+8,0]]
	par_ME_7_20_box1 = [mepsf_para[9*1+5,0],mepsf_para[9*1+6,0],mepsf_para[9*1+7,0],mepsf_para[9*1+8,0]]
	par_ME_7_20_box2 = [mepsf_para[9*2+5,0],mepsf_para[9*2+6,0],mepsf_para[9*2+7,0],mepsf_para[9*2+8,0]]
	;burst_dex = 1
	
	;energy_LE = '2-6'
	;energy_ME = '7-20'
	;burst_dex = 24
	;----------------------------------------------

	use_LE_only = use_LE and (Not use_ME)
	use_ME_only = use_ME and (Not use_LE)
	use_LE_ME = use_LE and use_ME
	;------------------------------------------------
	counts_le_0 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g0*.lc',1,header,/silent,/DSCALE)   
	counts_le_1 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g1*.lc',1,header,/silent,/DSCALE)
	counts_le_2 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g2*.lc',1,header,/silent,/DSCALE)
	time_lc_le = counts_le_0.time
	lc_le_0 = counts_le_0.COUNTS  &  err_lc_le_0 = counts_le_0.ERROR
	lc_le_1 = counts_le_1.COUNTS  &  err_lc_le_1 = counts_le_1.ERROR
	lc_le_2 = counts_le_2.COUNTS  &  err_lc_le_2 = counts_le_2.ERROR
	
	n_time_le = n_elements(time_lc_le)
	
	;------------------------
	GTIDesc_le_0 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g0*.lc',3,header,/silent)
	GTIDesc_le_1 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g1*.lc',3,header,/silent)
	GTIDesc_le_2 = mrdfits('LE_'+energy_LE+'keV/P021406400301_LE_'+energy_LE+'_g2*.lc',3,header,/silent)
	pixel_num_le_0 = total(GTIDesc_le_0.PIXEL,/INTEGER)
	pixel_num_le_1 = total(GTIDesc_le_1.PIXEL,/INTEGER)
	pixel_num_le_2 = total(GTIDesc_le_2.PIXEL,/INTEGER)
	
	print, 'Num of Good Pixel (LE):', pixel_num_le_0, pixel_num_le_1, pixel_num_le_2
	
	factor_le_0 = 20.d/pixel_num_le_0
	factor_le_1 = 20.d/pixel_num_le_1
	factor_le_2 = 20.d/pixel_num_le_2
	
	;------------------------------------------------
	counts_me_0 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g0*.lc',1,header,/silent,/DSCALE)   
	counts_me_1 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g1*.lc',1,header,/silent,/DSCALE)
	counts_me_2 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g2*.lc',1,header,/silent,/DSCALE)
	time_lc_me = counts_me_0.time
	lc_me_0 = counts_me_0.COUNTS  &  err_lc_me_0 = counts_me_0.ERROR
	lc_me_1 = counts_me_1.COUNTS  &  err_lc_me_1 = counts_me_1.ERROR
	lc_me_2 = counts_me_2.COUNTS  &  err_lc_me_2 = counts_me_2.ERROR
	
	n_time_me = n_elements(time_lc_me)
	
	;------------
	GTIDesc_me_0 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g0*.lc',3,header,/silent)
	GTIDesc_me_1 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g1*.lc',3,header,/silent)
	GTIDesc_me_2 = mrdfits('ME_'+energy_ME+'keV/P021406400301_ME_'+energy_ME+'_g2*.lc',3,header,/silent)
	pixel_num_me_0 = total(GTIDesc_me_0.PIXEL,/INTEGER)
	pixel_num_me_1 = total(GTIDesc_me_1.PIXEL,/INTEGER)
	pixel_num_me_2 = total(GTIDesc_me_2.PIXEL,/INTEGER)
	
	print, 'Num of Good Pixel (ME):', pixel_num_me_0, pixel_num_me_1, pixel_num_me_2
	
	factor_me_0 = 480.d/pixel_num_me_0
	factor_me_1 = 480.d/pixel_num_me_1
	factor_me_2 = 480.d/pixel_num_me_2
	
	;------------------------
	;------------------------
	data = read_data('burst_P021406400301_le_and_me.txt')
	n_row = n_elements(data[0,*])
	time_bkg_le_0_a = data[6,1:n_row-1];data[0,1:n_row-1]
	time_bkg_le_0_b = data[7,1:n_row-1];data[1,1:n_row-1]
	time_burst_le_a = data[8,1:n_row-1];data[2,1:n_row-1]
	time_burst_le_b = data[9,1:n_row-1];data[3,1:n_row-1]
	time_bkg_le_1_a = data[10,1:n_row-1];data[4,1:n_row-1]
	time_bkg_le_1_b = data[11,1:n_row-1];data[5,1:n_row-1]
	time_bkg_me_0_a = data[6,1:n_row-1]
	time_bkg_me_0_b = data[7,1:n_row-1]
	time_burst_me_a = data[8,1:n_row-1]
	time_burst_me_b = data[9,1:n_row-1]
	time_bkg_me_1_a = data[10,1:n_row-1]
	time_bkg_me_1_b = data[11,1:n_row-1]
	
	n_burst = n_row - 1
	
	n_burst_le = n_elements(time_burst_le_a)
	n_burst_me = n_elements(time_burst_me_a)
	
	;------
	for i=0, n_burst-1 do begin
		if time_bkg_le_0_a[i] lt time_bkg_le_0_b[i] - 60.d then time_bkg_le_0_a[i] = time_bkg_le_0_b[i] - 60.d
		if time_bkg_le_1_b[i] gt time_bkg_le_1_a[i] + 60.d then time_bkg_le_1_b[i] = time_bkg_le_1_a[i] + 60.d
		if time_bkg_me_0_a[i] lt time_bkg_me_0_b[i] - 60.d then time_bkg_me_0_a[i] = time_bkg_me_0_b[i] - 60.d
		if time_bkg_me_1_b[i] gt time_bkg_me_1_a[i] + 60.d then time_bkg_me_1_b[i] = time_bkg_me_1_a[i] + 60.d
	endfor
	;------
		
	if time_lc_le[0] eq data[0,0] then begin
		delta_time_le = time_lc_le[0]
	endif else begin
		print, 'LE data and burst list is not matching!'
	endelse
	
	if time_lc_me[0] eq data[6,0] then begin
		delta_time_me = time_lc_me[0]
	endif else begin
		print, 'ME data and burst list is not matching!'
	endelse
	;--------------------use ME start lc time------------
        delta_time_le = delta_time_me	
	;------------------------------------------------
	time_le = time_lc_le - delta_time_le
	
	dex_burst_le_0 = where(time_le ge time_burst_le_a[0] and time_le le time_burst_le_b[0])
	lc_burst_le_0 = lc_le_0[dex_burst_le_0]
	lc_burst_le_1 = lc_le_1[dex_burst_le_0]
	lc_burst_le_2 = lc_le_2[dex_burst_le_0]
	for i=1, n_burst_le-1 do begin
		dex_burst_le_tem = where(time_le ge time_burst_le_a[i] and time_le le time_burst_le_b[i])
		if dex_burst_le_tem[0] ne -1 then begin
			lc_burst_le_0 = [lc_burst_le_0, lc_le_0[dex_burst_le_tem]]
			lc_burst_le_1 = [lc_burst_le_1, lc_le_1[dex_burst_le_tem]]
			lc_burst_le_2 = [lc_burst_le_2, lc_le_2[dex_burst_le_tem]]
		endif
	endfor
	
	;------------------------
	time_me = time_lc_me - delta_time_me
	
	dex_burst_me_0 = where(time_me ge time_burst_me_a[0] and time_me le time_burst_me_b[0])
	lc_burst_me_0 = lc_me_0[dex_burst_me_0]
	lc_burst_me_1 = lc_me_1[dex_burst_me_0]
	lc_burst_me_2 = lc_me_2[dex_burst_me_0]
	for i=1, n_burst_me-1 do begin
		dex_burst_me_tem = where(time_me ge time_burst_me_a[i] and time_me le time_burst_me_b[i])
		if dex_burst_me_tem[0] ne -1 then begin
			lc_burst_me_0 = [lc_burst_me_0, lc_le_0[dex_burst_me_tem]]
			lc_burst_me_1 = [lc_burst_me_1, lc_le_1[dex_burst_me_tem]]
			lc_burst_me_2 = [lc_burst_me_2, lc_le_2[dex_burst_me_tem]]
		endif
	endfor
	
	;------------------------------------------------
    ; window,0,xsize=1200,ysize=460,xpos=0,ypos=420
    ; ;------------------------
    ; !P.MULTI = [0, 1, 2, 0, 1]
    ; !p.FONT=-1
	; !X.Margin=[6,6]
	; ;------------------------
	; plotsym,0,1,/FILL
	; ;------------------------
	; xr_min = min([time_lc_le,time_lc_me])
	; xr_max = max([time_lc_le,time_lc_me])
	; yr_min = 0
	; yr_max = max([lc_burst_le_0,lc_burst_le_1,lc_burst_le_2]) * 1.5
	; ;------------------------
	; multiplot,/doyaxis
	; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
	; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.72,0.96,0.96], $
	; 	  ytitle='!6Count Rate (cts/s)!N', /nodata	
	; plots, time_lc_le, lc_le_0, thick=1, color='0000ff'xl
	; plots, time_lc_le, lc_le_1, thick=1, color='00ff00'xl
	; plots, time_lc_le, lc_le_2, thick=1, color='ff0000'xl
	; for i=0, n_time_le-1 do begin
	; 	plots, [time_lc_le[i],time_lc_le[i]], [max([yr_min,lc_le_0[i]-err_lc_le_0[i]]),min([yr_max,lc_le_0[i]+err_lc_le_0[i]])], thick=1.6, color='0000ff'xl
	; 	plots, [time_lc_le[i],time_lc_le[i]], [max([yr_min,lc_le_1[i]-err_lc_le_1[i]]),min([yr_max,lc_le_1[i]+err_lc_le_1[i]])], thick=1.6, color='00ff00'xl
	; 	plots, [time_lc_le[i],time_lc_le[i]], [max([yr_min,lc_le_2[i]-err_lc_le_2[i]]),min([yr_max,lc_le_2[i]+err_lc_le_2[i]])], thick=1.6, color='ff0000'xl
	; endfor
	; ;------------------------
	; xr_min = min([time_lc_le,time_lc_me])
	; xr_max = max([time_lc_le,time_lc_me])
	; yr_min = 0
	; yr_max = max([lc_burst_me_0,lc_burst_me_1,lc_burst_me_2]) * 1.5
	; ;------------------------
	; multiplot,/doxaxis,/doyaxis
	; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
	; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.42,0.96,0.66], $
	; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata	
	; plots, time_lc_me, lc_me_0, thick=1, color='0000ff'xl
	; plots, time_lc_me, lc_me_1, thick=1, color='00ff00'xl
	; plots, time_lc_me, lc_me_2, thick=1, color='ff0000'xl
	; for i=0, n_time_le-1 do begin
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_0[i]-err_lc_me_0[i]]),min([yr_max,lc_me_0[i]+err_lc_me_0[i]])], thick=1.6, color='0000ff'xl
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_1[i]-err_lc_me_1[i]]),min([yr_max,lc_me_1[i]+err_lc_me_1[i]])], thick=1.6, color='00ff00'xl
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_2[i]-err_lc_me_2[i]]),min([yr_max,lc_me_2[i]+err_lc_me_2[i]])], thick=1.6, color='ff0000'xl
	; endfor
	; ;------------------------
	; multiplot,/doxaxis,/doyaxis
	; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
	; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.12,0.96,0.36], $
	; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata	
	; plots, time_lc_me, lc_me_0, thick=1, color='0000ff'xl
	; plots, time_lc_me, lc_me_1, thick=1, color='00ff00'xl
	; plots, time_lc_me, lc_me_2, thick=1, color='ff0000'xl
	; for i=0, n_time_le-1 do begin
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_0[i]-err_lc_me_0[i]]),min([yr_max,lc_me_0[i]+err_lc_me_0[i]])], thick=1.6, color='0000ff'xl
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_1[i]-err_lc_me_1[i]]),min([yr_max,lc_me_1[i]+err_lc_me_1[i]])], thick=1.6, color='00ff00'xl
	; 	plots, [time_lc_me[i],time_lc_me[i]], [max([yr_min,lc_me_2[i]-err_lc_me_2[i]]),min([yr_max,lc_me_2[i]+err_lc_me_2[i]])], thick=1.6, color='ff0000'xl
	; endfor
	
		
	;------------------------------------------------	
	time_le = time_lc_le - delta_time_le
	time_me = time_lc_me - delta_time_me
	
	dex_src_le = where(time_le ge time_burst_le_a[burst_dex] and time_le le time_burst_le_b[burst_dex])
	dex_bkg_le = where((time_le ge time_bkg_le_0_a[burst_dex] and time_le lt time_bkg_le_0_b[burst_dex]) or $
					   (time_le gt time_bkg_le_1_a[burst_dex] and time_le le time_bkg_le_1_b[burst_dex]))
	
	dex_src_me = where(time_me ge time_burst_me_a[burst_dex] and time_me le time_burst_me_b[burst_dex])
	dex_bkg_me = where((time_me ge time_bkg_me_0_a[burst_dex] and time_me lt time_bkg_me_0_b[burst_dex]) or $
					   (time_me gt time_bkg_me_1_a[burst_dex] and time_me le time_bkg_me_1_b[burst_dex]))
	
	if dex_src_le[0] ne -1 and dex_src_me[0] ne -1 then begin
		   
		time_burst_le = time_le[dex_src_le]
		time_burst_me = time_me[dex_src_me]
		
		lc_LE_src_0 = lc_le_0[dex_src_le]  &  err_lc_LE_src_0 = err_lc_le_0[dex_src_le]    &    lc_le_bkg_0 = lc_le_0[dex_bkg_le]  &  err_lc_le_bkg_0 = err_lc_le_0[dex_bkg_le]
		lc_LE_src_1 = lc_le_1[dex_src_le]  &  err_lc_LE_src_1 = err_lc_le_1[dex_src_le]    &    lc_le_bkg_1 = lc_le_1[dex_bkg_le]  &  err_lc_le_bkg_1 = err_lc_le_1[dex_bkg_le]
		lc_LE_src_2 = lc_le_2[dex_src_le]  &  err_lc_LE_src_2 = err_lc_le_2[dex_src_le]    &    lc_le_bkg_2 = lc_le_2[dex_bkg_le]  &  err_lc_le_bkg_2 = err_lc_le_2[dex_bkg_le]
		
		lc_me_src_0 = lc_me_0[dex_src_me]  &  err_lc_me_src_0 = err_lc_me_0[dex_src_me]    &    lc_me_bkg_0 = lc_me_0[dex_bkg_me]  &  err_lc_me_bkg_0 = err_lc_me_0[dex_bkg_me]
		lc_me_src_1 = lc_me_1[dex_src_me]  &  err_lc_me_src_1 = err_lc_me_1[dex_src_me]    &    lc_me_bkg_1 = lc_me_1[dex_bkg_me]  &  err_lc_me_bkg_1 = err_lc_me_1[dex_bkg_me]
		lc_me_src_2 = lc_me_2[dex_src_me]  &  err_lc_me_src_2 = err_lc_me_2[dex_src_me]    &    lc_me_bkg_2 = lc_me_2[dex_bkg_me]  &  err_lc_me_bkg_2 = err_lc_me_2[dex_bkg_me]
		
		n_time_src_le = n_elements(lc_le_src_0)  &  n_time_bkg_le = n_elements(lc_le_bkg_0)
		n_time_src_me = n_elements(lc_me_src_0)  &  n_time_bkg_me = n_elements(lc_me_bkg_0)
		
		;;------
		;burst_LE_0 = (lc_le_src_0 - mean(lc_le_bkg_0)) * factor_le_0  &  err_burst_LE_0 = (err_lc_le_src_0^2 + total(lc_le_bkg_0)/n_time_bkg_le^2)^0.5 * factor_le_0
		;burst_LE_1 = (lc_le_src_1 - mean(lc_le_bkg_1)) * factor_le_1  &  err_burst_LE_1 = (err_lc_le_src_1^2 + total(lc_le_bkg_1)/n_time_bkg_le^2)^0.5 * factor_le_1
		;burst_LE_2 = (lc_le_src_2 - mean(lc_le_bkg_2)) * factor_le_2  &  err_burst_LE_2 = (err_lc_le_src_2^2 + total(lc_le_bkg_2)/n_time_bkg_le^2)^0.5 * factor_le_2
		;
		;burst_ME_0 = (lc_me_src_0 - mean(lc_me_bkg_0)) * factor_me_0  &  err_burst_ME_0 = (err_lc_me_src_0^2 + total(lc_me_bkg_0)/n_time_bkg_me^2)^0.5 * factor_me_0
		;burst_ME_1 = (lc_me_src_1 - mean(lc_me_bkg_1)) * factor_me_1  &  err_burst_ME_1 = (err_lc_me_src_1^2 + total(lc_me_bkg_1)/n_time_bkg_me^2)^0.5 * factor_me_1
		;burst_ME_2 = (lc_me_src_2 - mean(lc_me_bkg_2)) * factor_me_2  &  err_burst_ME_2 = (err_lc_me_src_2^2 + total(lc_me_bkg_2)/n_time_bkg_me^2)^0.5 * factor_me_2
		
		;------
		burst_LE_0 = (total(lc_le_src_0)/n_time_src_le - total(lc_le_bkg_0)/n_time_bkg_le) * factor_le_0
		burst_LE_1 = (total(lc_le_src_1)/n_time_src_le - total(lc_le_bkg_1)/n_time_bkg_le) * factor_le_1
		burst_LE_2 = (total(lc_le_src_2)/n_time_src_le - total(lc_le_bkg_2)/n_time_bkg_le) * factor_le_2
		err_burst_LE_0 = (total(lc_le_src_0)/n_time_src_le^2 + total(lc_le_bkg_0)/n_time_bkg_le^2)^0.5 * factor_le_0
		err_burst_LE_1 = (total(lc_le_src_1)/n_time_src_le^2 + total(lc_le_bkg_1)/n_time_bkg_le^2)^0.5 * factor_le_1
		err_burst_LE_2 = (total(lc_le_src_2)/n_time_src_le^2 + total(lc_le_bkg_2)/n_time_bkg_le^2)^0.5 * factor_le_2
		
		burst_ME_0 = (total(lc_me_src_0)/n_time_src_me - total(lc_me_bkg_0)/n_time_bkg_me) * factor_me_0
		burst_ME_1 = (total(lc_me_src_1)/n_time_src_me - total(lc_me_bkg_1)/n_time_bkg_me) * factor_me_1
		burst_ME_2 = (total(lc_me_src_2)/n_time_src_me - total(lc_me_bkg_2)/n_time_bkg_me) * factor_me_2
		err_burst_ME_0 = (total(err_lc_me_src_0^2)/n_time_src_me^2 + total(err_lc_me_bkg_0^2)/n_time_bkg_me^2)^0.5 * factor_me_0
		err_burst_ME_1 = (total(err_lc_me_src_1^2)/n_time_src_me^2 + total(err_lc_me_bkg_1^2)/n_time_bkg_me^2)^0.5 * factor_me_1
		err_burst_ME_2 = (total(err_lc_me_src_2^2)/n_time_src_me^2 + total(err_lc_me_bkg_2^2)/n_time_bkg_me^2)^0.5 * factor_me_2
		
		print, 'Burst (LE):'
		print, burst_LE_0, err_burst_LE_0
		print, burst_LE_1, err_burst_LE_1
		print, burst_LE_2, err_burst_LE_2
		print, 'Burst (ME):'
		print, burst_ME_0, err_burst_ME_0
		print, burst_ME_1, err_burst_ME_1
		print, burst_ME_2, err_burst_ME_2
		
		;------------------------------------------------------------------------------------------------
		;------------------------------------------------------------------------------------------------
    	; window,1,xsize=600,ysize=460,xpos=0,ypos=320
    	; ;------------------------
    	; !P.MULTI = [0, 1, 2, 0, 1]
    	; !p.FONT=-1
		; !X.Margin=[6,6]
		; ;------------------------
		; plotsym,0,1,/FILL
		; ;------------------------
		; xr_min = min([time_burst_le]) - 120
		; xr_max = max([time_burst_le]) + 120
		; yr_min = 0
		; yr_max = max([lc_LE_src_0,lc_LE_src_1,lc_LE_src_2]) * 1.5
		; ;------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
		; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.60,0.96,0.96], $
		; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata
		
		; dex = where(time_le ge xr_min and time_le le xr_max)
		; time_p = time_le[dex]
		; lc_0_p = lc_le_0[dex]  &  err_lc_0_p = err_lc_le_0[dex]
		; lc_1_p = lc_le_1[dex]  &  err_lc_1_p = err_lc_le_1[dex]
		; lc_2_p = lc_le_2[dex]  &  err_lc_2_p = err_lc_le_2[dex]
		
		; plots, time_p, lc_0_p, thick=1, color='0000ff'xl
		; plots, time_p, lc_1_p, thick=1, color='00ff00'xl
		; plots, time_p, lc_2_p, thick=1, color='ff0000'xl
		; for i=0, n_elements(dex)-1 do begin
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_0_p[i]-err_lc_0_p[i]]),min([yr_max,lc_0_p[i]+err_lc_0_p[i]])], thick=1.6, color='0000ff'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_1_p[i]-err_lc_1_p[i]]),min([yr_max,lc_1_p[i]+err_lc_1_p[i]])], thick=1.6, color='00ff00'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_2_p[i]-err_lc_2_p[i]]),min([yr_max,lc_2_p[i]+err_lc_2_p[i]])], thick=1.6, color='ff0000'xl
		; endfor
		
		; plots, [time_burst_le_a[burst_dex],time_burst_le_a[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; plots, [time_burst_le_b[burst_dex],time_burst_le_b[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; if time_bkg_le_0_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_le_0_a[burst_dex],time_bkg_le_0_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_le_0_b[burst_dex],time_bkg_le_0_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		; if time_bkg_le_1_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_le_1_a[burst_dex],time_bkg_le_1_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_le_1_b[burst_dex],time_bkg_le_1_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		
		; ;------------------------
		; xr_min = min([time_burst_me]) - 120
		; xr_max = max([time_burst_me]) + 120
		; yr_min = 0
		; yr_max = max([lc_ME_src_0,lc_ME_src_1,lc_ME_src_2]) * 1.5
		; ;------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
		; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.12,0.96,0.48], $
		; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata
		
		; dex = where(time_me ge xr_min and time_me le xr_max)
		; time_p = time_me[dex]
		; lc_0_p = lc_me_0[dex]  &  err_lc_0_p = err_lc_me_0[dex]
		; lc_1_p = lc_me_1[dex]  &  err_lc_1_p = err_lc_me_1[dex]
		; lc_2_p = lc_me_2[dex]  &  err_lc_2_p = err_lc_me_2[dex]
		
		; plots, time_p, lc_0_p, thick=1, color='0000ff'xl
		; plots, time_p, lc_1_p, thick=1, color='00ff00'xl
		; plots, time_p, lc_2_p, thick=1, color='ff0000'xl
		; for i=0, n_elements(dex)-1 do begin
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_0_p[i]-err_lc_0_p[i]]),min([yr_max,lc_0_p[i]+err_lc_0_p[i]])], thick=1.6, color='0000ff'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_1_p[i]-err_lc_1_p[i]]),min([yr_max,lc_1_p[i]+err_lc_1_p[i]])], thick=1.6, color='00ff00'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_2_p[i]-err_lc_2_p[i]]),min([yr_max,lc_2_p[i]+err_lc_2_p[i]])], thick=1.6, color='ff0000'xl
		; endfor
		
		; plots, [time_burst_me_a[burst_dex],time_burst_me_a[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; plots, [time_burst_me_b[burst_dex],time_burst_me_b[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; if time_bkg_me_0_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_me_0_a[burst_dex],time_bkg_me_0_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_me_0_b[burst_dex],time_bkg_me_0_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		; if time_bkg_me_1_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_me_1_a[burst_dex],time_bkg_me_1_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_me_1_b[burst_dex],time_bkg_me_1_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		
		; ;------------------------------------------------------------------------------------------------
		; ;------------------------------------------------------------------------------------------------
    	; window,2,xsize=600,ysize=460,xpos=600,ypos=320
    	; ;------------------------
    	; !P.MULTI = [0, 1, 2, 0, 1]
    	; !p.FONT=-1
		; !X.Margin=[6,6]
		; ;------------------------
		; plotsym,0,1,/FILL
		; ;------------------------
		; xr_min = min([time_burst_le]) - 60
		; xr_max = max([time_burst_le]) + 60
		; yr_min = -10
		; yr_max = max([lc_LE_src_0,lc_LE_src_1,lc_LE_src_2]) * factor_le_1/1.6889557 * 1.2
		; ;------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
		; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.60,0.96,0.96], $
		; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata
		
		; dex = where(time_le ge xr_min and time_le le xr_max)
		; time_p = time_le[dex]
		; lc_0_p = (lc_le_0[dex] - mean(lc_le_bkg_0)) * factor_le_0/1.6944942  &  err_lc_0_p = (err_lc_le_0[dex]^2 + total(lc_le_bkg_0)/n_time_bkg_le^2)^0.5 * factor_le_0/1.6944942
		; lc_1_p = (lc_le_1[dex] - mean(lc_le_bkg_1)) * factor_le_1/1.6889557  &  err_lc_1_p = (err_lc_le_1[dex]^2 + total(lc_le_bkg_1)/n_time_bkg_le^2)^0.5 * factor_le_1/1.6889557
		; lc_2_p = (lc_le_2[dex] - mean(lc_le_bkg_2)) * factor_le_2/1.6251700  &  err_lc_2_p = (err_lc_le_2[dex]^2 + total(lc_le_bkg_2)/n_time_bkg_le^2)^0.5 * factor_le_2/1.6251700
		; print, mean(lc_le_bkg_0), mean(lc_le_0[dex_bkg_le]), lc_le_0[dex_bkg_le]
		; plots, time_p, lc_0_p, thick=1, color='0000ff'xl
		; plots, time_p, lc_1_p, thick=1, color='00ff00'xl
		; plots, time_p, lc_2_p, thick=1, color='ff0000'xl
		; for i=0, n_elements(dex)-1 do begin
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_0_p[i]-err_lc_0_p[i]]),min([yr_max,lc_0_p[i]+err_lc_0_p[i]])], thick=1.6, color='0000ff'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_1_p[i]-err_lc_1_p[i]]),min([yr_max,lc_1_p[i]+err_lc_1_p[i]])], thick=1.6, color='00ff00'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_2_p[i]-err_lc_2_p[i]]),min([yr_max,lc_2_p[i]+err_lc_2_p[i]])], thick=1.6, color='ff0000'xl
		; endfor
		
		; plots, [time_burst_le_a[burst_dex],time_burst_le_a[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; plots, [time_burst_le_b[burst_dex],time_burst_le_b[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; if time_bkg_le_0_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_le_0_a[burst_dex],time_bkg_le_0_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_le_0_b[burst_dex],time_bkg_le_0_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		; if time_bkg_le_1_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_le_1_a[burst_dex],time_bkg_le_1_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_le_1_b[burst_dex],time_bkg_le_1_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		
		; plots, [xr_min,xr_max], [0,0], linestyle=2
		
		; ;------------------------
		; xr_min = min([time_burst_me]) - 60
		; xr_max = max([time_burst_me]) + 60
		; yr_min = -10
		; yr_max = max([lc_ME_src_0,lc_ME_src_1,lc_ME_src_2]) * factor_me_1/0.818252 * 1.2
		; ;------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot, [0], [0], psym=3, symsize=1, xrange=[xr_min,xr_max], yrange=[yr_min,yr_max], xstyle=1, ystyle=1, $
		; 	  charthick=1.6, xthick=1.6, ythick=1.6, xcharsize=1.6, ycharsize=1.6, position=[0.12,0.12,0.96,0.48], $
		; 	  xtitle='!6Time (s)!N', ytitle='!6Count Rate (cts/s)!N', /nodata
		
		; dex = where(time_me ge xr_min and time_me le xr_max)
		; time_p = time_me[dex]
		; lc_0_p = (lc_me_0[dex] - mean(lc_me_bkg_0)) * factor_me_0/0.979319  &  err_lc_0_p = (err_lc_me_0[dex]^2 + total(lc_me_bkg_0)/n_time_bkg_me^2)^0.5 * factor_me_0/0.979319
		; lc_1_p = (lc_me_1[dex] - mean(lc_me_bkg_1)) * factor_me_1/0.818252  &  err_lc_1_p = (err_lc_me_1[dex]^2 + total(lc_me_bkg_1)/n_time_bkg_me^2)^0.5 * factor_me_1/0.818252
		; lc_2_p = (lc_me_2[dex] - mean(lc_me_bkg_2)) * factor_me_2/0.831341  &  err_lc_2_p = (err_lc_me_2[dex]^2 + total(lc_me_bkg_2)/n_time_bkg_me^2)^0.5 * factor_me_2/0.831341
		
		; plots, time_p, lc_0_p, thick=1, color='0000ff'xl
		; plots, time_p, lc_1_p, thick=1, color='00ff00'xl
		; plots, time_p, lc_2_p, thick=1, color='ff0000'xl
		; for i=0, n_elements(dex)-1 do begin
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_0_p[i]-err_lc_0_p[i]]),min([yr_max,lc_0_p[i]+err_lc_0_p[i]])], thick=1.6, color='0000ff'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_1_p[i]-err_lc_1_p[i]]),min([yr_max,lc_1_p[i]+err_lc_1_p[i]])], thick=1.6, color='00ff00'xl
		; 	plots, [time_p[i],time_p[i]], [max([yr_min,lc_2_p[i]-err_lc_2_p[i]]),min([yr_max,lc_2_p[i]+err_lc_2_p[i]])], thick=1.6, color='ff0000'xl
		; endfor
		
		; plots, [time_burst_me_a[burst_dex],time_burst_me_a[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; plots, [time_burst_me_b[burst_dex],time_burst_me_b[burst_dex]], [yr_min,yr_max], linestyle=2, thick=2
		; if time_bkg_me_0_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_me_0_a[burst_dex],time_bkg_me_0_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_me_0_b[burst_dex],time_bkg_me_0_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		; if time_bkg_me_1_a[burst_dex] ne -1 then begin
		; 	plots, [time_bkg_me_1_a[burst_dex],time_bkg_me_1_a[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; 	plots, [time_bkg_me_1_b[burst_dex],time_bkg_me_1_b[burst_dex]], [yr_min,yr_max], linestyle=1, thick=2
		; endif
		
		; plots, [xr_min,xr_max], [0,0], linestyle=2
		
		;------------------------------------------------------------------------------------------------
		;------------------------------------------------------------------------------------------------
		att = mrdfits(strcompress('ATT/HXMT_P0214064003_Att_FFFFFF_V2_L1P.FITS',/remove),3,header,/silent)
		time_att = double(att.time)
		Q1_att = att.Q1
		Q2_att = att.Q2
		Q3_att = att.Q3
		
		time_att_le = time_att - delta_time_le
		time_att_me = time_att - delta_time_me
		
		;------------------------
		Q1_le = INTERPOL(Q1_att, time_att_le, time_burst_le)
		Q2_le = INTERPOL(Q2_att, time_att_le, time_burst_le)
		Q3_le = INTERPOL(Q3_att, time_att_le, time_burst_le)
		
		n_time_le = n_elements(time_burst_le)
		
		;------------------------
		Q1_me = INTERPOL(Q1_att, time_att_me, time_burst_me)
		Q2_me = INTERPOL(Q2_att, time_att_me, time_burst_me)
		Q3_me = INTERPOL(Q3_att, time_att_me, time_burst_me)
		
		n_time_me = n_elements(time_burst_me)
		
		;;------------------------------------------------------------------------
		;;------------------------
		;DL_LE_1_6_box0_nangyi = [-0.175520176,       -5.25922d-05]
		;DL_LE_1_6_box1_nangyi = [-0.131411876,       -0.000280339]
		;DL_LE_1_6_box2_nangyi = [-0.143641291,       -0.002861802]
		;
		;euler_LE_1_6_box0_nangyi = [0.026018636,       0.066218707,      0.162763242]
		;euler_LE_1_6_box1_nangyi = [-0.001554774,       0.11566053,       0.112678894]
		;euler_LE_1_6_box2_nangyi = [-0.047662386,       0.081942721,      0.173199497]
		;
		;par_LE_1_6_box0_nangyi = [-0.001401647, -0.010797473, -0.051066905, 0.974555617]
		;par_LE_1_6_box1_nangyi = [-0.000922986, -0.012016397, -0.042864763, 1.014102451]
		;par_LE_1_6_box2_nangyi = [-0.001025824, -0.012060025, -0.052200055, 0.967361821]
		;
		;
		;------------------------
		DL_LE_1_6_box0 = [2.2721250,      0.40446761]
		DL_LE_1_6_box1 = [3.7387836,      0.69034524]
		DL_LE_1_6_box2 = [6.2906011,      0.39417481]
		
		euler_LE_1_6_box0 = [ 0.054455418,     0.077059705,    -0.010076819]
		euler_LE_1_6_box1 = [0.0075468207,     0.094741889,      0.35217648]
		euler_LE_1_6_box2 = [-0.052166667,     0.091302223,     -0.31558935]
		
		par_LE_1_6_box0 = [ 0.13719193,    -0.078033840,     -0.49036001,       1.6944942]
		par_LE_1_6_box1 = [0.064318660,     -0.12834982,     -0.61765350,       1.6889557]
		par_LE_1_6_box2 = [ 0.16017803,     -0.14397490,     -0.45583557,       1.6251700]
		
		;------------------------
		DL_LE_2_6_box0 = [3.922910,       0.399965]
		DL_LE_2_6_box1 = [2.399953,       0.602995]
		DL_LE_2_6_box2 = [2.929449,       0.238921]
		
		euler_LE_2_6_box0 = [0.016947,       0.054199,      -0.081497]
		euler_LE_2_6_box1 = [0.006143,       0.096848,       1.080478]
		euler_LE_2_6_box2 = [-0.074942,       0.107840,      -0.171996]
		
		par_LE_2_6_box0 = [0.006565, -0.042149, -0.239915, 0.891241]
		par_LE_2_6_box1 = [-0.005233, -0.050597, -0.305286, 0.899709]
		par_LE_2_6_box2 = [0.119697, -0.050229, -0.177212, 0.861934]
		
		;------------------------		
		DL_ME_7_12_box0 = [-0.21636257d,      0.77772640]
		DL_ME_7_12_box1 = [  6.5566438d,      0.94111247]
		DL_ME_7_12_box2 = [-0.39144551d,      0.72633634]
		
		euler_ME_7_12_box0 = [ 0.099142223d,    -0.048153759,     -0.50616208]
		euler_ME_7_12_box1 = [-0.030580797d,    -0.082359521,      0.61757356]
		euler_ME_7_12_box2 = [-0.058328317d,     -0.10376456,      0.55885838]
		
		par_ME_7_12_box0 = [-0.044484078d,   -0.0080666755,      -1.0661968,       1.0678338]
		par_ME_7_12_box1 = [ -0.51163237d,     -0.13021961,     -0.80195804,      0.81726050]
		par_ME_7_12_box2 = [ -0.19502559d,     0.081869619,     -0.94419083,      0.97263888]
		
		;------------------------
		; DL_ME_7_20_box0 = ;[-8.296641658550059084e-01,];[-0.69488468d,     0.031676293]
		; DL_ME_7_20_box1 = [ 2.01803950d,    -0.010488884]
		; DL_ME_7_20_box2 = [-0.80747338d,    -0.035860595]
		
		; euler_ME_7_20_box0 = [8.951564827392798074e-02,-3.764336043176106511e-02,-1.293231439483251810e+00];[ 0.092635436d,    -0.037681532,      -1.3641714]	
		; euler_ME_7_20_box1 = [-0.037142052d,    -0.086510975,       1.0368320]
		; euler_ME_7_20_box2 = [-0.092051502d,     -0.11820098,      0.47091772]
		
		; par_ME_7_20_box0 = [  0.23724572d,     0.028623222,    -0.046247147,      1.08326760]
		; par_ME_7_20_box1 = [-0.079325043d,    -0.093474350,     0.091077681,      0.88234257]
		; par_ME_7_20_box2 = [  0.18048410d,     0.099220646,      0.23518466,      0.96674312]
		
		;------------------------
		if energy_LE eq '1-6' then begin
			DL_LE_0 = DL_LE_1_6_box0  &  euler_LE_0 = euler_LE_1_6_box0  &  par_LE_0 = par_LE_1_6_box0
			DL_LE_1 = DL_LE_1_6_box1  &  euler_LE_1 = euler_LE_1_6_box1  &  par_LE_1 = par_LE_1_6_box1
			DL_LE_2 = DL_LE_1_6_box2  &  euler_LE_2 = euler_LE_1_6_box2  &  par_LE_2 = par_LE_1_6_box2
		endif else begin
			DL_LE_0 = DL_LE_2_6_box0  &  euler_LE_0 = euler_LE_2_6_box0  &  par_LE_0 = par_LE_2_6_box0
			DL_LE_1 = DL_LE_2_6_box1  &  euler_LE_1 = euler_LE_2_6_box1  &  par_LE_1 = par_LE_2_6_box1
			DL_LE_2 = DL_LE_2_6_box2  &  euler_LE_2 = euler_LE_2_6_box2  &  par_LE_2 = par_LE_2_6_box2
		endelse
		
		if energy_LE eq '7-12' then begin
			DL_ME_0 = DL_ME_7_12_box0  &  euler_ME_0 = euler_ME_7_12_box0  &  par_ME_0 = par_ME_7_12_box0
			DL_ME_1 = DL_ME_7_12_box1  &  euler_ME_1 = euler_ME_7_12_box1  &  par_ME_1 = par_ME_7_12_box1
			DL_ME_2 = DL_ME_7_12_box2  &  euler_ME_2 = euler_ME_7_12_box2  &  par_ME_2 = par_ME_7_12_box2
		endif else begin
			DL_ME_0 = DL_ME_7_20_box0  &  euler_ME_0 = euler_ME_7_20_box0  &  par_ME_0 = par_ME_7_20_box0
			DL_ME_1 = DL_ME_7_20_box1  &  euler_ME_1 = euler_ME_7_20_box1  &  par_ME_1 = par_ME_7_20_box1
			DL_ME_2 = DL_ME_7_20_box2  &  euler_ME_2 = euler_ME_7_20_box2  &  par_ME_2 = par_ME_7_20_box2
		endelse
    	
		
		;------------------------------------ Fit (MCMC) ------------------------------------
		;n_step=1d+06
		
		;RA_initial = 264.0
		;Dec_initial = -33.0
		RA_initial = 262.991
		Dec_initial = -33.834
		flux_le_initial = 100
		flux_me_initial = 100
		
		p_initial = [RA_initial, Dec_initial, flux_le_initial, flux_me_initial]
				
		p_new=dblarr(n_step+1,n_elements(p_initial))
		chi2=dblarr(n_step+1)
		p_new[0,*]=p_initial
		
		dof = n_elements([burst_le_0,burst_le_1,burst_le_2,burst_me_0,burst_me_1,burst_me_2])-n_elements(p_initial)
		
		chi2[0] = 1d+08
		
		p_random = double([0.1, 0.1,  2,  2])
		p_fix    = double([  1,   1,  1,  1])
		
		p_random = p_random * p_fix
		
		k=0L
		for ii=0L, n_step-1 do begin
			
			p_new_tem = p_new[k,*] + p_random * randomn(seed,n_elements(p_initial),/double)
			
			RA = p_new_tem[0]
			DEC= p_new_tem[1]
			Flux_LE = p_new_tem[2]
			Flux_ME = p_new_tem[3]
			
			burst_le_0_tem = modulate_curve_le_psf_modify(RA=RA, DEC=DEC, Flux=Flux_LE, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=0, DL=DL_LE_0, euler=euler_LE_0, par=par_LE_0)
			burst_le_1_tem = modulate_curve_le_psf_modify(RA=RA, DEC=DEC, Flux=Flux_LE, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=1, DL=DL_LE_1, euler=euler_LE_1, par=par_LE_1)
			burst_le_2_tem = modulate_curve_le_psf_modify(RA=RA, DEC=DEC, Flux=Flux_LE, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=2, DL=DL_LE_2, euler=euler_LE_2, par=par_LE_2)
			
			burst_me_0_tem = modulate_curve_me_psf_modify(RA=RA, DEC=DEC, Flux=Flux_ME, Qt1=Q1_me, Qt2=Q2_me, Qt3=Q3_me, BoxID=0);, DL=DL_ME_0, euler=euler_ME_0, par=par_ME_0)
			burst_me_1_tem = modulate_curve_me_psf_modify(RA=RA, DEC=DEC, Flux=Flux_ME, Qt1=Q1_me, Qt2=Q2_me, Qt3=Q3_me, BoxID=1);, DL=DL_ME_1, euler=euler_ME_1, par=par_ME_1)
			burst_me_2_tem = modulate_curve_me_psf_modify(RA=RA, DEC=DEC, Flux=Flux_ME, Qt1=Q1_me, Qt2=Q2_me, Qt3=Q3_me, BoxID=2);, DL=DL_ME_2, euler=euler_ME_2, par=par_ME_2)
					
			;------------------------------------		
			if use_LE_only then begin
				chi2_tem = total((total(burst_le_0_tem)/n_time_src_le - burst_le_0)^2/err_burst_le_0^2) + $
						   total((total(burst_le_1_tem)/n_time_src_le - burst_le_1)^2/err_burst_le_1^2) + $
						   total((total(burst_le_2_tem)/n_time_src_le - burst_le_2)^2/err_burst_le_2^2)
			endif else if use_ME_only then begin
				chi2_tem = total((total(burst_me_0_tem)/n_time_src_me - burst_me_0)^2/err_burst_me_0^2) + $
					   	total((total(burst_me_1_tem)/n_time_src_me - burst_me_1)^2/err_burst_me_1^2) + $
					   	total((total(burst_me_2_tem)/n_time_src_me - burst_me_2)^2/err_burst_me_2^2)
			endif else begin
				chi2_tem = total((total(burst_le_0_tem)/n_time_src_le - burst_le_0)^2/err_burst_le_0^2) + $
						   total((total(burst_le_1_tem)/n_time_src_le - burst_le_1)^2/err_burst_le_1^2) + $
						   total((total(burst_le_2_tem)/n_time_src_le - burst_le_2)^2/err_burst_le_2^2) + $
						   total((total(burst_me_0_tem)/n_time_src_me - burst_me_0)^2/err_burst_me_0^2) + $
						   total((total(burst_me_1_tem)/n_time_src_me - burst_me_1)^2/err_burst_me_1^2) + $
						   total((total(burst_me_2_tem)/n_time_src_me - burst_me_2)^2/err_burst_me_2^2)
			endelse	
			;------------------------
			like_ratio=exp(-(chi2_tem-chi2[k])/2)
			randomu_tem=randomu(seed,1)
			if (chi2_tem le chi2[k] or randomu_tem lt like_ratio) then begin
							
				p_new[k+1,*]=p_new_tem
				chi2[k+1]=chi2_tem
				
				k=k+1
				;print, '            ', string(ii,      format='(I12)'),    string(k,         format='(I12)'), $
				;	          	 	   string(chi2[k], format='(f16.2)'),  string(chi2[k-1], format='(f16.2)'),  string(dof, format='(I16)')
				;print, '    DL:     ', string(p_new[k,0:1], format='(2f14.6)')
				;print, '    R_0:    ', string(p_new[k,2:4], format='(3f14.6)')
				;print, '    R_1:    ', string(p_new[k,5:7], format='(3f14.6)')
				;print, '    R_2:    ', string(p_new[k,8:10], format='(3f14.6)')
				;print, '    Par:    ', string(p_new[k,11:14], format='(4f14.6)')
				;print, ''
				;------
				print, ' ', string(ii,      format='(I12)'),    string(k,         format='(I12)'), $
							string(p_new[k,0:n_elements(p_initial)-1], format='(99f9.2)'), $
							string(chi2[k], format='(f11.2)'),  string(chi2[k-1], format='(f11.2)'),  string(dof, format='(I9)')
				
			endif
    	
		endfor
		;n_burn = 1000
		RA_fit = p_new[n_burn:k,0]
		Dec_fit = p_new[n_burn:k,1]
		flux_le_fit = p_new[n_burn:k,2]
		flux_me_fit = p_new[n_burn:k,3]
		
		chi2 = chi2[n_burn:k]
		
		;------------
    	min_chi2=min(chi2,dex)
		
		RA = RA_fit[dex]
		Dec = Dec_fit[dex]
		flux_le = flux_le_fit[dex]
		flux_me = flux_me_fit[dex]
		
		err_RA = stddev(RA_fit)
		err_Dec = stddev(Dec_fit)
		err_flux_le = stddev(flux_le_fit)
		err_flux_me = stddev(flux_me_fit)
		
		deviation = dblarr(n_elements(RA_fit))
		for i=0L, n_elements(RA_fit)-1 do begin
			deviation[i] = angle_two_source([RA,Dec],[RA_fit[i],Dec_fit[i]])
		endfor
		error_position = total(deviation^2)^0.5/n_elements(deviation)^0.5
		
		;------------------------------------ 输出拟合值到文本 ------------------------------------
		if use_LE_only then begin
			mcmc_result_txtname = 'mcmc_LE_'+energy_LE+'_burst-'+string(burst_dex)+'.txt'
                        result_txtname = 'loc_LE_'+energy_LE+'_burst-'+string(burst_dex)+'.txt'
			plot_result_txtname = 'plot_LE_'+energy_LE+'_burst-'+string(burst_dex)+'.txt'
                endif else if use_ME_only then begin
                        mcmc_result_txtname = 'mcmc_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
                        result_txtname = 'loc_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
			plot_result_txtname = 'plot_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
                endif else begin
                        mcmc_result_txtname = 'mcmc_LE_'+energy_LE+'_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
                        result_txtname = 'loc_LE_'+energy_LE+'_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
			plot_result_txtname = 'plot_LE_'+energy_LE+'_ME_'+energy_ME+'_burst-'+string(burst_dex)+'.txt'
		endelse
		if Keyword_Set(print_result) then begin
			close,6
			;openw,6,strcompress('result_le_'+energy_LE+'_burst-'+string(burst_dex)+'.txt',/remove),width=4000
			openw,6,strcompress(mcmc_result_txtname,/remove),width=4000
			printf, 6, [transpose(p_new[n_burn:k,*]), transpose(chi2), replicate(dof,[1,n_elements(chi2)])]
			close,6
		endif
		
		;------------------------
		close,2
		openw,2,strcompress(result_txtname,/remove),width=4000
		printf,2, ''
		printf,2, '  chi2/dof:', string(min_chi2,format='(f10.2)'), ' / ', string(dof,format='(I5)')
		printf,2, ''
		printf,2, '              ', '(RA,Dec)', '      ', 'Flux'
		printf,2, '------------------------------------------------------------------------------------'
		printf,2, '  (RA,Dec):  ', string(RA, format='(f9.3)'), ' +/-', string(err_RA, format='(f9.3)'), '      ', $
								string(Dec, format='(f9.3)'), ' +/-', string(err_Dec, format='(f9.3)'), '      Error:', error_position
		printf,2, '------------------------------------------------------------------------------------'
		printf,2, '  Flux_LE:  ', string(flux_le, format='(f9.3)'), ' +/-', string(err_flux_le, format='(f9.3)')
		printf,2, '  Flux_ME:  ', string(flux_me, format='(f9.3)'), ' +/-', string(err_flux_me, format='(f9.3)')
		printf,2, '------------------------------------------------------------------------------------'
		
		printf,2, burst_dex
		printf,2, burst_le_0/burst_le_1, 1, burst_le_2/burst_le_1
		close,2
		close,3
		openw,3,strcompress(plot_result_txtname,/remove),width=4000
		printf,3,string([RA,err_RA,Dec,err_Dec,error_position,burst_dex])
		close,3
		
		;Off-axis
		burst_le_0_true = total(modulate_curve_le_psf_modify(RA=262.991, DEC=-33.834, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=0, DL=DL_LE_0, euler=euler_LE_0, par=par_LE_0))
		burst_le_1_true = total(modulate_curve_le_psf_modify(RA=262.991, DEC=-33.834, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=1, DL=DL_LE_1, euler=euler_LE_1, par=par_LE_1))
		burst_le_2_true = total(modulate_curve_le_psf_modify(RA=262.991, DEC=-33.834, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=2, DL=DL_LE_2, euler=euler_LE_2, par=par_LE_2))	
		print, burst_le_0_true/burst_le_1_true, 1, burst_le_2_true/burst_le_1_true
		
		;On-axis
		burst_le_0_true = total(modulate_curve_le_psf_modify(RA=263.353, DEC=-33.389, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=0, DL=DL_LE_0, euler=euler_LE_0, par=par_LE_0))
		burst_le_1_true = total(modulate_curve_le_psf_modify(RA=263.353, DEC=-33.389, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=1, DL=DL_LE_1, euler=euler_LE_1, par=par_LE_1))
		burst_le_2_true = total(modulate_curve_le_psf_modify(RA=263.353, DEC=-33.389, Flux=1, Qt1=Q1_le, Qt2=Q2_le, Qt3=Q3_le, BoxID=2, DL=DL_LE_2, euler=euler_LE_2, par=par_LE_2))	
		print, burst_le_0_true/burst_le_1_true, 1, burst_le_2_true/burst_le_1_true
		
		
		;------------------------------------------------------------------------
    	; window,3,xsize=720,ysize=760,xpos=720,ypos=280
    	; ;------------------------------------
    	; !P.MULTI = [0, 3, 4, 0, 1]
    	; !p.FONT=-1
		; ;------------------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot,findgen(k)+1, RA_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.80, 0.36, 0.96]
		; ;------------------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot,findgen(k)+1, Dec_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.80, 0.66, 0.96]
		; ;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, R_2_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.72, 0.80, 0.96, 0.96]	
		; ;------------------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot,findgen(k)+1, flux_le_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.60, 0.36, 0.76]	
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, flux_1_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.60, 0.66, 0.76]
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, flux_2_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.72, 0.60, 0.96, 0.76]
		; ;------------------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot,findgen(k)+1, flux_me_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.40, 0.36, 0.56]
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, flux_4_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.40, 0.66, 0.56]
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, flux_5_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.72, 0.40, 0.96, 0.56]
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, flux_6_fit, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.20, 0.36, 0.36]
		; ;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, Par_3_fit, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.20, 0.66, 0.36]
		; ;------------------------------------
		; multiplot,/doxaxis,/doyaxis
		; plot,findgen(k)+1, chi2, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.00, 0.36, 0.16]
		; ;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, Par_3_fit, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.00, 0.66, 0.16]
		; ;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, chi2, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.72, 0.00, 0.96, 0.16]
		; ;------------------------------------
		; multiplot, /reset
		; multiplot,[1,1],/init,/verbose
		; ;------------------------------------------------------------------------
		
		; ;;------------------------------------------------------------------------
    	; ;window,4,xsize=720,ysize=760,xpos=400,ypos=280
    	; ;;------------------------------------
    	; ;!P.MULTI = [0, 3, 4, 0, 1]
    	; ;!p.FONT=-1
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, RA_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.80, 0.36, 0.96]
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, Dec_fit,psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.42, 0.80, 0.66, 0.96]
		; ;;------------------------------------
		; ;for i=0, n_time_me*0+9-1 do begin
		; ;	;------------------------------------
		; ;	multiplot,/doxaxis,/doyaxis
		; ;	plot,findgen(k)+1, flux_me_fit[*,i],psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12+0.3*(i mod 3), 0.60-0.2*(i/3), 0.36+0.3*(i mod 3), 0.76-0.2*(i/3)]
		; ;endfor
		; ;;------------------------------------
		; ;multiplot,/doxaxis,/doyaxis
		; ;plot,findgen(k)+1, chi2, psym=3,symsize=0.6,xstyle=1,ystyle=1,position=[0.12, 0.00, 0.36, 0.16]
		; ;;------------------------------------
		; ;multiplot, /reset
		; ;multiplot,[1,1],/init,/verbose
		; ;;------------------------------------------------------------------------
		
	endif
	
	
end





