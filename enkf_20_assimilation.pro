PRO enkf_20_assimilation

  @set_time.inc
  COMMON SHARE_OBS_CH, nobs_ch, obsvalues_ch, obserrors_ch,   $
                       H_mat_ch, obs_i_ch, obs_j_ch
  COMMON SHARE_OBS_KR, nobs_kr, obsvalues_kr, obserrors_kr,   $
                       H_mat_kr, obs_i_kr, obs_j_kr

;; nx,ny,nz,nv,nobs were defined in enkf_cmaq_gen_obs.pro
nens = 40L
dx = 27;km
nx =144L & ny =105L & nz = 14L & nv = 1L
vec_len = nx*ny*nz*nv

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;comp = 'SO2'                                                                 ;;
;ecase = '08.enkf_pre7'                                                       ;;
comp = 'PM25'                                                                 ;;
ecase = '50.enkf_vs_3dvar'                                                    ;;
                                                                              ;;
itime = STRING(126               ,FORMAT='(I4.4)')                            ;;
kfint = 6                                                                     ;;
;ftime = '150' ;; for the 1st KF, otherwise block here                        ;;
ftime = STRING(    ftime         ,FORMAT='(I4.4)') ;; time for the xf         ;;
btime = STRING(FIX(ftime)-kfint  ,FORMAT='(I4.4)') ;; time for the xb         ;;
ptime = STRING(FIX(ftime)-kfint*2,FORMAT='(I4.4)') ;; time for RTPS           ;;
atime = ftime  ;; time for the xa, which is same as that of the forecast xf   ;;
rfile = 1 ; 1: read file at first, 0: don't read again                        ;;
wfile = 1 ;; 1: write or 0: not                                               ;;
nspc  = 63;; 1 for a quick check and 43 for all species                       ;;
                                                                              ;; 
nf = 6 ;; forecast periods (hour)                                             ;;
                                                                              ;;
sigma_obs = 0.05; Relative uncertainty in observation data                    ;;
obs_err   = 2   ; 1: observation error is a constant based on sigma_obs       ;;
                ; 2: using the calculated error by Elbern et al.(2007)        ;;
inf_pri   = 4   ; 0: turn off the inflation                                   ;;
                ; 1: additive inflation                                       ;;
                ; 2: multiplicative inflation (before filtering)              ;;
                ; 3: RTPS multiplicative inflation (before filtering)         ;;
                ; 4: RTPS multiplicative inflation based on NMC results       ;;
                ;   (RTPS, relaxation-to-prior-spread)                        ;;
inf_pos   = 0   ; 0: turn off the inflation                                   ;;
                ; 1: additive inflation                                       ;;
                ; 2: multiplicative inflation (after  filtering)              ;;
                ; 3: RTPS multiplicative inflation (after  filtering)         ;;
gam_2_pri= 4.0  ; inflation magnitude for prior     multiplicative inflation  ;;
gam_2_pos= 4.0  ; inflation magnitude for posterior multiplicative inflation  ;;
gam_3_pri= 1.0  ; inflation magnitude for prior     RTPS                      ;;
gam_3_pos= 1.0  ; inflation magnitude for posterior RTPS                      ;;
gam_4_pri= 1.0  ; inflation magnitude for prior     RTPS based on NMC         ;;
localization = 1; Gaspri and Cohn piecewise polynomial                        ;;
hwidth = 200 ; 025 =  2.5dx ( 25km)  2.5/49                                   ;;
             ; 050 =  5.0dx ( 50km)  5.0/49                                   ;;
             ; 100 = 10.0dx (100km) 10.0/49                                   ;;
             ; 200 = 20.0dx (200km) 20.0/49                                   ;;
vwidth = 2   ; km (about 7.5 th layer)                                        ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dir     ='/home/spark/IDLWorkspace82/Default_'+ecase+'/'
pwd_cctm='/home/spark/cmaqv5.2.1/DATA/'+ecase+'/cctm/'
pwd_icon='/home/spark/cmaqv5.2.1/DATA/'+ecase+'/icon/'

;; read ensemble data at the background time especially for the RTPS inflation
IF (inf_pri EQ 3) THEN BEGIN
  xi_b     = DBLARR(nx,ny,nz,nv,nens)
  FOR e=0,nens-1 DO BEGIN
    IF (ftime EQ itime) THEN BEGIN
    finpi1= pwd_icon+btime+'-'+btime+'_ie/' $
          +'ICON_'+btime+'-'+btime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
    finpi2= pwd_icon+btime+'-'+btime+'_ie/' $
          +'CCTM_PMDIAG_'+btime+'-'+btime+'.nc'
    ENDIF ELSE BEGIN
    finpi1= pwd_cctm+btime+'-'+btime+'_an/' $
          +'CCTM_CGRID_'+btime+'-'+btime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
    finpi2= pwd_cctm+btime+'-'+btime+'_an/' $
          +'CCTM_PMDIAG_'+btime+'-'+btime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
    ENDELSE
    pm25i_b = GET_PM25_COMBINE(finpi1,finpi2)
    xi_b[*,*,*,0,e] = pm25i_b[*,*,*]
    PRINT, finpi1
  ENDFOR
  xi_b_vec = REFORM(xi_b,vec_len,nens)
  x_b_vec= MEAN( xi_b_vec,DIM=2)
  x_b = REFORM(x_b_vec,nx,ny,nz,nv)
ENDIF 
PRINT

IF (rfile EQ 1) THEN BEGIN
  ;; read ensemble data at the forecast time
  xi_f      = DBLARR(nx,ny,nz,nv,nens)
  FOR e=0,nens-1 DO BEGIN
    finp1= pwd_cctm+btime+'-'+ftime+'_pr/CCTM_CGRID_'+btime+'-'+ftime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
    finp2= pwd_cctm+ftime+'-'+ftime+'_an/CCTM_PMDIAG_'+ftime+'-'+ftime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
    PRINT, finp1
    pm25i_f = GET_PM25_COMBINE(finp1,finp2)
    xi_f[*,*,*,0,e] = pm25i_f[*,*,*]
  ENDFOR
  xi_f_vec  = REFORM(xi_f,vec_len,nens)
  x_f_vec= MEAN( xi_f_vec,DIM=2)
  x_f = REFORM(x_f_vec,nx,ny,nz,nv)
  PRINT  
ENDIF

gam_sav2 = DBLARR(vec_len)

CASE inf_pri   OF
  2: BEGIN
    PRINT, 'MULTIPLICATIVE INFLATION (before filtering)' ; Constantinescu et al.(2007b)
    xi_nf_vec = DBLARR(vec_len,nens)
    FOR vec=0,vec_len-1 DO $
      xi_nf_vec[vec,*] = x_f_vec[vec] + gam_2_pri * (xi_f_vec[vec,*] - x_f_vec[vec] )
    xi_nf = REFORM(xi_nf_vec,nx,ny,nz,nv,nens)
    x_nf_vec = MEAN( xi_nf_vec,DIM=2)
    x_nf = REFORM(x_nf_vec,nx,ny,nz,nv)
  END  
  3: BEGIN    
    PRINT, 'RELAXATION-TO-PRIOR-SPREAD (before filtering)' ; Whitaker and Hamil (2012)
    sgm_b = DBLARR(vec_len)
    sgm_f = DBLARR(vec_len)
    gam_sav_pri = DBLARR(vec_len)
    xi_nf_vec = DBLARR(vec_len,nens)
    FOR vec=0,vec_len-1 DO BEGIN
      sgm_b[vec] = STDDEV(xi_b_vec[vec,*])
      sgm_f[vec] = STDDEV(xi_f_vec[vec,*])
;     IF (sgm_f[vec] GT sgm_b[vec]) THEN BEGIN
;       PRINT, vec, sgm_f[vec], sgm_b[vec], ARRAY_INDICES(xi_f,vec)
;     ENDIF
      gam_sav_pri[vec] = gam_3_pri * (sgm_b[vec] - sgm_f[vec]) / sgm_f[vec] + 1
      IF (sgm_b[vec] LE sgm_f[vec]) THEN gam_sav_pri[vec] = 1.0
    ENDFOR
    loc_inf = WHERE( FINITE(gam_sav_pri,/INFINITY) )
    loc_nan = WHERE( FINITE(gam_sav_pri,/NAN     ) )    
    gam_sav_pri[ WHERE( FINITE(gam_sav_pri,/INFINITY) ) ] = 1.0D
    gam_sav_pri[ WHERE( FINITE(gam_sav_pri,/NAN     ) ) ] = 1.0D    
    FOR vec=0,vec_len-1 DO $
      xi_nf_vec[vec,*] = x_f_vec[vec] + gam_sav_pri[vec] * (xi_f_vec[vec,*] - x_f_vec[vec] )
    xi_nf = REFORM(xi_nf_vec,nx,ny,nz,nv,nens)
    x_nf_vec = MEAN( xi_nf_vec,DIM=2)
    x_nf = REFORM(x_nf_vec,nx,ny,nz,nv)
  END
  4: BEGIN
    PRINT, 'RELAXATION-TO-PRIOR-SPREAD BASED ON NMC (before filtering)' ; Whitaker and Hamil (2012)
    xi_b_vec = DBLARR(vec_len,nens)  ;; dummy variable
    finp= pwd_cctm+'CMAQ_BERROR_NMC_PM25.NC'
    idinp=NCDF_OPEN(finp,/NOWRITE  )
    PRINT, finp
    NCDF_VARGET,idinp,'varce_pm2_5', var
    sgm_b_3d = SQRT(var)
    sgm_b = REFORM(sgm_b_3d,vec_len)
    sgm_f = DBLARR(vec_len)
    gam_sav_pri = DBLARR(vec_len)
    xi_nf_vec = DBLARR(vec_len,nens)
    FOR vec=0,vec_len-1 DO BEGIN
      sgm_f[vec] = STDDEV(xi_f_vec[vec,*])
;     IF (sgm_f[vec] GT sgm_b[vec]) THEN BEGIN
;       PRINT, vec, sgm_f[vec], sgm_b[vec], ARRAY_INDICES(xi_f,vec)
;     ENDIF
      gam_sav_pri[vec] = gam_4_pri * (sgm_b[vec] - sgm_f[vec]) / sgm_f[vec] + 1
      IF (sgm_b[vec] LE sgm_f[vec]) THEN gam_sav_pri[vec] = 1.0
    ENDFOR
    loc_inf = WHERE( FINITE(gam_sav_pri,/INFINITY) )
    loc_nan = WHERE( FINITE(gam_sav_pri,/NAN     ) )
    gam_sav_pri[ WHERE( FINITE(gam_sav_pri,/INFINITY) ) ] = 1.0D
    gam_sav_pri[ WHERE( FINITE(gam_sav_pri,/NAN     ) ) ] = 1.0D
    FOR vec=0,vec_len-1 DO $
      xi_nf_vec[vec,*] = x_f_vec[vec] + gam_sav_pri[vec] * (xi_f_vec[vec,*] - x_f_vec[vec])
    xi_nf = REFORM(xi_nf_vec,nx,ny,nz,nv,nens)
    x_nf_vec = MEAN( xi_nf_vec,DIM=2)
    x_nf = REFORM(x_nf_vec,nx,ny,nz,nv)
  END
  ELSE: BEGIN
    PRINT, 'NO INFLATION (before filtering)'
    xi_b_vec = DBLARR(vec_len,nens)  ;; dummy variable
    xi_nf_vec = xi_f_vec
    x_nf_vec = x_f_vec
    xi_nf = xi_f
    x_nf = x_f
  END  
ENDCASE

;; combine the obs. related variables in China and in Korea

nobs = nobs_ch + nobs_kr

CASE 1 OF

  (nobs LE 1): BEGIN
    PRINT, 'NO OBSERVATION DATA IS AVAILABLE AT THIS TIME!'
    PRINT, 'nobs= ',nobs
    PRINT, 'KF STEP WILL BE SKIPPED.'
    GOTO, NO_OBS
  END

  (nobs_ch EQ 0) AND (nobs_kr GE 2): BEGIN
    obsvalues = obsvalues_kr
    obserrors = obserrors_kr
    obs_i     = obs_i_kr
    obs_j     = obs_j_kr
    H_mat     = H_mat_kr
  END

  (nobs_kr EQ 0) AND (nobs_ch GE 2): BEGIN
    obsvalues = obsvalues_ch
    obserrors = obserrors_ch
    obs_i     = obs_i_ch
    obs_j     = obs_j_ch
    H_mat     = H_mat_ch
  END

  ELSE: BEGIN
    obsvalues = [obsvalues_ch,obsvalues_kr]
    obserrors = [obserrors_ch,obserrors_kr]
    obs_i     = [obs_i_ch    ,obs_i_kr    ]
    obs_j     = [obs_j_ch    ,obs_j_kr    ]
    H_mat     = [[H_mat_ch]  ,[H_mat_kr]  ]
  END

ENDCASE

;; FOR CACULATING SIGMA_b and SIGMA_f
hxi_b_vec = TRANSPOSE(H_mat ## TRANSPOSE(xi_b_vec))  ; hxi_b_vec[nobs,nens]
hxi_f_vec = TRANSPOSE(H_mat ## TRANSPOSE(xi_f_vec))  ; hxi_f_vec[nobs,nens]
hxi_nf_vec= TRANSPOSE(H_mat ## TRANSPOSE(xi_nf_vec))  ; hxi_nf_vec[nobs,nens]
 hx_nf_vec= MEAN(hxi_nf_vec,DIM=2)

PRINT, 'CALCULATE THE MODEL ERROR COVARIANCE'  
PH_sum  = DBLARR(vec_len,nobs)
HPH_sum = DBLARR(nobs,nobs)

FOR e=0,nens-1 DO BEGIN
  PH_sum  = PH_sum  + MATRIX_MULTIPLY(( xi_nf_vec[*,e]- x_nf_vec), $
                                      (hxi_nf_vec[*,e]-hx_nf_vec))
  HPH_sum = HPH_sum + MATRIX_MULTIPLY((hxi_nf_vec[*,e]-hx_nf_vec), $
                                      (hxi_nf_vec[*,e]-hx_nf_vec))
ENDFOR

PH = PH_sum / (nens-1)
HPH=HPH_sum / (nens-1)

IF (localization EQ 1) THEN BEGIN
PRINT, 'LOCALIZE THE MODEL ERROR COVARIANCE'
PRINT, 'LOCALIZATION SCALE IS ',hwidth,' KM'
  ; Distance and decorrelation function
  decorr_x2o = DBLARR(nx,ny,nz,nv,nobs)
  FOR o=0,nobs-1 DO BEGIN
    FOR j=0,ny-1 DO BEGIN
      FOR i=0,nx-1 DO BEGIN
        distance = SQRT( (obs_i[o]-i)^2.0 + (obs_j[o]-j)^2.0 ) * dx
        zoc = distance / hwidth
        IF ((distance GE 0) AND (distance LE hwidth)) THEN BEGIN
          decorr_x2o[i,j,*,0,o]=-(1/4.0)*zoc^5.0 + (1/2.0)*zoc^4.0 + (5/8.0)*zoc^3.0 $
                                -(5/3.0)*zoc^2.0 + 1.0
        ENDIF ELSE IF ((distance GE hwidth) AND (distance LE hwidth*2)) THEN BEGIN
          decorr_x2o[i,j,*,0,o]= (1/12.0)*zoc^5.0 - (1/2.0)*zoc^4.0 + (5/8.0)*zoc^3.0 $
                                +(5/3.0)*zoc^2.0 - (5/1.0)*zoc + 4.0 - (2/3.0)*(1/zoc)        
        ENDIF ELSE IF (distance GE hwidth*2) THEN BEGIN
          decorr_x2o[i,j,*,0,o]=0.0
        ENDIF
      ENDFOR    
    ENDFOR
  ENDFOR
  
  decorr_o2o = DBLARR(nobs,nobs)
  FOR o2=0,nobs-1 DO BEGIN
    FOR o1=0,nobs-1 DO BEGIN  
      distance = SQRT( (obs_i[o1]-obs_i[o2])^2.0 + (obs_j[o1]-obs_j[o2])^2.0 ) * dx
      zoc = distance / hwidth
      IF ((distance GE 0) AND (distance LE hwidth)) THEN BEGIN
        decorr_o2o[o1,o2]=-(1/4.0)*zoc^5.0 + (1/2.0)*zoc^4.0 + (5/8.0)*zoc^3.0 $
                          -(5/3.0)*zoc^2.0 + 1.0
      ENDIF ELSE IF ((distance GE hwidth) AND (distance LE hwidth*2)) THEN BEGIN
        decorr_o2o[o1,o2]= (1/12.0)*zoc^5.0 - (1/2.0)*zoc^4.0 + (5/8.0)*zoc^3.0 $
                          +(5/ 3.0)*zoc^2.0 - (5/1.0)*zoc + 4.0 - (2/3.0)*(1/zoc)
      ENDIF ELSE IF (distance GE hwidth*2) THEN BEGIN
        decorr_o2o[o1,o2]=0.0
      ENDIF
    ENDFOR
  ENDFOR

  zh = [   16.9,   55.8,  120.5,  244.5,  413.2 $
       ,  628.7,  983.5, 1492.0, 2076.2, 2857.4 $
       , 3945.7, 5585.5, 8490.8,13934.2]        ; elevation (m)

  FOR k=0,nz-1 DO BEGIN
    decorr_x2o[*,*,k,0,*] = decorr_x2o[*,*,k,0,*] $ 
                           * EXP(-((zh[k]-zh[0])/(vwidth*1000))^2.0)
  ENDFOR  
  
  PH  = REFORM(decorr_x2o,vec_len,nobs) * PH
  HPH = decorr_o2o * HPH 
ENDIF

; Observations at the end of the time interval
PRINT, 'CREATE OBSERVATION ERROR COVARIANCE AND PERTURBED OBSERATION (yi)'
y = obsvalues

;; fixed random number for perturbed y (yi)
seed = 2L

CASE obs_err   OF
  1: BEGIN
    ; Observation error and covariance
    obse = sigma_obs*y
    R = DIAG_MATRIX( (obse)^2.0 )
    ; Perturbed ensemble observation data set
    yi = REBIN(y,nobs,nens)
    std_obs = sigma_obs * yi
    FOR o=0,nobs-1 DO BEGIN
      yi[o,*] = yi[o,*] + RANDOMN(seed,nens)*std_obs[o,*]
    ENDFOR
  END
  2: BEGIN
    ; Observation error and covariance
    obse = obserrors
    R = DIAG_MATRIX( obse^2.0 )
    ; Perturbed ensemble observation data set
    yi = REBIN(y,nobs,nens)
    FOR o=0,nobs-1 DO BEGIN
      yi[o,*] = yi[o,*] + RANDOMN(seed,nens)*obse[o]
    ENDFOR  
  END
  ELSE: BEGIN
    PRINT, 'The obs_err should be 1 or 2.'
    EXIT
  END
ENDCASE

PRINT, 'CALCULATE KALMAN GAIN MATRIX'
; Kalman gain matrix
K = TRANSPOSE(PH) ## INVERT(HPH+R)
 
xi_a_vec = DBLARR(vec_len,nens)
FOR e=0,nens-1 DO $
  xi_a_vec[*,e] = xi_nf_vec[*,e] + K ## (yi[*,e] - (H_mat ## xi_nf_vec[*,e]))
xi_a = REFORM(xi_a_vec,nx,ny,nz,nv,nens)

x_a_vec = MEAN(xi_a_vec,DIM=2)     ; ensemble mean of analysis
x_a = REFORM(x_a_vec,nx,ny,nz,nv)

IF (nobs LE 1) THEN BEGIN
NO_OBS : PRINT, 'USE THE INFLATED CGRID DATA W/O FILTERING (XI_NF) FOR AN ANALYSIS.'
xi_a_vec = xi_nf_vec
 x_a_vec =  x_nf_vec
xi_a     = xi_nf
 x_a     =  x_nf
  inf_pos = 0  ; NO INFLATION IS NEEDED FOR NO OBS DATA, USE XI_NF.
ENDIF

CASE inf_pos   OF
  1: BEGIN
    PRINT, 'ADDITIVE INFLATION (after filterfing)' ; Constantinescu et al.(2007b)
    ;; Generate new background (nb) data by perturbing the analysis (a) mean 
    xi_nb_vec = REBIN(x_a_vec,vec_len,nens)
    std_mod = sigma_x0 * xi_nb_vec
    ;   PRINT, itime, std_mod[*,0]
    FOR vec=0,vec_len-1 DO $
      xi_nb_vec[vec,*] = xi_nb_vec[vec,*] + RANDOMN(seed,nens)*std_mod[vec,*]
    END
  2: BEGIN
    PRINT, 'MULTIPLICATIVE INFLATION (after filterfing)' ; Constantinescu et al.(2007b)
    xi_nb_vec = DBLARR(vec_len,nens)
    FOR vec=0,vec_len-1 DO $
      xi_nb_vec[vec,*] = x_a_vec[vec] + gam_2_pos * (xi_a_vec[vec,*] - x_a_vec[vec] )
    END
  3: BEGIN
    PRINT, 'RELAXATION-TO-PRIOR-SPREAD (after filtering)' ; Whitaker and Hamil (2012)
    sgm_nf = DBLARR(vec_len)
    sgm_a = DBLARR(vec_len)
    gam_sav_pos = DBLARR(vec_len)    
    xi_nb_vec = DBLARR(vec_len,nens)
    FOR vec=0,vec_len-1 DO BEGIN
      sgm_nf[vec] = STDDEV(xi_nf_vec[vec,*])
      sgm_a[vec] = STDDEV(xi_a_vec[vec,*])
      gam_sav_pos[vec] = gam_3_pos * (sgm_nf[vec] - sgm_a[vec]) / sgm_a[vec] + 1
    ENDFOR
    gam_sav_pos[ WHERE( FINITE(gam_sav_pos,/INFINITY) ) ] = 1.0D
    gam_sav_pos[ WHERE( FINITE(gam_sav_pos,/NAN     ) ) ] = 1.0D    
    FOR vec=0,vec_len-1 DO $
      xi_nb_vec[vec,*] = x_a_vec[vec] + gam_sav_pos[vec] * (xi_a_vec[vec,*] - x_a_vec[vec] )
    END
  ELSE: BEGIN    
    PRINT, 'NO INFLATION (after filterfing)'
    xi_nb_vec = xi_a_vec
    END
ENDCASE

  hxi_a_vec = TRANSPOSE(H_mat ## TRANSPOSE(xi_a_vec ))  ; hxi_a_vec[nobs,nens]
  hxi_nb_vec= TRANSPOSE(H_mat ## TRANSPOSE(xi_nb_vec))  ; hxi_nb_vec[nobs,nens]

;;; convert the negative con. to zero
;xi_nb_vec[WHERE(xi_nb_vec LE 0.0)] = 1.0D-30

xi_nb = REFORM(xi_nb_vec,nx,ny,nz,nv,nens)
x_nb = x_a ; basically x_a from the mean of xi_a is same as x_nb from the mean of xi_nb

;; plot the quick results at observation sites
IF (nobs GT 2) THEN BEGIN
PRINT,'PLOT QUCIK RESULTS OF THE ASIMILATION'

SET_PLOT, 'Z'
DEVICE, DECOMPOSED=0
DEVICE, SET_RESOLUTION=[2000,800], SET_PIXEL_DEPTH=24
loadct,39

;xi_f_plt = H_mat ## Transpose(xi_f_vec) ; xi_f_plt[nens,nobs]
; x_f_plt = H_mat ## Transpose( x_f_vec) ;  x_f_plt[     nobs]
;xi_a_plt = H_mat ## Transpose(xi_a_vec) ; xi_a_plt[nens,nobs]
; x_a_plt = H_mat ## Transpose( x_a_vec) ;  x_a_plt[     nobs]

;; mean of the propagated ensemble forecast and assimilated ensemble analysis
hx_b_vec= MEAN(hxi_b_vec,DIM=2) ;; [nobs,nens] -> [nobs]
hx_f_vec= MEAN(hxi_f_vec,DIM=2) ;; [nobs,nens] -> [nobs]
hx_a_vec= MEAN(hxi_a_vec,DIM=2) ;; [nobs,nens] -> [nobs]
hx_nb_vec=MEAN(hxi_nb_vec,DIM=2);; [nobs,nens] -> [nobs]

;; calculate the innovation vector and the analyais errors for ensemble mean
inno = y - hx_f_vec ;; [nobs]
anle = y - hx_a_vec ;; [nobs]

;; save the diagonal components in HPH
hph_sgm = dblarr(nobs)
for o=0,nobs-1 do hph_sgm[o]=SQRT(HPH[o,o])

for fig=0,5 do begin
; window, 1, xs=2000, ys=800
  plot, y, yrange=[0,120], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[0+(175*fig),174+(175*fig)]                                    $
      , ytitle='PM2.5',charsize=1.5,xtickinterval=40, yticklen=0.01, xminor=4
  if ((nobs_ch ge 0+(175*fig)) and (nobs_ch lt 174+(175*fig))) then begin
    oplot, intarr(10)+nobs_ch, indgen(10)*999, color=0, thick=2, linestyle=2
  endif

  for e=0,nens-1 do oplot,  hxi_nf_vec[*,e],color=0
  for e=0,nens-1 do oplot,  hxi_f_vec[*,e], color=50
  oplot,  hx_f_vec, color=100, thick=2

  for e=0,nens-1 do oplot,  yi[*,e],color=30
  oplot, y, color=150, thick=2.0 ;, xrange=[30,50]

  img = tvrd(/true)
  write_png, dir+'/assimilation1'+STRING(fig,FORMAT='(I1)')+'.png', img
endfor

for fig=0,5 do begin
; window, 1, xs=2000, ys=800
  plot, y, yrange=[0,120], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[0+(175*fig),174+(175*fig)]                                    $
      , ytitle='PM2.5',charsize=1.5,xtickinterval=40, yticklen=0.01, xminor=4
  if ((nobs_ch ge 0+(175*fig)) and (nobs_ch lt 174+(175*fig))) then begin
    oplot, intarr(10)+nobs_ch, indgen(10)*999, color=0, thick=2, linestyle=2
  endif
  
  for e=0,nens-1 do oplot,  hxi_nb_vec[*,e],color=0
  for e=0,nens-1 do oplot,  hxi_a_vec[*,e], color=200
  oplot,  hx_a_vec, color=240, thick=2
 
  for e=0,nens-1 do oplot,  yi[*,e],color=30
  oplot, y, color=150, thick=2.0 ;, xrange=[30,50]

  img = tvrd(/true)
  write_png, dir+'/assimilation2'+STRING(fig,FORMAT='(I1)')+'.png', img
endfor

for fig=0,5 do begin
; window, 1, xs=2000, ys=800
  plot, y, yrange=[-120,120], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[0+(175*fig),174+(175*fig)]                                    $
      , ytitle='PM2.5',charsize=1.5,xtickinterval=40, yticklen=0.01, xminor=4
  oplot, intarr(nobs),color=0  
  if ((nobs_ch ge 0+(175*fig)) and (nobs_ch lt 174+(175*fig))) then begin
    oplot, intarr(10)+nobs_ch, indgen(10)*900-450, color=0, thick=2, linestyle=2
  endif

  oplot,         y, color=150, thick=2 ;, xrange=[30,50]  oplot, 
  oplot,  hx_f_vec, color= 50, thick=2
  oplot,  hx_a_vec, color=200, thick=2
  oplot,      inno, color=100, thick=2
  oplot,      anle, color=240, thick=2

  img = tvrd(/true)
  write_png, dir+'/assimilation3'+STRING(fig,FORMAT='(I1)')+'.png', img
endfor

;; Calculate the RMSE (i.e. sigma) to measure the spread for each component
rmse_b = fltarr(nobs) 
rmse_f = fltarr(nobs) 
rmse_nf= fltarr(nobs) 
rmse_a = fltarr(nobs) 
rmse_nb= fltarr(nobs) 
rmse_y = fltarr(nobs) 
for o=0,nobs-1 do begin
  rmse_b[o] = sqrt(mean(( hxi_b_vec[o,*]- hx_b_vec[o])^2.0))
  rmse_f[o] = sqrt(mean(( hxi_f_vec[o,*]- hx_f_vec[o])^2.0))
  rmse_nf[o]= sqrt(mean((hxi_nf_vec[o,*]-hx_nf_vec[o])^2.0))
  rmse_a[o] = sqrt(mean(( hxi_a_vec[o,*]- hx_a_vec[o])^2.0))
  rmse_nb[o]= sqrt(mean((hxi_nb_vec[o,*]-hx_nb_vec[o])^2.0))
  rmse_y[o] = sqrt(mean((        yi[o,*]-        y[o])^2.0))
endfor

IF (inf_pri EQ 4) THEN BEGIN
  rmse_b =  TRANSPOSE(H_mat ## TRANSPOSE(sgm_b))
ENDIF

IF (nobs_kr GT 2) THEN BEGIN
  ;; normalized sigma only for Korean obs sites
  PRINT,'PLOT THE SPREAD INFORMATION ONLY AT KOREAN STATIONS'
  plot, y, yrange=[0,100], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[-10,nobs_kr+10], /nodata                                     $
      , ytitle='stddev (%)',charsize=1.5,xtickinterval=5, yticklen=0.01, xminor=4
  oplot, stddev( hxi_b_vec[nobs_ch:*,*],dim=2)/ hx_b_vec[nobs_ch:*]*100, color=0
  oplot, stddev( hxi_f_vec[nobs_ch:*,*],dim=2)/ hx_f_vec[nobs_ch:*]*100, color=50
  oplot, stddev(hxi_nf_vec[nobs_ch:*,*],dim=2)/hx_nf_vec[nobs_ch:*]*100, color=100
  oplot, stddev( hxi_a_vec[nobs_ch:*,*],dim=2)/ hx_a_vec[nobs_ch:*]*100, color=200
  oplot, stddev(hxi_nb_vec[nobs_ch:*,*],dim=2)/hx_nb_vec[nobs_ch:*]*100, color=240
  oplot, stddev(        yi[nobs_ch:*,*],dim=2)/        y[nobs_ch:*]*100, color=150
  img = tvrd(/true)
  write_png, dir+'/plot01_stddev_kr.png', img

  ;; rmse (i.e. sigma) only for Korean obs sites
  plot, y, yrange=[0,20], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[-10,nobs_kr+10], /nodata                                     $
      , ytitle='RMSE',charsize=1.5,xtickinterval=5, yticklen=0.01, xminor=4
  oplot,  rmse_b[nobs_ch:*] ,color=0
  oplot,  rmse_f[nobs_ch:*] ,color=50
  oplot, rmse_nf[nobs_ch:*] ,color=100
  oplot,  rmse_a[nobs_ch:*] ,color=200
  oplot, rmse_nb[nobs_ch:*] ,color=240
  oplot,  rmse_y[nobs_ch:*] ,color=150
  img = tvrd(/true)
  write_png, dir+'/plot02_rmse_kr.png', img
ENDIF

IF (nobs_ch GT 2) THEN BEGIN
  ;; normalized sigma only for Chinese obs sites
  PRINT,'PLOT THE SPREAD INFORMATION ONLY AT CHINESE STATIONS'
  plot, y, yrange=[0,100], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[-10,nobs_ch+10], /nodata                                     $
      , ytitle='stddev (%)',charsize=1.5,xtickinterval=20, yticklen=0.01, xminor=4
  oplot, stddev( hxi_b_vec[0:nobs_ch-1,*],dim=2)/ hx_b_vec[0:nobs_ch-1]*100, color=0
  oplot, stddev( hxi_f_vec[0:nobs_ch-1,*],dim=2)/ hx_f_vec[0:nobs_ch-1]*100, color=50
  oplot, stddev(hxi_nf_vec[0:nobs_ch-1,*],dim=2)/hx_nf_vec[0:nobs_ch-1]*100, color=100
  oplot, stddev( hxi_a_vec[0:nobs_ch-1,*],dim=2)/ hx_a_vec[0:nobs_ch-1]*100, color=200
  oplot, stddev(hxi_nb_vec[0:nobs_ch-1,*],dim=2)/hx_nb_vec[0:nobs_ch-1]*100, color=240
  oplot, stddev(        yi[0:nobs_ch-1,*],dim=2)/        y[0:nobs_ch-1]*100, color=150
  img = tvrd(/true)
  write_png, dir+'/plot01_stddev_ch.png', img
  ;; rmse (i.e. sigma) only for Chinese obs sites
  plot, y, yrange=[0,20], thick=2.0, color=0, back=255, xtitle='obs number'  $
      , xrange=[-10,nobs_ch+10], /nodata                                     $
      , ytitle='RMSE',charsize=1.5,xtickinterval=20, yticklen=0.01, xminor=4
  oplot,  rmse_b[0:nobs_ch-1], color=0
  oplot,  rmse_f[0:nobs_ch-1], color=50
  oplot, rmse_nf[0:nobs_ch-1],color=100
  oplot,  rmse_a[0:nobs_ch-1], color=200
  oplot, rmse_nb[0:nobs_ch-1],color=240
  oplot,  rmse_y[0:nobs_ch-1], color=150
  img = tvrd(/true)
  write_png, dir+'/plot02_rmse_ch.png', img
ENDIF

DEVICE, /CLOSE
SET_PLOT, 'X'
ENDIF

; color=0,   black,  inflated one (xi_nf and xi_nb)
; color=50,  blue,   xi_f (fcs)
; color=100, cyan,   x_f (emean fcs)
; color=150, green,  y (obs)
; color=200, yellow, xi_a (ana)
; color=240, red,    x_a (emean ana)

;; Calculate the ratio of analised pm25 to forecasted pm25, which will be applied to aerosol species
ratio_i = xi_nb / xi_f

;ratio_cm  = MEAN(ratio_i,DIM=5)
;ratio_mc = x_nb / x_f

scld_var_si = DBLARR(nx,ny,nz,nv,nspc,nens)
;scld_var_s  = DBLARR(nx,ny,nz,nv,nspc     )

; TOTAL AITKEN MODE
  spcs = ['ASO4I'  ,'ANO3I'  ,'ANH4I'  ,                    $ ; SNA
          'ALVPO1I','ASVPO1I','ASVPO2I',                    $ ; POM
          'ALVOO1I','ALVOO2I','ASVOO1I','ASVOO2I',          $ ; SOM
          'ANAI'   ,'ACLI'   ,                              $ ; SS
          'AECI'   ,                                        $ ; EC
          'AOTHRI' ,                                        $ ; OTHER
; TOTAL ACCUMULATION MODE
          'ASO4J'  ,'ANO3J'  ,'ANH4J'  ,                    $ ; SNA
          'ALVPO1J','ASVPO1J','ASVPO2J','ASVPO3J','AIVPO1J',$ ; POM
          'ALVOO1J','ALVOO2J','ASVOO1J','ASVOO2J',          $ ; SOM
          'ASVOO3J','APCSOJ' ,                              $
          'AXYL1J' ,'AXYL2J' ,'AXYL3J' ,                    $
          'ATOL1J' ,'ATOL2J' ,'ATOL3J' ,                    $
          'ABNZ1J' ,'ABNZ2J' ,'ABNZ3J' ,                    $
          'AISO1J' ,'AISO2J' ,'AISO3J' ,                    $
          'ATRP1J' ,'ATRP2J' ,'ASQTJ'  ,                    $
          'AALK1J' ,'AALK2J' ,                              $
          'APAH1J' ,'APAH2J' ,'APAH3J' ,                    $
          'AORGCJ' ,'AOLGBJ' ,'AOLGAJ' ,'ANAJ',             $
          'ACLJ'   ,                                        $ ; SS
          'AECJ'   ,                                        $ ; EC
          'AOTHRJ' ,                                        $ ; OTHER
          'AFEJ'   ,'ASIJ'   ,                              $ ; SOIL
          'ATIJ'   ,'ACAJ'   ,  'AMGJ' , 'AMNJ'  , 'AALJ'   ,'AKJ']

IF (wfile EQ 1) THEN BEGIN
PRINT,'WRITING THE SCALED AEROSOL VARIABLES IN CGRID FILES'
finp=STRARR(nens) & idinp=LONARR(nens)
fout=STRARR(nens) & idout=LONARR(nens)
PRINT
FOR e=0,nens-1 DO BEGIN
  finp[e]=pwd_cctm+btime+'-'+ftime+'_pr/CCTM_CGRID_'+btime+'-'+ftime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
  idinp[e]=NCDF_OPEN(finp[e],/NOWRITE  )
  fout[e]=pwd_cctm+atime+'-'+atime+'_an/CCTM_CGRID_'+atime+'-'+atime+'_e'+STRING(e+1,FORMAT='(I3.3)')+'.nc'
  idout[e]=NCDF_OPEN(fout[e],/WRITE  )
  PRINT, fout[e]
  ;; Save the newly analysed and inflated initial ensemble
  FOR s=0,nspc-1 DO BEGIN
    NCDF_VARGET,idinp[e], spcs[s], var
    scld_var_si[*,*,*,0,s,e] = var * ratio_i[*,*,*,0,e]
  ENDFOR

  scld_var_si[WHERE(scld_var_si LE 0.0D)] = 1.0D-05

  FOR s=0,nspc-1 DO BEGIN
;   PRINT, 'e=',e+1,', s=',s,' ',spcs[s]
    NCDF_VARPUT,idout[e], spcs[s], scld_var_si[*,*,*,0,s,e]
  ENDFOR
ENDFOR
FOR e=0,nens-1 DO NCDF_CLOSE,idinp[e]
FOR e=0,nens-1 DO NCDF_CLOSE,idout[e]

;; Save the mean of newly analysed initial esnsemble
;fout_emean=pwd_cctm+atime+'-'+atime+'_an/CCTM_CGRID_'+atime+'-'+atime+'_e000.nc'
;idout_emean=NCDF_OPEN(fout_emean,/WRITE  )
;scld_var_s = MEAN(scld_var_si,DIM=6)
;PRINT, fout_emean
;FOR s=0,nspc-1 DO BEGIN
;  PRINT, 'e=',0,', s=',s,' ',spcs[s]
;  NCDF_VARPUT,idout_emean, spcs[s], scld_var_s[*,*,*,0,s]
;ENDFOR
;NCDF_CLOSE,idout_emean

ENDIF

;; mapping the ensemble mean variables on ref_* variables to be simple
ref_b =hx_b_vec
ref_f =hx_f_vec
ref_nf=hx_nf_vec
ref_a =hx_a_vec
ref_nb=hx_nb_vec
ref_y =MEAN(yi,DIM=2)

;; Save the informations of ensemble spread and mean at all the stations
PRINT
PRINT,'SAVE THE SPREAD INFORMATIONS, ENS MEAN VALUES, AND OMB & OMA AT ALL STATIONS'
OPENW ,20, dir+'rmse_emean_omb_oma_all.txt'
PRINTF,20,'ATIME' ,'#OBS','CH_KR'                                                  $
         ,'RMSE_B','RMSE_F','RMSE_NF','HPH_sgm','RMSE_A','RMSE_NB','RMSE_Y','OBSE' $
         , 'REF_B', 'REF_F', 'REF_NF',           'REF_A', 'REF_NB', 'REF_Y',  'Y'  $
         ,   'OMB',   'OMA',FORMAT='(3(A5,X),17A9)'
for o=0,nobs-1 do begin
  if (o LT nobs_ch) then begin
    chkr = 'CHN'
  endif else begin
    chkr = 'KOR'
  endelse
PRINTF,20, atime   ,o+1      ,chkr                                                         $
         ,rmse_b[o],rmse_f[o],rmse_nf[o],hph_sgm[o],rmse_a[o],rmse_nb[o],rmse_y[o],obse[o] $
         , ref_b[o], ref_f[o], ref_nf[o],            ref_a[o], ref_nb[o], ref_y[o],   y[o] $
         ,  inno[o],  anle[o],FORMAT='(A5,X,I5,X,A5,X,17F9.3)'
endfor
CLOSE, 20

;; Save the informations of ensemble spread and mean averaged in Chinese stations
OPENW ,30, dir+'rmse_emean_omb_oma_ave_ch.dat', /APPEND
IF (ftime EQ itime) THEN BEGIN
  PRINTF,30,'ATIME'                                                                  $
           ,'RMSE_B','RMSE_F','RMSE_NF','HPH_sgm','RMSE_A','RMSE_NB','RMSE_Y','OBSE' $
           , 'REF_B', 'REF_F', 'REF_NF',           'REF_A', 'REF_NB', 'REF_Y',  'Y'  $
           ,   'OMB',   'OMA',FORMAT='(A5,X,17A9)'
ENDIF
IF (nobs_ch GT 2) THEN BEGIN
  PRINTF,30,atime                                                                            $
           ,MEAN( rmse_b[0:nobs_ch-1]),MEAN( rmse_f[0:nobs_ch-1]),MEAN(rmse_nf[0:nobs_ch-1]) $
                                                                 ,MEAN(hph_sgm[0:nobs_ch-1]) $
           ,MEAN( rmse_a[0:nobs_ch-1]),MEAN(rmse_nb[0:nobs_ch-1]),MEAN( rmse_y[0:nobs_ch-1]) $
                                                                 ,MEAN(   obse[0:nobs_ch-1]) $
           ,MEAN(  ref_b[0:nobs_ch-1]),MEAN(  ref_f[0:nobs_ch-1]),MEAN( ref_nf[0:nobs_ch-1]) $
           ,MEAN(  ref_a[0:nobs_ch-1]),MEAN( ref_nb[0:nobs_ch-1]),MEAN(  ref_y[0:nobs_ch-1]) $
                                                                 ,MEAN(      y[0:nobs_ch-1]) $
           ,SQRT(MEAN(inno[0:nobs_ch-1]^2.0)),SQRT(MEAN(anle[0:nobs_ch-1]^2.0)),FORMAT='(A5,X,17F9.3)'
  CLOSE, 30
ENDIF ELSE BEGIN
  PRINTF,30,atime, -999.0+FLTARR(17), FORMAT='(A5,X,17F9.3)'
  CLOSE, 30
ENDELSE

;; Save the informations of ensemble spread and mean averaged in Korean stations
OPENW ,40, dir+'rmse_emean_omb_oma_ave_kr.dat', /APPEND
IF (ftime EQ itime) THEN BEGIN
  PRINTF,40,'ATIME'                                                                  $
           ,'RMSE_B','RMSE_F','RMSE_NF','HPH_sgm','RMSE_A','RMSE_NB','RMSE_Y','OBSE' $
           , 'REF_B', 'REF_F', 'REF_NF',           'REF_A', 'REF_NB', 'REF_Y',  'Y'  $
           ,   'OMB',   'OMA',FORMAT='(A5,X,17A9)'
ENDIF
IF (nobs_kr GT 2) THEN BEGIN
  PRINTF,40,atime                                                                      $
           ,MEAN( rmse_b[nobs_ch:*]),MEAN( rmse_f[nobs_ch:*]),MEAN(rmse_nf[nobs_ch:*]) $
                                                             ,MEAN(hph_sgm[nobs_ch:*]) $
           ,MEAN( rmse_a[nobs_ch:*]),MEAN(rmse_nb[nobs_ch:*]),MEAN( rmse_y[nobs_ch:*]) $
                                                             ,MEAN(   obse[nobs_ch:*]) $
           ,MEAN(  ref_b[nobs_ch:*]),MEAN(  ref_f[nobs_ch:*]),MEAN( ref_nf[nobs_ch:*]) $
           ,MEAN(  ref_a[nobs_ch:*]),MEAN( ref_nb[nobs_ch:*]),MEAN(  ref_y[nobs_ch:*]) $
                                                             ,MEAN(      y[nobs_ch:*]) $
           ,SQRT(MEAN(inno[nobs_ch:*]^2.0)),SQRT(MEAN(anle[nobs_ch:*]^2.0)),FORMAT='(A5,X,17F9.3)'
  CLOSE, 40
ENDIF ELSE BEGIN
  PRINTF,40,atime, -999.0+FLTARR(17), FORMAT='(A5,X,17F9.3)'
  CLOSE, 40
ENDELSE

;; Save the informations of ensemble spread and mean averaged in all the stations
OPENW ,50, dir+'rmse_emean_omb_oma_ave.dat', /APPEND
IF (ftime EQ itime) THEN BEGIN
  PRINTF,50,'ATIME'                                                                  $
           ,'RMSE_B','RMSE_F','RMSE_NF','HPH_sgm','RMSE_A','RMSE_NB','RMSE_Y','OBSE' $
           , 'REF_B', 'REF_F', 'REF_NF',           'REF_A', 'REF_NB', 'REF_Y',  'Y'  $
           ,   'OMB',   'OMA',FORMAT='(A5,X,17A9)'
ENDIF
IF (nobs GT 2) THEN BEGIN
  PRINTF,50,atime                                     $
           ,MEAN( rmse_b),MEAN( rmse_f),MEAN(rmse_nf) $
                                       ,MEAN(hph_sgm) $
           ,MEAN( rmse_a),MEAN(rmse_nb),MEAN( rmse_y) $
                                       ,MEAN(   obse) $
           ,MEAN(  ref_b),MEAN(  ref_f),MEAN( ref_nf) $
           ,MEAN(  ref_a),MEAN( ref_nb),MEAN(  ref_y) $
                                       ,MEAN(      y) $
           ,SQRT(MEAN(inno^2.0)),SQRT(MEAN(anle^2.0)),FORMAT='(A5,X,17F9.3)'
  CLOSE, 50
ENDIF ELSE BEGIN
  PRINTF,50,atime, -999.0+FLTARR(17), FORMAT='(A5,X,17F9.3)'
  CLOSE, 50  
ENDELSE

;; Save PH for later analysis
PRINT
PRINT,'SAVE PH MATRIX AS A BINARY FORMAT FOR LATER ANALYSIS'
OPENW,  60, dir+'PH.bin'
WRITEU, 60, PH
CLOSE,  60

;; Create a flag file to check it's done
OPENW ,10, dir+'FINISH_FLAG'
PRINTF,10, "ftime       = '"+STRTRIM(STRING(ftime),2)+"'"
CLOSE ,10

END
