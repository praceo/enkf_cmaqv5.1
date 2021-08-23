PRO enkf_15_gen_obs_kr

  @set_time.inc
  COMMON SHARE_OBS_KR, nobs_kr, obsvalues_kr, obserrors_kr,   $
                       H_mat_kr, obs_i_kr, obs_j_kr


comp = 'PM25'
ecase = '50.enkf_vs_3dvar'

dir    ='/home/spark/IDLWorkspace82/Default_'+ecase+'/'
dir_obs='/home/spark/cmaqv5.2.1/DATA/'+ecase+'/obs/'
dir_grd='/home/spark/cmaqv5.2.1/DATA/'+ecase+'/mcip/'

fobs = dir_obs+'PM25_'+ymdh_utc+'_UTC.txt'
fgrd = dir_grd+'GRIDCRO2D_20160425_000000.nc'

idgrd = NCDF_OPEN(fgrd,/NOWRITE)
NCDF_VARGET, idgrd, 'LON', lon
NCDF_VARGET, idgrd, 'LAT', lat
NCDF_CLOSE, idgrd

nx =144L & ny =105L & nz = 14L
nv = 1 ;; number of state       variables, which are PM25 i.c. and PSO4 em. now
nvo= 1 ;; number of observation variables, which is  PM25 only now
vec_len = nx*ny*nz*nv

intp_H = 2 ;; interpolation methods for obs. operator (H)
           ;; 1 = 4-pnts weighting average
           ;; 2 = 4-pnts linear interpolation

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; SAVE THE OBSERVATION DATA EXCLUDING MISSING DATA
cnt=0
OPENR,20, fobs
lbyl = ''
WHILE ~EOF(20) DO BEGIN
  READF,20, lbyl
  splt_lbyl = STRSPLIT(lbyl,' ',/EXTRACT)
  IF ( ( splt_lbyl[0] EQ 'KOR'                         )  AND $
       ( splt_lbyl[3] GE 103.0                         )  AND $ ;DEFINE
       ( splt_lbyl[2] GE 21.0 AND splt_lbyl[2] LE 44.0 )  AND $ ;DEFINE
       ( splt_lbyl[4] NE '-999'                        )  AND $
       ( splt_lbyl[4] NE 0                             ) )THEN BEGIN
    cnt=cnt+1
    IF (cnt EQ 1) THEN BEGIN
      obs_cod=[splt_lbyl[1]] 
      obs_lon=[splt_lbyl[3]] & obs_lat=[splt_lbyl[2]]
      obs_val=[FLOAT(splt_lbyl[4])] & obs_err=[FLOAT(splt_lbyl[5])]
    ENDIF ELSE BEGIN        
      obs_cod=[obs_cod,splt_lbyl[1]]
      obs_lon=[obs_lon,splt_lbyl[3]] & obs_lat=[obs_lat,splt_lbyl[2]]
      obs_val=[obs_val,FLOAT(splt_lbyl[4])] & obs_err=[obs_err,FLOAT(splt_lbyl[5])]
    ENDELSE
  ENDIF
ENDWHILE
CLOSE,20

nsite=cnt

;; Nobs = number of observations = number of sites (nsite) * number of obs. variables (nvo)
nobs_kr = nsite * nvo

IF (nobs_kr EQ 0) THEN BEGIN
  PRINT, 'NO OBSERVATION DATA ARE AVAILABLE AT THIS TIME!'
  PRINT, 'nobs_kr= ',nobs_kr
  PRINT, 'GEN_OBS STEP WILL BE SKIPPED.'
  GOTO, NO_OBS
ENDIF

;; SORT SAVED OBS. DATA WITH INCREASING STATION NUMBER AND CODE
 obs_cod_srt_kr=obs_cod[SORT(obs_cod)]
 obs_lon_srt_kr=obs_lon[SORT(obs_cod)]
 obs_lat_srt_kr=obs_lat[SORT(obs_cod)]
 obs_val_srt_kr=obs_val[SORT(obs_cod)]
 obs_err_srt_kr=obs_err[SORT(obs_cod)]

;; SKIP SORTING
;obs_cod_srt_kr=obs_cod
;obs_lon_srt_kr=obs_lon
;obs_lat_srt_kr=obs_lat
;obs_val_srt_kr=obs_val
;obs_err_srt_kr=obs_err

obs_gri_kr    =INTARR(5,nobs_kr) ; i grid index [ [pt0 ~ pt4] , [0:nsite] ]
obs_grj_kr    =INTARR(5,nobs_kr) ; j grid index [ [pt0 ~ pt4] , [0:nsite] ]
obs_i_kr      =FLTARR(  nobs_kr)   ; x grid location of obs. relative to ij grid index 
obs_j_kr      =FLTARR(  nobs_kr)   ; y grid location of obs. relative to ij grid index

;; Observation operator, linear interpolation using 4 surrounding grid point in CMAQ
H = DBLARR(nx,ny,nz,nv,nobs_kr)

FOR o=0,nobs_kr-1 DO BEGIN
FOR v=0,nv-1 DO BEGIN

  dis = SQRT((lon-obs_lon_srt_kr[o])^2.0 + (lat-obs_lat_srt_kr[o])^2.0)
  pt0=array_indices(dis,WHERE(dis EQ MIN(dis)))

  IF ( (obs_lon_srt_kr[o] GE lon[pt0[0],pt0[1]]) AND $          ;;  pt1 -----> pt2
       (obs_lat_srt_kr[o] GT lat[pt0[0],pt0[1]]) ) THEN BEGIN   ;;              |
    pt1=[pt0[0]  ,pt0[1]+1] & pt2=[pt0[0]+1,pt0[1]+1]           ;;              V
    pt3=[pt0[0]+1,pt0[1]  ] & pt4=[pt0              ]           ;;  pt4 <----- pt3
  ENDIF
  IF ( (obs_lon_srt_kr[o] GT lon[pt0[0],pt0[1]]) AND $
       (obs_lat_srt_kr[o] LE lat[pt0[0],pt0[1]]) ) THEN BEGIN
    pt1=[pt0              ] & pt2=[pt0[0]+1,pt0[1] ]
    pt3=[pt0[0]+1,pt0[1]-1] & pt4=[pt0[0]  ,pt0[1]-1]
  ENDIF
  IF ( (obs_lon_srt_kr[o] LE lon[pt0[0],pt0[1]]) AND $
       (obs_lat_srt_kr[o] LT lat[pt0[0],pt0[1]]) ) THEN BEGIN
    pt1=[pt0[0]-1,pt0[1]   ] & pt2=[pt0              ]
    pt3=[pt0[0]  ,pt0[1]-1 ] & pt4=[pt0[0]-1,pt0[1]-1]
  ENDIF
  IF ( (obs_lon_srt_kr[o] LT lon[pt0[0],pt0[1]]) AND $
       (obs_lat_srt_kr[o] GE lat[pt0[0],pt0[1]]) ) THEN BEGIN
    pt1=[pt0[0]-1,pt0[1]+1] & pt2=[pt0[0]  ,pt0[1]+1]
    pt3=[pt0              ] & pt4=[pt0[0]-1,pt0[1]  ]
  ENDIF

  ;; Save the indexes for the closest grid and surrounding grids
  obs_gri_kr[0,o]=pt0[0] & obs_grj_kr[0,o]=pt0[1]
  obs_gri_kr[1,o]=pt1[0] & obs_grj_kr[1,o]=pt1[1]
  obs_gri_kr[2,o]=pt2[0] & obs_grj_kr[2,o]=pt2[1]
  obs_gri_kr[3,o]=pt3[0] & obs_grj_kr[3,o]=pt3[1]
  obs_gri_kr[4,o]=pt4[0] & obs_grj_kr[4,o]=pt4[1]

  ;; Calculate the relative grid location of observation site
  obs_i_kr[o] = pt4[0] + ( (obs_lon_srt_kr[o]  - lon[pt4[0],pt4[1]]) $
                         / (lon[pt3[0],pt3[1]] - lon[pt4[0],pt4[1]]) )
  obs_j_kr[o] = pt4[1] + ( (obs_lat_srt_kr[o]  - lat[pt4[0],pt4[1]]) $
                         / (lat[pt1[0],pt1[1]] - lat[pt4[0],pt4[1]]) )
  ENDFOR
ENDFOR


CASE intp_H   OF
  1: BEGIN
    FOR o=0,nobs_kr-1 DO BEGIN
    FOR v=0,nv-1 DO BEGIN
    ;; Distance between obs. site and 4 surrounding grids by means of latlon coordinate
    df = SQRT((lon[obs_gri_kr[4,o],obs_grj_kr[4,o]]-lon[obs_gri_kr[3,o],obs_grj_kr[3,o]])^2.0 $
             +(lat[obs_gri_kr[4,o],obs_grj_kr[4,o]]-lat[obs_gri_kr[3,o],obs_grj_kr[3,o]])^2.0)
    d1 = SQRT((lon[obs_gri_kr[1,o],obs_grj_kr[1,o]]-obs_lon_srt_kr[o])^2.0 $
             +(lat[obs_gri_kr[1,o],obs_grj_kr[1,o]]-obs_lat_srt_kr[o])^2.0)
    d2 = SQRT((lon[obs_gri_kr[2,o],obs_grj_kr[2,o]]-obs_lon_srt_kr[o])^2.0 $
             +(lat[obs_gri_kr[2,o],obs_grj_kr[2,o]]-obs_lat_srt_kr[o])^2.0)
    d3 = SQRT((lon[obs_gri_kr[3,o],obs_grj_kr[3,o]]-obs_lon_srt_kr[o])^2.0 $
             +(lat[obs_gri_kr[3,o],obs_grj_kr[3,o]]-obs_lat_srt_kr[o])^2.0)
    d4 = SQRT((lon[obs_gri_kr[4,o],obs_grj_kr[4,o]]-obs_lon_srt_kr[o])^2.0 $
             +(lat[obs_gri_kr[4,o],obs_grj_kr[4,o]]-obs_lat_srt_kr[o])^2.0)
    w1 = MAX([0.0D,(df-d1)]) & w2 = MAX([0.0D,(df-d2)])
    w3 = MAX([0.0D,(df-d3)]) & w4 = MAX([0.0D,(df-d4)])
    wf = w1 + w2 + w3 + w4
    ; w1/wf, w2/wf, w2/wf, w4/wf will be used for obs. operator
    H[obs_gri_kr[1,o],obs_grj_kr[1,o],0,v,o] = w1/wf
    H[obs_gri_kr[2,o],obs_grj_kr[2,o],0,v,o] = w2/wf
    H[obs_gri_kr[3,o],obs_grj_kr[3,o],0,v,o] = w3/wf
    H[obs_gri_kr[4,o],obs_grj_kr[4,o],0,v,o] = w4/wf
    ENDFOR
    ENDFOR
  END
  2: BEGIN
    FOR o=0,nobs_kr-1 DO BEGIN
    FOR v=0,nv-1 DO BEGIN
    delxp = obs_gri_kr[3,o] - obs_i_kr[o]
    delx  = 1.0 - delxp
    delyp = obs_grj_kr[1,o] - obs_j_kr[o]
    dely  = 1.0 - delyp
    H[obs_gri_kr[1,o],obs_grj_kr[1,o],0,v,o] = delxp * dely
    H[obs_gri_kr[2,o],obs_grj_kr[2,o],0,v,o] = delx  * dely
    H[obs_gri_kr[3,o],obs_grj_kr[3,o],0,v,o] = delx  * delyp
    H[obs_gri_kr[4,o],obs_grj_kr[4,o],0,v,o] = delxp * delyp
    ENDFOR
    ENDFOR
  END
  ELSE: BEGIN
    PRINT, 'The interpolation methods for H should be 1 or 2.'
    EXIT
  END
ENDCASE

H_mat_kr = REFORM(H,vec_len,nobs_kr)

obsvalues_kr = obs_val_srt_kr
obserrors_kr = obs_err_srt_kr

SET_PLOT, 'Z'
DEVICE, DECOMPOSED=0
DEVICE, SET_RESOLUTION=[600,400], SET_PIXEL_DEPTH=24

LOADCT,39

;WINDOW,10,xs=600,ys=400
PLOT, obsvalues_kr, psym=4
img1= tvrd(/true)
write_png, dir+'obs_kr.png', img1

DEVICE, SET_RESOLUTION=[1000,600], SET_PIXEL_DEPTH=24

;WINDOW,20,xs=1000,ys=600
PLOT, lon,lat, xrange=[min(lon),max(lon)] $
    , yrange=[min(lat),max(lat)],psym=3,/iso
OPLOT, obs_lon_srt_kr,obs_lat_srt_kr,psym=4, symsize=0.7,color=100
img2= tvrd(/true)
write_png, dir+'site_kr.png', img2

DEVICE, /CLOSE
SET_PLOT, 'X'

;; Save the assimilated observation data
OPENW,100, dir+'OBS_DATA_KR.txt'
PRINTF,100,'#','code','i','j','obs_i_kr','obs_j_kr'         $
              ,'lon','lat','obs[ug]y','err[ug]'             $
              ,'H_pt1','H_pt2','H_pt3','H_pt4'              $
              ,FORMAT='(A4,A7,2A5,10A10)'
FOR o=0,nobs_kr-1 DO PRINTF,100,o+1,obs_cod_srt_kr[o],obs_gri_kr[0,o],obs_grj_kr[0,o]      $
                               ,obs_i_kr[o],obs_j_kr[o],obs_lon_srt_kr[o],obs_lat_srt_kr[o]$
                               ,obsvalues_kr[o],obserrors_kr[o]                            $
                               ,H[obs_gri_kr[1,o],obs_grj_kr[1,o],0,0,o]                   $
                               ,H[obs_gri_kr[2,o],obs_grj_kr[2,o],0,0,o]                   $
                               ,H[obs_gri_kr[3,o],obs_grj_kr[3,o],0,0,o]                   $
                               ,H[obs_gri_kr[4,o],obs_grj_kr[4,o],0,0,o]                   $
                           , FORMAT='(I4,1X,A7,2I5,10F10.5)'
PRINT, 'In Korea,'
PRINT, 'Assimilated observation data is saved.'
CLOSE,100

IF (nobs_kr EQ 0) THEN BEGIN
  NO_OBS : PRINT,'enkf_15.gen_obs_kr is done'
ENDIF

END
