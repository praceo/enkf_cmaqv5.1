PRO enkf_00_set_time,ftime,yy,mm,dd,utc,kfint

print, 'start to find obs'

ecase = '50.enkf_vs_3dvar'

  dir    ='/home/spark/IDLWorkspace82/Default_'+ecase+'/'
  dmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  IF ( ((yy MOD 4) EQ 0) AND ((yy MOD 100) NE 0) OR ((yy MOD 400) EQ 0) ) THEN dmon[1]=29
 
  ymd_utc_en=STRING(yy ,FORMAT='(I4.4)')+STRING(mm ,FORMAT='(I2.2)')$
            +STRING(dd ,FORMAT='(I2.2)')+STRING(utc,FORMAT='(I2.2)')

 OPENW ,10, dir+'set_time.inc'
 PRINTF,10, "ymdh_utc    = '"+ymd_utc_en              +"'"
 PRINTF,10, "ftime       = '"+STRTRIM(STRING(ftime),2)+"'"
 CLOSE ,10

END
