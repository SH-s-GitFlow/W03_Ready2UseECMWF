:tada: Feat: "W03_Ready2UseECMWF"

## 내용
Make ECMWF forecast results ready to use. The code is applied to every subfolders in --folder.

## 실행 명령어
```sh
python ECMWF_API.py --save <savepath> --folder <saved grib path from ECMWF_API.py> --geojson <geojson file with abs path> --maxhour <maximum forecast hour>
python D:\Oilspill_Total\Code_Python\Wave_spectra\KCGSA_Codes\2_Ready2use_ECMWF.py 
                -s D:\Extract_spectra\ECMWF_forecast_ready -g "D:/TROPOMI/API/map.geojson" -f   # 예시D:\Extract_spectra\ECMWF_forecast -m 48
```


Resolves: #2
