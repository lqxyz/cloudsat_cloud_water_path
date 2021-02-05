## Process the CloudSat HDF files to get cloud water path netCDF dataset (T42)

CloudSat dataset is available at http://www.cloudsat.cira.colostate.edu/

You can download `2B-CWC-RO.P1_R05` dataset via the following bash script with your own `username` and `password`:
```bash
for year in {2012..2016}
do
    wget -r --user=username --password=password ftp://ftp.cloudsat.cira.colostate.edu//2B-CWC-RO.P1_R05/$year
done
```
The dataset should be located at `ftp.cloudsat.cira.colostate.edu/2B-CWC-RO.P1_R05/`, sorted by year.

### Required packages
Package [`pyhdf`](https://github.com/fhs/pyhdf) needs to be installed to read the HDF files, for example, you can install by `conda`:
```bash
conda install -c conda-forge pyhdf
```

### Run the scripts

For year 2012 to 2014, we only generate annual statiscal temperary datasets; while for year 2015 to 2016, we generate both monthly and yearly datasets.

For year 2012 to 2014, run:
```bash
./run_year_2012_2014.sh
python write_lwp_iwp_mon_year_2012_2014.py
```

For year 2015 to 2016, run:
```bash
./run_mon_2015_2016.sh
python write_lwp_iwp_mon_year_2015_2016.py
```

