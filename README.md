## Processing the CloudSat cloud water path dataset

CloudSat dataset is available at http://www.cloudsat.cira.colostate.edu/

You can download `2B-CWC-RO.P1_R05` dataset via the following bash script with your own `username` and `password`:
```bash
for year in {2012..2016}
do
    wget -r --user=username --password=password ftp://ftp.cloudsat.cira.colostate.edu//2B-CWC-RO.P1_R05/$year
done
```

