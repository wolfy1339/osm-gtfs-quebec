# osm-gtfs-quebec

A series of Python functions to visualize RTC GTFS data and to convert to JOSM/OSM format for importing. It's a work in progress.


## **DO NOT UPLOAD GENERATED DATA TO OSM, generated data is for personal uses only. The RTC's data is licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/), which isn't fully compatible with the [OBdL](https://opendatacommons.org/licenses/odbl/summary/) <small>[1]**
[1] https://blog.openstreetmap.org/2017/03/17/use-of-cc-by-data/

## Requirements
- Python requirements should be installed automatically by your IDE or use `pip3 install -r requirements.txt` to install manually
- You need to have GDAL installed on your system

## How to run
1. Download the [desired (latest) gtfs archive](https://cdn.rtcquebec.ca/Site_Internet/DonneesOuvertes/googletransit.zip) and pass it as an input. The script is not meant to be ran by a cli at this point.
2. Look at the functions and run the ones you need.


## Notes
- The script was previously specifically written having STL data in mind and currently works with data from the RTC.
- RTC's GTFS data is licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/)
