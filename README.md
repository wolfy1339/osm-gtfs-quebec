# osm-gtfs-quebec

A series of Python functions to visualize RTC GTFS data and to convert to JOSM/OSM format for importing. It's a work in progress.


# **There currently is no permission from the RTC to use it's data, _this is for personal uses only_, <u>DO NOT UPLOAD TO OSM</u>**

## Requirements
- Python requirements should be installed automatically by your IDE or use `pip3 install -r requirements.txt` to install manually
- You need to have GDAL installed on your system

## How to run
1. Download the [desired (latest) gtfs archive](https://cdn.rtcquebec.ca/Site_Internet/DonneesOuvertes/googletransit.zip) and pass it as an input. The script is not meant to be ran by a cli at this point.
2. Look at the functions and run the ones you need.


## Notes
The script was previously specifically written having STL data in mind and currently works with data from the RTC.
