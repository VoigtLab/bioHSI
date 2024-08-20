#!/bin/bash

mkdir -p 00_data
cd 00_data
curl -L -o uav_flight_hsi.zip "https://www.dropbox.com/scl/fi/clsudhx60z4vaqd9x0aoe/uav_flight_hsi.zip?rlkey=soluuzol3iihpap9emxy3w6yi&st=jkxlsazd&dl=1"
unzip uav_flight_hsi.zip
cd ..

