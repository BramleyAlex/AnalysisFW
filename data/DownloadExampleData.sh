#!/bin/bash

# Define remote path
remote="https://cernbox.cern.ch/remote.php/dav/public-files/FlfgWyzeyTu9MSC"

# Dataset folders
Data="user.dbaronmo.v26.data18_13TeV.periodI.physics_Main.M4.grp23_v01_p4238.sv1_Le"
Wplus="user.dbaronmo.v26.mc.361102.PoPy8_Wplustaunu.M4.e3601_s3126_r10724_p4512.sv1_Le"
Wminus="user.dbaronmo.v26.mc.361105.PoPy8_Wminustaunu.M4.e3601_s3126_r10724_p4512.sv1_Le"
Ztautau="user.dbaronmo.v26.mc.361108.PoPy8_Ztt.M4.e3601_s3126_r9364_p4512.sv1_Le"

# Download datasets
mkdir $Data && cd $Data
wget $remote/$Data/user.dbaronmo.25819190._000001.LepUniv_ttbar.root
wget $remote/$Data/user.dbaronmo.25819190._000002.LepUniv_ttbar.root
wget $remote/$Data/user.dbaronmo.25819190._000003.LepUniv_ttbar.root
wget $remote/$Data/user.dbaronmo.25819190._000004.LepUniv_ttbar.root
wget $remote/$Data/user.dbaronmo.25819190._000005.LepUniv_ttbar.root
wget $remote/$Data/user.dbaronmo.25819190._000006.LepUniv_ttbar.root
cd ..

mkdir $Wplus && cd $Wplus
wget $remote/$Wplus/user.dbaronmo.25819141._000002.LepUniv_ttbar.root
wget $remote/$Wplus/user.dbaronmo.25819141._000003.LepUniv_ttbar.root
cd ..

mkdir $Wminus && cd $Wminus
wget $remote/$Wminus/user.dbaronmo.25819146._000001.LepUniv_ttbar.root
wget $remote/$Wminus/user.dbaronmo.25819146._000002.LepUniv_ttbar.root
cd ..

mkdir $Ztautau && cd $Ztautau
wget $remote/$Ztautau/user.dbaronmo.25819000._000001.LepUniv_ttbar.root
wget $remote/$Ztautau/user.dbaronmo.25819000._000002.LepUniv_ttbar.root
wget $remote/$Ztautau/user.dbaronmo.25819000._000004.LepUniv_ttbar.root
cd ..

echo "Download complete!"





