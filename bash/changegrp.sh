#!/bin/bash

GRP="$1"
GRP=${GRP:-youngn}
export HOME=/home/$GRP/zhoux379
export home=$HOME
export data=$home/data
export genome=$home/data/genome
export misc1=$data/misc1
export misc2=$data/misc2
export misc3=$data/misc3
export misc4=$data/misc4
cd
pwd
newgrp $GRP
