#!/bin/csh

cd /d/bip3/ezbc/perseus/data/cfa

fits in=DHT21_Taurus_raw.fits op=xyin out=DHT21_Taurus_raw.mir

reorder in=DHT21_Taurus_raw.mir out=DHT21_Taurus_raw.reorder.mir \
    mode=231

imsub in=DHT21_Taurus_raw.reorder.mir out=perseus.cfa.cube.mir \
    region='boxes(200,5,300,120)'

fits in=perseus.cfa.cube.mir op=xyout out=perseus.cfa.cube.fits


