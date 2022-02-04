#!/bin/sh
echo new run name?
read newrun

mkdir $newrun

cp ~/MWF/install/* ./$newrun
rm ./$newrun/new32.job
cp ~/MWF/randIC.out ./$newrun

echo done
