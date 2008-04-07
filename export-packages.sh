#!/bin/sh

uname_a=`uname -a`
repos=$SSS
rm -rf /tmp/invearthquake_export
mkdir /tmp/invearthquake_export
cd /tmp/invearthquake_export
svn export $repos/kinherd/invearthquake/trunk invearthquake > exportlog || exit 1
tail -n 1 exportlog > invearthquake/VERSION
date >> invearthquake/VERSION
version=`cat exportlog | tail -n 1 | tr -d -c '[0-9]'`
svn log $repos/kinherd/invearthquake/trunk > invearthquake/CHANGES
mv invearthquake invearthquake-r$version
tar -czvf invearthquake-r$version.tar.gz invearthquake-r$version
mv invearthquake-r$version invearthquake
cd invearthquake
machine=`./hostinfo.pl --machine`
os=`./hostinfo.pl --os`

cat > Makefile.local <<EOF
FORTRANC = g95
CFLAGS += \$(CFAST)
EOF

make || exit 1

cat > bin/README <<EOF
This directory contains the Kinherd Invearthquake programs.

Compiled for:
  $os $machine

Documentation:
  http://kinherd.org

Compilation host:
  $uname_a

Revision:
  r$version 
EOF
cp VERSION bin/
cp CHANGES bin/

cp LICENSE bin/
cp NOTICE bin/

vom=r$version-$os-$machine
cp -R bin invearthquake-$vom-bin
tar -czvf invearthquake-$vom-bin.tar.gz invearthquake-$vom-bin
mv invearthquake-$vom-bin.tar.gz ..
cd ..

echo "Everything OK? Shall I copy the packages to the server? [y/N]"
read answer

if [ $answer == 'y' ]; then
    scp invearthquake-r$version.tar.gz invearthquake-$vom-bin.tar.gz emolch:www/kinherd/software || exit 1
    echo "removing build files..."
    rm -rf /tmp/invearthquake_export
else 
    echo "I will leave the files in /tmp/invearthquake_export"
fi
