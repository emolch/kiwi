#!/bin/sh

uname_a=`uname -a`
repos=$SSS
rm -rf /tmp/kiwi_export
mkdir /tmp/kiwi_export
cd /tmp/kiwi_export
svn export $repos/kiwi/trunk kiwi > exportlog || exit 1
tail -n 1 exportlog > kiwi/VERSION
date >> kiwi/VERSION
version=`cat exportlog | tail -n 1 | tr -d -c '[0-9]'`
svn log $repos/kiwi/trunk > kiwi/CHANGES
mv kiwi kiwi-r$version
tar -czvf kiwi-r$version.tar.gz kiwi-r$version

echo "Everything OK? Shall I copy the packages to the server? [y/N]"
read answer

if [ $answer == 'y' ]; then
    scp kiwi-r$version.tar.gz emolch:www/kinherd/software/ || exit 1
    echo "removing exported files..."
    rm -rf /tmp/kiwi_export
else 
    echo "I will leave the files in /tmp/kiwi_export"
fi
