#!/bin/bash

ip_address=`hostname -I`
ips=($ip_address)
ip=${ips[0]}
datestamp=`date +%Y%m%d-%T`
log=versions-run-${ip}-${datestamp}.txt

echo "# Git" > $log
git log|head -n 3 >> $log

echo "# FastQC" >> $log
./FastQC/fastqc --version >> $log

echo "# PEAR" >> $log
ls pear-* >> $log

echo "# BWA" >> $log
ls | grep bwa >> $log

echo "# Samtools" >> $log
samtools 2>&1 >/dev/null | grep 'Version' >> $log

echo "# Picard tools" >> $log
ls |grep picard >> $log

echo "# Varscan" >> $log
ls VarScan* >> $log

echo "# Python" >> $log
python --version 2>&1 >/dev/null | cat >> $log

echo "# SQLite" >> $log
sqlite3 -version >> $log
