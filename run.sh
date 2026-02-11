#!/bin/bash

#fofn-checker

CASES=$1
#2025_cases.csv
DATES=$2
#../2025_UNIX-start-end.csv
OUTFILE=$3

# generate random prefix for all tmp files
RAND_1=`echo $((1 + RANDOM % 100))`
RAND_2=`echo $((100 + RANDOM % 200))`
RAND_3=`echo $((200 + RANDOM % 300))`
RAND=`echo "${RAND_1}${RAND_2}${RAND_3}"`

echo "START_unix_time,END_unix_time,START_date,END_date,Exposure day (4.5 months before window center),Cases in treatment zone,Cases in control zone,Control minus treatment,TOTAL cases,IRR,95% CI LOWER,95% CI UPPER,PVAL" > ${OUTFILE}


for LINE in $(seq 1 152); do

	START=`head -${LINE} ${DATES} | tail -1 | cut -f 1 -d ','`
	END=`head -${LINE} ${DATES} | tail -1 | cut -f 2 -d ','`
	START_DATE=`head -${LINE} ${DATES} | tail -1 | cut -f 3 -d ','`
	END_DATE=`head -${LINE} ${DATES} | tail -1 | cut -f 4 -d ','`
	EXP_DAY=`head -${LINE} ${DATES} | tail -1 | cut -f 5 -d ','`


	
	python 650m_zone_counts.py ${CASES} Treatment_lat_lon.csv Control_lat_lon.csv ${START} ${END} > ${RAND}_tmp.csv
	IN_T=`grep "Cases inside treatment zone" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	IN_C=`grep "Cases inside control zone" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	DIFF=`grep "Control-minus-treatment case count difference:" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	TOTAL=`grep "Total unique cases in window" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	
	echo ${IN_T},${IN_C},${TOTAL}

	python PRT.py ${IN_T} ${IN_C} > ${RAND}_tmp.csv
	LRT=`grep "LRT statistic" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	PVAL=`grep "Two-sided p-value" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	IRR=`grep "IRR (treat/control)" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' '`
	LOWER=`grep "95% CI for IRR (score)" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' ' | tr -d '(' | tr -d ')' | tr ',' '-' | cut -f 1 -d '-'`
	UPPER=`grep "95% CI for IRR (score)" ${RAND}_tmp.csv | cut -f 2 -d ':' | tr -d ' ' | tr -d '(' | tr -d ')' | tr ',' '-' | cut -f 2 -d '-'`

	echo "${START},${END},${START_DATE},${END_DATE},${EXP_DAY},${IN_T},${IN_C},${DIFF},${TOTAL},${IRR},${LOWER},${UPPER},${PVAL}" >> ${OUTFILE}

done

		