#!/usr/bin/env python

# script to parse a gmrt online flagging file and apply the flags using the
#flagdata2 task in casa.

# To use edit the base name of the visibility file and the flagfilename in
#this file.  Then run in casa with execfile('gmrtolflags.py')

# Note that I have not trapped the situation if day2 in the time range goes
#across a month boundary.  

months = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,
          'Oct':10,'Nov':11,'Dec':12}

# set parameters 
base = '16jan19_007'      # base name of the visibility measurement set (the part before .ms)
ms = base+'.ms'
flagfilename = '19_007.16jan.flag'    # name of the gmrt online flag file

flagfile = open(flagfilename,'r')

# parse the flag file
allLines = flagfile.readlines()

year = allLines[0].split()[7]
month = allLines[0].split()[4]
day = allLines[0].split()[5]

mon = str(months[month])

date = year + '/' + mon + '/' + day

print "\n Observing date is %s %s %s = %s" % (year, month, day, date)
for i in range(len(allLines)):
    test1 = allLines[i][0:3]
    if(test1 == 'ANT'):
#        print allLines[i]
        ant = allLines[i][9:11]
        d1 = allLines[i][22:24]
        t0h = allLines[i][25:27]
        t0m = allLines[i][28:30]
        t0s = allLines[i][31:33]
        d2 = allLines[i][34:36]
        t1h = allLines[i][37:39]
        t1m = allLines[i][40:42]
        t1s = allLines[i][43:45]
        trange = date+'/'+t0h+':'+t0m+':'+t0s+' ~ '+date+'/'+t1h+':'+t1m+':'+t1s
        if (d2 > d1):
            day2 = str(int(day)+1)
            date2 = year + '/' + mon + '/' + day2
            trange = date+'/'+t0h+':'+t0m+':'+t0s+' ~ '+date2+'/'+t1h+':'+t1m+':'+t1s           
        print "Flagging antenna %s: timerange %s" % (ant,trange)
        flagdata2(vis=ms, manualflag = true, mf_spw = '', mf_antenna = ant,
                  mf_timerange = trange, flagbackup = false, async = false)

