import numpy as np

class peak:
	def __init__(self):
		self.nrun = None
		self.detnum = None
		self.h = None
		self.k = None
		self.l = None
		self.ipk = None
		self.peakID = None

def getPeaks(inFile):
	peakList = list()
	with open(inFile,'r') as f:
		r = f.readline()
		breakFlag = False
		while True:
			if r[0] == '1': #new detector
				nrun = int(r.split()[1])
				detnum = int(r.split()[2])
				r = f.readline() #2 - headers
				r = f.readline() #3 - first guy
				
				while r[0] == '3':
					thisPeak = peak()
					thisPeak.nrun = nrun
					thisPeak.detnum = detnum
					thisPeak.peakID = int(r.split()[1])
					thisPeak.h = int(r.split()[2])
					thisPeak.k = int(r.split()[3])
					thisPeak.l = int(r.split()[4])
					peakList.append(thisPeak)
					r = f.readline()
					if r == '':
						breakFlag = True
						break
				if breakFlag:
					break
			else:
				r = f.readline()
	return peakList
'''
inFile2016 = '/SNS/users/vel/Dropbox (ORNL)/first62.peaks'
peakList2016 = getPeaks(inFile2016)
inFile2017 = '/SNS/TOPAZ/IPTS-18474/shared/2017A/Si/TOPAZ_DetScan_2-2mmBN_use_DetCal_mantid_0p75/Si2mm_Cubic_F.integrate'
peakList2017 = getPeaks(inFile2017)
#Check for equality do this the lazy way
for p1 in peakList2017:
	for p2 in peakList2016:
		#if p1.h == p2.h and p1.k == p2.k and p1.l == p2.l and p1.detnum == p2.detnum:
		if (p1.h**2 + p1.k**2 + p1.l**2) == (p2.h**2 + p2.k**2 + p2.l**2) and (p1.detnum == p2.detnum):	
			if p1.peakID < 200:
				print p1.peakID, p2.peakID
#Write out easier to read version of these
with open('peakList2017.csv','w') as f:
	for p in peakList2017:
		r = '%i %i %i %i %i\n'%(p.peakID, p.detnum, p.h, p.k, p.l)
		f.write(r)
		
with open('peakList2016.csv','w') as f:
	for p in peakList2016:
		r = '%i %i %i %i %i\n'%(p.peakID, p.detnum, p.h, p.k, p.l)
		f.write(r)
print 'wrote files'	
'''
