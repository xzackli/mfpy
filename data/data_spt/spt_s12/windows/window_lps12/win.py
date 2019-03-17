from numpy import *

l,win1 = loadtxt('window_1',unpack=True)
l,win2 = loadtxt('window_2',unpack=True)
l,win3 = loadtxt('window_3',unpack=True)
l,win4 = loadtxt('window_4',unpack=True)
l,win5 = loadtxt('window_5',unpack=True)
l,win6 = loadtxt('window_6',unpack=True)
l,win7 = loadtxt('window_7',unpack=True)
l,win8 = loadtxt('window_8',unpack=True)
l,win9 = loadtxt('window_9',unpack=True)
l,win10 = loadtxt('window_10',unpack=True)
l,win11 = loadtxt('window_11',unpack=True)
l,win12 = loadtxt('window_12',unpack=True)
l,win13 = loadtxt('window_13',unpack=True)
l,win14 = loadtxt('window_14',unpack=True)
l,win15 = loadtxt('window_15',unpack=True)
l,win16 = loadtxt('window_16',unpack=True)
l,win17 = loadtxt('window_17',unpack=True)
l,win18 = loadtxt('window_18',unpack=True)
l,win19 = loadtxt('window_19',unpack=True)
l,win20 = loadtxt('window_20',unpack=True)
l,win21 = loadtxt('window_21',unpack=True)
l,win22 = loadtxt('window_22',unpack=True)
l,win23 = loadtxt('window_23',unpack=True)
l,win24 = loadtxt('window_24',unpack=True)
l,win25 = loadtxt('window_25',unpack=True)
l,win26 = loadtxt('window_26',unpack=True)
l,win27 = loadtxt('window_27',unpack=True)
l,win28 = loadtxt('window_28',unpack=True)
l,win29 = loadtxt('window_29',unpack=True)
l,win30 = loadtxt('window_30',unpack=True)
l,win31 = loadtxt('window_31',unpack=True)
l,win32 = loadtxt('window_32',unpack=True)
l,win33 = loadtxt('window_33',unpack=True)
l,win34 = loadtxt('window_34',unpack=True)
l,win35 = loadtxt('window_35',unpack=True)
l,win36 = loadtxt('window_36',unpack=True)
l,win37 = loadtxt('window_37',unpack=True)
l,win38 = loadtxt('window_38',unpack=True)
l,win39 = loadtxt('window_39',unpack=True)
l,win40 = loadtxt('window_40',unpack=True)
l,win41 = loadtxt('window_41',unpack=True)
l,win42 = loadtxt('window_42',unpack=True)
l,win43 = loadtxt('window_43',unpack=True)
l,win44 = loadtxt('window_44',unpack=True)
l,win45 = loadtxt('window_45',unpack=True)
l,win46 = loadtxt('window_46',unpack=True)
l,win47 = loadtxt('window_47',unpack=True)






win=matrix(zeros([3251,48]))
for i in range(3251):
	win[i,0]=l[i]
	win[i,1]=win1[i]
	win[i,2]=win2[i]
	win[i,3]=win3[i]
	win[i,4]=win4[i]
	win[i,5]=win5[i]
	win[i,6]=win6[i]
	win[i,7]=win7[i]
	win[i,8]=win8[i]
	win[i,9]=win9[i]
	win[i,10]=win10[i]
	win[i,11]=win11[i]
	win[i,12]=win12[i]
	win[i,13]=win13[i]
	win[i,14]=win14[i]
	win[i,15]=win15[i]
        win[i,16]=win16[i]
        win[i,17]=win17[i]
        win[i,18]=win18[i]
        win[i,19]=win19[i]
        win[i,20]=win20[i]
        win[i,21]=win21[i]
        win[i,22]=win22[i]
        win[i,23]=win23[i]
        win[i,24]=win24[i]
        win[i,25]=win25[i]
        win[i,26]=win26[i]
        win[i,27]=win27[i]
        win[i,28]=win28[i]
        win[i,29]=win29[i]
        win[i,30]=win30[i]
        win[i,31]=win31[i]
        win[i,32]=win32[i]
        win[i,33]=win33[i]
        win[i,34]=win34[i]
        win[i,35]=win35[i]
        win[i,36]=win36[i]
        win[i,37]=win37[i]
        win[i,38]=win38[i]
        win[i,39]=win39[i]
        win[i,40]=win40[i]
        win[i,41]=win41[i]
        win[i,42]=win42[i]
        win[i,43]=win43[i]
        win[i,44]=win44[i]
        win[i,45]=win45[i]
        win[i,46]=win46[i]
        win[i,47]=win47[i]

savetxt('BblMean_150x150.dat', win)	
