import os
import subprocess

abcs = ['111','311','331','611','631','661']
dust_szds = ['size1','size2','size3','size4','size5','size6','size7','size8','size9','size10',]
dust_ns = ['MgS-be','SiC-Pg']

"""
for dust_szd in dust_szds:
    for abc in abcs:
	for dust_n in dust_ns:
	    a, b, c = abc
	    os.system('mv {0}t{1}{2}{3}s{4}.rfi {0}{1}{2}{3}s{4}.rfi'.format(dust_n,a,b,c,dust_szd[-1]))
		
"""
"""
with open('make_opc.in','a') as f:
    f.write('compile grains "graphite.rfi" "small2.szd" 10\n')
    f.write('compile grains "be1-amcarb.rfi" "big2.szd" 10\n')
    for dust_szd in dust_szds:
	for abc in abcs:
	    for dust_n in dust_ns:
		a, b, c = abc
		f.write('compile grains "{}{}{}{}s{}.rfi" "{}.szd" 01\n'.format(dust_n,a,b,c,dust_szd[-1],dust_szd))
"""
with open('make_opc.in','r') as f:
    for command in f.readlines():
	subprocess.call(["echo '{}'| /usr/local/Cloudy/c13.03/source/cloudy.exe".format(command)], shell=True)
