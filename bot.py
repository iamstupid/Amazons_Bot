# An botzone-flavored python frontend for bot.cpp
import numpy as np
from ctypes import *
t=np.zeros((8,8),dtype=np.uint8)
uct=CDLL("./uct.so")
inits=[(0,2,1),(2,0,1),(5,0,1),(7,2,1),(0,5,2),(2,7,2),(5,7,2),(7,5,2)]
for i in inits:
	a,b,c=i
	t[a][b]=c
tr=int(input())
for i in range(1,2*tr):
	r=list(map(int,input().strip().split()))
	if r[0]==-1:
		continue
	a,b,c,d,e,f=r
	t[c][d]=t[a][b]
	t[a][b]=0
	t[e][f]=255
pos=[[],[]]
def decord(a):
	return a>>3,a&7
def encord(x,y):
	return x*8+y
def nzer(x):
	if x==0:
		return 0
	return 1
for i in range(8):
	for j in range(8):
		if t[i][j]>0 and t[i][j]<255:
			pos[t[i][j]-1].append(encord(i,j))
rs=np.array(pos,dtype=np.uint8).flatten()
ts=t.flatten()
tply=0
tbrd=0
for i in range(8):
	tply=tply|(rs[i]<<(8*i))
for i in range(64):
	tbrd=tbrd|(nzer(ts[i])<<i)
uct.mov.restype=c_uint
def int_to_bytes(value, length):
	result = []
	for i in range(0, length):
		result.append(value >> (i * 8) & 0xff)
		result.reverse()
	return bytes(result)
import random
prt=uct.mov(c_ulonglong(tbrd),c_ulonglong(tply),c_char_p(int_to_bytes(random.getrandbits(256),32)))
ar=prt&255
prt=prt>>8
to=prt&255
prt=prt>>8
fr=prt&255
a,b,c,d,e,f=np.array([decord(fr),decord(to),decord(ar)],dtype=np.uint8).flatten()
print(a,b,c,d,e,f)
