#!/usr/bin/python
'''Need to rename HACKS into a constant framework, so we can drag and drop
in to srclists for future versions of '''
import subprocess

class hack_data():
	def __init__(self):
		self.cat_name = ''
		self.cat_index = ''
		self.hack_name = ''
		self.source_name = ''
		self.original_name = ''


simple_cats = "/home/jline/Documents/cataloguing/RTS/extended/simple_cats"
mrc_data = open('%s/mrc_simple.txt' %simple_cats).read().split('\n')
#sumss_data = open('%s/sumss_simple.txt' %simple_cats).read().split('\n')
#vlssr_data = open('%s/vlssr_simple.txt' %simple_cats).read().split('\n')
#nvss_data = open('%s/nvss_simple.txt' %simple_cats).read().split('\n')
#A154_data = open('%s/A154_simple.txt' %simple_cats).read().split('\n')
#A182_data = open('%s/A182_simple.txt' %simple_cats).read().split('\n')
mwacs_data = open('%s/mwacs_simple.txt' %simple_cats).read().split('\n')

mwacs_names = [line.split()[1] for line in mwacs_data]
mrc_names = [line.split()[1] for line in mrc_data]

hack_info = open('puma-v1_hacklog.txt','r').read().split('\n')

HACKS = []

for line in hack_info:
	if '(' in line and ')' in line and 'HACK' in line:
		HACK = hack_data()
		info = line.split()
		source_name = info[2]
		HACK.original_name = source_name
		if source_name in mwacs_names:
			HACK.cat_index = mwacs_names.index(source_name)
			HACK.cat_name = 'mwacs'
		elif source_name in mrc_names:
			HACK.cat_index = mrc_names.index(source_name)
			HACK.cat_name = 'mrc'
		
		for piece in info:
			if 'HACK' in piece:
				HACK.hack_name = piece
		HACKS.append(HACK)

hack_sources = open('srclist_hack_v1.txt','r').read().split('ENDSOURCE')
hack_sources = hack_sources[:-1]

#for HACK in HACKS:
	#subprocess.call('cp /home/jline/Documents/cataloguing/RTS/extended/puma_v1/%s_full.png EXT%s%s_full.png' %(HACK.hack_name,HACK.cat_name,HACK.cat_index),shell=True)
	#subprocess.call('cp /home/jline/Documents/cataloguing/RTS/extended/puma_v1/%s.png EXT%s%s.png' %(HACK.hack_name,HACK.cat_name,HACK.cat_index),shell=True)
	
new_srclist = open("srclist_hack_v1_renamed.txt",'w+')
	
for src in hack_sources:
	src_lines = src.split('\n')
	if src_lines[0]=='': del src_lines[0]
	if src_lines[-1]=='': del src_lines[-1]
	src_line1 = src_lines[0].split()
	src_name = src_line1[1]
	
	for HACK in HACKS:
		if HACK.hack_name==src_name: 
			print HACK.hack_name,src_name,HACK.cat_name,HACK.cat_index,HACK.original_name
			src_line1[1] = 'EXT%s%s' %(HACK.cat_name,HACK.cat_index)
	
	new_srclist.write("\n%s %s %s %s" %(src_line1[0],src_line1[1],src_line1[2],src_line1[3]))
		
	for line in src_lines[1:]:
		new_srclist.write('\n'+line)
	new_srclist.write("\nENDSOURCE")

new_srclist.close()	
	
	




